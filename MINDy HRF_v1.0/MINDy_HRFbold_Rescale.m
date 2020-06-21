function[Out]=MINDy_HRFbold_Rescale(X,Pre,ParStr,Res,ReScale,varargin)
MINDy_type='HRF'; %#ok<NASGU>
%% Res Allows adding custom resolution to HRFpoly estimation
%% Multiplies Outgoing HRF by Rescale
%% 
dH=.0001;
%% Remark: mu is just a placeholder
mu=1;
doRobust='y';

derivKern=[1 -1];

if isempty(varargin)
    tLength=30;
else
    tLength=varargin{1};
end

HRFout=MINDy_MakeHRFpoly_Res(X,Pre,ParStr,'n',Res,tLength);
Out.HRFout=HRFout;
DeconvLength=HRFout.DeconvLength; %#ok<NASGU>
if iscell(X)
    dY0=cell(1,numel(X));
else
    dY0=cell(1);
end

xHRFpoly=HRFout.Poly; 
DxHRFpoly=HRFout.dPoly;
PolyVec=HRFout.PolyMat;

nX=size(HRFout.Poly,1);


tVec=[0 Pre.TR*(1:HRFout.DeconvLength)];

a2H=16;
b2H=1;
cH=1/6;
%h=@(h1,h2,t)((t.^(h1-1).*exp(-h2.*t).*(h2.^h1))./gamma(h1)-cH*((t.^(a2H-1).*exp(-b2H.*t).*(b2H.^a2H))./gamma(a2H)));

Hfix=-ReScale*cH*((tVec.^(a2H-1).*exp(-b2H.*tVec).*(b2H.^a2H))./gamma(a2H));

h=@(h1,h2)(ReScale*(tVec.^(h1-1).*exp(-h2.*tVec).*((h2.^h1)./gamma(h1))));

%lnT=log(tVec);lnT(tVec==0)=0;

Y0=HRFout.ReMean;
for iY=1:numel(HRFout.ReMean)
    tTmp=size(HRFout.OrigTime{iY},2);    
    Y0{iY}=[Y0{iY} zeros(nX,tTmp-size(Y0{iY},2))];
    dY0{iY}=[convn(Y0{iY},derivKern,'valid') zeros(nX,1)];
end
Y0=[Y0{:}];
dY0=[dY0{:}];


MINDy_Initialize


doRecord=(strcmpi(ParStr.RecCorr,'y')+strcmpi(ParStr.RecW,'y')+strcmpi(ParStr.RecDdiff,'y')+strcmpi(ParStr.RecA,'y'))~=0;



C=zeros(nX,1); %#ok<NASGU>

ParTmp.nStep=numel(tVec)+BatchSz-2;
StartTimes=MINDy_MakeStarts(HRFout.OrigTime,ParTmp);
nT=numel(StartTimes);
BatchInd=StartTimes(1,randi(nT,1,NBatch))'+(0:ParTmp.nStep);

%H2=repmat(1.1,nX,1);

if BatchSz<numel(tVec)
    error('Batch Size should be greater than HRF kernel length')
else
zeroMAT=zeros(nX,BatchSz-1);
end

Decay=Dmin+D.^2; %#ok<NODEF,NASGU>

WK=Wk1*Wk2; %#ok<NODEF>
W0=W1+WK; %#ok<NODEF>
nTV=numel(tVec);
for iBatch=1:NBatch
    
    if mod(iBatch,50)==1
        disp([iBatch NBatch])
        Out.RecH1(:,ceil(iBatch/50))=H1;
        Out.RecH2(:,ceil(iBatch/50))=H2;
        if doRecord
            iRep=1; %#ok<NASGU>
            ModularRecord_00
        end
    end
    
HRFcomb=permute((H1.^PolyVec(1,:)).*(H2.^PolyVec(2,:)),[1 3 2]);

bInd=BatchInd(iBatch,:);

tX=xHRFpoly(:,bInd,:);
dYx=dY0(:,bInd(1:BatchSz));
Yx=Y0(:,bInd(1:BatchSz));

X0=sum(HRFcomb.*tX,3);

Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;
bbTemp=X0+B1(:,2); %#ok<NODEF>

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    

NonLin=W0*PolyZ;

hTmp=h(H1,H2);
HRFkern=hTmp+Hfix;

dHRFdH1=(h(H1+dH,H2)-hTmp)/dH;
%dHRFdH2=(h(H1,H2+dH)-HRFkern)/dH;

%dHRFdH1=(((1./log(H2))-psi(H1))+lnT).*hTmp;
dHRFdH2=ReScale*((H1./H2)-tVec).*hTmp;

FFTnonLin=fft(NonLin,[],2);

%NonLinConv=MINDy_MatConv(NonLin,Hstack,2);
    
tmpPred=ifft(fft([HRFkern zeroMAT],[],2).*FFTnonLin,[],2);
E=dYx+(Decay.*Yx)-tmpPred(:,nTV:size(NonLin,2));


%dEdH=mean(E.*NonLinConv(:,:,[2 3]),2);
tmp1h=ifft(fft([dHRFdH1 zeroMAT],[],2).*FFTnonLin,[],2);
dEdH1=mean(E.*tmp1h(:,nTV:size(NonLin,2)),2);
tmp2h=ifft(fft([dHRFdH2 zeroMAT],[],2).*FFTnonLin,[],2);
dEdH2=mean(E.*tmp2h(:,nTV:size(NonLin,2)),2);
%dEdH1=mean(E.*MINDy_MatConv(NonLin,dHRFdH1,2),2);
%dEdH2=mean(E.*MINDy_MatConv(NonLin,dHRFdH2,2),2);


EXXH=MINDy_LongMatConv(E,HRFkern(:,end:-1:1),2);

W1E=W0'*EXXH;

gradD=-2*D.*mean(E.*Yx,2);
%gradC=sum(EXXH,2)/BatchSz;

gradW1=(EXXH*PolyZ')/BatchSz;

if L2SpPlsEN~=0
%gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1)-L2SpPlsEN*WK*Wk2'/sqrt(sum(WK(:).^2));
%gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2)-L2SpPlsEN*Wk1'*WK/sqrt(sum(WK(:).^2));
WcTemp=gradW1-L2SpPlsEN*WK;
gradWk1=(WcTemp)*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*WcTemp-SpPls2*sign(Wk2);
else
gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2);
end

mP=mean(abs(PolyZ),2)'; %J=mean(abs(W1),1);%./(mP.^2);
wCost=((SpScale)*sign(W1)+(2*ENL2)*W1).*(mP/mean(mP));
gradW1=gradW1-wCost-diag(SpDiag*sign(diag(W1)));


DpTemp=(B1(:,1).^-1).*W1E.*((1./P1)-(1./P2));
dZdX=((bbTemp+(B1(:,1))/2).*DpTemp+(W1E./P2));

%% Really fitting grad (A.^-.5)
gradA=-A.^(1.5).*(mean(DpTemp,2)/2);

dH1=permute(PolyVec(1,:).*(H1.^(PolyVec(1,:)-1)).*(H2.^PolyVec(2,:)),[1 3 2]);
dH2=permute(PolyVec(2,:).*(H1.^(PolyVec(1,:))).*(H2.^(PolyVec(2,:)-1)),[1 3 2]);

H1X=sum(dH1.*tX,3);
H2X=sum(dH2.*tX,3);

gradH1=dEdH1+(sum(dZdX.*H1X,2)/BatchSz);
gradH2=dEdH2+(sum(dZdX.*H2X,2)/BatchSz);

    if mod(iBatch,50)==1
        Out.E(:,ceil(iBatch/50))=mean(E.^2,2);
        Out.gradH1(:,ceil(iBatch/50))=gradH1;
        Out.gradH2(:,ceil(iBatch/50))=gradH2;
    end

mH1=mu*mH1+(1-mu)*gradH1;
mH2=mu*mH2+(1-mu)*gradH2;
mW1=mu*mW1+(1-mu)*gradW1;
mD=mu*mD+(1-mu)*gradD;
mA=mu*mA+(1-mu)*gradA;
mWk1=mu*mWk1+(1-mu)*gradWk1;
mWk2=mu*mWk2+(1-mu)*gradWk2;

nH1=v*nH1+(1-v)*gradH1.^2;
nH2=v*nH2+(1-v)*gradH2.^2;
nW1=v*nW1+(1-v)*gradW1.^2;
nD=v*nD+(1-v)*gradD.^2;
nA=v*nA+(1-v)*gradA.^2;
nWk1=v*nWk1+(1-v)*gradWk1.^2;
nWk2=v*nWk2+(1-v)*gradWk2.^2;


gHat=(1-mu)/(1-mu^iBatch);
mHat=mu/(1-mu^(iBatch+1));
v0=1/(1-v^iBatch);


Rtv0=sqrt(v0);

%if iBatch<1000
Wk1=Wk1+(((wk1Rate*mHat/Rtv0)*mWk1)+(wk1Rate*gHat/Rtv0)*gradWk1)./(sqrt(nWk1)+(Reg/Rtv0));
Wk2=Wk2+(((wk2Rate*mHat/Rtv0)*mWk2)+(wk2Rate*gHat/Rtv0)*gradWk2)./(sqrt(nWk2)+(Reg/Rtv0));
W1=W1+(((W1Rate*mHat/Rtv0)*mW1)+(W1Rate*gHat/Rtv0)*gradW1)./(sqrt(nW1)+(Reg/Rtv0));

A=((A.^-.5)+(((ARate*mHat/Rtv0)*mA)+(ARate*gHat/Rtv0)*gradA)./(sqrt(nA)+(Reg/Rtv0))).^-2;

if DRegScale>0 
D=abs(D+DRegScale*DRate*((mHat*mD)+gHat*gradD)./(sqrt(v0*nD)+DRegScale));
else    
D=abs(D+DRate*((mHat*mD)+gHat*gradD));
end
%end

%if iBatch>1000
H1=H1+(((H1Rate*mHat/Rtv0)*mH1)+(H1Rate*gHat/Rtv0)*gradH1)./(sqrt(nH1)+(Reg/Rtv0));
H2=H2+(((H2Rate*mHat/Rtv0)*mH2)+(H2Rate*gHat/Rtv0)*gradH2)./(sqrt(nH2)+(Reg/Rtv0));
%end

%H1=(H1+H1Rate*((mHat*mH1)+gHat*gradH1));
%H2=(H2+H2Rate*((mHat*mH2)+gHat*gradH2));

%W1(~Mask)=0;
A=min(max(A,.25*(B1(:,1).^2)+AdiffMin),.25*(B1(:,1).^2)+AdiffMax);

H1=max(H1min,min(H1max,H1));
H2=max(H2min,min(H2max,H2));
end
%% Remove Rescaling from HRF functions
Hfix=-cH*((tVec.^(a2H-1).*exp(-b2H.*tVec).*(b2H.^a2H))./gamma(a2H)); %#ok<NASGU>
h=@(h1,h2)((tVec.^(h1-1).*exp(-h2.*tVec).*((h2.^h1)./gamma(h1)))); %#ok<NASGU>
%% Apply Rescaling to W's
WK=ReScale*WK;W0=ReScale*W0; %#ok<NASGU>
W1=ReScale*W1;
Wk1=sqrt(ReScale)*Wk1;
Wk2=sqrt(ReScale)*Wk2;

HRFcomb=permute((H1.^PolyVec(1,:)).*(H2.^PolyVec(2,:)),[1 3 2]);
X0=sum(HRFcomb.*xHRFpoly,3);
dX=sum(HRFcomb.*DxHRFpoly,3);

%D=D.*std(dX,[],2);
%W1=W1.*std(dX,[],2);
%Wk1=Wk1.*std(dX,[],2);


Out.AutoCorr=DiagCorr(X0',dX');
Out.HRF={H1,H2};

Decay=Dmin+D.^2; 

WK=Wk1*Wk2;
W0=W1+WK;
bbTemp=X0+B1(:,2);

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    

A0=sqrt((A./(B1(:,1).^2)-.25));
B10(:,1)=B1(:,1).^-1;   B10(:,2)=B1(:,2)./B1(:,1);
A=A0;   B1=B10; %#ok<NASGU>

MINDy_StoreVal;
Out.Param={Val.W1,Val.A,Val.B1,Val.C,Val.W1+Val.Wk1*Val.Wk2,Val.D}; 
Out.Training_Rescale=ReScale;

%% Adjust for region-dependent normalization in deconvolved domain
DeconvDat=MINDy_CrossHRF_Deconv(X,Pre,Out,tLength);

%madD=median(abs([DeconvDat{:}]-median([DeconvDat{:}],2)),2);
%madX=median(abs(X0-median(X0,2)),2);
if iscell(DeconvDat)
stdD=std([DeconvDat{:}],[],2);
else
    stdD=std(DeconvDat,[],2);
end
    
stdX=std(X0,[],2);

%Decay=Out.Param{6}.*(stdD./stdX);
%Val.D=Decay;
%Out.Param{6}=Decay;
Out.Corr=DiagCorr(dX',(W0*PolyZ-(Decay.*X0))');
%% Ratio of un-normalized and normalized std in deconvolved domain
Out.HRFscale=stdD./stdX;
%% Adjust for global regularization bias:
Out=MINDy_Inflate_HRF(Out,HRFout,derivKern,doRobust);

%Out=MINDy_Inflate(Out,X0,dX,doRobust);
end



