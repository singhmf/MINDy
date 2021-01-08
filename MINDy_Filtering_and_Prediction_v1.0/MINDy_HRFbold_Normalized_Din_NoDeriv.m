function[Out]=MINDy_HRFbold_Normalized_Din_NoDeriv(X,Pre,ParStr,Res,ReScale,varargin)
%% Res Allows adding custom resolution for HRFpoly estimation
%% Multiplies Outgoing HRF by Rescale (recommended to use a ReScale=1 -> no further rescaling)
%% varargin can contain (in order): tLength=HRF length (in TRs), startDrop, endDrop
%% startDrop and endDrop denote how much of timeseries to clip at ends due to deconvolution

MINDy_type='HRF'; %#ok<NASGU>
dH=.0001;

%% Remark: mu here is just a placeholder to prevent interference with matlab's mu function
mu=[];
doRobust='n';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This line indicates errors are computed in terms of predicting BOLD(t+1) instead of BOLD(t+1)-BOLD(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
derivKern=[1 0];

if isempty(varargin)
    tLength=30;
    HRFin={tLength};
else
    HRFin=varargin;
end

HRFout=MINDy_MakeHRFpoly_Normalize_PrePost_Silent(X,Pre,ParStr,'n',Res,HRFin{:});
Out.HRFout=HRFout;
DeconvLength=HRFout.DeconvLength; 
if iscell(X)
    dY0=cell(1,numel(X));
else
    dY0=cell(1);
end

SDpoly=HRFout.SDpoly;

xHRFpoly=HRFout.Poly; 
DxHRFpoly=HRFout.dPoly;
PolyVec=HRFout.PolyMat;

nX=size(HRFout.Poly,1);

%% time vector for HRF length
tVec=[0 Pre.TR*(1:HRFout.DeconvLength)];

a2H=16;
b2H=1;
cH=1/6;

%% Fixed part of HRF (the undershoot)
Hfix=-cH*((tVec.^(a2H-1).*exp(-b2H.*tVec).*(b2H.^a2H))./gamma(a2H));

%% Parameterized part of HRF
h=@(h1,h2)((tVec.^(h1-1).*exp(-h2.*tVec).*((h2.^h1)./gamma(h1))));

Y0=cell(1,numel(HRFout.OrigTime));
for ii=1:numel(HRFout.OrigTime)
    Y0{ii}=(X{ii}(:,HRFout.OrigTime{ii}(1,(1+ceil((DeconvLength-1)/2)):(end-(floor((DeconvLength-1)/2)))))')';
end

for iY=1:numel(HRFout.ReMean)
    tTmp=size(HRFout.OrigTime{iY},2);    
    Y0{iY}=[Y0{iY} zeros(nX,tTmp-size(Y0{iY},2))];
    dY0{iY}=[convn(Y0{iY},derivKern,'valid') zeros(nX,1)];
end
Y0=[Y0{:}]; %#ok<NASGU>
dY0=[dY0{:}];

MINDy_Initialize

doRecord=(strcmpi(ParStr.RecCorr,'y')+strcmpi(ParStr.RecW,'y')+strcmpi(ParStr.RecDdiff,'y')+strcmpi(ParStr.RecA,'y'))~=0;

%% Constant input not used in this model
C=zeros(nX,1); %#ok<NASGU>

ParTmp.nStep=numel(tVec)+BatchSz-2;
StartTimes=MINDy_MakeStarts(HRFout.OrigTime,ParTmp);
nT=numel(StartTimes);
BatchInd=StartTimes(1,randi(nT,1,NBatch))'+(0:ParTmp.nStep);

if BatchSz<numel(tVec)
    error('Batch Size should be greater than HRF kernel length')
else
zeroMAT=zeros(nX,BatchSz-1);
end

Decay=-1+Dmin+D.^2; %#ok<NODEF,NASGU>

WK=Wk1*Wk2; %#ok<NODEF>
W0=W1+WK; %#ok<NASGU,NODEF>
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
%Yx=Y0(:,bInd(1:BatchSz));

X0=sum(HRFcomb.*tX,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Changed this line so becomes (1-D)X since predicting BOLD(t+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Decay=-1+Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;
bbTemp=X0+B1(:,2); %#ok<NODEF>

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    

NonLin=ReScale*(W0*PolyZ)-Decay.*X0;

hTmp=h(H1,H2);
HRFkern0=hTmp+Hfix;

%% Numeric gradients of the HRF w.r.t. H1, H2
dHRFdH1=(h(H1+dH,H2)-hTmp)/dH;
dHRFdH2=(h(H1,H2+dH)-hTmp)/dH;

%% Analytic gradients are actually slightly slower and require a toolbox for the psi function
%dHRFdH1=(((1./log(H2))-psi(H1))+lnT).*hTmp;
%dHRFdH2=ReScale*((H1./H2)-tVec).*hTmp;

dH1=permute(PolyVec(1,:).*(H1.^(PolyVec(1,:)-1)).*(H2.^PolyVec(2,:)),[1 3 2]);
dH2=permute(PolyVec(2,:).*(H1.^(PolyVec(1,:))).*(H2.^(PolyVec(2,:)-1)),[1 3 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalization re-scaling %%
NormScale=sum(squeeze(HRFcomb).*SDpoly,2);
dHRFdH1=NormScale.*dHRFdH1+sum(squeeze(dH1).*SDpoly,2).*HRFkern0;
dHRFdH2=NormScale.*dHRFdH2+sum(squeeze(dH2).*SDpoly,2).*HRFkern0;

HRFkern=NormScale.*HRFkern0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFTnonLin=fft(NonLin,[],2);

%% Predicting BOLD(t+1) by reconvolving predicted neural with HRF
tmpPred=ifft(fft([HRFkern zeroMAT],[],2).*FFTnonLin,[],2);

E=dYx-tmpPred(:,nTV:size(NonLin,2));


tmp1h=ifft(fft([dHRFdH1 zeroMAT],[],2).*FFTnonLin,[],2);
dEdH1=mean(E.*tmp1h(:,nTV:size(NonLin,2)),2);
tmp2h=ifft(fft([dHRFdH2 zeroMAT],[],2).*FFTnonLin,[],2);
dEdH2=mean(E.*tmp2h(:,nTV:size(NonLin,2)),2);

%% Function that accelerates vectorized convolution (errors backprop'd through the HRF)
EXXH=MINDy_LongMatConv(E,HRFkern(:,end:-1:1),2);

W1E=W0'*EXXH*ReScale;

gradD=-2*D.*mean(EXXH.*X0,2);

gradW1=(EXXH*PolyZ')*(ReScale/BatchSz);

%% Regularization terms

if L2SpPlsEN~=0
WcTemp=gradW1-L2SpPlsEN*WK;
gradWk1=(WcTemp)*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*WcTemp-SpPls2*sign(Wk2);
else
gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2);
end
wCost=((SpScale)*sign(W1)+(2*ENL2)*W1);
gradW1=gradW1-wCost-diag(SpDiag*sign(diag(W1)));

DpTemp=(B1(:,1).^-1).*W1E.*((1./P1)-(1./P2));
dZdX=((bbTemp+(B1(:,1))/2).*DpTemp+(W1E./P2))-Decay.*EXXH;

%% Really fitting grad (A.^-.5)
gradA=-A.^(1.5).*(mean(DpTemp,2)/2);

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

Wk1=Wk1+(((wk1Rate*mHat/Rtv0)*mWk1)+(wk1Rate*gHat/Rtv0)*gradWk1)./(sqrt(nWk1)+(Reg/Rtv0));
Wk2=Wk2+(((wk2Rate*mHat/Rtv0)*mWk2)+(wk2Rate*gHat/Rtv0)*gradWk2)./(sqrt(nWk2)+(Reg/Rtv0));
W1=W1+(((W1Rate*mHat/Rtv0)*mW1)+(W1Rate*gHat/Rtv0)*gradW1)./(sqrt(nW1)+(Reg/Rtv0));

A=((A.^-.5)+(((ARate*mHat/Rtv0)*mA)+(ARate*gHat/Rtv0)*gradA)./(sqrt(nA)+(Reg/Rtv0))).^-2;

if DRegScale>0 
D=abs(D+DRegScale*DRate*((mHat*mD)+gHat*gradD)./(sqrt(v0*nD)+DRegScale));
else    
D=abs(D+DRate*((mHat*mD)+gHat*gradD));
end

H1=H1+(((H1Rate*mHat/Rtv0)*mH1)+(H1Rate*gHat/Rtv0)*gradH1)./(sqrt(nH1)+(Reg/Rtv0));
H2=H2+(((H2Rate*mHat/Rtv0)*mH2)+(H2Rate*gHat/Rtv0)*gradH2)./(sqrt(nH2)+(Reg/Rtv0));

A=min(max(A,.25*(B1(:,1).^2)+AdiffMin),.25*(B1(:,1).^2)+AdiffMax);

H1=max(H1min,min(H1max,H1));
H2=max(H2min,min(H2max,H2));
end


HRFcomb=permute((H1.^PolyVec(1,:)).*(H2.^PolyVec(2,:)),[1 3 2]);
NormScale=sum(squeeze(HRFcomb).*SDpoly,2);
Out.NormScalePre=sum(squeeze(HRFcomb).*HRFout.SDpolyPre,2);
Out.NormScale=NormScale;

%% Remove Rescaling from HRF functions
Hfix=-cH*((tVec.^(a2H-1).*exp(-b2H.*tVec).*(b2H.^a2H))./gamma(a2H)); 
h=@(h1,h2)((tVec.^(h1-1).*exp(-h2.*tVec).*((h2.^h1)./gamma(h1))));
%% Apply Rescaling to W's
Wk1=sqrt(ReScale)*Wk1;
Wk2=sqrt(ReScale)*Wk2;
W1=ReScale*W1;
WK=Wk1*Wk2;
W0=W1+WK;
HRFcomb=permute((H1.^PolyVec(1,:)).*(H2.^PolyVec(2,:)),[1 3 2]);
X0=sum(HRFcomb.*xHRFpoly,3);
dX=sum(HRFcomb.*DxHRFpoly,3);
Out.deconvX=X0;


hTmp=h(H1,H2);
HRFkern0=hTmp+Hfix;
Out.HRF_kern=HRFkern0;
Out.normHRF_kern=NormScale.*HRFkern0;
Out.normHRF_kernPre=Out.NormScalePre.*HRFkern0;

Out.AutoCorr=DiagCorr(X0',dX');
Out.HRF={H1,H2};

Decay=(Dmin+D.^2); 

bbTemp=X0+B1(:,2);

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    

A0=sqrt((A./(B1(:,1).^2)-.25));
B10(:,1)=B1(:,1).^-1;   B10(:,2)=B1(:,2)./B1(:,1);
A=A0;   B1=B10; %#ok<NASGU>

MINDy_StoreVal;
Out.tVec=tVec;
Out.Param={Val.W1,Val.A,Val.B1,Val.C,Val.W1+Val.Wk1*Val.Wk2,Val.D}; 
Out.Training_Rescale=ReScale;
Out.Mat={W1,Wk1,Wk2};
Out.Corr=DiagCorr(dX',(W0*PolyZ-(Decay.*X0))');

%% Adjust for global regularization bias:
Out=MINDy_Inflate_HRF_PrePost_NoDeriv(X,Out,Pre,doRobust);

if exist('Rec','var')
    Out.Rec=Rec;
end
end