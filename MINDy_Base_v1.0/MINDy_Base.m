function[Out]=MINDy_Base(X0,dX0,Pre,ParStr,Wmask)  %#ok<INUSL>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input arguments: X0=``Neural" (i.e. deconvolved-fMRI) time series 
%%                       n x time or a cell of n x time_i time series
%%                  dX0=Derivatives of X (e.g. X(t+1)-X0 or (X(t+2)-X)/2 etc.
%%                      NOTE: derivatives must be same format/size as X0 [so trim last values of X]
%%                  Pre=Preprocessing structure (generate from ChosenPARSTR)
%%                  ParStr=Parameter structure  (generate from ChosenPARSTR)
%%                  Optional Wmask=Mask matrix for connectivity (binary nxn)

MINDy_type='Base'; %#ok<NASGU>
X=X0;
dX=dX0;
%% To overrule functional definition of mu
mu=1;
MINDy_Initialize

iRep=1; %#ok<NASGU>
doRecord=(strcmpi(ParStr.RecCorr,'y')+strcmpi(ParStr.RecW,'y')+strcmpi(ParStr.RecDdiff,'y')+strcmpi(ParStr.RecA,'y'))~=0;

nX=size(X0,1);

C=zeros(nX,1);

if iscell(X)
    X0=[X{:}];
else
    X0=X;
end
if iscell(dX)
    dX0=[dX{:}];
else
    dX0=dX;
end

nT=size(X0,2);
BatchInd=randi(nT,NBatch,BatchSz);

if nargin<5
    Wmask=[];
elseif (mean(size(Wmask)==nX)~=1)||~isempty(find(~ismember(Wmask(:),[0 1]),1))
    error('Wmask should be n x n binary or empty/unspecified[default= all connections]')
end


%% Optimization for curvature is in terms of the linearized variable:
%%  Xi=b/sqrt(a^2+.25)=max_(x)(d/dx[Psi_a(x)])
for iBatch=1:NBatch

Decay=Dmin+D.^2;
WK=Wk1*Wk2;
W0=W1+WK;

%% Apply optional mask that excludes certain connections
if ~isempty(Wmask)
    W0=W0.*Wmask;
    W1=W1.*Wmask;
    WK=WK.*Wmask;
end
    

    if mod(iBatch,50)==1
        disp([iBatch NBatch])
        if doRecord
            ModularRecord
        end
    end

bInd=BatchInd(iBatch,:);

X=X0(:,bInd);
dX=dX0(:,bInd);

bbTemp=X+B1(:,2); %#ok<NODEF>

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    

if strcmpi(ParStr.FixC,'n')    
E=dX-(W0*PolyZ-(Decay.*X)+C);
else
E=dX-(W0*PolyZ-(Decay.*X));
end
W1E=W0'*E;

gradD=-2*D.*mean(E.*X,2);
gradW1=(E*PolyZ')/BatchSz;

%% Backprop gradients through optional connectivity mask
if ~isempty(Wmask)
    gradW1=Wmask.*gradW1;
end

if L2SpPlsEN~=0
WcTemp=gradW1-L2SpPlsEN*WK;
gradWk1=(WcTemp)*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*WcTemp-SpPls2*sign(Wk2);
else
gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2);
end

%% Commented out is a version that normalizes regularization for W_(i,j) by E[|Psi_j|}
%%      to prevent bias by curvature influencing weights (we found this did not matter).
% mP=mean(abs(PolyZ),2)';
% wCost=((SpScale)*sign(W1)+(2*ENL2)*W1);%.*(mP/mean(mP));

wCost=((SpScale)*sign(W1)+(2*ENL2)*W1);
gradW1=gradW1-wCost-diag(SpDiag*sign(diag(W1)));

DpTemp=(B1(:,1).^-1).*W1E.*((1./P1)-(1./P2));

%% Really fitting for (A.^-.5)
gradA=-A.^(1.5).*mean(DpTemp,2);

%% Record evolution of error (MSE across samples for each parcel)
    if mod(iBatch,50)==1
        Out.E(:,ceil(iBatch/50))=mean(E.^2,2);
    end

%% If want to include a non-zero forcing term (will be recorded in out.Param{4})
if strcmpi(ParStr.FixC,'n')
gradC=mean(E,2);
mC=mu*mC+(1-mu)*gradC;
nC=v*n+(1-v)*gradC.^2;
C=C+(((W1Rate*mHat/Rtv0)*mC)+(W1Rate*gHat/Rtv0)*gradC)./(sqrt(nC)+(Reg/Rtv0));
end    
    
%% note: gradA mA and nA 

mW1=mu*mW1+(1-mu)*gradW1;
mD=mu*mD+(1-mu)*gradD;
mA=mu*mA+(1-mu)*gradA;
mWk1=mu*mWk1+(1-mu)*gradWk1;
mWk2=mu*mWk2+(1-mu)*gradWk2;

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

if ~isempty(Wmask)
    W1=Wmask.*W1;
end

if DRegScale>0 
D=abs(D+DRegScale*DRate*((mHat*mD)+gHat*gradD)./(sqrt(v0*nD)+DRegScale));
else    
D=abs(D+DRate*((mHat*mD)+gHat*gradD));
end

%% Admissible values of Xi for A bounds [Xi=b/sqrt(a^2+.25)]
A=min(max(A,.25*(B1(:,1).^2)+AdiffMin),.25*(B1(:,1).^2)+AdiffMax);

end



Out.AutoCorr=DiagCorr(X',dX');

Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;

if ~isempty(Wmask)
    W0=W0.*Wmask;
    W1=W1.*Wmask;
    WK=WK.*Wmask; %#ok<NASGU>
end
    


bbTemp=X+B1(:,2);

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    
Out.Corr=DiagCorr(dX',(W0*PolyZ-(Decay.*X))');

%% Convert from Xi back to A=Xi/(b^2-.25)
A0=sqrt((A./(B1(:,1).^2)-.25));
B10(:,1)=B1(:,1).^-1;   B10(:,2)=B1(:,2)./B1(:,1);
A=A0;   B1=B10; %#ok<NASGU>

if doRecord
    gField=fields(Rec);
    for kk=1:numel(gField)
        Out.(gField{kk})=Rec.(gField{kk});
    end
end
MINDy_StoreVal;
Out.Param={Val.W1,Val.A,Val.B1,Val.C,W0,Val.D};
Out.Wmask=Wmask;
Out.Mat={W1,Wk1,Wk2};
end



