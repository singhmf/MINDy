function[Out]=MINDy_D2D(Xorig,Pre,ParStr,Res,reScale,tLength,varargin)
%% Used for re-estimating models in which decay is estimated inside convolution (good for rs Weights)
%% to models in which decay is estimated outside of convolution (better for filtering)
%% This step minimizes bias as the portions that commute with convolution (univariate linear) are 
%% allowed to remain outside of the convolution-deconvolution sequence

%% Xorig = BOLD data (Parc x Time matrix per scan) cell array for multiple scans
%% Pre=Preprocessing structure (make sure TR is right)
%% ParStr=Hyperparameter structure
%% Res=Sampling resolution for surrogate surface (more is better/slower, recommended 15-25)
%% reScale=additional rescaling in deconvolved domain (recommended: 1 [no rescaling])
%% tLength=HRF length


commNSR=Pre.ConvLevel;
%% Original model with linear terms inside convolution
ooP=MINDy_HRFbold_Normalized_Din_NoDeriv(Xorig,Pre,ParStr,Res,reScale,tLength,varargin{:});
ooP=MakeMINDyFunction(ooP);
wDat=cell(size(Xorig));linDat=cell(size(Xorig));
dDat=cell(size(Xorig));

%% Extract region-specific HRF's
fTrue=MINDy_MakeHRF_H1H2(ooP.HRF{1},ooP.HRF{2});

for ii=1:numel(Xorig)

%% Fourier transform of HRF
fftHRF=fft(fTrue(Pre.TR*(0:(size(Xorig{ii},2)-1))),[],2);

%% Wiener deconvolution using normalized HRF
DeconvDat=real(ifft(fft(Xorig{ii},[],2).*(conj(fftHRF)./(commNSR+abs(fftHRF).^2)),[],2))...
    ./ooP.NormScalePre;
%% Prediction of nonlinear/multivariate part inside convolution
wBrainPred=ooP.Param{5}*ooP.Tran(DeconvDat);
%% Reconvolution with normalized HRF
tmp=ooP.NormScale.*real(ifft(fft(wBrainPred,[],2).*fftHRF,[],2));
    wDat{ii}=tmp(:,1:(end-1));
%% BOLD(t)
linDat{ii}=Xorig{ii}(:,1:(end-1));
%% BOLD(t+1)
dDat{ii}=Xorig{ii}(:,2:end);
end

nX=size(dDat{1},1);

wDat=[wDat{:}];
linDat=[linDat{:}];
dDat=[dDat{:}];

%% Solve: BOLD(t+1)=(1-D)BOLD(t)+c*NonLin(t)
wL2=sum(sum(wDat.^2));
wNorm=wDat/wL2;
Ypart=dDat-sum(sum(wNorm.*dDat))*wNorm;
newD=zeros(nX,1);
for ii=1:nX
    newD(ii)=Ypart(ii,:)/linDat(ii,:);
end
Yleft=dDat-newD.*linDat;
wScale=Yleft(:)'/wDat(:)';

%% Update parameter structure
Out=ooP;
Out.wScale=wScale;
Out.DinParam=ooP.Param;
Out.Param{5}=wScale*ooP.Param{5};
Out.Param{1}=wScale*ooP.Param{1};
Out.Param{6}=1-newD;
end