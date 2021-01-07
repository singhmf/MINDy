function[PredDat]=MINDy_Predict_NoDeriv_PrePost_Dout(Xorig,ooP,Pre,preDrop,endDrop)
%% Predicts according to a hybrid MINDy model:
%% Deconvolution for nonlinear part, but no convolution for linear decay part 
%%      (since it commutes) to reduce bias
%% Xorig = BOLD data
%% ooP=MINDy output model (recommended: Fit with MINDy_Din_2_Douter)
%% Pre=Preprocessing structure (contains TR and "ConvLevel" (NSR) for Wiener deconvolution)
%% preDrop=number of time points to drop from start
%% endDrop=number of time points to drop from end
%%  e.g. if want to use BOLD(:,1:(end-1)) to predict BOLD(:,2:end) use preDrop=0;endDrop=1.

commNSR=Pre.ConvLevel;
disp('Remember: Already factor normalization into transfer scale')
if ~iscell(Xorig)
    Xorig={Xorig};
    wasCell='n';
else
    wasCell='y';
end
PredDat=cell(size(Xorig));
ooP=MakeMINDyFunction(ooP);
for ii=1:numel(Xorig)
    
fTrue=MINDy_MakeHRF_H1H2(ooP.HRF{1},ooP.HRF{2});
fftHRF=fft(fTrue(Pre.TR*(0:(size(Xorig{ii},2)-1))),[],2);

DeconvDat=real(ifft(fft(Xorig{ii},[],2).*(conj(fftHRF)./(commNSR+abs(fftHRF).^2)),[],2))...
    ./ooP.NormScalePre;
BrainPred=ooP.Param{5}*ooP.Tran(DeconvDat);
    PrTmp=ooP.NormScale.*real(ifft(fft(BrainPred,[],2).*fftHRF,[],2));
    PrTmp=PrTmp+(1-ooP.Param{6}).*Xorig{ii};
    PredDat{ii}=PrTmp(:,(1+preDrop):(end-endDrop));
end
if strcmpi(wasCell,'n')
    PredDat=PredDat{1};
end
end