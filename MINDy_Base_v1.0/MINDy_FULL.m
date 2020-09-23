function[Out]=MINDy_FULL(Dat,Pre,ParStr,isDnSampDeriv,doPreProc,doInflate,doSmooth)
%% This is the master-function for fully customizable MINDy (without HRF estimation).

%% Dat is input as a cell of timeseries (region x time) or a single time-series
%% Pre is a structure of parameters for doing preprocessing (includes TR, deconvolution filtering etc.)
%% ParStr is a structure of hyperparameters for model fitting
%% isDnSampDeriv denotes whether to use   dx=x(t+1)-x(t) (='n') or 
%%                                        dx=(x(t+2)-x(t))/2 (='y')
%% doPreProc=whether to do deconvolution/spike removal
%% doInflate=whether to do global readjustment of weights vs. decay post-fitting
%% doSmooth=whether to do nearest neighbor smoothing: new x(t+1)<-(x(t)+x(t+1))/2

if ~iscell(Dat)
    Dat={Dat};
end

%% Preprocessing (if desired)
if strcmpi(doPreProc,'y')
for ii=1:numel(Dat)
Dat{ii}=Dat{ii}(:,~isnan(sum(Dat{ii},1)));
end
Dat=cellfun(@(xx)(zscore(xx')'),Dat,'UniformOutput',0);
Dat=MINDy_RestingPreProcInterp(Dat,Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR);
Dat=cellfun(@(xx)(zscore(xx(:,50:(end-50))')'),Dat,'UniformOutput',0);
end


if strcmpi(doSmooth(1),'y')
    for ii=1:numel(Dat)
        Dat{ii}=convn(Dat{ii},[1 1]/2,'valid');
    end
end
%% Calculating derivatives
if ischar(isDnSampDeriv)
    if strcmpi(isDnSampDeriv(1),'y')
        DerivDat=Uncellfun(@(xx)(convn(xx,[1 0 -1]/2,'valid')),Dat);
        Dat=Uncellfun(@(xx)(xx(:,1:end-2)),Dat);
    elseif strcmpi(isDnSampDeriv(1),'n')        
        DerivDat=Uncellfun(@(xx)(convn(xx,[1 -1],'valid')),Dat);
        Dat=Uncellfun(@(xx)(xx(:,1:end-1)),Dat);
    else
        error('isDnSampDeriv should be y/n or numeric')
    end
else
        DerivDat=Uncellfun(@(xx)(convn(xx,[1 -1]*isDnSampDeriv,'valid')),Dat);
        Dat=Uncellfun(@(xx)(xx(:,1:end-1)),Dat);
end 
Out=MINDy_Base(Dat,DerivDat,Pre,ParStr);
[X1,dX1]=ModularCensor(Out,Dat,DerivDat);
if strcmpi(doInflate(1),'y')
    if ~isfield(ParStr,'doRobust')        
    Out=MINDy_Inflate(Out,X1,dX1,'y');
    else
    Out=MINDy_Inflate(Out,X1,dX1,ParStr.doRobust);
    end
end
Out=MakeMINDyFunction(Out);
Out.Corr=DiagCorr(Out.FastFun([X1{:}])',[dX1{:}]');
end