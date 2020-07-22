function[Out,Wfin,Dfin]=MINDy_Simple(Dat,isDnSampDeriv,varargin)
%% Just Input Data as region x time (a single matrix or cell of matrices)
%% IsDnSampDeriv denotes whether to take the one step derivative ('n') or...
%% 2-step derivative: x(t+2)-x(t)/2 ('y'/default)

if isempty(isDnSampDeriv)
    isDnSampDeriv='y';
end


%% Output.Param={Wsparse,A,b,c,Wfull,D}
%% optional second and third outputs give Wfull and D

%% Optional input argument is whether to do preprocessing (filtering and deconvolution)
if isempty(varargin)
    doPreProc='y';
else
    doPreProc=varargin{1};
end

ChosenPARSTR;
ParStr.BatchSz=300;ParStr.NBatch=5000;
doRobust='n';
%% If you are using a different TR insert it here (in seconds) as:
%% Pre.TR=(your TR)

ParStr.H1min=5;ParStr.H1max=7;ParStr.H2min=.7;ParStr.H2max=1.3;ParStr.H1Rate=.1;ParStr.H2Rate=.1;
ParStr.L2SpPlsEN=0;

if ~iscell(Dat)
Dat={Dat};
end
for i=1:numel(Dat)
Dat{i}=Dat{i}(:,~isnan(sum(Dat{i},1)));
end


Dat=cellfun(@(xx)(zscore(xx')'),Dat,'UniformOutput',0);

if strcmpi(doPreProc,'y')
Dat=MINDy_RestingPreProcInterp(Dat,Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR);
Dat=cellfun(@(xx)(zscore(xx(:,20:(end-20))')'),Dat,'UniformOutput',0);
end

%% For HCP TRs down-sample by a factor of two
if strcmpi(isDnSampDeriv(1),'y')
dDat=cellfun(@(xx)(convn(xx,[1 0 -1]/2,'valid')),Dat,'UniformOutput',0);
Dat=cellfun(@(xx)(xx(:,1:end-2)),Dat,'UniformOutput',0);
elseif strcmpi(isDnSampDeriv(1),'n')
dDat=cellfun(@(xx)(convn(xx,[1 -1],'valid')),Dat,'UniformOutput',0);
Dat=cellfun(@(xx)(xx(:,1:end-1)),Dat,'UniformOutput',0);
end

Out=MINDy_Base(Dat,dDat,Pre,ParStr);

[X1,dX1]=MINDy_Censor(Out,Dat,dDat);

%% Global regression on W and D to get rid of global regularization bias
Out=MINDy_Inflate(Out,X1,dX1,doRobust);
Out=MakeMINDyFunction(Out);

%% Recalculate goodness-of-fit
Out.Corr=DiagCorr(Out.FastFun([X1{:}])',[dX1{:}]');
Out=rmfield(Out,{'FastFun','Tran'});
Out.Pre=Pre;Out.ParStr=ParStr;
Wfin=Out.Param{5};Dfin=Out.Param{6};
end