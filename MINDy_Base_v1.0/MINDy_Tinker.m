function[Out,Wfin,Dfin]=MINDy_Tinker(Dat,TR,isDnSampDeriv,doPreProc,lam1,lam2,SpDiag,MaxRank,BatchSz,NBatch,ConvNSR,Wmask,doRecord)
%% Just Input Data as region x time (a single matrix or cell of matrices)
%% TR is the sampling TR for fMRI in seconds
%% IsDnSampDeriv denotes whether to take the one step derivative ('n') or...
%% 2-step derivative: x(t+2)-x(t)/2 ('y'/default)
%% doPreProc denotes whether to do pre-processing (normalization and deconvolution)

%% Parameter Tinkering (these are the main parameters)
%% lam1: The L1 regularization for the sparse component (original=.075)
%% lam2: The L1 regularization for the low-rank components (original=.05)
%% SpDiag: The additional L1 regularization added to diagonal (recurrent connections) above Lam1 (original=.2)
%%          note: negative values are admissible (less regularization applied to recurrents)
%% MaxRank: The maximum dimension of the low-rank component (original = 150 for 419 parcels)
%% BatchSz: Number of Time-Points sampled each iteration (original=300)
%% NBatch: The number of iterations (original=5000)
%% ConvNSR: The Noise-Signal-Ratio regularizer used to determine filtering level (original=.02)
%%          (smaller values->less filtering [recommended to be small])
%% Wmask: Optional connectivity mask: (binary nxn matrix)
%% doRecord: Whether to record the evolution of parameter estimates
%%              'y'=record all parameters, 'n'=don't record any (default)
%%              'W'=weights, 'D'=decay','A'=Curvature, can specify combinations like 'DA'
%% Empty values [] indicate default

%% Output.Param={Wsparse,A,b,c,Wfull,D}
%% optional second and third outputs give Wfull and D

if isempty(isDnSampDeriv)
    isDnSampDeriv='y';
end


ChosenPARSTR;
ParStr.BatchSz=300;ParStr.NBatch=5000;

%% Optional robust regression in performing post-processing regression to rescale W and D
doRobust='n';

Pre.TR=TR;

ParStr.H1min=5;ParStr.H1max=7;ParStr.H2min=.7;ParStr.H2max=1.3;ParStr.H1Rate=.1;ParStr.H2Rate=.1;
ParStr.L2SpPlsEN=0;

Tinker_Names={'SpScale','SpPls1','SpPls2','SpDiag','wPC','BatchSz','NBatch'};
Tinker_Vars={lam1,lam2,lam2,SpDiag,MaxRank,BatchSz,NBatch};
%% Set new values for non-empty hyperparameters
for ii=1:numel(Tinker_Names)
    if ~isempty(Tinker_Vars{ii})
        ParStr.(Tinker_Names{ii})=Tinker_Vars{ii};
    end
end
if ~isempty(ConvNSR)
    Pre.ConvLevel=ConvNSR;
end

if strcmpi(doRecord,'y')
    ParStr.RecW='y';ParStr.RecA='y';ParStr.RecDdiff='y';
elseif ~isempty(doRecord)
    if mean(lower(doRecord)=='w')~=0
        ParStr.RecW='y';
    end
    if mean(lower(doRecord)=='d')~=0
        ParStr.RecDdiff='y';
    end
    if mean(lower(doRecord)=='a')~=0
        ParStr.RecA='y';
    end
end


if ~iscell(Dat)
Dat={Dat};
end
for i=1:numel(Dat)
Dat{i}=Dat{i}(:,~isnan(sum(Dat{i},1)));
end


Dat=cellfun(@(xx)(zscore(xx')'),Dat,'UniformOutput',0);

if strcmpi(doPreProc,'y')
Dat=MINDy_RestingPreProcInterp(Dat,Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR);
%% Trimming first 19 and last 20 TRs after deconvolution
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

Out=MINDy_Base(Dat,dDat,Pre,ParStr,Wmask);

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