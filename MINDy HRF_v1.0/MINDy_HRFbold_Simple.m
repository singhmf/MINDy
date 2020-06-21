function[Out,HRF,W,D,DeconvDat]=MINDy_HRFbold_Simple(X,TR,doSmooth)
%% X is data: cell array of scans: region x time
%% TR is the TR in seconds: Default =.72s
%% doSmooth indicates whether to do additional smoothing; Default='n'

%% Outputs:
%% Out=MINDy model output structure
%% HRF=function handle that produces the HRFs for all nodes evaluated at user-selected input times
%%      input arguments are 1xtime(s) and outputs are region x time
%% W and D give the weight and decay parameters
%% DeconvDat gives the deconvolved data.

ChosenPARSTR_HRF;
if nargin==1
    TR=.72;
    doSmooth='n';
    disp('Assuming HCP TR: 720ms')
end

if nargin==2
    doSmooth='n';
    disp('Assuming no additional smoothing')
end
Pre.TR=TR;
ParRes=10;
ParStr.H1Rate=10;ParStr.H2Rate=1;ParStr.H2min=.5;ParStr.H2max=1.5;
ParStr.BatchSz=250;
ParStr.NBatch=6000;
ConvBase=Pre.ConvLevel;
Pre.ConvLevel=.002;


if ~iscell(X)
    X={X};
end

tHRF=30;

%% Canonical HRF to calculate rescaling
a1=6;
a2=16;
b1=1;
b2=1;
c=1/6;
h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));

dX0=cell(1,numel(X));
X0=cell(1,numel(X));
HdX0=cell(1,numel(X));
HX0=cell(1,numel(X));
%hSD=cell(1,numel(X));
BadFrames=cell(1,numel(X));

%% Calculate change in SD under deconvolution in order to rescale priors.
%%  i.e. rescale the problem so that original hyperparameters will work

for ii=1:numel(X)
    [tmp,BadFrames{ii}]=MyRestingPreProcInterpNoHRF(X{ii}(:,3:(end-3)),Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR,'n');
    if strcmpi(doSmooth(1),'y')
    X{ii}=zscore(convn(tmp{1},[1 1]/2,'valid')')';
    else
        X{ii}=zscore(tmp{1}')';
    end
    dX0{ii}=convn(X{ii},[1 -1],'valid');
    X0{ii}=X{ii}(:,1:end-1);
    Dtmp=MINDy_DeconvHRF(X{ii},Pre.TR,ConvBase,ceil(tHRF/Pre.TR));
    HdX0{ii}=convn(zscore(Dtmp')',[1 -1],'valid');
    HX0{ii}=zscore(Dtmp(:,1:(end-1))')';
   % hSD{ii}=std(tmp,[],2);
end

%% Default parameters for transfer fun:
defPsi=@(x)(sqrt(25+((20*x/3)+.5).^2)-sqrt(25+((20*x/3)-.5).^2));
%% Estimate SD due to non-autocorrelation for BOLD measurements
boldSD=std([dX0{:}]-mean([X0{:}].*[dX0{:}],2).*[X0{:}],[],2);
    %% =h*psi(zscore(deconv(x,h)))
boldPsi=std(convn(defPsi([HX0{:}]),h(0:Pre.TR:ceil(tHRF/Pre.TR)),'valid'),[],2);
%% Estimate SD due to non-autocorrelation for deconvolved measurements
deconvSD=std([HdX0{:}]-mean([HX0{:}].*[HdX0{:}],2).*[HX0{:}],[],2);
    %% =psi(zscore(deconv(x,h)))
deconvPsi=std(defPsi([HX0{:}]),[],2);


%% Estimate rescaling needed to use original MINDy parameters for BOLD model fitting
ReScale=median((deconvPsi.*boldSD)./(boldPsi.*deconvSD));

%% For non-HCP data, recommended to also factor in different TR (scale relative to HCP 720ms)
ReScale=ReScale*TR/.72;
%Out=MINDy_HRFbold_Rescale(Uncellfun(@single,X),Pre,ParStr,10,2);
Out=MINDy_HRFbold_Rescale(Uncellfun(@single,X),Pre,ParStr,ParRes,ReScale);
Out.ReScale=ReScale;
HRF=MINDy_MakeHRF_H1H2(Out.HRF{1},Out.HRF{2});


W=Out.Param{5};
D=Out.Param{6};
Out.BadFrames=BadFrames;

DeconvDat=MINDy_CrossHRF_Deconv(X,Pre,Out,ceil(tHRF/Pre.TR));
end