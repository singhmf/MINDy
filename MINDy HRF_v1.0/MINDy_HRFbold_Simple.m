function[Out,HRF,W,D,DeconvDat]=MINDy_HRFbold_Simple(X,TR,doSmooth)
%% X is data: cell array of scans: region x time
%% TR is the TR in seconds
%% doSmooth indicates whether to do additional smoothing ([.5 .5] kernel)

ChosenPARSTR_HRF;
if nargin==1
    TR=.72;
    doSmooth='n';
    disp('By default: TR=.72s and no smoothing (smoothing is recommended for this TR)')
end

if nargin==2
    doSmooth='n';    
    disp('By default: no smoothing (smoothing is recommended for short TRs)')
end
Pre.TR=TR;
ParRes=30;
ParStr.H1Rate=10;ParStr.H2Rate=1;ParStr.H2min=.5;ParStr.H2max=1.5;
%% BatchSz determines the number of sequential time points selected on each iteration--
%% and is limited by the shortest scanning run
ParStr.BatchSz=100;
%% NBatch determines the number of iterations
ParStr.NBatch=6000;
%% The Weiner Noise-to-Signal Ratio (you may need to tweak this depending upon your scanner)
Pre.ConvLevel=.002;
ConvBase=Pre.ConvLevel;



if ~iscell(X)
    X={X};
end

%% Rescale dimensions of low-rank component proportionately to the number of parcels
nX=size(X{1},1);
%% 419 parcels in the original MINDy papers
ParStr.wPC=ceil(ParStr.wPC*nX/419);

%% HRF kernel length in seconds
tHRF=20;

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
%%  i.e. rescale the original hyperparameters to the HRF model

for ii=1:numel(X)
    [tmp,BadFrames{ii}]=MyRestingPreProcInterpNoHRF(X{ii},Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR,'n');
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
Out=MINDy_HRFbold_OrigX(X,Pre,ParStr,ParRes,ReScale,20,20,20);
Out.ReScale=ReScale;
HRF=MINDy_MakeHRF_H1H2(Out.HRF{1},Out.HRF{2});


W=Out.Param{5};
D=Out.Param{6};
Out.BadFrames=BadFrames;

DeconvDat=MINDy_CrossHRF_Deconv(X,Pre,Out,ceil(tHRF/Pre.TR));
end