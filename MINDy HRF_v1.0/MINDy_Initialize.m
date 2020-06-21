%% Initialize the MINDy Software
%% For non-HRF clean data

%% Set Defaults
ParStr=MINDy_SetDefaults(ParStr);

if strcmpi(MINDy_type,'HRF')
SpecificRates={'ARate','DRate','wk1Rate','wk2Rate','H1Rate','H2Rate'};
else
SpecificRates={'ARate','DRate','wk1Rate','wk2Rate'};
end


for iSR=1:numel(SpecificRates)
    ParStr.(SpecificRates{iSR})=ParStr.(SpecificRates{iSR})*ParStr.Rate;
end
ParStr.W1Rate=ParStr.Rate;
mu=ParStr.Dec;v=ParStr.Dec2;


%% Extract
MINDy_ExtractParStr


%% Clean Data for non-HRF
if (~strcmpi(MINDy_type,'HRF'))&&(~strcmpi(MINDy_type,'Rec'))
    [X,dX,Out]=MINDy_Filter(X,dX,ParStr);
    nX=size(X,1);
    nT=size(X,2);
elseif strcmpi(MINDy_type,'Rec')
    RecStarts=MINDy_MakeStarts(X,ParStr);
    if iscell(X)
        X=[X{:}];
    end
    nX=size(X,1);
    nT=size(X,2);
else
    CleanStr=[];
    nX=size(xHRFpoly,1);
    nT=size(xHRFpoly,2);
end

%% I have nParc to generalize for the case when using BIG
if exist('numSource')==1
    nX=numSource;
end
nParc=nX;
%% Initialize Parameter Values
B1(:,1)=Bmin;
A=.5+(1+abs(randn(nX,1))).^(.5)*RandScale+AdiffMin+B1(:,1).^2/4;
B1(:,2)=0;
W1=RandScale*randn(nParc);
D=Dstart+5*sqrt(RandScale*abs(randn(nX,1)));
Wk1=RandScale*randn(nParc,wPC);Wk2=RandScale*randn(wPC,nParc);


%% Make NADAM Variables
mA=zeros(size(A));         mW1=zeros(size(W1));    mD=zeros(size(D));     
mWk1=zeros(size(Wk1));  mWk2=mWk1';
%%%%
nA=zeros(size(A));      nW1=zeros(size(W1));    nD=zeros(size(D));
nWk1=zeros(size(Wk1));  nWk2=nWk1';
%%%%

if strcmpi(MINDy_type,'HRF')
    H1=repmat(H1start,nX,1);
    H2=repmat(H2start,nX,1);
    mH1=zeros(size(H1));nH1=mH1;
    mH2=zeros(size(H2));nH2=mH2;
end

