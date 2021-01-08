function[OutStruct]=MINDy_MakeHRFpoly_Normalize_PrePost_Silent(X0,Pre,ParStr,NormDeriv,ParamRes,varargin)
%% Adds tLength (same as HRFpoly_Mod)
%% Accounts for change in SD of h*(h*^-1 x)

if isempty(varargin)
    tLength=30;
else
    tLength=varargin{1};
end

OutStruct.ParamRes=ParamRes;
if numel(varargin)>=2
StartDrop=varargin{2};
    if numel(varargin)==2
        EndDrop=StartDrop;
    else
        EndDrop=varargin{3};
    end
else
StartDrop=39;
EndDrop=40;
end

OutStruct.DeconvLength=tLength;

if ~iscell(X0)
    X0={X0};
end

FiltAmp=Pre.FiltAmp;

PolyC=cell(1,numel(X0));
dPolyC=cell(1,numel(X0));
OrigTime=cell(1,numel(X0));

%% Polynomial basis terms in surrogate surface estimation
PolyTerms={[0 0],[1 0],[0 1],[2 0],[1 1],[0 2],[3 0],[2 1],[1 2],[0 3]};

%SDc=cell(1,numel(X0));

nX=size(X0{1},1);
SDc=nan(ParamRes,ParamRes,nX,numel(X0));
SDpost=nan(ParamRes,ParamRes,nX,numel(X0));
nPt=cellfun(@numel,X0)/nX;
%% Count-based normalization
%% Remark: Using total count here before clipping
nProp=reshape(nPt(:),1,1,1,numel(X0))/sum(nPt(:));


ConvRes=Pre.TR;

ChangeVec=linspace(ParStr.H1min,ParStr.H1max,ParamRes);
ChangeVec2=linspace(ParStr.H2min,ParStr.H2max,ParamRes);

Xfix=cell(1,numel(X0));
%% bad frames;
BF=cell(1,numel(X0));

for iXcell=1:numel(X0)
X=zscore(X0{iXcell}')';

%find(isnan(X))    
for iParc=1:size(X,1)
        gg=find(abs(X(iParc,:))>FiltAmp);
%        BadFrames{i,iParc}=gg;
%% Remove bad end points (can't interpolate)
gg(ismember(gg,[1 size(X,2)]))=[];
        goodFrames=setdiff(1:size(X,2),gg);
        X(iParc,gg)=interp1(goodFrames,X(iParc,goodFrames),gg);
        BF{iXcell}{iParc}=gg;
end
Xfix{iXcell}=X;


for i=1:numel(ChangeVec)
    for j=1:numel(ChangeVec)
%a1=6;
a2=16;
%b1=1;
b2=1;
c=1/6;
b1=ChangeVec2(j);
a1=ChangeVec(i);
h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));
ConvVec=[0 h(ConvRes*(1:tLength))];
tmp=deconvwnr(X,ConvVec,Pre.ConvLevel);

SDc(i,j,:,iXcell)=nanstd(tmp,[],2);
%% This is the changed part!!!!!
SDpost(i,j,:,iXcell)=nanstd(tmp,[],2)./nanstd(convn(tmp,ConvVec,'valid'),[],2);
    end
end
end

%% Create Polynomial for appropriate rescaling during deconvolution
SDmean=sum(SDc.*nProp,4);
%% Create Polynomial for appropriate rescaling during reconvolution
SDmeanPost=sum(SDpost.*nProp,4);


OutStruct.nProp=nProp;


    OutStruct.BadFrames=BF;
for iXcell=1:numel(X0)
    X=Xfix{iXcell};

%% ConvRes=TR
a1=6;
a2=16;
b1=1;
b2=1;
c=1/6;
h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));
ConvVec=[0 h(ConvRes*(1:tLength))];
Outbase=deconvwnr(X,ConvVec,Pre.ConvLevel);


kLength=ceil((tLength)/2);

%% For the Autocorrelation
for iParc=1:size(X,1)
tmpXcorr=xcov(Outbase(iParc,(1+StartDrop):(end-EndDrop)),kLength,'coeff');
OutStruct.ChanAutoCorr{iXcell}(iParc,:)=tmpXcorr;
end
tmpX=mean(OutStruct.ChanAutoCorr{iXcell},1);
tmpMx=convmtx(tmpX,numel(tmpX));
OutStruct.AutoCorr{iXcell}=tmpMx(1:(tLength+1),(1+kLength):(kLength+tLength+1));


OutStruct.ReMean{iXcell}=convn(Outbase(:,(1+StartDrop):(end-EndDrop)),ConvVec,'valid');

Out=nan(numel(ChangeVec),numel(ChangeVec),nX,size(Outbase,2));
nC=[numel(ChangeVec2) numel(ChangeVec2)]; 
Coeff1=nan(nC);Coeff2=nan(nC);
%ConvAll=cell(nC);
for i=1:numel(ChangeVec)
    for j=1:numel(ChangeVec)
%a1=6;
a2=16;
%b1=1;
b2=1;
c=1/6;
b1=ChangeVec2(j);
a1=ChangeVec(i);
h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));

ConvVec=[0 h(ConvRes*(1:tLength))];

Out(i,j,:,:)=(deconvwnr(X,ConvVec,Pre.ConvLevel));
%disp([i j numel(ChangeVec)])
Coeff1(i,j)=a1;
Coeff2(i,j)=b1;
    end
end
Out=Out./SDmean;

for i=1:numel(PolyTerms)
    PolyReg(i,:)=((Coeff1(:)).^PolyTerms{i}(1)).*((Coeff2(:)).^PolyTerms{i}(2)); %#ok<AGROW>
end

Out=Out(:,:,:,(1+StartDrop):end-EndDrop);

Poly=nan(nX,size(Out,4)-1,numel(PolyTerms));
dPoly=nan(nX,size(Out,4)-1,numel(PolyTerms));
for iP=1:nX
%    disp([iP nX])
%% Do Original
Out2=permute(squeeze(Out(:,:,iP,1:end-1)),[3 1 2]);
Poly(iP,:,:)=(PolyReg'\(Out2(:,:))')';
%% Do Derivatives
Out3=permute(squeeze(Out(:,:,iP,2:end)-Out(:,:,iP,1:end-1)),[3 1 2]);

if strcmpi(NormDeriv,'y')
Out3=Out3./std(Out3,[],1);
end
dPoly(iP,:,:)=(PolyReg'\(Out3(:,:))')';
end
dPolyC{iXcell}=dPoly;
PolyC{iXcell}=Poly;
OrigTime{iXcell}(1,:)=(1+StartDrop):(size(X,2)-(1+EndDrop));
OrigTime{iXcell}(2,:)=iXcell;
end

SDmeanPre=permute(SDmean,[3 1 2]);
OutStruct.SDpolyPre=(PolyReg'\(SDmeanPre(:,:))')';


SDmeanPost=permute(SDmeanPost,[3 1 2]);
OutStruct.SDpoly=(PolyReg'\(SDmeanPost(:,:))')';



Poly=[PolyC{:}];
dPoly=[dPolyC{:}];
PolyTerms=Uncellfun(@(xx)(xx'),PolyTerms);
PolyMat=[PolyTerms{:}];

OutStruct.OrigTime=OrigTime;
OutStruct.StartDrop=StartDrop;
OutStruct.EndDrop=EndDrop;

OutStruct.Poly=Poly;
OutStruct.dPoly=dPoly;
OutStruct.PolyMat=PolyMat;


%% Make Combined Auto-corr across sessions
tmpCount=0;
tmpMat=zeros(size(OutStruct.AutoCorr{1}));
for iXcell=1:numel(X0)
    nT=size(OutStruct.ReMean{iXcell},2);
    tmpMat=tmpMat+OutStruct.AutoCorr{iXcell}*nT;
    tmpCount=tmpCount+nT;
end
OutStruct.CombAutoCorr=tmpMat/tmpCount;
OutStruct.TR=Pre.TR;
end
