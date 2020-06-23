function[OutStruct]=MINDy_MakeHRFpoly_Res_NoZ_Silent(X0,Pre,ParStr,NormDeriv,ParamRes,varargin)
%% Adds tLength (same as HRFpoly_Mod)

%% Reminder: This version automatically z-scores deconvolved data!!

if isempty(varargin)
    tLength=30;
else
    tLength=varargin{1};
end

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


PolyTerms={[0 0],[1 0],[0 1],[2 0],[1 1],[0 2],[3 0],[2 1],[1 2],[0 3]};

for iXcell=1:numel(X0)
X=X0{iXcell};
%X=zscore(X')';
%find(isnan(X))    
for iParc=1:size(X,1)
        gg=find(abs(X(iParc,:))>FiltAmp);
%        BadFrames{i,iParc}=gg;
%% Remove bad end points (can't interpolate)
gg(ismember(gg,[1 size(X,2)]))=[];
        goodFrames=setdiff(1:size(X,2),gg);
        X(iParc,gg)=interp1(goodFrames,X(iParc,goodFrames),gg);
        OutStruct.BadFrames{iParc}=gg;
end
%find(isnan(X))

%% ConvRes=TR
a1=6;
a2=16;
b1=1;
b2=1;
c=1/6;
nX=size(X,1);

ConvRes=Pre.TR;

h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));
ConvVec=[0 h(ConvRes*(1:tLength))];

ChangeVec=linspace(ParStr.H1min,ParStr.H1max,ParamRes);
ChangeVec2=linspace(ParStr.H2min,ParStr.H2max,ParamRes);

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


OutStruct.Base{iXcell}=Outbase(:,(1+StartDrop):(end-EndDrop));
OutStruct.ReMean{iXcell}=(convn(Outbase(:,(1+StartDrop):(end-EndDrop)),ConvVec,'valid'));

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

    
%ConvAll{i,j}=(ConvVec(1:EndDrop)/std(ConvVec(1:EndDrop)))';

%Out(i,j,:,:)=(deconvwnr(X,ConvVec,Pre.ConvLevel));
Out(i,j,:,:)=(deconvwnr(X,ConvVec,Pre.ConvLevel));
if max(ParamRes)>2
%disp([i j numel(ChangeVec)])
end
Coeff1(i,j)=a1;
Coeff2(i,j)=b1;
    end
end
Out=(Out-mean(Out,4))./std(Out,[],4);

for i=1:numel(PolyTerms)
    PolyReg(i,:)=((Coeff1(:)).^PolyTerms{i}(1)).*((Coeff2(:)).^PolyTerms{i}(2)); %#ok<AGROW>
end

Out=Out(:,:,:,(1+StartDrop):end-EndDrop);

Poly=nan(nX,size(Out,4)-1,numel(PolyTerms));
dPoly=nan(nX,size(Out,4)-1,numel(PolyTerms));
for iP=1:nX
    if max(ParamRes)>2
%    disp([iP nX])
    end
%% Do Original
Out2=permute(squeeze(Out(:,:,iP,1:end-1)),[3 1 2]);
Poly(iP,:,:)=(PolyReg'\(Out2(:,:))')';
%% Do Derivatives
Out3=permute(squeeze(Out(:,:,iP,2:end)-Out(:,:,iP,1:end-1)),[3 1 2]);
%Out3=RegressOut(Out3,Out2,1);
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