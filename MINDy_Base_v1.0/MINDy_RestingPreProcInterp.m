function [Y,BadFrames] = MINDy_RestingPreProcInterp(XDat,FiltAmp,ConvLevel,DownSamp,TR,varargin)

if isempty(varargin)
    tLength=30;
else
    tLength=varargin{1};
end


if iscell(XDat)
    nSub=numel(XDat);
else
nSub=size(XDat,3);
yy=XDat;XDat=cell(1,nSub);
for k=1:nSub
    XDat{k}=yy(:,:,k);
end
end
Y=cell(1,nSub);
for i=1:nSub
    disp([i nSub])
    sdX=nanstd(XDat{i},[],2);
    meanX=nanmean(XDat{i},2);
    X=zscore(XDat{i}')';
    for iParc=1:size(X,1)
        gg=(abs(X(iParc,:))>FiltAmp);
        BadFrames{i,iParc}=find(gg); %#ok<AGROW>
        X(iParc,gg)=interp1(find(~gg),X(iParc,~gg),BadFrames{i,iParc});
    end
    %% Filter nan's
    mm=mean(X,1);
    X(:,isnan(mm))=[];
    
X=(X.*sdX)+meanX;

ConvRes=TR;WeinerNoise=ConvLevel;
%% Deconvolves HRF using canonical form in SPM
%% ConvLength gives spacing between points (not end-time).

a1=6;
a2=16;
b1=1;
b2=1;
c=1/6;

h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));

ConvVec=h(ConvRes*(0:tLength));

X=deconvwnr(X,ConvVec,WeinerNoise);

if DownSamp~=1
X=PCinterp(X,DownSamp,'s','n');
end
    Y{i}=X;
end
end

