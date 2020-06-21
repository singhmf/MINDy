function [Y,BadFrames] = MyRestingPreProcInterpNoHRF(XDat,FiltAmp,ConvLevel,DownSamp,TR,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Filters by Point rather than by Frame!!

%% Default:
%% OLD: XC=MyRestingPreProc(TS,5,.5,.5,1)
%% NEW: XC=MyRestingPreProc(TS,5,.3,.5,.72)
%HRFconv=[0    0.0037    0.0433    0.1210    0.1875    0.2105    0.1926 0.1526    0.1081    0.0690    0.0385    0.0162    0.0008   -0.0093 -0.0153   -0.0182   -0.0187   -0.0175   -0.0154   -0.0129   -0.0103 -0.0079   -0.0058   -0.0042   -0.0029   -0.0020   -0.0013   -0.0008 -0.0005   -0.0003   -0.0002   -0.0001   -0.0001];

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
%    X=zscore(XDat{i}')';

    X=(XDat{i}-nanmean(XDat{i},2))./nanstd(XDat{i},[],2);
%    [~,tFilt]=FilterAmplitude(X,FiltAmp);
%    gg=unique(tFilt);
%    XFilt=X;
%    XFilt(:,gg)=[];
%    goodFrames=setdiff(1:size(X,2),gg);
%    goodFrames
%    gg
    for iParc=1:size(X,1)
        gg=find(abs(X(iParc,:))>FiltAmp);
        BadFrames{i,iParc}=gg;
        goodFrames=setdiff(1:size(X,2),gg);
        X(iParc,gg)=interp1(goodFrames,X(iParc,goodFrames),gg);
    end
    %% Filter nan's
    mm=mean(X,1);
    X(:,isnan(mm))=[];
%    X(:,unique(tFilt))=[];
%    X=DeconvHRF(X,TR,ConvLevel);
%    X=deconvwnr(X,HRFconv,ConvLevel);

if DownSamp~=1
X=PCinterp(X,DownSamp,'s','n');
end
%% REMEMBER THAT I INCLUDED THE SMOOTHING IN HERE!


if isempty(varargin)||~strcmpi(varargin{1},'n')
    disp('Smoothing Data by Default')
X=convn(X,[1 1]/2,'valid');
end

    Y{i}=X;
end
end

