function [Out,Kern] = MINDy_DeconvHRF(X,ConvRes,WeinerNoise,HRFlength,varargin)
%% Use this version instead of DeconvHRF (Problem with HRF length)
if nargin<4
    HRFlength=30;
end

a1=6;
a2=16;
b1=1;
b2=1;
c=1/6;

h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));

%ConvVec=repmat(h(ConvRes*(1:(size(X,2)))),size(X,1),1);
ConvVec=h(ConvRes*(0:(HRFlength)));
%ConvVec=h(0:.72:25);
%% NOTE THAT HAVE FIXED LENGTH OF CONVVEC BUT TIME DEPENDS ON TR
%%--for the deconvblind code which is now off
if ~isempty(varargin)&&strcmpi(varargin{1},'y')
for i=1:size(X,1)
    [~,yy]=deconvblind(15+zscore(X(i,:)),ConvVec+.3,20);
    cVec=[yy-yy(end)];
    Kern(i,:)=cVec; %#ok<AGROW>
    if i==1
        Out2=deconvwnr(zscore(X(i,:)),cVec,WeinerNoise); 
        Out=zeros(size(X,1),numel(Out2));Out(1,:)=Out2;
    else
        Out(i,:)=deconvwnr(zscore(X(i,:)),cVec,WeinerNoise); 
    end
end
else
    Kern=[];
Out=deconvwnr(X,ConvVec,WeinerNoise);
end
%X=[X X(:,end:-1:1)];
%ConvVec=[ConvVec ConvVec(:,end:-1:1)];


%%
%L1=log(fft(ConvVec,[],2));

%L1=log((L1+L1')/2);

%LX=log(fft(X,[],2));

%LX=log((LX+LX')/2);

%\Out=ifft(exp(LX-L1),[],2);
%Out=ifft(fft([ConvVec],[],2)./fft([X],[],2),[],2);
%Out=Out(:,1:end);
%Out=Out(:,1:size(Out,2)/2);
end

