function[Out]=MINDy_CrossHRF_Conv(X,Pre,ooP,tLength)
%% Performs Convolution with model's HRF kernels
ConvRes=Pre.TR;
H1=ooP.HRF{1};
H2=ooP.HRF{2};
if ~iscell(X)
    wasCell='n';
    X={X};
else
    wasCell='y';
end
nX=size(X{1},1);
Out=cell(1,numel(X));
for iScan=1:numel(X)
%a1=6;
a2=16;
%b1=1;
b2=1;
c=1/6;
X0=X{iScan};
nT=size(X0,2);
Out{iScan}=nan(size(X0));
for iChan=1:nX
    a1=H1(iChan);
    b1=H2(iChan);
    h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));
    ConvVec=[0 h(ConvRes*(1:(tLength-1)))];%(nT-1)))];
    Out{iScan}(iChan,:)=(convn(X0(iChan,:),ConvVec,'same'));
end
end
if strcmpi(wasCell,'n')
    Out=Out{1};
end
end    