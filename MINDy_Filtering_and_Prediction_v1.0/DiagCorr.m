function[Out,pOut]=DiagCorr(Mat1,Mat2,varargin)
%% Just gives diagonal elements of correlation (vectorized correlation)
%% Output is nx1
%% Changed to only calculate pOut if requested 9.12.18

Out=zeros(size(Mat1,2),1);
pOut=Out;
if sum(size(Mat1)~=size(Mat2))~=0
    error('Matrices should be the same size')
end
if nargout==2
for i=1:size(Mat1,2)
    [Out(i), pOut(i)]=corr(Mat1(:,i),Mat2(:,i),varargin{:});
end
else
    pOut=[];
    for i=1:size(Mat1,2)
        Out(i)=corr(Mat1(:,i),Mat2(:,i),varargin{:});
    end
end
