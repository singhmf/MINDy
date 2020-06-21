function[Out]=PCinterp(X,ResFactor,type,varargin)

%% Performs spline interpolation in PCA space for upsampling data
SampVec=1:ResFactor:size(X,2);

O=zeros(size(X,1),numel(SampVec));
nanFrame=or(isnan(sum(X,1)),isinf(sum(X,1)));
[u,~]=eig(cov(X(:,~nanFrame)'));
if numel(varargin)==1&&strcmpi(varargin(1),'n')
    oo=X;
else
oo=u'*X;
end
if strcmpi(type(1),'s')
for i=1:size(X,1)

    nanMark=find(~isnan(oo(i,:)));
    O(i,:)=spline(nanMark,oo(i,nanMark),SampVec);
end
elseif strcmpi(type(1),'p')
    for i=1:size(X,1)
   %     O(i,:)=pchip(1:size(X,2),oo(i,:),SampVec);
   
    nanMark=find(~isnan(oo(i,:)));
    O(i,:)=pchip(nanMark,oo(i,nanMark),SampVec);
    end
else
    error('Type should be s(pline) or p(chip)')
end
if numel(varargin)==1&&strcmpi(varargin(1),'n')
    Out=O;
else
Out=u*O;
end
end