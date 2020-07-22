function[Out,b]=FilterAmplitude(X,nSD)
% Filter's X to remove greater than SD number of SD's

X=zscore(X')';
[~,b]=find(abs(X)>nSD);
Out=X;%(:,setdiff();
Out(:,unique(b))=[];
end