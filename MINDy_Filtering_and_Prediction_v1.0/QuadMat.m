function [Out] = QuadMat(X,P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% computes x^T P x for multiple X's
Out=zeros(1,size(X,2));
PX=P*X;
for i=1:size(X,2)
    Out(i)=X(:,i)'*PX(:,i);
end
    
    
end

