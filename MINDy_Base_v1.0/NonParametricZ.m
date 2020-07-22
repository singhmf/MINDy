function[Out]=NonParametricZ(X)
if isvector(X)
    X=X(:)';
end
Out=(X-nanmedian(X,2));
Out=Out./sqrt(nanmedian(Out.^2,2));
end