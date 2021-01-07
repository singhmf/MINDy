function[Out]=MINDy_MakeStarts(X,ParStr)
if ~isfield(ParStr,'nBack')
    nBack=0;
else
nBack=ParStr.nBack;
end
nStep=ParStr.nStep;
if ~iscell(X)
    X={X};
end
CellLength=Uncellfun(@(xx)(size(xx,2)),X);
CellLength=[0 [CellLength{:}]];
nC=numel(X);
Out=cell(1,nC);
for i=1:nC
Out{i}=sum(CellLength(1:i))+((nBack+1):(CellLength(i+1)-nStep));
end
Out=[Out{:}];
end
