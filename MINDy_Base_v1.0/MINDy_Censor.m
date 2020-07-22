function[X1,dX1]=MINDy_Censor(ooP,X,dX)
if ~iscell(ooP)
    ooP={ooP};
end
tempEmpt=(cellfun(@isempty,ooP,'UniformOutput',0));
ooP=ooP(~[tempEmpt{:}]);
nSet=cellfun(@(xx)(numel(xx.GoodFrames)),ooP,'UniformOutput',0);
nSet=cumsum([0 [nSet{:}]]);
IndX=cell(1,numel(ooP));
for i=1:numel(ooP)
IndX{i}=(nSet(i)+1):nSet(i+1);
end
dX1=cell(1,numel(ooP));
X1=cell(1,numel(ooP));
for i=1:numel(ooP)
    for j=1:numel(IndX{i})
        X{IndX{i}(j)}=X{IndX{i}(j)}(:,ooP{i}.GoodFrames{j});        
        dX{IndX{i}(j)}=dX{IndX{i}(j)}(:,ooP{i}.GoodFrames{j});
    end
    dX1{i}=[dX{[IndX{i}]}];
    X1{i}=[X{[IndX{i}]}];
end
end