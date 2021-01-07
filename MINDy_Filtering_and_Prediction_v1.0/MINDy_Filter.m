function[X0,dY,Out]=MINDy_Filter(X0,dY,ParStr)
if ~iscell(X0)
    X0={X0};
end
if ~iscell(dY)
    dY={dY};
end
DerivFiltAmp=ParStr.DerivFiltAmp;
GoodFrames=cell(1,numel(X0));
BadFrames=cell(1,numel(X0));
for i=1:numel(X0)
GoodFrames{i}(1,:)=1:size(X0{i},2);
GoodFrames{i}(2,:)=i;
Out.oLength{i}=size(dY{i},2);
    if ParStr.DerivFiltAmp~=0
    [~,BadFrames{i}]=FilterAmplitude(dY{i},ParStr.DerivFiltAmp);
    dY{i}(:,BadFrames{i})=[];
    X0{i}(:,BadFrames{i})=[];
    GoodFrames{i}(:,BadFrames{i})=[];
    end
end
    if DerivFiltAmp~=0
    Out.BadFrames=BadFrames;
    else
        Out.BadFrames=cell(1,numel(X0));
    end
    X0=[X0{:}];
    dY=[dY{:}];
    GoodFrames=[GoodFrames{:}];
    
if (ParStr.NormDerivFiltAmp>0)||(ParStr.MahalDerivFiltAmp>0)
Out.BadDerivs=ModularDerivFilt(dY,ParStr);
else
    Out.BadDerivs=cell(1,2);
end
if ~isempty([Out.BadDerivs{:}])
dY(:,unique([Out.BadDerivs{:}]))=[];
X0(:,unique([Out.BadDerivs{:}]))=[];
GoodFrames(:,unique([Out.BadDerivs{:}]))=[];
end

for i=unique(GoodFrames(2,:))
Out.GoodFrames{i}=GoodFrames(1,GoodFrames(2,:)==i);
end
end