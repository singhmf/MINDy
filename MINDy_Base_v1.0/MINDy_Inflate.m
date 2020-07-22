function [Out] = MINDy_Inflate(ooP,X,dX,doRobust)
Out=ooP;
if ~isfield(ooP,'Tran')
    ooP=MakeMINDyFunction(ooP);
end

if iscell(dX)
    dX=[dX{:}];
end
if iscell(X)
    X=[X{:}];
end
if isempty(ooP.Param{4})
   ooP.Param{4}=0;
end

DecayTerm=(-ooP.Param{6}).*X;
NonLinTerm=ooP.Param{5}*ooP.Tran(X);
%% Remove the C component
dX=dX-ooP.Param{4};

if strcmpi(doRobust(1),'y')
    GLMweights(1,:)=robustfit([NonLinTerm(:) DecayTerm(:)],dX(:),[],[],'off');
    Out.isRobustInflate='y';
else
    GLMweights(1,:)=[NonLinTerm(:) DecayTerm(:)]\[dX(:)];
    Out.isRobustInflate='n';
end
Out.GLMweights=GLMweights;
Out.oldDecay=ooP.Param{6};
Out.oldW1=ooP.Param{1};
if GLMweights(2)>0
Out.Param{5}=GLMweights(1)*ooP.Param{5};
Out.Param{6}=GLMweights(2)*ooP.Param{6};
else
    Out.Warning='Didnt use GLMweights b/c coefficient for D is negative';
end
end