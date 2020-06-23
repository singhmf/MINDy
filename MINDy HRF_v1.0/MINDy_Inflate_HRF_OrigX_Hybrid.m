function[Out]=MINDy_Inflate_HRF_OrigX_Hybrid(ooP,Xorig,Pre,derivKern,doRobust,doNormSD,dropStart,dropEnd,tLength)
%% doNormSD indicates whether to recalculate scaling factor $b$ to take into account 
%% normalization in deconvolved domain.
if ~iscell(Xorig)
    Xorig={Xorig};
end
DecayTerm=cell(1,numel(Xorig));
NonLinTerm=cell(1,numel(Xorig));
dX=cell(1,numel(Xorig));
Deconv=MINDy_CrossHRF_Deconv(Xorig,Pre,ooP,tLength);
if strcmpi(doNormSD(1),'y')
    NewSD=std([Deconv{:}],[],2);
    if size(ooP.Param{3},1)==1
        ooP.Param{3}=repmat(ooP.Param{3},numel(NewSD),1);
    end
    ooP.Param{3}(:,1)=ooP.Param{3}(:,1)./NewSD;
end
Out=ooP;
ooP=MakeMINDyFunction(ooP);
for ii=1:numel(Deconv)
    tmp=MINDy_CrossHRF_Conv(ooP.Param{5}*ooP.Tran(Deconv{ii}),Pre,ooP,tLength);
    NonLinTerm{ii}=tmp(:,(1+dropStart):(end-dropEnd-numel(derivKern)+1));
    DecayTerm{ii}=-ooP.Param{6}.*Xorig{ii}(:,(1+dropStart):(end-dropEnd-numel(derivKern)+1));
    %% Also subtracting off Cz
    dX{ii}=convn(Xorig{ii}(:,(1+dropStart):(end-dropEnd)),derivKern,'valid')-ooP.Param{4};
end

dX=[dX{:}];
DecayTerm=[DecayTerm{:}];
NonLinTerm=[NonLinTerm{:}];
if strcmpi(doRobust(1),'y')
    GLMweights(1,:)=robustfit([NonLinTerm(:) DecayTerm(:)],dX(:),[],[],'off');
    Out.isRobustInflate='y';
else
    GLMweights(1,:)=[NonLinTerm(:) DecayTerm(:)]\(dX(:));
    Out.isRobustInflate='n';
end
Out.GLMweights=GLMweights;
Out.oldDecay=ooP.Param{6};
Out.oldW1=ooP.Param{1};
if mean(GLMweights>0)==1
Out.Param{5}=GLMweights(1)*ooP.Param{5};
Out.Param{6}=GLMweights(2)*ooP.Param{6};
else
    Out.Warning='Didnt use GLMweights b/c coefficient for D is negative';
end

Out.Acorr=DiagCorr(dX',DecayTerm');
Out.Corr=DiagCorr(dX',(GLMweights(1)*NonLinTerm+GLMweights(2)*DecayTerm)');
end