function [Out] = MINDy_Inflate_HRF_PrePost_NoDeriv(X,ooP,Pre,doRobust)
%% DeConv gives deconvolved time-series
Out=ooP;
if ~isfield(ooP,'Tran')
    ooP=MakeMINDyFunction(ooP);
end

if isempty(ooP.Param{4})
   ooP.Param{4}=0;
end

DecayTerm=cell(1,numel(X));
NonLinTerm=cell(1,numel(X));
Xshift=cell(1,numel(X));
DeconvLength=ooP.HRFout.DeconvLength;
reconvX=cell(1,numel(X));
for ii=1:numel(X)
    deconvX=zeros(size(X{ii},1),size(ooP.HRFout.OrigTime{ii},2));
    Xtmp=X{ii}(:,ooP.HRFout.OrigTime{ii}(1,(2+ceil((DeconvLength-1)/2)):(end-(floor((DeconvLength-1)/2)))));    
    for jj=1:size(X{1},1)
    tmp=deconvwnr(X{ii}(jj,:),ooP.normHRF_kernPre(jj,:),Pre.ConvLevel);    
    tmp=tmp((1+ooP.HRFout.StartDrop):(end-1-ooP.HRFout.EndDrop));
    deconvX(jj,:)=tmp;
    reconvX{ii}(jj,:)=convn(tmp,ooP.normHRF_kern(jj,:),'valid');
    end
    preNL=ooP.Param{5}*ooP.Tran(deconvX);
    for jj=1:size(X{1})
        NonLinTerm{ii}(jj,:)=convn(preNL(jj,:),ooP.normHRF_kern(jj,:),'valid');
    end
    DecayTerm{ii}=-ooP.Param{6}.*reconvX{ii};
    %% Pre-subtract the difference
    Xshift{ii}=Xtmp-reconvX{ii}-ooP.Param{4}.*sum(ooP.normHRF_kern,2);
end


dX=[Xshift{:}];
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
if GLMweights(2)>0
Out.Param{5}=GLMweights(1)*ooP.Param{5};
Out.Param{6}=GLMweights(2)*ooP.Param{6};
else
    Out.Warning='Didnt use GLMweights b/c coefficient for D is negative';
end
Out.Corr=DiagCorr(dX',(GLMweights(1)*NonLinTerm+GLMweights(2)*DecayTerm)');
end