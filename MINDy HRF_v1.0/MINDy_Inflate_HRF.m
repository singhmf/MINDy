function [Out] = MINDy_Inflate_HRF(ooP,HRFout,dXkern,doRobust)
%% DeConv gives deconvolved time-series
Out=ooP;
if ~isfield(ooP,'Tran')
    ooP=MakeMINDyFunction(ooP);
end

if isempty(ooP.Param{4})
   ooP.Param{4}=0;
end

HRF=MINDy_MakeHRF_H1H2(Out.HRF{1},Out.HRF{2});

tCount=0;
%% Line up 
H1=ooP.HRF{1};
H2=ooP.HRF{2};
tVec=HRFout.TR*(0:HRFout.DeconvLength);
HRFkern=HRF(tVec);
DecayTerm=cell(1,numel(HRFout.OrigTime));
NonLinTerm=DecayTerm;
dX=DecayTerm;
for ii=1:numel(HRFout.ReMean)
    dX{ii}=convn(HRFout.ReMean{ii}(:,1:(size(HRFout.OrigTime{ii},2)-HRFout.DeconvLength)),dXkern,'valid');
%% Remove the part due to C
    dX{ii}=dX{ii}-ooP.Param{4}.*sum(HRFkern,2);
    HRFcomb=permute((H1.^HRFout.PolyMat(1,:)).*(H2.^HRFout.PolyMat(2,:)),[1 3 2]);
    X0=sum(HRFcomb.*HRFout.Poly(:,(1+tCount):(tCount+size(dX{ii},2)+numel(tVec)-1),:),3);
    tCount=tCount+size(HRFout.OrigTime{ii},2);
    DecayTerm{ii}=MINDy_LongMatConv(X0,HRFkern,2);
    DecayTerm{ii}=-ooP.Param{6}.*DecayTerm{ii}(:,numel(tVec):(end+1-numel(tVec)));
    NonLinTerm{ii}=MINDy_LongMatConv(ooP.Param{5}*ooP.Tran(X0),HRFkern,2);    
    NonLinTerm{ii}=NonLinTerm{ii}(:,numel(tVec):(end+1-numel(tVec)));
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
if GLMweights(2)>0
Out.Param{5}=GLMweights(1)*ooP.Param{5};
Out.Param{6}=GLMweights(2)*ooP.Param{6};
else
    Out.Warning='Didnt use GLMweights b/c coefficient for D is negative';
end
Out.Corr=DiagCorr(dX',(GLMweights(1)*NonLinTerm+GLMweights(2)*DecayTerm)');
end