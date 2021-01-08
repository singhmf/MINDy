function[Out,Pred,Resid,R2]=MINDy_RAW_CV(Dat1,Dat2,varargin)
%% Dat1=training data (normalized/z-scored) [region x time]: use cells for multiple scans
%% Dat2=testing data (normalized/z-scored) [region x time]: use cells for multiple scans
%% varargin is {ParStr,Pre}

ChosenPARSTR_HRF;
ParStr.NBatch=15000;ParStr.BatchSz=500;
%% Number of Parcels
if ~iscell(Dat1)
    nX=size(Dat1,1);
else
    nX=size(Dat1{1},1);
end

if ~isempty(varargin)
    ParStr=varargin{1};
else    
%% Naive rescaling of hyperparameters based upon parcel size
ParStr=MINDy_Naive_Hyper_Scale(nX,ParStr);
end
if numel(varargin)>1
   Pre=varargin{2};
end

if ~iscell(Dat1)
    X0=Dat1(:,1:end-1);
    dX0=convn(Dat1,[1 -1],'valid');
else
    X0=cellfun(@(xx)(xx(:,1:end-1)),Dat1,'UniformOutput',0);
    dX0=cellfun(@(xx)(convn(xx,[1 -1],'valid')),Dat1,'UniformOutput',0);
end
Out=MakeMINDyFunction(MINDy_Base(X0,dX0,Pre,ParStr));

%% Predict Testing Data: X(t+1)
if iscell(Dat2)
    Pred=cell(size(Dat2));
    for ii=1:numel(Dat2)
        %% X(t+1)=X(t)+dX
        Pred{ii}=Dat2{ii}(:,1:end-1)+Out.FastFun(Dat2{ii}(:,1:end-1));
    end
    Pred0=[Pred{:}];
    xTrue=cellfun(@(xx)(xx(:,2:end)),Dat2,'UniformOutput',false);    
    Resid=cellfun(@(xx,yy)(xx-yy),xTrue,Pred,'UniformOutput',false);
    xTrue=[xTrue{:}];
else
    Pred=Dat2(:,1:end-1)+Out.FastFun(Dat2(:,1:end-1));
    Pred0=Pred;
    xTrue=Dat2(:,2:end);
    Resid=xTrue-Pred;
end
%% Evaluate goodness of fit (R2)
R2=1-sum((xTrue-Pred0).^2,2)./sum((xTrue-mean(xTrue,2)).^2,2);
end
