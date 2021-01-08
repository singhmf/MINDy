function[Out,Pred,Resid,R2]=MINDy_HRF_CV(Dat1,Dat2,TR,varargin)
%% Dat1=training data (normalized/z-scored) [region x time]: use cells for multiple scans
%% Dat2=testing data (normalized/z-scored) [region x time]: use cells for multiple scans
%% TR is in seconds
%% varargin is {ParStr,Pre}

ChosenPARSTR_HRF;
ParStr.NBatch=15000;ParStr.BatchSz=500;
Pre.TR=TR;
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

%% Resolution in forming the HRF surrogate surface (more is always better, but slower)
Res=15;
%% HRF length
tLength=30; 
%% Periods for clipping the HRF around beginning/end
drop1=10;drop2=10;
Out=MINDy_D2D(Dat1,Pre,ParStr,Res,1,tLength,drop1,drop2);
%% Predict Testing Data: BOLD(t+1)
Pred=MINDy_Predict_NoDeriv_PrePost_Dout(Dat2,Out,Pre,0,1);
if iscell(Pred)
    Pred0=[Pred{:}];
    xTrue=cellfun(@(xx)(xx(:,2:end)),Dat2,'UniformOutput',false);    
    Resid=cellfun(@(xx,yy)(xx-yy),xTrue,Pred,'UniformOutput',false);
    xTrue=[xTrue{:}];
else
    Pred0=Pred;
    xTrue=Dat2(:,2:end);
    Resid=xTrue-Pred;
end
%% Evaluate goodness of fit (R2)
R2=1-sum((xTrue-Pred0).^2,2)./sum((xTrue-mean(xTrue,2)).^2,2);
end
