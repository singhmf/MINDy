function[ParStr]=MINDy_Naive_Hyper_Scale(nX,varargin)
%% nX=number of parcels
%% varargin is old ParStr...otherwise will generate a new one
%% Rescales MINDy parameters according to Singh et al. (2020)--
%%  -- i.e. rescales priors according to the new expected variation per weight

%% Remark: It's worth looking into whether diagonal's (Recurrents) rescale the same as other weights
%% this would depend upon the space constant: whether they come from all over a parcel (scale with nX) 
%% or are very local (within voxel would mean invariant to parcel size like decay is)

if isempty(varargin)
    ChosenPARSTR_HRF
else
    ParStr=varargin{1};
end
%% Original MINDy paper used Schaefer 400 + 19 subcorticals
SizeRatio=nX/419;
ParStr.Dstart=max(.4,min(1,SizeRatio)*ParStr.Dstart);
%% Low rank dimension (scales with n)
ParStr.wPC=ceil(ParStr.wPC*SizeRatio);
%% L1 sparse components (scales with n)
ParStr.SpScale=ParStr.SpScale*SizeRatio;
ParStr.SpDiag=ParStr.SpDiag*SizeRatio;
%% L1 low-rank components (scales with sqrt(n))
ParStr.SpPls1=ParStr.SpPls1*sqrt(SizeRatio);
ParStr.SpPls2=ParStr.SpPls2*sqrt(SizeRatio);
%% L2 component (scales with n^2)
%% (or you could just set equal 0, this one's not too important)
ParStr.L2SpPlsEN=ParStr.L2SpPlsEN*(SizeRatio)^2;
end
