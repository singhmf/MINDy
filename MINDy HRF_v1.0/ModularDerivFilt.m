function[BadFrames]=ModularDerivFilt(dY,ParStr)
if ParStr.NormDerivFiltAmp~=0
    if strcmpi(ParStr.ParametricFilt(1),'y')
    BadNorm=find(abs(zscore(sqrt(sum(dY.^2,1))))>ParStr.NormDerivFiltAmp);
    else        
    BadNorm=find(abs(NonParametricZ(sqrt(sum(dY.^2,1))))>ParStr.NormDerivFiltAmp);
    end
else
    BadNorm=[];
end

if ParStr.MahalDerivFiltAmp~=0
    tmp=sqrt(QuadMat(dY,pinv(cov(dY(:,~or(isnan(sum(dY,1)),isinf(sum(dY,1))))'))));
    
    if strcmpi(ParStr.ParametricFilt(1),'y')        
    BadMahal=find(abs(zscore(tmp))>ParStr.MahalDerivFiltAmp);
    else
    BadMahal=find(abs(NonParametricZ(tmp))>ParStr.MahalDerivFiltAmp);
    end
else
    BadMahal=[];
end
BadFrames={BadNorm,BadMahal};
end
