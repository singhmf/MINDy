function [ParStr] = MINDy_SetDefaults(ParStr)
%% Sets default values for ParStr
FieldVals={ 'nConv',    2;
            'LinDec',   'n';
            'DerivFiltAmp', 0;
            'doConvo', 'n';
            'SubInd', 0;
            'DerivType','FD';
            'FixD','n';
            'FixC','n';
            'FixB1','n';
            'FixB2','n';
            'doDisp', 'y';
            'NormSp','n';
            'RecW', 'n';
            'RecB2','n';
            'RecE','n';
            'RecA','n';
            'RecCost','n';
            'SpSlope', 0;
            'SpDiag', 0;
            'DiagTran','n';
            'FixDB1', 'n';
            'FixDB2','n';
            'DBrate',1;
            'DiagPC','n';
            'wk1Rate',1;
            'wk2Rate',1;
            'L2SpPls1',0;
            'L2SpPls2',0;
            'SpPls1',0;
            'SpPls2',0;
            'W1Rate',1;
            'DBRate',1;
            'jRate',1;
            'wPC',0;
            'RandScale',.01;
            'NormDeriv','n';
            'NormReg','n';
            'DoDrop','n';
            'pDrop',0;
            'SymmW','n';
            'ConstA','n';
            'ConstD','n';
            'ConstB1','n';
            'ParcMask',[];
            'doTrue','n';
            'doBigger','n';
            'AReg',ParStr.Reg;
            'PruneTime',0;
            'PruneLim',0;
            'PruneSP',0;
            'NormDerivFiltAmp',0;            
            'MahalDerivFiltAmp',0;
            'nRec',0;
            'Dsp',0;
            'Alpha',1;
            'LassoLambda',.5;
            'DRegScale',0;
            'Bstart',ParStr.Bmin;
            'ParametricFilt','n';
            'ENL2',0;
            'L2SpPlsEN',0;
            'RecCorr','n';
            'German',0;
            'AdScale','n';
            'PostInflate','n';
            'InflateParc','n';
            'InflateRobust','y';
            'InflateD','y';
            'H1start',6;
            'H2start',1;
            'H1min',5;
            'H1max',7;
            'H2min',5/6
            'H2max',7/6;
            'RandScale',.001};
        
for i=1:size(FieldVals,1)
    if ~isfield(ParStr,FieldVals{i,1})
        ParStr.(FieldVals{i,1})=FieldVals{i,2};
    end
end
end            


