%% FINAL VERSION FOR THE PAPER!

%% Changes: Cut Arate in half

ParStr.ParametricFilt='n';
ParStr.L2SpPlsEN=.05;
ParStr.ENL2=0;
ParStr.MahalDerivFiltAmp=0;
ParStr.RandScale=.01;
ParStr.RecCorr='n';ParStr.RecW='n';ParStr.RecA='n';ParStr.RecDdiff='n';
Pre.FiltAmp=5;Pre.ConvLevel=.02;Pre.DownSamp=1; Pre.TR=.72;
Pre.SelfKern='n';
Pre.doSmooth='n';
ParStr.Dmin=.1; ParStr.Scale=Pre.DownSamp;
ParStr.Rate=.000025;
ParStr.Dec=.9;ParStr.Dec2=.95;
ParStr.Reg=.15;%% =1/b1

ParStr.AdiffMin=.1;ParStr.AdiffMax=8;ParStr.ARate=2.5;
%% Remember that Bmin here is really 1/B1 in the equation
ParStr.Bmin=.15; ParStr.Bmax=.15;ParStr.AdiffMin=.1;ParStr.AdiffMax=8;
ParStr.DRate=3.5;ParStr.Dstart=1.75; ParStr.AReg=.2;

ParStr.RecDdiff='n';
ParStr.RecB2='n';ParStr.RecW='n';ParStr.RecA='n';ParStr.RecE='n';
ParStr.SpScale=.075;
ParStr.FixC='y';
ParStr.RecWrate=50;
ParStr.NormSp='y';

%%  MASK
%ParStr.Mask=true(nX);
ParStr.SpDiag=.2;

ParStr.wPC=150;
ParStr.SpPls1=.05;%6;
ParStr.SpPls2=.05;%6;
ParStr.wk1Rate=2.5;%3;
ParStr.wk2Rate=2.5;%3;

ParStr.NormDerivFiltAmp=7;
ParStr.DerivFiltAmp=0;
ParStr.PreDeriv='y';
ParStr.PostW1Scale=1;

ParStr.DRegScale=200;