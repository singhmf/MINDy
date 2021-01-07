function[ooP]=MakeMINDyFunction(ooP)

%% Note:
%% Factors sqrt(A.^2+(B(:,1).*X+.5).^2)-sqrt(A.^2+(B(:,1).*X-.5).^2) into 
%% (B(:,1))*(sqrt(A2+(B5P+X).^2)-sqrt(A2+(B5N+X).^2))
%% For FastFun it factors the B(:,1) back into W
if ~isempty(ooP.Param{5})
    W=ooP.Param{5};
else
    W=ooP.Param{1};
end
A=ooP.Param{2};
D=ooP.Param{6};
B=ooP.Param{3};
A2=A.^2;

A2=A2./(B(:,1).^2);
B5P=(B(:,2)+.5)./B(:,1);
B5N=(B(:,2)-.5)./B(:,1);
Wfast=W.*(B(:,1)');
%% Transfer function Psi
ooP.Tran=@(xx)(((B(:,1)).*(sqrt(A2+(xx+B5P).^2)-sqrt(A2+(xx+B5N).^2))));
%% Derivative of Transfer Function w.r.t. X
ooP.dTran=@(xx)((B(:,1).*(((xx+B5P)./sqrt(A2+(xx+B5P).^2))-((xx+B5N)./sqrt(A2+(xx+B5N).^2)))));
%% Full function dx=Wpsi(x)-Dx+c
if (~isempty(ooP.Param{4}))&&(mean((ooP.Param{4})==0)~=1)
    c=ooP.Param{4};
    ooP.FastFun=@(yy)(Wfast*(sqrt(A2+(yy+B5P).^2)-sqrt(A2+(yy+B5N).^2))-D.*yy+c);
else
    ooP.FastFun=@(yy)(Wfast*(sqrt(A2+(yy+B5P).^2)-sqrt(A2+(yy+B5N).^2))-D.*yy);
end
end