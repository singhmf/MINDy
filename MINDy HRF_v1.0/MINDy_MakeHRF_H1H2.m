function[Out]=MINDy_MakeHRF_H1H2(H1,H2,varargin)
%% Make HRF function from H1H2 vectors or make array using time vector in varargin

a1=H1;
a2=16;
b1=H2;
b2=1;
c=1/6;

h=@(t)((t.^(a1-1).*exp(-b1.*t).*(b1.^a1))./gamma(a1)-c*((t.^(a2-1).*exp(-b2.*t).*(b2.^a2))./gamma(a2)));
if isempty(varargin)
    Out=h;
else
    Out=h(varargin{:});
end