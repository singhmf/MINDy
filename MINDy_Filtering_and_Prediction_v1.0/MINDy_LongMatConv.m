function[Out]=MINDy_LongMatConv(A,B,dim)
%% For vectorized 1d convolution along a given dimension (extended convolution)
szA=size(A,dim);
szB=size(B,dim);
szMax=szA+szB-1;%max(szA,szB);
szMin=1;%min(szA,szB);
%% Zero-Padding
zInd=repmat({':'},1,ndims(A));
zIndA=zInd;zIndA{dim}=(szA+1):szMax;
zIndB=zInd;zIndB{dim}=(szB+1):szMax;
A(zIndA{:})=0;
B(zIndB{:})=0;
%% FFT-based convolution
FF0=ifft(fft(A,[],dim).*fft(B,[],dim),[],dim);
%% Trimming
oInd=zInd;oInd{dim}=szMin:szMax;
Out=FF0(oInd{:});
end