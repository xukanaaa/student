function [Xf,f]=SpectrumViewer(xt,Fs,freqRange,df)
if nargin==2
    M=0;
    freqRange='onesided';
elseif nargin==3
    M=0
else
    M=Fs/df;
end
L=length(xt);
N=2^(max(nextpow2(M),nextpow2(L)));
X=fft(xt,N);
if strcmp(freqRange,'onesided')
    Xf=zeros(1,floor(N/2-0.5)+1)/L;
    Xf(1)=X(1)/L;
    Xf(2:end)=2*X(floor(1:N/2-0.5)+1)/L;
    f=floor(0:N/2-0.5)*Fs/N;
elseif strcmp(freqRange,'twosided')
    Xf=fftshift(X)/L;
    f=floor(-N/2+0.5:N/2-0.5)*Fs/N;
end
end