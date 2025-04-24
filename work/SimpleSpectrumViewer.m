function [Xf,f]=SimpleSpectrumViewer(xt,Fs)
N=length(xt);
X=fft(xt,N);
Xf=zeros(1,floor(N/2-0.5)+1);
Xf(1)=X(1)/N;
Xf(2:end)=2*X(floor(1:N/2-0.5)+1)/N;
f=floor(0:N/2-0.5)*Fs/N;
end