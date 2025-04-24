function x=idft(Xk)
N=length(Xk);
xn=zeros(1,N);
x=zeros(1,N);
for n=0:N-1
    for k=0:N-1
        xn(k+1)=Xk(k+1)*exp(1i*2*pi*n*k/N);
    end
    x(n+1)=sumx(n);
end
x=(1/N)*x;