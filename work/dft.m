function Xk=dft(x)
N=length(x);
Xk=zeros(1,N);
X=zeros(1,N);
for k=0:N-1
    for n=0:N-1
        X(n+1)=x(n+1)*exp(-1i*2*pi*k*n/N);
    end
    Xk(k+1)=sum(X);
end
    