function [tr,xx,xr]=sincinterp(x,Ts,Nr)
dT=1/Nr;
N=length(x);
t=0:dT:N-dT;
xr=zeros(1,N*Nr);
kk=2;
for k=1:N
    
    xr=xr+x(k)*sinc(t-(k-1)).*(heaviside(t-(k-1)+kk/2*Nr)-heaviside(t-(k-1)-kk/2*Nr));
end
xx=(1:Nr:N*Nr)=x(1:N);
xx=[xx zeros(1,Nr-1)];
NN=length(xx);
tr=0:Ts/Nr:NN-1;
end