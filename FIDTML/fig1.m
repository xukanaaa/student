clear all
MSE=zeros(6,100000);
for L=1:100000
N=10;
thetaAB=6*rand-3;
yitaAB=1+((20e-6)*rand-(10e-6));
taoAB=2*rand-1;
v=600*rand-300;

%B矩阵
B=zeros(2*N,2);
for k=1:2:2*N-1
    B(k,1)=1;
    B(k,2)=1;
    B(k+1,1)=1;
    B(k+1,2)=-1;
end

%不同thegemaf下的F矩阵
F1=zeros(2,N);
for i=1:N
    F1(1,i)=v/(3e8)+(1/yitaAB-1)+(0.04e-6)*randn;
    F1(2,i)=v/(3e8)-(1/yitaAB-1)+(0.04e-6)*randn;
end
F2=zeros(2,N);
for i=1:N
    F2(1,i)=v/(3e8)+(1/yitaAB-1)+(0.16e-6)*randn;
    F2(2,i)=v/(3e8)-(1/yitaAB-1)+(0.16e-6)*randn;
end
F3=zeros(2,N);
for i=1:N
   F3(1,i)=v/(3e8)+(1/yitaAB-1)+(0.32e-6)*randn;
   F3(2,i)=v/(3e8)-(1/yitaAB-1)+(0.32e-6)*randn;
end

%不同的thegemaf下的f矩阵
f1=zeros(2*N,1);
for i=1:2:2*N-1
    f1(i,1)=F1(1,(i+1)/2)+1;
    f1(i+1,1)=F1(2,(i+1)/2)-1;
end
f2=zeros(2*N,1);
for i=1:2:2*N-1
    f2(i,1)=F2(1,(i+1)/2)+1;
    f2(i+1,1)=F2(2,(i+1)/2)-1;
end
f3=zeros(2*N,1);
for i=1:2:2*N-1
    f3(i,1)=F3(1,(i+1)/2)+1;
    f3(i+1,1)=F3(2,(i+1)/2)-1;
end

guji1=zeros(2,1);guji2=zeros(2,1);guji3=zeros(2,1);
guji1=inv(B'*B)*B'*f1;
guji2=inv(B'*B)*B'*f2;
guji3=inv(B'*B)*B'*f3;
yitahat1=1/guji1(2,1);
yitahat2=1/guji2(2,1);
yitahat3=1/guji3(2,1);
MSE(1,L)=(yitahat1-yitaAB)^2;
MSE(2,L)=(yitahat2-yitaAB)^2;
MSE(3,L)=(yitahat3-yitaAB)^2;
MSE(4,L)=(yitaAB^4)*(0.04e-6)^2/(2*N);
MSE(5,L)=(yitaAB^4)*(0.16e-6)^2/(2*N);
MSE(6,L)=(yitaAB^4)*(0.32e-6)^2/(2*N);

end

%100000次实验后的均值
mse1=(mean(MSE(1,:)));
mse2=(mean(MSE(2,:)));
mse3=(mean(MSE(3,:)));
CRLB1=(mean(MSE(4,:)));
CRLB2=(mean(MSE(5,:)));
CRLB3=(mean(MSE(6,:)));
twlas=[1e-14 1.6e-14 2.5e-14 3e-14 4e-14 5.5e-14 6.7e-14 8e-14 1e-13 1.2e-13 1.4e-13 1.5e-13 1.8e-13 2e-13 2.2e-13 2.5e-13 2.8e-13];
h=semilogy(4:20,mse1*ones(1,17),'r*',4:20,mse2*ones(1,17),'ks',4:20,...
    mse3*ones(1,17),'b>',4:20,CRLB1*ones(1,17),'r-',4:20,...
    CRLB2*ones(1,17),'k-',4:20,CRLB3*ones(1,17),'b-',4:20,twlas,'m>');
ylim([1e-17 1e-10])
xlabel('消息时延segemat(μs)');ylabel('频率偏移MSE');
title('频率偏移MSE')
grid on
set(h,'linewidth',1.5)
legend('segemaf=0.04','segemaf=0.16','segemaf=0.32','CRLB 0.04',...
    'CRLB 0.16','CRLB 0.32','TWLAS')


