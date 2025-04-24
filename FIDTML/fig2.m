clear all
N=10;

%C矩阵
C=zeros(2*N,2);
for i=1:2:2*N-1
    C(i,1)=-1;C(i,2)=-1;C(i+1,1)=1;C(i+1,2)=-1;
end
%B矩阵
B=zeros(2*N,2);
for k=1:2:2*N-1
    B(k,1)=1;
    B(k,2)=1;
    B(k+1,1)=1;
    B(k+1,2)=-1;
end
thegemaf=[0.04e-6;0.16e-6;0.32e-6];
thegemat=4e-6:1e-6:20e-6;
T=zeros(4,N+1);%时间戳矩阵
thetaMSE=zeros(3,17);
thetaCRLB=zeros(3,17);
MSE=zeros(1,100000);%每一个点的MSE和CRLB的100000仿真
CRLB=zeros(1,100000);

for n=1:3
    for j=1:17
        for L=1:100000 
            thetaAB=6*rand-3;
            yitaAB=1+((20e-6)*rand-(10e-6));
            taoBA=2*rand-1;
            v=600*rand-300;
%F矩阵
           F=zeros(2,N);
           for i=1:N
               F(1,i)=v/(3e8)+(1/yitaAB-1)+thegemaf(n,1)*randn;
               F(2,i)=v/(3e8)-(1/yitaAB-1)+thegemaf(n,1)*randn;
           end
%f矩阵
           f=zeros(2*N,1);
           for i=1:2:2*N-1
               f(i,1)=F(1,(i+1)/2)+1;
               f(i+1,1)=F(2,(i+1)/2)-1;
           end
%速度和yitaAB估计
           guji=inv(B'*B)*B'*f;
           yitaABhat=1/guji(2,1);
           vhat=(3e8)*guji(1,1);
%时间戳T矩阵
           for i=1:10
               T(2,i)=yitaAB*(T(1,i)+taoBA+v*(T(1,i)-...
                   T(1,1))/(3e8-v)+thegemat(1,j)*randn)+thetaAB;
               T(3,i)=T(2,i)+0.2e-3;
               T(4,i)=T(3,i)/yitaAB+taoBA+v*(T(3,i)-...
                   T(2,1))/(3e8*yitaAB)-thetaAB/yitaAB+thegemat(1,j)*randn;
               T(1,i+1)=T(1,i)+2;
           end
%J矩阵
           J=zeros(2*N,1);
           for k=1:2:2*N-1
               J(k,1)=T(1,(k+1)/2)-T(2,(k+1)/2)/yitaABhat+...
                   vhat*(T(1,(k+1)/2)-T(1,1))/(3e8-vhat);
               J(k+1,1)=T(4,(k+1)/2)-T(3,(k+1)/2)/yitaABhat-...
                   vhat*(T(3,(k+1)/2)-T(2,1))/(3e8*yitaABhat);
           end
%单次仿真的theta的MSE
           thetahat=inv(C'*C)*C'*J;
           thetaABhat=yitaABhat*thetahat(2,1);
           mse=(thetaABhat-thetaAB)^2;
           MSE(1,L)=mse;
%CRLB仿真中的元素
           a=0;b=0;
           for i=1:N
               a=a-T(1,i)-T(2,i)+T(3,i)-T(2,1)+T(3,i)-2*thetaAB;
               b=b-T(1,i)-T(2,i)+T(3,i)-T(2,1)-T(3,i)+2*thetaAB;
           end
%单次仿真的theta的CRLB
           crlb=((thegemaf(n,1)^2)/(4*N)*(a*a)+(thegemaf(n,1)^2)/(4*N)*(b*b)+...
               2*N*thegemat(1,j)^2)/(4*N*N*(1/yitaAB+...
               (sum((thegemaf(n,1))*randn(1,10))-sum((thegemaf(n,1))*...
               randn(1,10)))/(2*N)));
           CRLB(1,L)=crlb;
        end
%100000次仿真的均值  
        x=mean(MSE(1,:));
        y=mean(CRLB(1,:));
        thetaMSE(n,j)=x;
        thetaCRLB(n,j)=y;
    end
end
twlas=[3e-12 5e-12 7e-12 1e-11 1.4e-11 1.6e-11 2e-11 2.5e-11 3e-11 3.5e-11 4e-11 4.5e-11 5e-11 5.8e-11 6.5e-11 7.3e-11 8e-11];  
    
h=semilogy(4:20,thetaMSE(1,:),'r*',4:20,thetaMSE(2,:),'k<',4:20,thetaMSE(3,:),'b>',...
    4:20,thetaCRLB(1,:),'r-',4:20,thetaCRLB(2,:),'k-',4:20,thetaCRLB(3,:),'b-',4:20,twlas,'m<')
set(h,'linewidth',0.5)
ylim([5e-13 1e-10]);
grid on
xlabel('消息延时segemat(μs)');ylabel('相位偏移MSE');title('相位偏移MSE')
legend('segemaf=0.04 FIDTML','segemaf=0.16 FIDTML','segemaf=0.32 FIDTML','segemaf=0.04 CRLB','segemaf=0.16 CRLB','segemaf=0.32 CRLB','TWLAS')

