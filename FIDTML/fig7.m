clear all
yitaMSE_N=zeros(1,100000);
yitaCRLB_N=zeros(1,100000);
yitamse_N=zeros(3,27);
yitacrlb_N=zeros(3,27);
v=[1e2 1e4 1e6];
for n=1:3
for N=1:3:80
    for L=1:100000
        thetaAB=6*rand-3;
        yitaAB=1+((20e-6)*rand-(10e-6));
        taoBA=2*rand-1;
 %B矩阵
        B=zeros(2*N,2);
        for k=1:2:2*N-1
            B(k,1)=1;
            B(k,2)=1;
            B(k+1,1)=1;
            B(k+1,2)=-1;
        end
 %F矩阵
        F=zeros(2,N);
        for i=1:N
            F(1,i)=v(1,n)/(3e8)+(1/yitaAB-1)+(0.04e-6)*randn;
            F(2,i)=v(1,n)/(3e8)-(1/yitaAB-1)+(0.04e-6)*randn;
        end
 %f矩阵
        f=zeros(2*N,1);
        for i=1:2:2*N-1
            f(i,1)=F(1,(i+1)/2)+1;
            f(i+1,1)=F(2,(i+1)/2)-1;
        end
        guji=zeros(2,1);
        guji=inv(B'*B)*B'*f;
        yitahat=1/guji(2,1);
        yitaMSE_N(1,L)=(yitahat-yitaAB)^2;
        yitaCRLB_N(1,L)=yitaAB^4*(0.04e-6)^2/(2*N);
    end
    yitamse_N(n,(N+2)/3)=mean(yitaMSE_N(1,:));
    yitacrlb_N(n,(N+2)/3)=mean(yitaCRLB_N(1,:));
end
end
h=semilogy(1:3:80,yitamse_N(1,:),'r*',1:3:80,yitamse_N(2,:),'k<',...
    1:3:80,yitamse_N(3,:),'bs',1:3:80,yitacrlb_N(1,:),'r-',...
    1:3:80,yitacrlb_N(2,:),'k-',1:3:80,yitacrlb_N(3,:),'b-');
set(h,'linewidth',1);
xlabel('消息交换轮次');ylabel('频率偏移MSE');title('频率偏移MSE');
grid on;
legend('v=1e2 FIDTML','v=1e4 FIDTML','v=1e6 FIDTML','v=1e2 CRLB','v=1e4 CRLB','v=1e6 CRLB');
ylim([1e-18 1e-13]);