clear all
yitaMSE_N=zeros(1,100000);
yitaCRLB_N=zeros(1,100000);
yitamse_N=zeros(3,27);
yitacrlb_N=zeros(3,27);
segemaf=[0.04e-6 0.16e-6 0.32e-6];
for n=1:3
for N=1:3:80
    for L=1:100000
        thetaAB=6*rand-3;
        yitaAB=1+((20e-6)*rand-(10e-6));
        taoBA=2*rand-1;
        v=600*rand-300;

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
            F(1,i)=v/(3e8)+(1/yitaAB-1)+segemaf(1,n)*randn;
            F(2,i)=v/(3e8)-(1/yitaAB-1)+segemaf(1,n)*randn;
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
        yitaCRLB_N(1,L)=yitaAB^4*segemaf(1,n)^2/(2*N);
    end
    yitamse_N(n,(N+2)/3)=mean(yitaMSE_N(1,:));
    yitacrlb_N(n,(N+2)/3)=mean(yitaCRLB_N(1,:));
end
end
h=semilogy(1:3:80,yitamse_N(1,:),'r*',1:3:80,yitamse_N(2,:),'k<',...
    1:3:80,yitamse_N(3,:),'bs',1:3:80,yitacrlb_N(1,:),'r-',...
    1:3:80,yitacrlb_N(2,:),'k-',1:3:80,yitacrlb_N(3,:),'b-');
set(h,'linewidth',1.5)
grid on
xlabel('消息交换轮次');ylabel('频率偏移MSE');title('频率偏移MSE')
legend('segemaf=0.04 FIDTML','segemaf=0.16 FIDTML','segemaf=0.32 FIDTML','segemaf=0.04 CRLB','segemaf=0.16 CRLB','segemaf=0.32 CRLB');
ylim([1e-18 1e-12]);