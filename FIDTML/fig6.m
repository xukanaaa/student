clear all
thetaMSE_N=zeros(1,100000);
thetaCRLB_N=zeros(1,100000);
thetamse_N=zeros(3,27);
thetacrlb_N=zeros(3,27);
thegemaf=[0.04e-6 0.16e-6 0.32e-6];

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
            F(1,i)=v/(3e8)+(1/yitaAB-1)+thegemaf(1,n)*randn;
            F(2,i)=v/(3e8)-(1/yitaAB-1)+thegemaf(1,n)*randn;
        end
 %f矩阵
        f=zeros(2*N,1);
        for i=1:2:2*N-1
            f(i,1)=F(1,(i+1)/2)+1;
            f(i+1,1)=F(2,(i+1)/2)-1;
        end
 %估计v和yitaAB
        guji=inv(B'*B)*B'*f;
        yitaABhat=1/guji(2,1);
        vhat=(3e8)*guji(1,1);
 %C矩阵
        C=zeros(2*N,2);
        for i=1:2:2*N-1
            C(i,1)=-1;C(i,2)=-1;C(i+1,1)=1;C(i+1,2)=-1;
        end
 %T矩阵
        T=zeros(4,N+1);
        for i=1:N
            T(2,i)=yitaAB*(T(1,i)+taoBA+v*(T(1,i)-...
                   T(1,1))/(3e8-v)+(10e-6)*randn)+thetaAB;
            T(3,i)=T(2,i)+0.2e-3;
            T(4,i)=T(3,i)/yitaAB+taoBA+v*(T(3,i)-...
                   T(2,1))/(3e8*yitaAB)-thetaAB/yitaAB+(10e-6)*randn;
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
        thetahat=inv(C'*C)*C'*J;
        thetaABhat=yitaABhat*thetahat(2,1);
        mse=(thetaABhat-thetaAB)^2;
        thetaMSE_N(1,L)=mse;
 %CRLB仿真中的元素
        a=0;b=0;
        for i=1:N
            a=a-T(1,i)-T(2,i)+T(3,i)-T(2,1)+T(3,i)-2*thetaAB;
            b=b-T(1,i)-T(2,i)+T(3,i)-T(2,1)-T(3,i)+2*thetaAB;
        end
%单次仿真的theta的CRLB
        crlb=((thegemaf(1,n)^2)/(4*N)*(a*a)+(thegemaf(1,n)^2)/(4*N)*(b*b)+...
               2*N*(10e-6)^2)/(4*N*N*(1/yitaAB+...
               (sum((thegemaf(1,n))*randn(1,10))-sum((thegemaf(1,n))*...
               randn(1,10)))/(2*N)));
         thetaCRLB_N(1,L)=crlb;
    end
    x=mean(thetaMSE_N(1,:));
    y=mean(thetaCRLB_N(1,:));
    thetamse_N(n,(N+2)/3)=x;
    thetacrlb_N(n,(N+2)/3)=y;
end
    
end
h=semilogy(1:3:80,thetamse_N(1,:),'r*',1:3:80,thetamse_N(2,:),'k<',...
    1:3:80,thetamse_N(3,:),'bs',1:3:80,thetacrlb_N(1,:),'r-',1:3:80,...
    thetacrlb_N(2,:),'k-',1:3:80,thetacrlb_N(3,:),'b-');
set(h,'linewidth',1)
grid on
xlabel('消息交换轮次');ylabel('相位偏移MSE');title('相位偏移MSE')
legend('segemaf=0.04 FIDTML','segemaf=0.16 FIDTML','segemaf=0.32 FIDTML','segemaf=0.04 CRLB','segemaf=0.16 CRLB','segemaf=0.32 CRLB');