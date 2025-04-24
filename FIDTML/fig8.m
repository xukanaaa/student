clear all
thetaMSE_N=zeros(1,100000);
thetaCRLB_N=zeros(1,100000);
thetamse_N=zeros(3,27);
thetacrlb_N=zeros(3,27);
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
            T(2,i)=yitaAB*(T(1,i)+taoBA+v(1,n)*(T(1,i)-...
                   T(1,1))/(3e8-v(1,n))+(10e-6)*randn)+thetaAB;
            T(3,i)=T(2,i)+0.2e-3;
            T(4,i)=T(3,i)/yitaAB+taoBA+v(1,n)*(T(3,i)-...
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
        crlb=(((0.04e-6)^2)/(4*N)*(a*a)+((0.04e-6)^2)/(4*N)*(b*b)+...
               2*N*(10e-6)^2)/(4*N*N*(1/yitaAB+...
               (sum(((0.04e-6))*randn(1,10))-sum(((0.04e-6))*...
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
set(h,'linewidth',1);
grid on
xlabel('消息交换轮次');ylabel('相位偏移MSE');title('相位偏移MSE');
legend('v=1e2 FIDTML','v=1e4 FIDTML','v=1e6 FIDTML','v=1e4 CRLB','v=1e4 CRLB','v=1e6 CRLB');
ylim([1e-13 1e-10]);