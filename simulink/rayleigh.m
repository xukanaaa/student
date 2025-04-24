function [h]=rayleigh(fd,t)
%fd：信道的最大多普勒频移，单位赫兹
%t:信号的抽样时间序列
%h:为输出的瑞利信道函数，是一个时间函数复序列
N=40;%假设入射波数目
wm=2*pi*fd;
N0=N/4;%每个象限入射波数目，即振荡器数目
Tc=zeros(1,length(t));%信道函数的实部
Ts=zeros(1,length(t));%信道函数的虚部
P_nor=sqrt(1/N0);%归一化功率系数
theta=2*pi*rand(1,1)-pi;%区别每个路径的均匀分布的随机相位
for ii=1:N0
alfa(ii)=(2*pi*ii-pi+theta)/N;%第i条入射波的入射角
fi_tc=2*pi*rand(1,1)-pi;
fi_ts=2*pi*rand(1,1)-pi;
Tc=Tc+cos(cos(alfa(ii))*wm*t+fi_tc);
Ts=Ts+cos(sin(alfa(ii))*wm*t+fi_ts);%计算冲激响应函数
end;
h=P_nor*(Tc+1i*Ts);
end