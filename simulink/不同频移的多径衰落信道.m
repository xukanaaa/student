clear all
ts=0.0025;%信号采样时间间隔
t=0:ts:10-ts;%时间矢量
fs=1/ts;df=fs/length(t);%fft的频率分辨率
msg=randi([-3,3],100,1);%生成消息序列（100行1列）
msg1=msg*ones(1,fs/10);%把生成的随机数拓展成100行40列的矩阵，10/ts/100=fs/10
msg2=reshape(msg1.',1,length(t));%将矩阵重构,变成抽样信号形式，最后形式为-3-3-3-3,...,-2-2-2,...,-1-1-1-1,...
Pm=fft(msg2)/fs;
f=-fs/2:df:fs/2-df;%对抽样信号做傅里叶变换
subplot(3,1,1);plot(f,fftshift(abs(Pm)));title('消息信号频谱');
A=4;fc=100;Sam=(A+msg2).*cos(2*pi*fc*t);%以调信号表达式
Pam=fft(Sam)/fs;%一条信号频谱
subplot(3,1,2);plot(f,fftshift(abs(Pam)));title('已调信号频谱');
axis([-200 200 0 23]);
Pc=sum(abs(Sam).^2)/length(Sam);%已调信号功率
Ps=Pc-A^2/2;%消息信号功率
eta=Ps/Pc;%调制效率
y=awgn(Sam,20,'measured');%信号通过AWGN信道
dems2=abs(hilbert(y))-A;%包络检波，并且去掉直流分量
subplot(3,1,3);plot(t,dems2);title('解调信号');