%读取信号的matlab程序
[x,Fs,bits]=wavread('mail.wav');%读出信号，采样率和采样位数。
%x=x(:,1);                      %这里语音信号是双声道，如只取单声道作分析，可加x=x(:,1)或x=x(:,2)
sound(x,Fs,bits);               %对原始信号的声音进行回放
X=fft(x,4096);                  %对信号函数进行快速离散傅里叶变换分析
magX=abs(X);                    %求幅值
angX=angle(X);                  %求相位
figure(1)                       %画图
subplot(221);plot(x);title('原始信号波形');
subplot(222);plot(X); title('原始信号频谱');
subplot(223);plot(magX);title('原始信号幅值');
subplot(224);plot(angX);title('原始信号相位');
%加噪matlab程序
fs=4096;                        %采样频率取值
f=fs*(0:511)/1024;              %对采样频率进行取值
t=0:1/1400:(size(x)-1)/1400;    %所加噪声信号的点数调整到与原始信号相同
d=0.004*randn(size(x));         %用randn产生均值为0方差为1的正态分布白噪声
x1=x+d;                         %叠加噪声
sound(x1,4096);                 %播放加噪声后的语音信号
y1=fft(x,4096);
y2=fft(x1,4096);                %对加噪信号函数进行快速离散傅里叶变换分析
figure(2)                       %画图
plot(t,x1)
title('加噪后的信号');
xlabel('time n');
ylabel('幅值 n');
figure(3)
subplot(2,1,1);plot(f,abs(y1(1:512)));title('原始信号频谱');
xlabel('Hz');ylabel('幅值');
subplot(2,1,2);plot(f,abs(y2(1:512)));title('加噪后的信号频谱');
xlabel('Hz');ylabel('幅值'); 
%去噪matlab程序
N=10;wc=0.3;
[b,a]=butter(N,wc);             %设计N为10阶的低通滤波器，wc为它的0.3dB边缘频率，向量b和a分别表示系统函数的分子、分母多项式的系数
X=fft(x);
figure(4)                       %画图
subplot(321);plot(x);title('滤波前信号的波形');
subplot(322);plot(X);title('滤波前信号的频谱');
y=filter(b,a,x1);               %采用数字滤波器对数据进行滤波
Y=fft(y);
subplot(323);plot(y);title('IIR滤波后信号的波形');
subplot(324);plot(Y);title('IIR滤波后信号的频谱');
sound(y)                        %IIR滤波后的恢复信号声音回放
z=fftfilt(b,x1);                %基于FFT的重叠相加法对数据进行滤波
Z=fft(z);
subplot(325);plot(z);title('FIR滤波后信号的波形');
subplot(326);plot(Z);title('FIR滤波后信号的频谱');
sound(z)                        %IIR滤波后的恢复信号声音回放
