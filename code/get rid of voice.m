%��ȡ�źŵ�matlab����
[x,Fs,bits]=wavread('mail.wav');%�����źţ������ʺͲ���λ����
%x=x(:,1);                      %���������ź���˫��������ֻȡ���������������ɼ�x=x(:,1)��x=x(:,2)
sound(x,Fs,bits);               %��ԭʼ�źŵ��������лط�
X=fft(x,4096);                  %���źź������п�����ɢ����Ҷ�任����
magX=abs(X);                    %���ֵ
angX=angle(X);                  %����λ
figure(1)                       %��ͼ
subplot(221);plot(x);title('ԭʼ�źŲ���');
subplot(222);plot(X); title('ԭʼ�ź�Ƶ��');
subplot(223);plot(magX);title('ԭʼ�źŷ�ֵ');
subplot(224);plot(angX);title('ԭʼ�ź���λ');
%����matlab����
fs=4096;                        %����Ƶ��ȡֵ
f=fs*(0:511)/1024;              %�Բ���Ƶ�ʽ���ȡֵ
t=0:1/1400:(size(x)-1)/1400;    %���������źŵĵ�����������ԭʼ�ź���ͬ
d=0.004*randn(size(x));         %��randn������ֵΪ0����Ϊ1����̬�ֲ�������
x1=x+d;                         %��������
sound(x1,4096);                 %���ż�������������ź�
y1=fft(x,4096);
y2=fft(x1,4096);                %�Լ����źź������п�����ɢ����Ҷ�任����
figure(2)                       %��ͼ
plot(t,x1)
title('�������ź�');
xlabel('time n');
ylabel('��ֵ n');
figure(3)
subplot(2,1,1);plot(f,abs(y1(1:512)));title('ԭʼ�ź�Ƶ��');
xlabel('Hz');ylabel('��ֵ');
subplot(2,1,2);plot(f,abs(y2(1:512)));title('�������ź�Ƶ��');
xlabel('Hz');ylabel('��ֵ'); 
%ȥ��matlab����
N=10;wc=0.3;
[b,a]=butter(N,wc);             %���NΪ10�׵ĵ�ͨ�˲�����wcΪ����0.3dB��ԵƵ�ʣ�����b��a�ֱ��ʾϵͳ�����ķ��ӡ���ĸ����ʽ��ϵ��
X=fft(x);
figure(4)                       %��ͼ
subplot(321);plot(x);title('�˲�ǰ�źŵĲ���');
subplot(322);plot(X);title('�˲�ǰ�źŵ�Ƶ��');
y=filter(b,a,x1);               %���������˲��������ݽ����˲�
Y=fft(y);
subplot(323);plot(y);title('IIR�˲����źŵĲ���');
subplot(324);plot(Y);title('IIR�˲����źŵ�Ƶ��');
sound(y)                        %IIR�˲���Ļָ��ź������ط�
z=fftfilt(b,x1);                %����FFT���ص���ӷ������ݽ����˲�
Z=fft(z);
subplot(325);plot(z);title('FIR�˲����źŵĲ���');
subplot(326);plot(Z);title('FIR�˲����źŵ�Ƶ��');
sound(z)                        %IIR�˲���Ļָ��ź������ط�
