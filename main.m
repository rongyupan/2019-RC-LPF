clc
clear 
close all
A=4;        % A 幅度值
fs=5000000; % fs 采样率//最后转化为频域的最远距离最大频率//采样的频率
F=10000;    % F=10KHz频率,小于采样率的一半（奈奎斯特）//信号频率
 
N=5000;     % N 采样个数//取了多少个点//
            % 采样频率与采样个数没有什么必然关系，采样个数乘上采样频率的倒数（时间间隔）即为总的时间长度
dt=1/fs;    % 时间间隔//决定了时间点与点之间的距离，当然随便指定一个数也行，只不过使信号的时域上的长度不同
t=0:dt:(N-1)*dt;    % 时间向量

% 通带
F1=5.0*10^3;
F2=10.0*10^3;
F3=14.0*10^3;

% 过渡带
F4=20*10^3;
F5=30*10^3;
F6=50*10^3;

% 阻带
F7=107.67*10^3;
F8=218.54*10^3;
F9=329.55*10^3;

% 截止频率
F10=15.6*10^3;

%------------------------信源产生-------------------------------------------
dataSourceType=2;   % 波形类型 0表示正弦波，1表示三角波，2表示方波
switch dataSourceType
    case 0
  %      y=A*sin(2*pi*F1*t);
  %      y=A*sin(2*pi*F2*t);
  %      y=A*sin(2*pi*F3*t);
  %      y=A*sin(2*pi*F4*t);
  %      y=A*sin(2*pi*F5*t);
  %      y=A*sin(2*pi*F6*t);
  %      y=A*sin(2*pi*F7*t);
  %      y=A*sin(2*pi*F8*t);
        y=A*sin(2*pi*F9*t);
  %      y=A*sin(2*pi*F10*t);         
    case 1
        y=A*sawtooth(2*pi*F*t,0.5);   %三角波
    case 2
        y=A*square(2*pi*F*t,50);      %方波
       
    otherwise
end

%--------------------画出原始输入信号的时域与频域图像----------------------
figure(1);
subplot(2,1,1);
plot(t,y);
axis([-inf,inf,-2,2]);
title('信号的时域波形');%对原始图像进行时域画图
xlabel('时间/s');
ylabel('电压/v');
h=fft(y,5000);       %快速傅里叶变换
h_d=abs(fftshift(h));%使频域图像中间为零
f=(-N/2:N/2-1)*fs/N; %将取得时间上的点转化为频率上的点
subplot(2,1,2);
plot(f,h_d/N);    %画原始图像频域上图
%axis([-4*10^4 4*10^4 0 1]);
title('信号的频域波形');
xlabel('频率/hz');
ylabel('电压/v');

%---------------------------滤波器部分------------------------------------
%-----------------------滤波器本身特性（波特图）----------------------------
r=10.2*10^3;         %电阻阻值（ome）
c=1000*10^(-12);     %电容（f）
w=2*pi*f;            %角频率
Para=r*c*1i;

% Preallocae space for S
S=zeros(1,length(f));
P=zeros(1,length(f));
for n=1:length(f)
   S(n)=abs(1/(1+Para*w(n)));           %滤波器对于不同频率幅值衰减系数
   P(n)=angle(1/(1+Para*w(n)))*180/pi;  %滤波器对于不同频率相移系数,*180/pi表示为角度制
end
figure(2);
subplot(2,1,1);
plot(f,S,'r');      %幅值曲线
title('幅值衰减特性');
xlabel('频率/hz');
ylabel('增益倍数');
subplot(2,1,2);
plot(f,P,'blue');   %相位曲线
title('相位特性');
xlabel('频率/hz');
ylabel('相位');
Func=tf(1,[r*c,1]); %系统的传递函数
%调整波特图显示
p=bodeoptions;
p.Grid='on';        %网格
p.FreqUnits = 'Hz'; %横坐标以Hz单位显示
figure(3);
bode(Func,p);         %系统的波特图
title('幅频特性');

%--------------------信号通过滤波器----------------------------------
[yout,tout] = lsim(Func,y,t);%以方波y为输入通过低通滤波函数后输出信号图像
%-----------------------时域图像------------------------------------
figure(4);
subplot(2,1,1);
plot(t,y);
title('原始信号');
axis([0,0.0005,-1.2,1.2]);
xlabel('时间/s');
ylabel('电压/v');
subplot(2,1,2);
plot(tout,yout);
title('滤波后的时域波形');
axis([0,0.0005,-1.2,1.2]);
xlabel('时间/s');
ylabel('电压/v');
%------------------------------频域图像---------------------------------
q=fft(yout,N);
q_d=abs(fftshift(q));
figure(5);
subplot(2,1,1);
plot(f,h_d/N);
title('输入信号的频谱');
axis([-2*10.^5,2*10.^5,0,0.65]);
xlabel('频率/Hz');
ylabel('电压/v');
pks=findpeaks(h_d);

subplot(2,1,2);
plot(f,q_d/N);
title('滤波后的频谱');
axis([-2*10.^5,2*10.^5,0,0.6]);
xlabel('频率/Hz');
ylabel('电压/v');

%-------------------------------自相关函数--------------------------------
figure(6)
subplot(2,1,1);
[Rx,maxlags]=xcorr(y,'unbiased');  %信号的自相关
if fs>10000  %调整时间轴单位及标签,便于观测波形
    plot(maxlags/fs,Rx);
else
    plot(maxlags/fs,Rx);
end
title('输入信号自相关函数');
axis([-0.0005,0.0005,-1.2,1.2]);
xlabel('时间/s');
ylabel('R(t)');
subplot(2,1,2);
[Rx1,maxlags1]=xcorr(yout,'unbiased');  %信号的自相关
if fs>10000  %调整时间轴单位及标签,便于观测波形
    plot(maxlags1/fs,Rx1);
else
    plot(maxlags1/fs,Rx1);
end
title('输出信号自相关函数');
axis([-0.0005,0.0005,-1,1]);
xlabel('时间/s');
ylabel('R(t)');
 

% %-------------------------------功率谱------------------------------------
figure(7);
subplot(2,1,1);
ypsd=(h_d/N).*(conj(h_d)/N);
plot(f,ypsd);
title('输入信号功率谱');
axis([-2*10.^5,2*10.^5,0,0.5]);
xlabel('频率/Hz');
ylabel('W/Hz');
subplot(2,1,2);
youtpsd=q_d.*conj(q_d)/(N*N);%conj共轭
plot(f,youtpsd);
title('输出信号功率谱');
axis([-2*10.^5,2*10.^5,0,0.32]);
xlabel('频率/Hz');
ylabel('W/Hz');

