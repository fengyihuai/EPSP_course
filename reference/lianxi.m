clc;
clear all;

%小波去燥后的信号与傅里叶变换去燥后的信号比较

snr = 4;
[xref, x] = wnoise(1, 11, snr);
xref = xref(1:2000);
x = x(1:2000);


%用全局默认阈值进行去燥处理
[thr, sorh, keepapp] = ddencmp('den', 'wv', 'x'); %获取全局默认阈值
xd = wdencmp('gbl', 'x', 'sym8', 3, thr, sorh, keepapp); %利用全局默认阈值对信号进行去燥

figure;
subplot(231);plot(xref); % 画出方波信号的图
title('原始信号');
subplot(234); plot(x); %画出含有噪声信号的图
title('含有噪声的原始信号')


%用傅里叶变换对方波信号和噪声信号进行频谱分析
dt = 1/(2^11);  % 时域分辨率
fs = 1/dt;    %频域分辨率
df = fs/2000;

xxref = fft(xref);   %对原始信号做快速傅里叶变换
xxref = fftshift(xxref);
xxref = abs(xxref);

xx = fft(x);         %对噪声信号做快速傅里叶变换
xx = fftshift(xx);
xx = abs(xx);

ff = -1000*df:df:1000*df - df;  %设置频率轴
subplot(232);plot(ff, xxref);   
title('原始信号的频谱图');
subplot(235);plot(ff, absxx);
title('含噪声信号的频谱图');

%进行低通滤波，滤波频率为0~200的相对频率
indd2 = 1:800;
xx(indd2) = zeros(size(indd2));
indd2 = 1201:2000;
xx(indd2) = zeros(size(indd2));

xden = ifft(xx);
xden = abs(xden);
subplot(233);plot(xd);
title('小波分析去噪后的信号');
subplot(233);plot(xden);
title('傅里叶分析去燥后的信号')





