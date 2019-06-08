clc;
clear all;

%С��ȥ�����ź��븵��Ҷ�任ȥ�����źűȽ�

snr = 4;
[xref, x] = wnoise(1, 11, snr);
xref = xref(1:2000);
x = x(1:2000);


%��ȫ��Ĭ����ֵ����ȥ�ﴦ��
[thr, sorh, keepapp] = ddencmp('den', 'wv', 'x'); %��ȡȫ��Ĭ����ֵ
xd = wdencmp('gbl', 'x', 'sym8', 3, thr, sorh, keepapp); %����ȫ��Ĭ����ֵ���źŽ���ȥ��

figure;
subplot(231);plot(xref); % ���������źŵ�ͼ
title('ԭʼ�ź�');
subplot(234); plot(x); %�������������źŵ�ͼ
title('����������ԭʼ�ź�')


%�ø���Ҷ�任�Է����źź������źŽ���Ƶ�׷���
dt = 1/(2^11);  % ʱ��ֱ���
fs = 1/dt;    %Ƶ��ֱ���
df = fs/2000;

xxref = fft(xref);   %��ԭʼ�ź������ٸ���Ҷ�任
xxref = fftshift(xxref);
xxref = abs(xxref);

xx = fft(x);         %�������ź������ٸ���Ҷ�任
xx = fftshift(xx);
xx = abs(xx);

ff = -1000*df:df:1000*df - df;  %����Ƶ����
subplot(232);plot(ff, xxref);   
title('ԭʼ�źŵ�Ƶ��ͼ');
subplot(235);plot(ff, absxx);
title('�������źŵ�Ƶ��ͼ');

%���е�ͨ�˲����˲�Ƶ��Ϊ0~200�����Ƶ��
indd2 = 1:800;
xx(indd2) = zeros(size(indd2));
indd2 = 1201:2000;
xx(indd2) = zeros(size(indd2));

xden = ifft(xx);
xden = abs(xden);
subplot(233);plot(xd);
title('С������ȥ�����ź�');
subplot(233);plot(xden);
title('����Ҷ����ȥ�����ź�')





