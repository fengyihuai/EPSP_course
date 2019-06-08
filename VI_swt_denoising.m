% (6)选择4个以上上述信号，做离散小波变换，并通过合适地处理小波系数，滤除噪声
clc; clear; close all;
%%  Basics of time series data
Fs = 1000; % sampling rate of 1000Hz
t = 0:1/Fs:15;
n = length(t);

%% Simulation signal
% Sine waves
a = [1 1 1];
f_sw = [20 70 50];
% non-overlapping time chunks
tchunks = round(linspace(1, n, length(a)+1));
swave = 0;
for i = 1:length(a)-1
    swave = cat(2, swave, a(i)*sin(2*pi*f_sw(i)*t(tchunks(i):tchunks(i+1)-1)));
end
swave = cat(2, swave, a(1)*sin(2*pi*f_sw(1)*t(tchunks(i):tchunks(i+1)-1))+...
    a(3)*sin(2*pi*f_sw(3)*t(tchunks(i):tchunks(i+1)-1)));

% Box wave, freq = 20Hz, duty = 50%
f_sq = 20;
boxwave = double(sin(2*pi*f_sq*t)>0);

% Morlet wavelet
f_mw = 1;
sinewave = cos(2*pi*f_mw*t);
w  = 2*( 5/(2*pi*f_mw) )^2;
gaussian = exp( (-(t-5).^2)/w );
mwavelet = sinewave .* gaussian;

% White noise
unoise = rand(n, 1); 
nnoise = randn(n, 1);

% Mixed signals - normal noise
% mix = [0.2 2.0];
% sinnoise1 = mix(1)*swave' + mix(2)*nnoise;
% mwnoise1 = mix(1)*mwavelet' + mix(2)*nnoise;
% bwnoise1 = mix(1)*boxwave' + mix(2)*nnoise;
mix = [2 4];
sinnoise = mix(1)*swave' + mix(2)*nnoise;
mwnoise = mix(1)*mwavelet' + mix(2)*nnoise;
bwnoise = mix(1)*boxwave' + mix(2)*nnoise;

% Mixed signals - uniform noise
sinnoise = mix(1)*swave' + mix(2)*unoise;
mwnoise = mix(1)*mwavelet' + mix(2)*unoise;
bwnoise = mix(1)*boxwave' + mix(2)*unoise;

%% Discrete wavelet transform - signal + normal noise
nlevel = 6;
wname = 'db3';
sig_end = floor(length(t)/ 2^nlevel) * 2^nlevel;
t_swt = t(1:sig_end);

% Sine waves + noise
sinnoiseSWC = swt(sinnoise(1:sig_end), nlevel, wname);
figure, subplot(size(sinnoiseSWC, 1)+1, 1, 1), plot(t, sinnoise),  % xlim([0, 1]) ,
ylabel('original'), title('sinewave + normal noise');
for i = 2:size(sinnoiseSWC, 1)+1
    subplot(size(sinnoiseSWC, 1)+1, 1, i), plot(t_swt, sinnoiseSWC(i-1, :)),  % xlim([0, 1]);
    if i-1<=nlevel
        ylabel(['d', num2str(i-1)]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    else
        ylabel(['a', num2str(i-2), ]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    end
end
xlabel('time/s')

% Box waves + noise
bwnoiseSWC = swt(bwnoise(1:sig_end), nlevel, wname);
figure, subplot(size(bwnoiseSWC, 1)+1, 1, 1), plot(t, bwnoise),  % xlim([0, 1]) ,
ylabel('original'), title('sinewave + normal noise');
for i = 2:size(bwnoiseSWC, 1)+1
    subplot(size(bwnoiseSWC, 1)+1, 1, i), plot(t_swt, bwnoiseSWC(i-1, :)),  % xlim([0, 1]);
    if i-1<=nlevel
        ylabel(['d', num2str(i-1)]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    else
        ylabel(['a', num2str(i-2), ]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    end
end

% Mwavelet waves + noise


%% denoise: Universal threshold - statitic threshold
% reference: https://blog.csdn.net/EbowTang/article/details/40481393

% sine waves + noise
% threshold calculation
% [SWC_o SWD_o] = swt(sinnoise(1:sig_end), nlevel, wname);
% sigma_o = wnoisest(SWD_o);
sigma = median(abs(sinnoiseSWC(1:end-1, :)), 2)/0.6745;
N = numel(sinnoise(1:sig_end));
thr = sigma.*sqrt(2.*log(N));
% thr = sigma .*thselect(SWD', 'sqtwolog')
for ilevel = 1:size(sinnoiseSWC, 1)-1
    indices = find(abs(sinnoiseSWC(ilevel, :))<thr(ilevel));
    sinnoiseSWC_de(ilevel, indices) = 0;
end

figure, subplot(size(sinnoiseSWC_de, 1)+1, 1, 1), plot(t, sinnoise),  % xlim([0, 1]) ,
ylabel('original'), title('sinewave + normal noise');
for i = 2:size(sinnoiseSWC_de, 1)+1
    subplot(size(sinnoiseSWC_de, 1)+1, 1, i), plot(t_swt, sinnoiseSWC_de(i-1, :)),  % xlim([0, 1]);
    if i-1<=nlevel
        ylabel(['d', num2str(i-1)]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    else
        ylabel(['a', num2str(i-2), ]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    end
end
xlabel('time/s')

% box waves + noise
sigma = median(abs(bwnoiseSWC(1:end-1, :)), 2)/0.6745;
N = numel(bwnoise(1:sig_end));
thr = sigma.*sqrt(2.*log(N));
% thr = sigma .*thselect(SWD', 'sqtwolog')
for ilevel = 1:size(bwnoiseSWC, 1)-1
    indices = find(abs(bwnoiseSWC(ilevel, :))<thr(ilevel));
    bwnoiseSWC_de(ilevel, indices) = 0;
end

figure, subplot(size(bwnoiseSWC_de, 1)+1, 1, 1), plot(t, bwnoise),  % xlim([0, 1]) ,
ylabel('original'), title('sinewave + normal noise');
for i = 2:size(bwnoiseSWC_de, 1)+1
    subplot(size(bwnoiseSWC_de, 1)+1, 1, i), plot(t_swt, bwnoiseSWC_de(i-1, :)),  % xlim([0, 1]);
    if i-1<=nlevel
        ylabel(['d', num2str(i-1)]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    else
        ylabel(['a', num2str(i-2), ]), title([num2str(Fs/2^i), 'Hz - ', num2str(Fs/2^(i-1)), 'Hz']);
    end
end
xlabel('time/s');

