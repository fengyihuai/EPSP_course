% (5)选择4个以上上述信号，做连续小波变换，观察时间频率特性
clc, clear, close all
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

%% cwt - matlabR2016b or later
% Old version cwt is not recommended
% reference: https://ww2.mathworks.cn/help/wavelet/ref/cwt.html
path = "\cwt_screenshots\";
prefix = ["swave", "boxwave", "mwavelet", "unoise", "nnoise"];
signal = [swave(:), boxwave(:), mwavelet(:), unoise(:), nnoise(:)];

% cwt, Morse wavelet
wname = 'morse';
for i = 1:size(signal, 2)
    figure, cwt(signal(:, i), wname, Fs);
    title(strcat(prefix(i), " - Analytic CWT using ", wname, " Wavelet"));
    saveas(gcf, strcat(pwd, path, wname, "_", prefix(i), '.png'));
end

% cwt, Morlet wavelet
wname = 'amor';
for i = 1:size(signal, 2)
    figure, cwt(signal(:, i), wname, Fs);
    title(strcat(prefix(i), " - Analytic CWT using ", wname, " Wavelet"));
    saveas(gcf, strcat(pwd, path, wname, "_", prefix(i), '.png'));
end

% cwt, bump wavelet
wname = 'bump';
for i = 1:size(signal, 2)
    figure, cwt(signal(:, i), wname, Fs);
    title(strcat(prefix(i), " - Analytic CWT using ", wname, " Wavelet"));
    saveas(gcf, strcat(pwd, path, wname, "_", prefix(i), '.png'));
end
