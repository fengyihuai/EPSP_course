% （3）采用不同的混合系数，混合上述两种信号，观察结果
clc, clear, close all;
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

%% Mixed signals
mix = [0.2 2.0];
sinnoise1 = mix(1)*swave' + mix(2)*nnoise;
mwnoise1 = mix(1)*mwavelet' + mix(2)*nnoise;
bwnoise1 = mix(1)*boxwave' + mix(2)*nnoise;
figure;
subplot(311), plot(t, sinnoise1), xlabel('Time(s)'), ylabel('Amplitude'), title('swave/normal-noise = 1/10');
subplot(312), plot(t, mwnoise1), xlabel('Time(s)'), ylabel('Amplitude'), title('mwavelet/normal-noise = 1/10');
subplot(313), plot(t, bwnoise1), xlabel('Time(s)'), ylabel('Amplitude'), title('boxwave/normal-noise = 1/10')

mix = [2 0.2];
sinnoise2 = mix(1)*swave' + mix(2)*nnoise;
mwnoise2 = mix(1)*mwavelet' + mix(2)*nnoise;
bwnoise2 = mix(1)*boxwave' + mix(2)*nnoise;
figure;
subplot(311), plot(t, sinnoise2), xlabel('Time(s)'), ylabel('Amplitude'), title('swave/normal-noise = 10');
subplot(312), plot(t, mwnoise2), xlabel('Time(s)'), ylabel('Amplitude'), title('mwavelet/normal-noise = 10');
subplot(313), plot(t, bwnoise2), xlabel('Time(s)'), ylabel('Amplitude'), title('boxwave/normal-noise = 10');
