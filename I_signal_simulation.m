clc; clear; close all;
%%  Basics of time series data
Fs = 1000; % sampling rate of 1000Hz
t = 0:1/Fs:15;
n = length(t);

%% Sine waves
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
figure;
plot(t, swave), ylabel('Amplitude'), xlabel('Time(s)'), set(gca, 'ylim', [-max(a)*2.3, max(a)*2.3]), grid minor;

%% Morlet wavelet
f_mw = 1;
sinewave = cos(2*pi*f_mw*t);
w  = 2*( 5/(2*pi*f_mw) )^2;
gaussian = exp( (-(t-5).^2)/w );
mwavelet = sinewave .* gaussian;
figure;
plot(t, mwavelet), ylabel('Amplitude'), xlabel('Time(s)'), grid minor; % set(gca, 'xlim', [-1 10])

%% Box wave, freq = 20Hz, duty = 50%
f_sq = 20;
boxwave = double(sin(2*pi*f_sq*t)>0);
figure;
plot(t, boxwave), grid minor, set(gca, 'xlim', [0, 1], 'ylim', [-.5, 1.5]), ylabel('Amplitude'), xlabel('Time(s)');
