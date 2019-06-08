% (4)上述三类信号对应的傅里叶变换，观察结果
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

% Mixed signals
mix = [0.2 2.0];
sinnoise1 = mix(1)*swave' + mix(2)*nnoise;
mwnoise1 = mix(1)*mwavelet' + mix(2)*nnoise;
bwnoise1 = mix(1)*boxwave' + mix(2)*nnoise;
mix = [2 0.2];
sinnoise2 = mix(1)*swave' + mix(2)*nnoise;
mwnoise2 = mix(1)*mwavelet' + mix(2)*nnoise;
bwnoise2 = mix(1)*boxwave' + mix(2)*nnoise;

%% Fourier spectrum
hz = linspace(-Fs/2, Fs/2, n);

% The spectrum of sine waves, mwavelet, boxwave
swaveFFT = fftshift(fft(swave)/n);
mwaveletFFT = fftshift(fft(mwavelet)/n);
boxwaveFFT = fftshift(fft(boxwave)/n);
figure;
subplot(311), plot(hz, 2*abs(swaveFFT)), set(gca, 'xlim', [-max(f_sw)*1.3 max(f_sw)*1.3], 'ylim', [0 max(a)*1.3]), grid on,
title('Sine waves - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');
subplot(312), plot(hz, 2*abs(mwaveletFFT)), set(gca, 'xlim', [-max(f_mw)*3 max(f_mw)*3]), grid on,
title('Mwavelet - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');
subplot(313), plot(hz, 2*abs(boxwaveFFT), '-o'), set(gca, 'xlim', [-max(f_sq)*4 max(f_sq)*4]), grid on,
title('Boxwave - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');

% The spectrum of white noise 
% unoiseFFT = fftshift(fft(unoise)/n);
unoiseFFT = fftshift(fft(detrend(unoise))/n); % Remove linear trend of uniformlly noise
nnoiseFFT = fftshift(fft(nnoise)/n);
figure;
subplot(211), plot(hz, 2*abs(unoiseFFT(1:length(hz)))), grid on, title('Uniformly noise spectrum'),
xlabel('Frequencies (Hz)'), ylabel('Amplitude');
subplot(212), plot(hz, 2*abs(nnoiseFFT(1:length(hz)))), grid on, title('Normally noise spectrum'),
xlabel('Frequencies (Hz)'), ylabel('Amplitude');

% The spectrum of mixed signal
sinnoise1FFT = fftshift(fft(sinnoise1)/n);
mwnoise1FFT = fftshift(fft(mwnoise1)/n);
bwnoise1FFT = fftshift(fft(bwnoise1)/n);
figure;
subplot(311), plot(hz, 2*abs(sinnoise1FFT)), set(gca, 'xlim', [-max(f_sw)*1.3 max(f_sw)*1.3], 'ylim', [0 max(a)*1.3]), grid on,
title('Sine+noise waves[amp=0.2:2] - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');
subplot(312), plot(hz, 2*abs(mwnoise1FFT)), set(gca, 'xlim', [-max(f_mw)*3 max(f_mw)*3]), grid on,
title('Mwavelet+noise[amp=0.2:2] - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');
subplot(313), plot(hz, 2*abs(bwnoise1FFT), 'linewidth', 1.2), set(gca, 'xlim', [-max(f_sq)*4 max(f_sq)*4]), grid on,
title('Boxwave+noise[amp=0.2:2] - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');

sinnoise2FFT = fftshift(fft(sinnoise2)/n);
mwnoise2FFT = fftshift(fft(mwnoise2)/n);
bwnoise2FFT = fftshift(fft(bwnoise2)/n);
figure;
subplot(311), plot(hz, 2*abs(sinnoise2FFT)), set(gca, 'xlim', [-max(f_sw)*1.3 max(f_sw)*1.3], 'ylim', [0 max(a)*1.3]), grid on,
title('Sine+noise waves[amp=2:0.2] - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');
subplot(312), plot(hz, 2*abs(mwnoise2FFT)), set(gca, 'xlim', [-max(f_mw)*3 max(f_mw)*3]), grid on,
title('Mwavelet+noise[amp=2:0.2] - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');
subplot(313), plot(hz, 2*abs(bwnoise2FFT)), set(gca, 'xlim', [-max(f_sq)*4 max(f_sq)*4]), grid on,
title('Boxwave+noise[amp=2:0.2] - spectrum'), xlabel('Frequencies (Hz)'), ylabel('Amplitude');