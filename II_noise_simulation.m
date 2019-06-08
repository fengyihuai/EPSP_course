clc; clear; close all;
%% white noise (normally [Gaussian] distributed or uniformly distributed)
Fs = 1000; % sampling rate
t = 0:1/Fs:5;
n = length(t);

unoise = rand(n, 1); nnoise = randn(n, 1);
figure, subplot(211), plot(t, nnoise, t, unoise), title('random noise over time'),
legend('Uniform', 'Normal'), xlabel('Time(s)'), ylabel('Amplitude');
subplot(223), hist(unoise, 200), title('Distribution:uniform'), xlabel('Amplitude bins'), ylabel('Counts');
subplot(224), hist(nnoise, 200), title('Distribution:normal'), xlabel('Amplitude bins'), ylabel('Counts');

%% high-frequency noise

noise = randn(n, 1);
cut_f = 30; % Hz

f = Fs*(-n/2:(n/2-1))/n;

Noise = fftshift(fft(noise));
f = Fs*(-n/2:(n/2-1))/n;
shiftF = floor(n*cut_f/Fs);

% Delete the higher frequency components
% Hp_noise = Noise;
% Hp_noise(floor(n/2-shiftF):floor(n/2+shiftF)) = 0;

% High-pass filtering based on self-constructed gauss kernel
f_filter = (abs(f)>(cut_f-5)).*f;
g_filter = 1 - exp(-abs(f_filter)./(cut_f));  
Hp_noise = Noise.*g_filter';
hp_noise = abs(ifft(fftshift(Hp_noise)));

% Frequency band
n_freq = [floor(n/2)+1:floor(n/2+100*n/Fs)]; 

figure, subplot(2,2,1), plot(t, noise), title('white[normal] noise'), xlabel('Time(s)'), ylabel('Amplitude'), grid minor;
subplot(2,2,2), plot(Fs/n*[1:length(n_freq)],abs(Noise(n_freq))), title('noise spectrum'), grid minor;
subplot(2,2,3), plot(t, hp_noise), title('white hp-noise'), xlabel('Time(s)'), ylabel('Amplitude'), grid minor;
subplot(2,2,4), plot(Fs/n*[1:length(n_freq)], abs(Hp_noise(n_freq))), title('hp-noise spectrum'), grid minor;
