% (6)选择4个以上上述信号，做离散小波变换，并通过合适地处理小波系数，滤除噪声
% reference: https://blog.csdn.net/EbowTang/article/details/40481393
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
for ilevel = 1:length(a)-1
    swave = cat(2, swave, a(ilevel)*sin(2*pi*f_sw(ilevel)*t(tchunks(ilevel):tchunks(ilevel+1)-1)));
end
swave = cat(2, swave, a(1)*sin(2*pi*f_sw(1)*t(tchunks(ilevel):tchunks(ilevel+1)-1))+...
    a(3)*sin(2*pi*f_sw(3)*t(tchunks(ilevel):tchunks(ilevel+1)-1)));

% Box wave, freq = 20Hz, duty = 50%
f_sq = 20;
boxwave = double(sin(2*pi*f_sq*t)>0);

% Morlet wavelet
f_mw = 5;
sinewave = cos(2*pi*f_mw*t);
w  = 2*( 5/(2*pi*f_mw) )^2;
gaussian = exp( (-(t-5).^2)/w );
mwavelet = sinewave .* gaussian;

% White noise
unoise = rand(1, n); 
nnoise = randn(1, n);

% Mixed signals - normal noise
mix = [2 0.5];
sinnnoise = mix(1)*swave + mix(2)*nnoise;
mwnnoise = mix(1)*mwavelet + mix(2)*nnoise;
bwnnoise = mix(1)*boxwave + 0.2*mix(2)*nnoise;

%% Denoise [sine waves] added with normal noise

nlevel = 4; % Decompositionof the signal at level 6
wname = 'db3';

% wavelet function plot
[phi,psi,xval] = wavefun(wname);
figure, subplot(211), plot(xval,phi), title(strcat(wname, ' Scaling Function')), grid on;
subplot(212), plot(xval,psi), title(strcat(wname, ' Wavelet')), grid on;

% Multilevel 1-D wavelet decomposition - Mallat
[c, l] = wavedec(sinnnoise, nlevel, wname);

% Original wavelet coefficients
figure, subplot(nlevel+2, 1, 1), plot(t, sinnnoise), ylabel('original'),title('DWT coefficients- sine waves + normal noise');
for ilevel = 1:nlevel
    subplot(nlevel+2, 1, ilevel+1), plot(detcoef(c, l, ilevel)), ylabel(['d', num2str(ilevel)]),
    title([num2str(Fs/2^(ilevel+1)), 'Hz - ', num2str(Fs/2^ilevel), 'Hz']);
end
subplot(nlevel+2, 1, nlevel+2), plot(appcoef(c, l, wname, nlevel)), ylabel(['a', num2str(nlevel)]),
title([num2str(0), 'Hz - ', num2str(Fs/2^(nlevel+1)), 'Hz']);

% Threshold calculation - global for each detail coefficients
glb = 1;
% sigma = wnoisest(c, l, glb);
N = numel(sinnnoise);
% thr = sigma*sqrt(2*log(N)); 

denoise_c = c;
denoise_level = 2;
for ilevel = 1:denoise_level
    sigma = wnoisest(c, l, ilevel);
    thr = sigma.*sqrt(2.*log(N)); % threshold
    ilevel_indices = sum(l(1:nlevel+1-ilevel))+1:sum(l(1:nlevel+2-ilevel));
    denoise_c(ilevel_indices) = wthresh(detcoef(c, l, ilevel), 'h', thr); % hard thresholding
end
% a = appcoef(c, l, wname, nlevel);
% denoise_c = [a, denoise_c];

denoise_sine = waverec(denoise_c, l, wname);

% denoised wavelet coefficients
figure, subplot(nlevel+2, 1, 1), plot(t, sinnnoise), ylabel('original'), title('Denoised DWT coefficients- sine waves + normal noise');
for ilevel = 1:nlevel
    subplot(nlevel+2, 1, ilevel+1), plot(detcoef(denoise_c, l, ilevel)), ylabel(['d', num2str(ilevel)]),
    title([num2str(Fs/2^(ilevel+1)), 'Hz - ', num2str(Fs/2^ilevel), 'Hz']);
end
subplot(nlevel+2, 1, nlevel+2), plot(appcoef(denoise_c, l, wname, nlevel)), ylabel(['a', num2str(nlevel)]),
title([num2str(0), 'Hz - ', num2str(Fs/2^(nlevel+1)), 'Hz']);

% Result of denoising
figure, subplot(211), plot(t, sinnnoise), xlim([0 1]), ylabel('amp'), title('Sinewave + normal noise'),  grid minor;
subplot(212), plot(t, denoise_sine), xlim([0 1]), xlabel('time/s'), ylabel('amp'), title('Denoised sinewave'),  grid minor;

%% Denoise [morlet wavelet wave] added with normal noise
nlevel = 6; % Decompositionof the signal at level 4
wname = 'coif2'; % similar

% wavelet function plot
[phi,psi,xval] = wavefun(wname);
figure, subplot(211), plot(xval,phi), title(strcat(wname, ' Scaling Function')), grid on;
subplot(212), plot(xval,psi), title(strcat(wname, ' Wavelet')), grid on;

% Multilevel 1-D wavelet decomposition - Mallat
[c, l] = wavedec(mwnnoise, nlevel, wname);

% Original wavelet coefficients
figure, subplot(nlevel+2, 1, 1), plot(t, mwnnoise), ylabel('original'),title('DWT coefficients- morlet wave + normal noise');
for ilevel = 1:nlevel
    subplot(nlevel+2, 1, ilevel+1), plot(detcoef(c, l, ilevel)), ylabel(['d', num2str(ilevel)]),
    title([num2str(Fs/2^(ilevel+1)), 'Hz - ', num2str(Fs/2^ilevel), 'Hz']);
end
subplot(nlevel+2, 1, nlevel+2), plot(appcoef(c, l, wname, nlevel)), ylabel(['a', num2str(nlevel)]),
title([num2str(0), 'Hz - ', num2str(Fs/2^(nlevel+1)), 'Hz']);

% Threshold calculation - global for each detail coefficients
% glb = 1;
% sigma = wnoisest(c, l, glb);
N = numel(mwnnoise);
% thr = sigma*sqrt(2*log(N)); 

denoise_c = c;
denoise_level = 5;
for ilevel = 1:denoise_level
    sigma = wnoisest(c, l, ilevel);
    thr = sigma.*sqrt(2.*log(N)); % threshold
    ilevel_indices = sum(l(1:nlevel+1-ilevel))+1:sum(l(1:nlevel+2-ilevel));
    denoise_c(ilevel_indices) = wthresh(detcoef(c, l, ilevel), 'h', thr); % hard thresholding
end

denoise_mw = waverec(denoise_c, l, wname);

figure, subplot(nlevel+2, 1, 1), plot(t, mwnnoise), ylabel('original'),title('Denoised DWT coefficients- morlet waves + normal noise');
for ilevel = 1:nlevel
    subplot(nlevel+2, 1, ilevel+1), plot(detcoef(denoise_c, l, ilevel)), ylabel(['d', num2str(ilevel)]),
    title([num2str(Fs/2^(ilevel+1)), 'Hz - ', num2str(Fs/2^ilevel), 'Hz']);
end
subplot(nlevel+2, 1, nlevel+2), plot(appcoef(denoise_c, l, wname, nlevel)), ylabel(['a', num2str(nlevel)]),
title([num2str(0), 'Hz - ', num2str(Fs/2^(nlevel+1)), 'Hz']);

% Result of denoising
figure, subplot(211), plot(t, mwnnoise), xlim([0 9]), ylabel('amp'), title('Morlet wave + normal noise');
subplot(212), plot(t, denoise_mw), xlim([0 9]), xlabel('time/s'), ylabel('amp'), title('Denoised morlet wave');

%% Denoise [box wave] added with normal noise
nlevel = 6; % Decompositionof the signal at level 4
wname = 'haar'; % similar

% wavelet function plot
[phi,psi,xval] = wavefun(wname);
figure, subplot(211), plot(xval,phi), title(strcat(wname, ' Scaling Function')), grid minor;
subplot(212), plot(xval,psi), title(strcat(wname, ' Wavelet')), grid minor;

% Multilevel 1-D wavelet decomposition - Mallat
[c, l] = wavedec(bwnnoise, nlevel, wname);

% Original wavelet coefficients
figure, subplot(nlevel+2, 1, 1), plot(t, bwnnoise), ylabel('original'),title('DWT coefficients- box wave + normal noise');
for ilevel = 1:nlevel
    subplot(nlevel+2, 1, ilevel+1), plot(detcoef(c, l, ilevel)), ylabel(['d', num2str(ilevel)]),
    title([num2str(Fs/2^(ilevel+1)), 'Hz - ', num2str(Fs/2^ilevel), 'Hz']);
end
subplot(nlevel+2, 1, nlevel+2), plot(appcoef(c, l, wname, nlevel)), ylabel(['a', num2str(nlevel)]),
title([num2str(0), 'Hz - ', num2str(Fs/2^(nlevel+1)), 'Hz']);

% Threshold calculation - global for each detail coefficients
glb = 1;
sigma = wnoisest(c, l, glb);
N = numel(bwnnoise);
thr = sigma*sqrt(2*log(N)); 

denoise_c = c;
denoise_level = 3;
for ilevel = 1:denoise_level
     % sigma = wnoisest(c, l, ilevel);
     % thr = sigma.*sqrt(2.*log(N)); % threshold
    ilevel_indices = sum(l(1:nlevel+1-ilevel))+1:sum(l(1:nlevel+2-ilevel));
    denoise_c(ilevel_indices) = wthresh(detcoef(c, l, ilevel), 'h', thr); % hard thresholding
end

denoise_bw = waverec(denoise_c, l, wname);

figure, subplot(nlevel+2, 1, 1), plot(t, bwnnoise), ylabel('original'),title('Denoised DWT coefficients- sine waves + normal noise');
for ilevel = 1:nlevel
    subplot(nlevel+2, 1, ilevel+1), plot(detcoef(denoise_c, l, ilevel)), ylabel(['d', num2str(ilevel)]),
    title([num2str(Fs/2^(ilevel+1)), 'Hz - ', num2str(Fs/2^ilevel), 'Hz']);
end
subplot(nlevel+2, 1, nlevel+2), plot(appcoef(denoise_c, l, wname, nlevel)), ylabel(['a', num2str(nlevel)]),
title([num2str(0), 'Hz - ', num2str(Fs/2^(nlevel+1)), 'Hz']);

% Result of denoising
figure, subplot(211), plot(t, bwnnoise), xlim([0 1]), ylabel('amp'), title('Boxwave + normal noise');
subplot(212), plot(t, denoise_bw), xlim([0 1]), xlabel('time/s'), ylabel('amp'), title('Denoised boxwave');
