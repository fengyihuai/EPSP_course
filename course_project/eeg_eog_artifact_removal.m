% 1.	脑电信号中的眼电去除。
% （1）	数据： EEG_FPZ_PZ.mat，两个通道，采样率为2000 Hz。
% （2）	参考文献举例：10.30 R10b Comparative Study of Wavelet-Based.pdf。
% （3）	数据读取与数据片段显示程序如下。
clc, close all, clear;
%% EEG example data read and display
%
% EEG_FPZ_PZ.mat
 
load EEG_FPZ_PZ.mat; % to get two channels of EEG data and sampling rate Fs
 
%% Basic parameters, observation of raw EEG 
% n = Fs * duration2Show;
% tTick = [1:n]/Fs;
raw_eeg = [LSW_FPZ_DS2_200sec; LSW_PZ_DS2_200sec];
num_channels = 2;
n = length(raw_eeg(1, :));
% t = 0:1/Fs:n/Fs; % time label
t = linspace(0, n/Fs, n); % time label
show_dt = [0 10];
show_ft = [-60 60];
 
figure, subplot(211), plot(t, LSW_FPZ_DS2_200sec(1:n), t, LSW_PZ_DS2_200sec(1:n)), xlabel('Time(s)'), ylabel('uv')
xlim(show_dt), grid minor, legend('FPZ', 'PZ'), title('Raw EEG with OA artifact');
% Fourier spectrum
subplot(212), plot(linspace(-Fs/2, Fs/2, length(raw_eeg(1, :))), fftshift(2* abs(fft(raw_eeg(1, :)))/n)), 
hold on, plot(linspace(-Fs/2, Fs/2, length(raw_eeg(2, :))), fftshift(2* abs(fft(raw_eeg(2, :)))/n)), legend('Fpz', 'Pz'),
xlabel('frequency(Hz)'), ylabel('Amp'), xlim(show_ft), grid minor,  title('Raw EEG spectrum')

%% Eog threshold denoising - 4 methods
rs_Fs = 500;
wname = 'sym3'; % the first time: sym3
nlevel = 8;
denoised_levels = [5 6 7 8];

% Downsampling - Fs is converted to resampled Fs
rs_t = linspace(1/Fs, n/Fs, floor(n*rs_Fs/Fs));

for ichannel = 1:num_channels
    rs_eeg(ichannel, :) = resample(raw_eeg(ichannel, :), rs_Fs, Fs);
end


methods = 4;
% uninersal + hard thresholding
denoised_eeg{1} = dwtUtEogDenoise(rs_eeg, rs_Fs, wname, nlevel, denoised_levels, 'h');
% uninersal + soft thresholding
denoised_eeg{2} = dwtUtEogDenoise(rs_eeg, rs_Fs, wname, nlevel, denoised_levels, 's');
% statistical+ hard thresholding
denoised_eeg{3} = dwtStEogDenoise(rs_eeg, rs_Fs, wname, nlevel, denoised_levels, 'h');
% statistical+ soft thresholding
denoised_eeg{4} = dwtStEogDenoise(rs_eeg, rs_Fs, wname, nlevel, denoised_levels, 's');

%% Raw eeg(downsampled) vs. Denoised eeg
figure;
y_shift = 1:2:8;
subplot(211), plot(rs_t, rs_eeg(1, :)), hold on,
for imethod = 1:methods
    plot(rs_t, denoised_eeg{imethod}(1, :) - y_shift(imethod), 'linewidth', 1.2),
end
xlim([0 5]), ylabel('uv'), title('Fpz - raw EEG vs. denoised eeg'), xlabel('time(s)'), grid minor, legend('raw-eeg(resampled)', 'UT-hard', 'UT-soft', 'ST-hard', 'ST-soft');

% figure, plot(rs_t, rs_eeg(1, :), rs_t, denoised_eeg{3}(1, :)), xlim([0 20]), xlabel([0 30]), ylabel('uv'), title('Fpz - ST-hard- raw EEG vs. denoised eeg'),  grid minor;

subplot(212), plot(rs_t, rs_eeg(2, :)), hold on,
for imethod = 1:methods
    plot(rs_t, rs_eeg(2, :), rs_t, denoised_eeg{imethod}(2, :) + y_shift(imethod), 'linewidth', 1.2),
end
 xlim([1 5]), ylabel('uv'), title('Pz - raw EEG vs. denoised eeg'), xlabel('time(s)'), grid minor, legend('raw-eeg(resampled)', 'UT-hard', 'UT-soft', 'ST-hard', 'ST-soft');
%%  eye-region detection, channel 1-2, Fpz, Pz

eog_expand = [0.2 0.2];
[eog_start_point eog_end_point] = eog_detection (rs_t, rs_eeg, rs_Fs, num_channels, eog_expand);

%% Splice of eeg non-eog region - artifact-free segment
% close all;
% The splice of non-eog region in resampled eeg
for ichannel = 1:num_channels
    rs_neog_splice{ichannel} = rs_eeg(ichannel, :);
    for j = length(eog_start_point{ichannel}):-1:1
        indices = eog_start_point{ichannel}(j):eog_end_point{ichannel}(j);
%         raw_eog_section = [raw_eog_section; raw_neog_splice{ichannel}(indices)];
%         raw_neog_splice = [raw_neog_splice; raw_neog_splice{ichannel}(indices)];
        rs_neog_splice{ichannel}(indices) = [];
    end
end

% The splice of non-eog region in denoised eeg(resampled)
for imethod = 1:methods
    for ichannel = 1:num_channels
        denoised_neog_splice{imethod, ichannel} = denoised_eeg{imethod}(ichannel, :);
        for j = length(eog_start_point{ichannel}):-1:1
            indices = eog_start_point{ichannel}(j):eog_end_point{ichannel}(j);
            denoised_neog_splice{imethod, ichannel}(indices) = [];
        end
    end
end

%% NMSE and IM calculation of raw eeg and denoised eeg in the non-eog region
for imethod = 1:methods
    for ichannel = 1:num_channels
    
    % Normalized mean squared error
    NMSE(imethod, ichannel) = nmserr(rs_neog_splice{ichannel}, denoised_neog_splice{imethod, ichannel});
    % Mutual information
    MI(imethod, ichannel) = minfo(rs_neog_splice{ichannel}, denoised_neog_splice{imethod, ichannel});
    
    end
end

figure, bar(NMSE), legend('Fpz', 'Pz'), set(gca,'xticklabel',{'UT-hard', 'UT-soft', 'ST-hard', 'ST-soft'}), ylabel('Normalized mean squared error', 'FontSize',14),
title('Normalized mean squared error', 'FontSize',14), set(gca,'FontSize',20);

figure, bar(MI), legend('Fpz', 'Pz'), set(gca,'xticklabel',{'UT-hard', 'UT-soft', 'ST-hard', 'ST-soft'}), ylabel('Mutual Information', 'FontSize',14),
title('Mutual information', 'FontSize',14), set(gca,'FontSize',20);

% %% Magnitude squared coherence measure
% nfft = 1024;
% noverlap = nfft/2;
% [cxy, f] = mscohere(rs_eeg(1, :), denoised_eeg{imethod}(1, :), hann(nfft), noverlap, nfft, rs_Fs);
% 
% figure, plot(f, cxy)
% 
