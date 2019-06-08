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
duration = 200; % in seconds
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
hold on, plot(linspace(-Fs/2, Fs/2, length(raw_eeg(2, :))), fftshift(2* abs(fft(raw_eeg(2, :)))/n)), 
xlabel('frequency(Hz)'), ylabel('Amp'), xlim(show_ft), grid minor,  title('Raw EEG spectrum')

%% Block test
% rs_Fs = 500;
% wname = 'sym5'; % the first time: sym3
% nlevel = 6;
% denoised_levels = [5 6];
% thr_type = 'h';
% 
% denoised_eeg = dwtStEogDenoise(raw_eeg, Fs, wname, nlevel, denoised_levels, thr_type, rs_Fs);


%% Downsampling - 2000Hz is converted to 500Hz
% Attention:resample applies an antialiasing FIR lowpass filter to x and compensates for 
% the delay introduced by the filter.

rs_Fs = 500;
rs_t = linspace(0, duration, floor(n*rs_Fs/Fs));

for ichannel = 1:num_channels
    rs_eeg(ichannel, :) = resample(raw_eeg(ichannel, :), rs_Fs, Fs);
end

%% EOG artifact removal based on wavelet tranform coefficients thresholding
wname = 'sym5'; % the first time: sym3
nlevel = 6;
denoised_levels = [5 6];
thr_type = 'h';

% wavelet function plot
[phi,psi,xval] = wavefun(wname);
figure, subplot(211), plot(xval,phi), title(strcat(wname, ' Scaling Function')), grid on;
subplot(212), plot(xval,psi), title(strcat(wname, ' Wavelet')), grid on;

for ichannel = 1:num_channels
    % Mallat wavelet transform - decomposition and reconstruction of each frequency band
    
    s = rs_eeg(ichannel, :);
    % wavelet decomposition
    [c{ichannel}, l{ichannel}] = wavedec(s, nlevel, wname);
    
    % EOG removal - statistical thresholding calculation
    % remove eog (0-16Hz)
    denoise_c{ichannel} = c{ichannel};
    for ilevel = denoised_levels
        hk = detcoef(c{ichannel}, l{ichannel}, ilevel);
        thr = 1.5*std(hk); % threshold
        ilevel_indices = sum(l{ichannel}(1:nlevel+1-ilevel))+1:sum(l{ichannel}(1:nlevel+2-ilevel));
        denoise_c{ichannel}(ilevel_indices) = iwthresh(hk, thr_type, thr); % hard/soft thresholding
    end
    
    % addition: denoising the approximation coefficients
    hk_a = appcoef(c{ichannel}, l{ichannel}, wname);
    thr = 1.5*std(hk_a); % threshold
    denoise_c{ichannel}(1:l{ichannel}(1)) = iwthresh(hk_a, thr_type, thr); % hard/soft thresholding
    
    denoised_eeg(ichannel, :) = waverec(denoise_c{ichannel}, l{ichannel}, wname);
    
end

%% Original wavelet coefficients  vs. Denoised wavelet coefficients
% 

for ichannel = 1:num_channels
    
    l_l = length(denoised_levels);
    % Original wavelet coefficients
    figure, subplot(l_l+2, 1, 1), plot(rs_t, rs_eeg(ichannel, :)), ylabel('original'),title('Original DWT coefficients');
    for ilevel = denoised_levels
        subplot(l_l+2, 1, find(ilevel==denoised_levels)+1), plot(detcoef(c{ichannel}, l{ichannel}, ilevel)), ylabel(['d', num2str(ilevel)]),
        title([num2str(rs_Fs/2^(ilevel+1)), 'Hz - ', num2str(rs_Fs/2^ilevel), 'Hz']);
    end
    subplot(l_l+2, 1, l_l+2), plot(appcoef(c{ichannel}, l{ichannel}, wname, l_l)), ylabel(['a', num2str(l_l)]),
    title([num2str(0), 'Hz - ', num2str(rs_Fs/2^(nlevel+1)), 'Hz']);
    
    % Denoised wavelet coefficients
    figure, subplot(l_l+2, 1, 1), plot(rs_t, rs_eeg(ichannel, :)), ylabel('original'),title('Denoised DWT coefficients');
    for ilevel = denoised_levels
        subplot(l_l+2, 1, find(ilevel==denoised_levels)+1), plot(detcoef(denoise_c{ichannel}, l{ichannel}, ilevel)), ylabel(['d', num2str(ilevel)]),
        title([num2str(rs_Fs/2^(ilevel+1)), 'Hz - ', num2str(rs_Fs/2^ilevel), 'Hz']);
    end
    subplot(l_l+2, 1, l_l+2), plot(appcoef(denoise_c{ichannel}, l{ichannel}, wname, l_l)), ylabel(['a', num2str(l_l)]),
    title([num2str(0), 'Hz - ', num2str(rs_Fs/2^(nlevel+1)), 'Hz']);
    
end

%% raw eeg vs. denoised eeg
figure, subplot(211), plot(t, raw_eeg(1, :), rs_t, denoised_eeg(1, :)), xlim([0 10]), ylabel('uv'), title('raw EEG vs. denoised eeg'),  grid minor, legend('', '');
subplot(212), plot(t, raw_eeg(2, :), rs_t, denoised_eeg(2, :)), xlim([0 10]), ylabel('uv'), title('raw EEG vs. denoised eeg'),  grid minor;

%% eye-region detection

%% channel 1
eog_detect_l = 4;

[c_d, l_d] = wavedec(rs_eeg(1, :), eog_detect_l, wname);
c_da = appcoef(c_d, l_d, wname);

thr_line = repmat(mean(c_da) + 1.5*std(c_da), 1, length(c_da));

% % dwt长度的反演化因子
% dwt_Ltrans = 2^eog_detect_l - 30*eog_detect_l+2^eog_detect_l - 1 - 16;

t_a = dwt_Ltrans(1:length(c_da), eog_detect_l)/rs_Fs;
figure, plot(t_a, c_da), title('a4'),
hold on, plot(t_a, thr_line), xlim([0 100]);

eog_thr = c_da>thr_line;
eog_expand = [0.5 0.3];
expand_points = floor(eog_expand ./ 2^eog_detect_l .* rs_Fs);
% close all;

eog_thrmerge = eog_thr;
% expand the eog starting point's location
eog_start_point = find([0 diff(eog_thr)]==1);
start_expand = eog_start_point - expand_points(1);
for i = 1:length(eog_start_point)
    if start_expand(i)<0
        start_expand(i) = 1;
    end
    eog_thrmerge(start_expand(i):eog_start_point(i)) =  true;
end

figure, plot(eog_thr, 'Linewidth', 1.5),hold on, plot(eog_thrmerge)
xlim([0 1000]), ylim([0 2])

% expand the eog ending point's location
eog_end_point = find([0 diff(eog_thr)]==-1);
end_expand = eog_end_point + expand_points(2);
for i = 1:length(eog_end_point)
    if end_expand(i)>length(eog_thrmerge)
        end_expand(i) = length(eog_thrmerge);
    end
    eog_thrmerge(eog_end_point(i):end_expand(i)) =  true;
end

figure, plot(eog_thr, 'Linewidth', 1.5), hold on, plot(eog_thrmerge)
xlim([0 1000]), ylim([0 2])

% update the eog starting and ending points
eog_start_point = find([0 diff(eog_thrmerge)]==1);
eog_start = dwt_Ltrans(eog_start_point, eog_detect_l)/rs_Fs;
eog_end_point = find([0 diff(eog_thrmerge)]==-1);
eog_end = dwt_Ltrans(eog_end_point, eog_detect_l)/rs_Fs;

% label the eog regiion
figure, plot(rs_t, rs_eeg(1, :)), title('eog detection'), hold on
label_range= [-100 500]; % uv
for istart = eog_start
    plot([istart istart], label_range, 'r')
end
for iend= eog_end
    plot([iend iend], label_range, 'g')
end
grid minor, ylim([-200 600]), xlim([0 10]),xlabel('time(s)'), ylabel('uv');

% close all;
% %% channel 2 -- lots of differences with channel 1
% eog_detect_l = 4;
% [c_d, l_d] = wavedec(rs_eeg(2, :), eog_detect_l, wname);
% c_da = appcoef(c_d, l_d, wname);
% 
% thr_line = repmat(mean(c_da) - 1.5*std(c_da), 1, length(c_da));
% 
% t_a = linspace(0, length(c_da)*2^eog_detect_l/rs_Fs, length(c_da));
% figure, plot(t_a, c_da), title('a4'),
% hold on, plot(t_a, thr_line), xlim([0 100]);
% 
% eog_thr = c_da<thr_line;
% % 0.2 for channel 2
% eog_expand = [0.4 0.2];
% expand_points = floor(eog_expand ./ 2^eog_detect_l .* rs_Fs);
% 
% eog_thrmerge = eog_thr;
% % expand the eog starting point's location
% eog_start_point = find(diff(eog_thr)==1);
% start_expand = eog_start_point - expand_points(1);
% for i = 1:length(eog_start_point)
%     if start_expand(i)<0
%         start_expand(i) = 0;
%     end
%     eog_thrmerge(start_expand(i):eog_start_point(i)) =  true;
% end
% 
% figure, plot(eog_thr, 'Linewidth', 1.5),hold on, plot(eog_thrmerge)
% xlim([0 1000]), ylim([0 2])
% 
% % expand the eog ending point's location
% eog_end_point = find(diff(eog_thr)==-1);
% end_expand = eog_end_point + expand_points(2);
% for i = 1:length(eog_end_point)
%     if end_expand(i)>length(eog_thrmerge)
%         end_expand(i) = length(eog_thrmerge);
%     end
%     eog_thrmerge(eog_end_point(i):end_expand(i)) =  true;
% end
% 
% figure, plot(eog_thr, 'Linewidth', 1.5),hold on, plot(eog_thrmerge)
% xlim([0 1000]), ylim([0 2])
% 
% % update the eog starting and ending points
% eog_start_point = find(diff(eog_thrmerge)==1);
% eog_start = eog_start_point*2^eog_detect_l/rs_Fs;
% eog_end_point = find(diff(eog_thrmerge)==-1);
% eog_end = eog_end_point*2^eog_detect_l/rs_Fs;
% 
% % label the eog regiion
% figure, plot(rs_t, rs_eeg(2, :)), title('eog detection'), hold on
% label_range= [-10 5]; % uv
% for istart = eog_start
%     plot([istart istart], label_range, 'r')
% end
% for iend= eog_end
%     plot([iend iend], label_range, 'g')
% end
% % grid minor, ylim([-200 600]), xlim([0 10]),xlabel('time(s)'), ylabel('uv');
