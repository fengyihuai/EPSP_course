%% eeg_preprocess.m
% 
% Artifact removal method:
% Mallat wavelet transform - WTCoefficients threshold setting method
% 


clc; clear all; close all;

Fs = 2000; % Sampling rate

num_channelsArray = 1:6; % selected channels for strooptest 
num_channels = 6;
nblocks =5;

%% Define the data folder
[fileName, pathName] = uigetfile('*.cnt', 'Select the NeuroScan data file');
% fileName = '1_zouwenbin.cnt';
% fileName = subjectName{iname};
% pathName = 'D:\Program Files\MATLAB\WorkDir\StroopTest_EEG\EEG_Wave_PSD_w\';
dataFile = [pathName fileName];
fprintf([fileName ' is selected!!!\n']);

%% Read the raw data and extract parameters
cnt = loadcnt(dataFile);
event = load_event(dataFile); % event
raw_eeg = double(cnt.data(num_channelsArray, :));

% Segregate the raw eeg data into experiment blocks
block_type = {'Arrows'; 'noConflict1'; 'noConflict2';...
    'Conflict1'; 'Conflict2'}';
channel_type ={};
for ilab=1:num_channels
    channel_type = [channel_type; cnt.electloc(1, ilab).lab];
end
raw_eeg_seg{1, 1} = raw_eeg(:, event(1,2):event(40,2));

for iblock= 1:4
    raw_eeg_seg{iblock+1, 1} = raw_eeg(:, event((iblock-1)*100+41, 2):event(iblock*100+40, 2));
end


%% Downsampling - 2000Hz is converted to 500Hz
% Attention:resample applies an antialiasing FIR lowpass filter to x and compensates for 
% the delay introduced by the filter.

rs_Fs = 500;

for iblock =1:nblocks
    for ichannel = 1:num_channels
        rs_eeg{iblock, 1}(ichannel, :) = resample(raw_eeg_seg{iblock, 1}(ichannel, :), rs_Fs, Fs);
    end
end

rs_event = floor(event(:, 2)/(Fs/500));

%% Mallat wavelet transform - WTCoefficients threshold setting method
wavename = 'sym5'; % the first time: sym3
nlevel = 6;

for iblock =1:nblocks
    for ichannel = 1:num_channels
        %         iblock = 1; % Arrows
        %         ichannel = 1; %Fp1
        %         iblock = 2; % Arrows
        %         ichannel = 2; %Fp1
        % Mallat wavelet transform - decomposition and reconstruction of each frequency band
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % based on Power 2
        
        s = rs_eeg{iblock, 1}(ichannel, :);
        l_s = length(s);
        % wavelet decomposition
        [C{iblock, ichannel}, L{iblock, ichannel}] = wavedec(s, nlevel, wavename);
        cA{iblock, ichannel} = appcoef(C{iblock, ichannel}, L{iblock, ichannel}, wavename, nlevel);
        for ilevel = 1:nlevel
            cD{iblock, ichannel, ilevel} = detcoef(C{iblock, ichannel},L{iblock, ichannel}, ilevel);
        end
        % wavelet reconstruction
        A{iblock, ichannel} = wrcoef('a', C{iblock, ichannel}, L{iblock, ichannel}, wavename, nlevel);
        for ilevel = 1:nlevel
            D{iblock, ichannel, ilevel} = wrcoef('d', C{iblock, ichannel}, L{iblock, ichannel}, wavename, ilevel);
        end
        
        
        % EOG and EMG removal - Wavelet coefficient threshold setting
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % remove eog (0-16Hz) and emg(16-65Hz:principally)
        % using the module maximun-mean_std_threshold
        % Attention:eyeblink arifact may exist in 4-65Hz frequency band
        
        % !!!only Fp1, Fp2, Fpz need to be remove eog and emg
        % !!!Fz, Fcz, Pz doesn't need
        if ichannel<4;
            m = [0 0 2.0 2.0 1.5 1.5]; % recommed:0.5-2
            n = [0 0 0.1 0.1 0.1 0.1]; % recommed:0.01-0.2
        else (ichannel>3)&(ichannel<7);
            m = [0 0 0 0 0 0];
            n = [0 0 0 0 0 0];
        end
        
        for ireset = 1:nlevel;
            th_wave = cD{iblock, ichannel, ireset};
            wpeak = waveMmax(th_wave, length(th_wave));
            threshold{iblock, ichannel, ireset} = mean(wpeak)+m(ireset)*std(wpeak);
            
            if (ireset>2);
                % reset detail 3-6 coefficients
                for isn = 1:length(cD{iblock, ichannel, ireset})
                    if (abs(cD{iblock, ichannel, ireset}(1, isn))> threshold{iblock, ichannel, ireset})&(m(ireset)*n(ireset)~=0)
                        reset_cD{iblock, ichannel, ireset}(1, isn) = n(ireset)*cD{iblock, ichannel, ireset}(1, isn); %reset approximation wavelet co
                    else
                        reset_cD{iblock, ichannel, ireset}(1, isn) = cD{iblock, ichannel, ireset}(1, isn);
                    end
                end
            else 
                %set detail 1-2 coefficients to 
                reset_cD{iblock, ichannel, ireset} = zeros(1, length(cD{iblock, ichannel, ireset}));
            end
            
        end
        
        % set ca6 to 0 and update the changes to reset_C
        reset_cA{iblock, ichannel} = zeros(1, length(cA{iblock, ichannel}));
        reset_Carray = reset_cA{iblock, ichannel};
        for jlevel = nlevel:-1:1
            cD_j = reset_cD{iblock, ichannel, jlevel};
            reset_Carray = [reset_Carray cD_j];
        end
        reset_C{iblock, ichannel} = reset_Carray;
        
        % Reconstructed wavelet coefficients which are reset
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % wavelet reconstruction
        reset_A{iblock, ichannel} = wrcoef('a', reset_C{iblock, ichannel}, L{iblock, ichannel}, wavename, nlevel);
        for ilevel = 1:nlevel
            reset_D{iblock, ichannel, ilevel} = wrcoef('d', reset_C{iblock, ichannel}, L{iblock, ichannel}, wavename, ilevel);
        end
        
        % 40Hz low-pass filter just for d3
        reset_D{iblock, ichannel, 3} = lowpassfilt(reset_D{iblock, ichannel, 3}, rs_Fs, 40);
        
        % Reconstruct the de-noise/corrected eeg
        corrected_eeg{iblock, ichannel} = reset_D{iblock, ichannel, 3} + reset_D{iblock, ichannel, 4} +...
            reset_D{iblock, ichannel, 5} + reset_D{iblock, ichannel, 6};
        
    end
end

%% Neural osicillation of EEG
for iblock = 1:nblocks
    for ichannel = 1:num_channels
        osc{iblock, ichannel}.theta = reset_D{iblock, ichannel, 6};
        osc{iblock, ichannel}.alpha = reset_D{iblock, ichannel, 5};
        osc{iblock, ichannel}.beta = reset_D{iblock, ichannel, 4};
    end
end

%% Data saving for further analysis

% extract subject's number
str_serial = isstrprop(fileName, 'digit');
num_sub = fileName(str_serial);

save(['corrected_eeg_' num_sub '.mat'], 'corrected_eeg', 'osc', 'block_type', 'channel_type', ...
    'rs_Fs', 'rs_event', 'nblocks', 'num_channels');
fprintf(['corrected_eeg_' num_sub '.mat' 'is saved!!!\n']);

%% Display the preprocess-eeg and plot the curve
 
iblock = 2; % Block1
ichannel = 2; %Fp1
% ichannel =1; %Fpz
% iblock = 2; % Block1
% ichannel = 4; %Fz
% plot the WT-Curve
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
s = rs_eeg{iblock, 1}(ichannel, :);
raw_s = raw_eeg_seg{iblock, 1}(ichannel, :);

l_s = length(s);
figure;
subplot(nlevel/2, 1, 1);
plot((1:l_s)/rs_Fs, s, 'b'); xlabel('t/s'); ylabel('Amplitude/uv');
title([block_type{iblock} ' \color{red}' cnt.electloc(1,ichannel).lab ' \color{black}channel']);
for iplot = 1:nlevel/2-1;
    subplot(nlevel/2, 1, iplot+1);
    plot((1:l_s)/rs_Fs, D{iblock, ichannel, iplot}, 'b');
    hold on;plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, iplot}, 'Color', 'r', 'Linewidth', 1.5); 
%     hold on;plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, iplot}, 'r', 'linespec', '.-'); 
    xlabel('t/s'); ylabel(['d_' num2str(iplot) '/uv']);
    title(['Detail-frequency band: ' num2str(rs_Fs/(2^(iplot+1))) 'Hz -' num2str(rs_Fs/(2^iplot)) 'Hz']);
    legend('original wave', 'corrected wave');
end
figure;
for iplot = 1+nlevel/2-1:nlevel-1
    subplot(nlevel/2, 1, iplot-(nlevel/2-1));
    plot((1:l_s)/rs_Fs, D{iblock, ichannel, iplot}, 'b'); 
    hold on;plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, iplot}, 'Color', 'r', 'Linewidth', 1.5); 
    xlabel('t/s'); ylabel(['d_' num2str(iplot) '/uv']);
    title(['Detail-frequency band: ' num2str(rs_Fs/(2^(iplot+1))) 'Hz -' num2str(rs_Fs/(2^iplot)) 'Hz']);
    legend('original wave', 'corrected wave');
end
figure;
subplot(nlevel/2, 1, 1);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, nlevel}, 'b'); 
hold on; plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, nlevel}, 'Color', 'r', 'Linewidth', 1.5); 
xlabel('t/s'); ylabel(['d_' num2str(nlevel) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(nlevel+1))) 'Hz -' num2str(rs_Fs/(2^nlevel)) 'Hz']);
legend('original wave', 'corrected wave');
subplot(nlevel/2, 1, 2);
plot((1:l_s)/rs_Fs, A{iblock, ichannel}, 'b'); 
hold on; plot((1:l_s)/rs_Fs, reset_A{iblock, ichannel}, 'Color', 'r', 'Linewidth', 1.5); 
xlabel('t/s'); ylabel(['A_' num2str(nlevel) '/uv']);
title(['Approximation-' 'frequency band:' '0.05Hz -' num2str(rs_Fs/(2^(nlevel+1))) 'Hz']);
legend('original wave', 'corrected wave');


% plot the corrected eeg
figure;
plot((1:l_s)/rs_Fs, s, 'b'); 
hold on;
plot((1:l_s)/rs_Fs, corrected_eeg{iblock, ichannel}, 'Color', 'r', 'Linewidth', 1.5);
xlabel('t/s'); ylabel('Amplitude/uv');
title(['Raw EEG with noise' ' \color{red}VS' ' \color{black}Corrected EEG']);
legend('Raw EEG with noise', 'Corrected EEG');

%% Single curve plot - baseline showing

figure;
plot((1:l_s)/rs_Fs, s, 'b'); 
hold on;
xlabel('t/s'); ylabel('Amplitude/uv');
title('Raw EEG with noise');

timeWindow = 0.10;       %second
windowLeng = timeWindow*rs_Fs+1;

s_wander = medfilt1(s, windowLeng);
% s_lp = lowpassfilt(s, rs_Fs, 1);

plot((1:l_s)/rs_Fs, s_wander, 'r', 'linewidth', 2); 
legend('original wave', 'baseline');

%% Single curve plot - downsample showing

figure;

subplot(2, 1, 1)
plot((1:length(raw_s))/Fs, raw_s, 'b'); 
hold on;
xlabel('t/s'); ylabel('Amplitude/uv');
title('Raw EEG');

subplot(2, 1, 2)
plot((1:l_s)/rs_Fs, s, 'b'); 
hold on;
xlabel('t/s'); ylabel('Amplitude/uv');
title('Downsampled EEG');

%% Single curve plot - 50Hz interference

% mscohere function

figure;
plot((1:l_s)/rs_Fs, s, 'b'); 
hold on;
xlabel('t/s'); ylabel('Amplitude/uv');
title('Raw EEG with noise');

figure;

N_s = length(s);
s_fft = fft(s, N_s); % Fast fourier transform
s_fftMag = abs(s_fft);
f = (0:length(s_fft)-1)'*rs_Fs/length(s_fft);
plot(f(1:N_s/2),s_fftMag(1:N_s/2));  % Plot the FFT
xlabel('f/Hz');ylabel('Magnitude');
title('Spectrogram');

% s_lp = lowpassfilt(s, rs_Fs, 1);
% 
% plot((1:l_s)/rs_Fs, s_lp, 'r', 'linewidth', 2); 

%% pwelch psd
w_length = 1000; % time length-2s 1000
noverlap = 500; % 500
nfft = 2000; % frequency resolution-0.5Hz while fs=500Hz 2000
range='onesided';

 window = blackman(w_length);

 [Pxx_s f_s] = pwelch(s, window, noverlap, nfft, rs_Fs, range);
 [Pxx_corrected f_corrected]= pwelch(corrected_eeg{iblock, ichannel}, window, noverlap, nfft, rs_Fs, range);

 figure;
 plot(f_s, Pxx_s);
 hold on;
 plot(f_corrected, Pxx_corrected); 
 xlabel('f/Hz'); ylabel('W/Hz');
 upper_limit = 1.5*max(Pxx_s);
 axis([0 40 0 upper_limit]);
 title(['\color{blue}Absolute ' 'Power Spectra Density ' '\color{red}' channel_type{ichannel}]);
 grid on;

%%  Magnitude squared coherence measure plot for one subject


[Cxy, F] = mscohere(s, corrected_eeg{iblock, ichannel}, window, noverlap, nfft, rs_Fs);

figure;
stem(F, Cxy);
xlabel('Frequency/Hz'); ylabel('Coherence');
title('Magnitude squared coherence for raw and corrected eeg')

% xlabel('t/s'); ylabel('Amplitude/uv');

% title('Downsampled EEG');

%% Display the preprocessing-eeg and plot the curve
 
iblock = 2; % Block1
ichannel = 2; %Fp1
% ichannel =1; %Fpz
% iblock = 2; % Block1
ichannel = 4; %Fz
% plot the WT-Curve
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
s = rs_eeg{iblock, 1}(ichannel, :);
raw_s = raw_eeg_seg{iblock, 1}(ichannel, :);

l_s = length(s);
figure;

subplot(4, 1, 1);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 3}, 'b');
hold on;plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, 3}, 'Color', 'r', 'Linewidth', 1.5);
xlabel('t/s'); ylabel(['d_' num2str(3) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(3+1))) 'Hz -' num2str(rs_Fs/(2^3)) 'Hz']);
legend('original wave', 'corrected wave');

subplot(4, 1, 2);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 4}, 'b');
hold on;plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, 4}, 'Color', 'r', 'Linewidth', 1.5);
xlabel('t/s'); ylabel(['d_' num2str(4) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(4+1))) 'Hz -' num2str(rs_Fs/(2^4)) 'Hz']);
legend('original wave', 'corrected wave');

subplot(4, 1, 3);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 5}, 'b');
hold on;plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, 5}, 'Color', 'r', 'Linewidth', 1.5);
xlabel('t/s'); ylabel(['d_' num2str(5) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(5+1))) 'Hz -' num2str(rs_Fs/(2^5)) 'Hz']);
legend('original wave', 'corrected wave');

subplot(4, 1, 4);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 6}, 'b');
hold on;plot((1:l_s)/rs_Fs, reset_D{iblock, ichannel, 6}, 'Color', 'r', 'Linewidth', 1.5);
xlabel('t/s'); ylabel(['d_' num2str(6) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(6+1))) 'Hz -' num2str(rs_Fs/(2^6)) 'Hz']);
legend('original wave', 'corrected wave');

%% Display the Mallat-Reconstruction wave
 
iblock = 2; % Block1
ichannel = 2; %Fp1
% ichannel =1; %Fpz
% iblock = 2; % Block1
% ichannel = 4; %Fz
% plot the WT-Curve
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
s = rs_eeg{iblock, 1}(ichannel, :);
raw_s = raw_eeg_seg{iblock, 1}(ichannel, :);

l_s = length(s);
figure;

subplot(4, 1, 1);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 1}, 'b');
xlabel('t/s'); ylabel(['d_' num2str(1) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(1+1))) 'Hz -' num2str(rs_Fs/(2^1)) 'Hz']);

subplot(4, 1, 2);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 2}, 'b');
xlabel('t/s'); ylabel(['d_' num2str(2) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(2+1))) 'Hz -' num2str(rs_Fs/(2^2)) 'Hz']);

subplot(4, 1, 3);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 3}, 'b');
xlabel('t/s'); ylabel(['d_' num2str(3) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(3+1))) 'Hz -' num2str(rs_Fs/(2^3)) 'Hz']);

subplot(4, 1, 4);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 4}, 'b');
xlabel('t/s'); ylabel(['d_' num2str(4) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(4+1))) 'Hz -' num2str(rs_Fs/(2^4)) 'Hz']);

figure;
subplot(4, 1, 1);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 5}, 'b');
xlabel('t/s'); ylabel(['d_' num2str(5) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(5+1))) 'Hz -' num2str(rs_Fs/(2^5)) 'Hz']);

subplot(4, 1, 2);
plot((1:l_s)/rs_Fs, D{iblock, ichannel, 6}, 'b');
xlabel('t/s'); ylabel(['d_' num2str(6) '/uv']);
title(['Detail-frequency band: ' num2str(rs_Fs/(2^(6+1))) 'Hz -' num2str(rs_Fs/(2^6)) 'Hz']);

subplot(4, 1, 3);
plot((1:l_s)/rs_Fs, A{iblock, ichannel}, 'b');
xlabel('t/s'); ylabel(['A_' num2str(6) '/uv']);
title(['Approximation-' 'frequency band:' '0.05Hz -' num2str(rs_Fs/(2^(6+1))) 'Hz']);

%% Display the Mallat-Decomposition wave
 
% iblock = 2; % Block1
% ichannel = 2; %Fp1
% % ichannel =1; %Fpz
% % iblock = 2; % Block1
% % ichannel = 4; %Fz
% % plot the WT-Curve
% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% s = rs_eeg{iblock, 1}(ichannel, :);
% raw_s = raw_eeg_seg{iblock, 1}(ichannel, :);
% 
% l_s = length(s);
% figure;
% 
% subplot(4, 1, 1);
% % plot((1:l_s)/rs_Fs, cD{iblock, ichannel, 1}, 'b');
% plot(cD{iblock, ichannel, 1}, 'b');
% xlabel('t/s'); ylabel(['d_' num2str(1) '/uv']);
% title(['Detail-frequency band: ' num2str(rs_Fs/(2^(1+1))) 'Hz -' num2str(rs_Fs/(2^1)) 'Hz']);
% 
% subplot(4, 1, 2);
% plot((1:l_s)/rs_Fs, cD{iblock, ichannel, 2}, 'b');
% xlabel('t/s'); ylabel(['d_' num2str(2) '/uv']);
% title(['Detail-frequency band: ' num2str(rs_Fs/(2^(2+1))) 'Hz -' num2str(rs_Fs/(2^2)) 'Hz']);
% 
% subplot(4, 1, 3);
% plot((1:l_s)/rs_Fs, cD{iblock, ichannel, 3}, 'b');
% xlabel('t/s'); ylabel(['d_' num2str(3) '/uv']);
% title(['Detail-frequency band: ' num2str(rs_Fs/(2^(3+1))) 'Hz -' num2str(rs_Fs/(2^3)) 'Hz']);
% 
% subplot(4, 1, 4);
% plot((1:l_s)/rs_Fs, cD{iblock, ichannel, 4}, 'b');
% xlabel('t/s'); ylabel(['d_' num2str(4) '/uv']);
% title(['Detail-frequency band: ' num2str(rs_Fs/(2^(4+1))) 'Hz -' num2str(rs_Fs/(2^4)) 'Hz']);
% 
% figure;
% subplot(4, 1, 1);
% plot((1:l_s)/rs_Fs, cD{iblock, ichannel, 5}, 'b');
% xlabel('t/s'); ylabel(['d_' num2str(5) '/uv']);
% title(['Detail-frequency band: ' num2str(rs_Fs/(2^(5+1))) 'Hz -' num2str(rs_Fs/(2^5)) 'Hz']);
% 
% subplot(4, 1, 2);
% plot((1:l_s)/rs_Fs, cD{iblock, ichannel, 6}, 'b');
% xlabel('t/s'); ylabel(['d_' num2str(6) '/uv']);
% title(['Detail-frequency band: ' num2str(rs_Fs/(2^(6+1))) 'Hz -' num2str(rs_Fs/(2^6)) 'Hz']);
% 
% subplot(4, 1, 3);
% plot((1:l_s)/rs_Fs, cA{iblock, ichannel}, 'b');
% xlabel('t/s'); ylabel(['A_' num2str(6) '/uv']);
% title(['Approximation-' 'frequency band:' '0.05Hz -' num2str(rs_Fs/(2^(6+1))) 'Hz']);
