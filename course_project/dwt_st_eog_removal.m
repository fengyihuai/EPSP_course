% 1.	脑电信号中的眼电去除。
% 测试
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
hold on, plot(linspace(-Fs/2, Fs/2, length(raw_eeg(2, :))), fftshift(2* abs(fft(raw_eeg(2, :)))/n)), 
xlabel('frequency(Hz)'), ylabel('Amp'), xlim(show_ft), grid minor,  title('Raw EEG spectrum')

%% dwt-st-eog threshold denoising
rs_Fs = 500;
wname = 'sym5'; % the first time: sym3
nlevel = 6;
denoised_levels = [5 6];
thr_type = 'h';

rs_t = linspace(1/Fs, n/Fs, floor(n*rs_Fs/Fs));

denoised_eeg = dwtStEogDenoise(rs_eeg, Fs, wname, nlevel, denoised_levels, thr_type);

%% raw eeg vs. denoised eeg
figure, subplot(211), plot(t, raw_eeg(1, :), rs_t, denoised_eeg(1, :)), xlim([0 10]), ylabel('uv'), title('raw EEG vs. denoised eeg'),  grid minor, legend('', '');
subplot(212), plot(t, raw_eeg(2, :), rs_t, denoised_eeg(2, :)), xlim([0 10]), ylabel('uv'), title('raw EEG vs. denoised eeg'),  grid minor;

%%  eye-region detection, channel 1-2, Fpz, Pz

for ichannel = 1:num_channels
    s = rs_eeg(ichannel, :);
    
%     [b,a] = butter(6, 20/(rs_Fs/2)); % 20Hz, 6-pole lowpass
%     freqs(b, a)
%     s = filter(b,a,s);
%     [b,a] = butter(6, 0.5/(rs_Fs/2), 'high'); %0.5 Hz, 6-pole highpass
%     s = filter(b,a,s);    
    if abs(max(s))>abs(min(s))
        thr_line = repmat(mean(s) + 1.5*std(s), 1, length(s));  % for Fpz
        eog_thr = s>thr_line;
    else
        thr_line = repmat(mean(s) - 1.5*std(s), 1, length(s));  % for Pz
        eog_thr = s<thr_line;
    end
    
    figure, plot(rs_t, s), title('a4'),
    hold on, plot(rs_t, thr_line), xlim([0 100]);
    
    eog_expand = [0.2 0.2];
    expand_points = floor(eog_expand .* rs_Fs);
    % close all;
    
    figure, plot(rs_t, s, rs_t, max(s)*eog_thr), xlabel('time(s)'), xlim([190 200])
    
    eog_thrmerge = eog_thr;
    % expand the eog starting point's location
    peak_start_point = find([0 diff(eog_thr)]==1);
    peak_end_point = find([0 diff(eog_thr)]==-1);
%     % Prevention of mis-detection, just detect the peak
%     for element = length(peak_end_point)
%         del_label = find(peak_start_point==element);
%         peak_start_point(del_label) = [];
%         peak_end_point(find(peak_end_point==element)) = [];
%     end
    
    start_expand = peak_start_point - expand_points(1);
    for i = 1:length(peak_start_point)
        if start_expand(i)<0
            start_expand(i) = 1;
        end
        eog_thrmerge(start_expand(i):peak_start_point(i)) =  true;
    end
    
    figure, plot(eog_thr, 'Linewidth', 1.5),hold on, plot(eog_thrmerge)
    xlim([0 1000]), ylim([0 2])
    figure, plot(rs_t, s, rs_t, eog_thrmerge), xlim([190 200]), % ylim([ 2])
    
    % expand the eog ending point's location

    end_expand = peak_end_point + expand_points(2);
    for i = 1:length(peak_end_point)
        if end_expand(i)>length(eog_thrmerge)
            end_expand(i) = length(eog_thrmerge)-1;
        end
        eog_thrmerge(peak_end_point(i):end_expand(i)) =  true;
    end
    
    figure, plot(eog_thr, 'Linewidth', 1.5), hold on, plot(eog_thrmerge)
    xlim([0 1000]), ylim([0 2]);
    
    % update the eog starting and ending points
    eog_start_point{ichannel} = find([0 diff(eog_thrmerge)]==1);
    eog_end_point{ichannel} = find([0 diff(eog_thrmerge)]==-1);
    
    % time label calculation
    eog_start = eog_start_point{ichannel}/rs_Fs;
    eog_end = eog_end_point{ichannel}/rs_Fs;
    
    % label the eog regiion
    figure, plot(rs_t, s), title('eog detection'), hold on
    label_range= [0.5*min(s) 0.5*max(s)]; % uv
    for istart = eog_start % eog_start
        plot([istart istart], label_range, 'r')
    end
    for iend= eog_end % eog_end
        plot([iend iend], label_range, 'g')
    end
    grid minor, xlim([190 200]),xlabel('time(s)'), ylabel('uv');

end

%% Splice of eeg non-eog region - artifact-free segment
% close all;
for ichannel = 1:num_channels
    raw_neog_splice{ichannel} = rs_eeg(ichannel, :);
    for j = length(eog_start_point{ichannel}):-1:1
        indices = eog_start_point{ichannel}(j):eog_end_point{ichannel}(j);
%         raw_eog_section = [raw_eog_section; raw_neog_splice{ichannel}(indices)];
%         raw_neog_splice = [raw_neog_splice; raw_neog_splice{ichannel}(indices)];
        raw_neog_splice{ichannel}(indices) = [];
    end
end

for ichannel = 1:num_channels
    denoised_neog_splice{ichannel} = denoised_eeg(ichannel, :);
    for j = length(eog_start_point{ichannel}):-1:1
        indices = eog_start_point{ichannel}(j):eog_end_point{ichannel}(j);
        denoised_neog_splice{ichannel}(indices) = [];
    end
end

%% NMSE and IM calculation of raw eeg and denoised eeg in the non-eog region
for ichannel = 1:num_channels
    % Normalized mean squared error
    NMSE(ichannel) = nmserr(raw_neog_splice{ichannel}, denoised_neog_splice{ichannel});
    %     figure, bar([NMSE_sine, NMSE_mw, NMSE_bw]), set(gca,'xticklabel',{'sine', 'morletW', 'boxwave'}),ylabel('Normalized mean squared error', 'FontSize',14)
    %     title('Normalized mean squared error -signals + normal white noise', 'FontSize',14), set(gca,'FontSize',20);
    
    % Mutual information
    MI(ichannel) = minfo(raw_neog_splice{ichannel}, denoised_neog_splice{ichannel});
    % MI_sine = minfo(rs_eeg(ichannel, :), denoised_eeg(ichannel, :));
    %     figure, bar([MI_sine, MI_mw, MI_bw]), set(gca,'xticklabel',{'sine'; 'morletW'; 'boxwave'}),ylabel('Mutual Information', 'FontSize',14)
    %     title('Mutual information -signals + normal white noise', 'FontSize',14), set(gca,'FontSize',20);
end
