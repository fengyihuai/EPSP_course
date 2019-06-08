function [eog_start_point eog_end_point] = eog_detection (t, eeg, Fs, num_channels, eog_expand)
% 


n = size(eeg, 2);
t = linspace(1/Fs, n/Fs, n);

for ichannel = 1:num_channels
    s = eeg(ichannel, :);
    
%     [b,a] = butter(6, 20/(Fs/2)); % 20Hz, 6-pole lowpass
%     freqs(b, a)
%     s = filter(b,a,s);
%     [b,a] = butter(6, 0.5/(Fs/2), 'high'); %0.5 Hz, 6-pole highpass
%     s = filter(b,a,s);    
    if abs(max(s))>abs(min(s))
        thr_line = repmat(mean(s) + 1.5*std(s), 1, length(s));  % for Fpz etc.- positive eog
        eog_thr = s>thr_line;
    else
        thr_line = repmat(mean(s) - 1.5*std(s), 1, length(s));  % for Pz etc.- negative eog
        eog_thr = s<thr_line;
    end
    
    figure, plot(t, s), title('a4'),
    hold on, plot(t, thr_line), xlim([0 100]);
    
%     eog_expand = [0.2 0.2];
    expand_points = floor(eog_expand .* Fs);
    % close all;
    
    figure, plot(t, s, t, max(s)*eog_thr), xlabel('time(s)'), xlim([190 200])
    
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
    figure, plot(t, s, t, eog_thrmerge), xlim([190 200]), % ylim([ 2])
    
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
    eog_start = eog_start_point{ichannel}/Fs;
    eog_end = eog_end_point{ichannel}/Fs;
    
    % label the eog regiion
    figure, plot(t, s), title('Eog detection'), hold on
    label_range= [0.5*min(s) 0.5*max(s)]; % uv
    for istart = eog_start % eog_start
        plot([istart istart], label_range, 'r')
    end
    for iend= eog_end % eog_end
        plot([iend iend], label_range, 'g')
    end
    grid minor, xlim([190 200]),xlabel('time(s)'), ylabel('uv');

    % ouput picture - ÅÐ±ð´°¿ÚÍØÕ¹
    figure, plot(t, s), xlabel('time(s)'), ylabel('uv'), xlim([0 5]), ylim([-50 150]), grid minor, title('Eog detection - window expanding'),
    hold on , plot(t, 0.8*max(s)*eog_thr, 'linewidth', 1.2), plot(t, 0.8*max(s)*eog_thrmerge, 'linewidth', 1.2),
    legend('raw eeg', 'Initial window', 'Expanded window');
    
end