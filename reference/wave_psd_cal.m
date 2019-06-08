% wave_psd_cal.m
% EEG wave power spectr calculation
% Energy calculation method:
% 

clc; clear; close all;

%% selected subjects and define the eeg wave band
subject_num = [1 2 4 5 6 7 9 10 11 14 15 16];
% subject_num = [1];
beseline_number = 1;

choice = {'APSD' ; 'RPSD'; 'RRPSD'};
% PSDType = choice{1}; % choose absolute PSD
PSDType = choice{2}; % choose relative PSD
% PSDType = choice{3}; % choose relative-relative PSD

thetaBand = [4  7];
alphaBand = [8 13];
betaBand = [13 30];
slowA = [8 11];
fastA = [11 13];

for isub = 1:length(subject_num)
    %% Define the data folder and read the preprocessed eeg data
    
%     [fileName, pathName] = uigetfile('*.mat', 'Select the preprocessed eeg data file');   
    fileName = ['corrected_eeg_' num2str(subject_num(isub)) '.mat'];
    
    % % fileName = 'corrected_eeg_4_chenkai.mat';
    pathName = 'D:\Program Files\MATLAB\WorkDir\StroopTest_EEG\EEG_Wave_PSD\';
    
    dataFile = [pathName fileName];
    fprintf([fileName ' is selected!!!\n']);
    load(dataFile);
    
    %% Wave Energy calculation for each trial
    
    w_length = 1000; % time length-2s 1000
    noverlap = 500; % 500
    nfft = 2000; % frequency resolution-0.5Hz while fs=500Hz 2000
    range='onesided';
    
    % frequency index serials
    cut_offId = (nfft/rs_Fs)*[1 40]+1; % for 1-40Hz cut off frequency band
    thetaBId = (nfft/rs_Fs)*thetaBand+1;
    alphaBId = (nfft/rs_Fs)*alphaBand;
    betaBId = (nfft/rs_Fs)*betaBand;
    
    window = blackman(w_length);
    
    figure(subject_num(isub));

    for ichannel = 1:num_channels
        for iblock = 1:nblocks
            
            [Pxx{isub, ichannel, iblock}, f] = pwelch(corrected_eeg{iblock, ichannel}, window, noverlap, nfft, rs_Fs, range);
            
        end
        
        % Absolute Power Spectra Density
        if strcmp(PSDType, 'APSD')
            
            max_tmp = 0;
            for iblock = 1:nblocks
                subplot(3, 2, ichannel);
                plot(f, Pxx{isub, ichannel, iblock}); xlabel('f/Hz'); ylabel('W/Hz');
                hold on;
                upper_limit = 1.5*max([max(Pxx{isub, ichannel, iblock}), max_tmp]);
                axis([0 40 0 upper_limit]);
                title(['\color{blue}Absolute ' 'Power Spectra Density ' '\color{red}' channel_type{ichannel}]);
                grid on;
                
                max_tmp = max(Pxx{isub, ichannel, iblock});
                
            end
            legend('BL', 'NCF1', 'NCF2', 'CF1', 'CF2');
            plot([thetaBand(1) thetaBand(1)], [0 max(Pxx{isub, ichannel, iblock})], 'k');
            plot([alphaBand(1) alphaBand(1)], [0 max(Pxx{isub, ichannel, iblock})], 'k');
            plot([betaBand(1) betaBand(1)], [0 max(Pxx{isub, ichannel, iblock})], 'k');
            plot([betaBand(2) betaBand(2)], [0 max(Pxx{isub, ichannel, iblock})], 'k');
            
        end
        
        % Relative Power Spectra Density
        if strcmp(PSDType, 'RPSD')
            
            baseline(ichannel) = sum(Pxx{isub, ichannel, beseline_number}([cut_offId(1):cut_offId(2)]'));
            
            for iblock = 1:nblocks
                
                relative_Pxx{isub, ichannel, iblock} = Pxx{isub, ichannel, iblock}./baseline(ichannel);
                
                subplot(3, 2, ichannel);
                plot(f, relative_Pxx{isub, ichannel, iblock}); xlabel('f/Hz'); ylabel('%');
                hold on;
                axis([0 40 0 0.2]);
                title(['\color{green}Realtive ' '\color{black}Power Spectra Density ' '\color{red}' channel_type{ichannel}]);
                grid on;
            end
            legend('BL', 'NCF1', 'NCF2', 'CF1', 'CF2');
            plot([thetaBand(1) thetaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            plot([alphaBand(1) alphaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            plot([betaBand(1) betaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            plot([betaBand(2) betaBand(2)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            
        end
        
          % Relative-relative Power Spectra Density
        if strcmp(PSDType, 'RRPSD')
            
            
            for iblock = 1:nblocks
                
                baseline(ichannel) = sum(Pxx{isub, ichannel, iblock}([cut_offId(1):cut_offId(2)]'));
                relative_Pxx{isub, ichannel, iblock} = Pxx{isub, ichannel, iblock}./baseline(ichannel);
                
                subplot(3, 2, ichannel);
                plot(f, relative_Pxx{isub, ichannel, iblock}); xlabel('f/Hz'); ylabel('%');
                hold on;
                axis([0 40 0 0.2]);
                title(['\color{green}Realtive ' '\color{black}Power Spectra Density ' '\color{red}' channel_type{ichannel}]);
                grid on;
            end
            legend('BL', 'NCF1', 'NCF2', 'CF1', 'CF2');
            plot([thetaBand(1) thetaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            plot([alphaBand(1) alphaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            plot([betaBand(1) betaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            plot([betaBand(2) betaBand(2)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
            
        end
        
    end

    %% Mean relative power calculation in theta(4-8Hz), alpha(8-16Hz) and 
    % beta(16-32Hz) band
    % Attention: resort the cell's dimensions for the subject number, block
    % number and channel number
    
    if strcmp(PSDType, 'RPSD')|strcmp(PSDType, 'RRPSD')
        for ichannel = 1:num_channels
            for iblock = 1:nblocks
                thetaBR_Pxx{isub, iblock, ichannel} = max(relative_Pxx{isub, ichannel, iblock}([thetaBId(1):thetaBId(2)]'));
                alphaBR_Pxx{isub, iblock, ichannel} =  max(relative_Pxx{isub, ichannel, iblock}([alphaBId(1):alphaBId(2)]'));
                betaBR_Pxx{isub, iblock, ichannel} =  max(relative_Pxx{isub, ichannel, iblock}([betaBId(1):betaBId(2)]'));
                                
            end
        end
    end
    
end

%% Export the mean relative power to excel
warning off;

if strcmp(PSDType, 'RPSD')
    for ichannel = 1:num_channels
        theta_sheetTmp = [thetaBR_Pxx{:, :, ichannel}];
        alpha_sheetTmp = [alphaBR_Pxx{:, :, ichannel}];
        beta_sheetTmp = [betaBR_Pxx{:, :, ichannel}];
        
        theta_sheet = [];
        alpha_sheet = [];
        beta_sheet = [];
                
        for iblock = 1:nblocks
            startId = (iblock-1)*length(subject_num)+1;
            endId = iblock*length(subject_num);
            theta_sheet = [theta_sheet; theta_sheetTmp(startId:endId)];
            alpha_sheet = [alpha_sheet; alpha_sheetTmp(startId:endId)];
            beta_sheet = [beta_sheet; beta_sheetTmp(startId:endId)];
        end
        
        sheet_name = channel_type{ichannel};
        
        xlswrite('thetaBandRPSD.xlsx', theta_sheet, sheet_name);        
        xlswrite('alphaBandRPSD.xlsx', alpha_sheet, sheet_name);
        xlswrite('betaBandRPSD.xlsx', beta_sheet, sheet_name);
        
    end
    fprintf('\nthetaWave.xlsx is saved!!!\n');
    fprintf('alphaWave.xlsx is saved!!!\n');
    fprintf('betaWave.xlsx is saved!!!\n');
    
else
    fprintf('\nNo xlsx data file output\n');
end


%% Single RPSD Plot

isub = 12;
ichannel = 1;

figure;
for iblock =1:nblocks

    plot(f, relative_Pxx{isub, ichannel, iblock}, 'linewidth', 2); 
    xlabel('f/Hz', 'FontWeight','bold'); ylabel('%', 'FontWeight','bold');
    hold on;
    axis([0 40 0 0.2]);
    title(['\color{green}Realtive ' '\color{black}Power Spectra Density ' '\color{red}' channel_type{ichannel}], ...
        'FontWeight','bold');
    grid on;
    
    
    legend('BL', 'NCF1', 'NCF2', 'CF1', 'CF2');
end
    plot([thetaBand(1) thetaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k' );
    plot([alphaBand(1) alphaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
    plot([betaBand(1) betaBand(1)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');
    plot([betaBand(2) betaBand(2)], [0 1.5*max(relative_Pxx{isub, ichannel, iblock})], 'k');

%% Plot 3 different sections
    
num_section = 3;

for isub = 1:length(subject_num)
    
    figure(1000+subject_num(isub));
    for ichannel = 1:num_channels
        
        RRPxx_section{isub, ichannel, 1} = relative_Pxx{isub, ichannel, 1};
        RRPxx_section{isub, ichannel, 2} = (relative_Pxx{isub, ichannel, 2}+relative_Pxx{isub, ichannel, 3})/2;
        RRPxx_section{isub, ichannel, 3} = (relative_Pxx{isub, ichannel, 4}+relative_Pxx{isub, ichannel, 5})/2;
        
        for isec = 1:num_section
            
            subplot(3, 2, ichannel);
            plot(f, RRPxx_section{isub, ichannel, isec}); xlabel('f/Hz'); ylabel('%');
            hold on;
            axis([0 40 0 0.2]);
            title(['\color{green}Realtive ' '\color{black}Power Spectra Density ' '\color{red}' channel_type{ichannel}]);
            grid on;
        end
        legend('BL', 'NCF', 'CF');
        plot([thetaBand(1) thetaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        plot([alphaBand(1) alphaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        plot([betaBand(1) betaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        plot([betaBand(2) betaBand(2)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        
    end
    
     fprintf([num2str(subject_num(isub)) ' subject''s RRPxx_section is calculated.\n']);
end

%% Single RPSD Plot

isub = 12; % Zouwenbin
ichannel = 3; % Fp2

figure;
for isec =1:num_section

    plot(f, RRPxx_section{isub, ichannel, isec}, 'linewidth', 2); 
    xlabel('f/Hz', 'FontWeight','bold'); ylabel('%', 'FontWeight','bold');
    hold on;
    axis([0 40 0 0.2]);
    title(['\color{green}Realtive ' '\color{black}Power Spectra Density ' '\color{red}' channel_type{ichannel}], ...
        'FontWeight','bold');
    grid on;
    
    
    legend('BL', 'NCF1', 'NCF2', 'CF1', 'CF2');
end
    plot([thetaBand(1) thetaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k' );
    plot([alphaBand(1) alphaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
    plot([betaBand(1) betaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
    plot([betaBand(2) betaBand(2)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');

%% Mean RPSD Plot of all subjects

for isub = 1:length(subject_num)
    
    figure(1000+subject_num(isub));
    for ichannel = 1:num_channels
        
        RRPxx_section{isub, ichannel, 1} = relative_Pxx{isub, ichannel, 1};
        RRPxx_section{isub, ichannel, 2} = (relative_Pxx{isub, ichannel, 2}+relative_Pxx{isub, ichannel, 3})/2;
        RRPxx_section{isub, ichannel, 3} = (relative_Pxx{isub, ichannel, 4}+relative_Pxx{isub, ichannel, 5})/2;
        
        for isec = 1:num_section
            
            subplot(3, 2, ichannel);
            plot(f, RRPxx_section{isub, ichannel, isec}); xlabel('f/Hz'); ylabel('%');
            hold on;
            axis([0 40 0 0.2]);
            title(['\color{green}Realtive ' '\color{black}Power Spectra Density ' '\color{red}' channel_type{ichannel}]);
            grid on;
        end
        legend('BL', 'NCF', 'CF');
        plot([thetaBand(1) thetaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        plot([alphaBand(1) alphaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        plot([betaBand(1) betaBand(1)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        plot([betaBand(2) betaBand(2)], [0 1.5*max(RRPxx_section{isub, ichannel, isec})], 'k');
        
    end
    
     fprintf([num2str(subject_num(isub)) ' subject''s RRPxx_section is calculated.\n']);
end

