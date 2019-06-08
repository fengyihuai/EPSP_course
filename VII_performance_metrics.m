% (7)使用今天讨论的算法性能评价指标中的2种，评价第（6）项信号处理算法的性能
% The VIxxx.m must have been runned

%% Simulated signals contaminated normal white noise
if exist('sinnnoise')
    % Normalized mean squared error 
    NMSE_sine = nmserr(swave, denoise_sine)
    NMSE_mw = nmserr(mwavelet, denoise_mw)
    NMSE_bw = nmserr(boxwave, denoise_bw)   
    figure, bar([NMSE_sine, NMSE_mw, NMSE_bw]), set(gca,'xticklabel',{'sine', 'morletW', 'boxwave'}),ylabel('Normalized mean squared error', 'FontSize',14)
    title('Normalized mean squared error -signals + normal white noise', 'FontSize',14), set(gca,'FontSize',20);
    
    % Mutual information
    MI_sine = minfo(swave, denoise_sine)
    MI_mw = minfo(mwavelet, denoise_mw)
    MI_bw = minfo(boxwave, denoise_bw)
    figure, bar([MI_sine, MI_mw, MI_bw]), set(gca,'xticklabel',{'sine'; 'morletW'; 'boxwave'}),ylabel('Mutual Information', 'FontSize',14)
    title('Mutual information -signals + normal white noise', 'FontSize',14), set(gca,'FontSize',20);
end

%% Simulated signals contaminated uniform white noise
if exist('sinunoise')
    % Normalized mean squared error 
    NMSE_sine = nmserr(swave, denoise_sine)
    NMSE_mw = nmserr(mwavelet, denoise_mw)
    NMSE_bw = nmserr(boxwave, denoise_bw)
    figure, bar([NMSE_sine, NMSE_mw, NMSE_bw]), set(gca,'xticklabel',{'sine', 'morletW', 'boxwave'}),ylabel('Normalized mean squared error', 'FontSize',14),
    title('Normalized mean squared error -signals + uniform white noise', 'FontSize',14), set(gca,'FontSize',20);
    
    % Mutual information
    MI_sine = minfo(swave, denoise_sine)
    MI_mw = minfo(mwavelet, denoise_mw)
    MI_bw = minfo(boxwave, denoise_bw)
    figure, bar([MI_sine, MI_mw, MI_bw]), set(gca,'xticklabel',{'sine'; 'morletW'; 'boxwave'}),ylabel('Mutual Information', 'FontSize',14),
    title('Mutual information -signals + uniform white noise', 'FontSize',14), set(gca,'FontSize',20);
end
