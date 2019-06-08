function denoised_eeg = dwtUtEogDenoise(eeg, Fs, wname, nlevel, denoised_levels, thr_type)

n = size(eeg, 2);
num_channels = size(eeg, 1);
t = linspace(1/Fs, n/Fs, n);

%% EOG artifact removal based on wavelet tranform coefficients thresholding

% % wavelet function plot
% [phi,psi,xval] = wavefun(wname);
% figure, subplot(211), plot(xval,phi), title(strcat(wname, ' Scaling Function')), grid on;
% subplot(212), plot(xval,psi), title(strcat(wname, ' Wavelet')), grid on;

for ichannel = 1:num_channels
    % Mallat wavelet transform - decomposition and reconstruction of each frequency band
    
    s = eeg(ichannel, :);
    % wavelet decomposition
    [c{ichannel}, l{ichannel}] = wavedec(s, nlevel, wname);
    
    % EOG removal - universal thresholding calculation
    % remove eog (0-16Hz)
    denoise_c{ichannel} = c{ichannel};
    for ilevel = denoised_levels
        hk = detcoef(c{ichannel}, l{ichannel}, ilevel);
        sigma = wnoisest(c{ichannel}, l{ichannel}, ilevel);
        thr = sigma*sqrt(2*log(n)); % threshold  
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
    figure, subplot(l_l+2, 1, 1), plot(t, eeg(ichannel, :)), ylabel('original'),title(['UT-', thr_type, '-Original & denoised DWT coefficients - channel', num2str(ichannel)]);
    for ilevel = denoised_levels
        subplot(l_l+2, 1, find(ilevel==denoised_levels)+1), plot(detcoef(c{ichannel}, l{ichannel}, ilevel)), hold on,
        plot(detcoef(denoise_c{ichannel}, l{ichannel}, ilevel)), legend('original coif', 'denoised coif.'),
        ylabel(['d', num2str(ilevel)]), title([num2str(Fs/2^(ilevel+1)), 'Hz - ', num2str(Fs/2^ilevel), 'Hz']);
    end
    subplot(l_l+2, 1, l_l+2), plot(appcoef(c{ichannel}, l{ichannel}, wname, l_l)), hold on,
    plot(appcoef(denoise_c{ichannel}, l{ichannel}, wname, l_l)), legend('original coif', 'denoised coif.'),
    ylabel(['a', num2str(nlevel)]), title([num2str(0), 'Hz - ', num2str(Fs/2^(nlevel+1)), 'Hz']);
    
end
