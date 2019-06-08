%% generate a square wave, fundaFreq = 20 Hz
% Fs = 1000 Hz; Amp = 1; simDuration  = 2 Seconds;

% hardcoded parameters
Fs = 200.; % Hz; 
fundaFreq = 20; % Hz
Amp = 3;  
sigDuration  = 2; % Seconds;
dutyCycle = 50; % percentage

% derived variant values
tTick = [1/Fs:1/Fs:sigDuration];
sigLength = length(tTick);
fTick = [0:sigLength-1]/sigLength*Fs;

sqrWave = Amp * sin(2*pi* fundaFreq*tTick); %, dutyCycle);

sqrWaveFFT = fft(sqrWave);

figure;
subplot(2, 1, 1); plot(tTick, sqrWave); ylim([-Amp-0.2, Amp+0.2]);
subplot(2, 1, 2); plot(fTick, abs(sqrWaveFFT), '*'); xlim([0, 200]);
hold on; stem(fTick, abs(sqrWaveFFT))



%% white noise contaminated square wave and wavelet denoising
nLevels = 5;

sqrWave_with_Noise = sqrWave + randn(1, sigLength);
sigLengh4swt = floor(sigLength/ 2^nLevels) * 2^nLevels;

SWC = swt(sqrWave_with_Noise(1:sigLengh4swt), nLevels, 'db3');

figure; 
subplot(nLevels + 2, 1, 1); plot(sqrWave_with_Noise(1:sigLengh4swt));
for iLevel = 1:nLevels+1
    subplot(nLevels + 2, 1, 1+iLevel); plot(SWC(iLevel, :));
end

SWC_filtered = SWC; 
% wavelet filtering - the simplest version

SWC_filtered(1, :) = 0;
SWC_filtered(2, :) = 0;
SWC_filtered(4, :) = 0; 
SWC_filtered(5, :) = 0;
SWC_filtered(6, :) = 0; 

sqrWave_wFiltered = iswt(SWC_filtered, 'db3');


figure;
subplot(2, 1, 1); plot(sqrWave_with_Noise);
subplot(2, 1, 2); plot(sqrWave_wFiltered);

figure; 
plot(sqrWave_with_Noise, '-r'); hold on;
plot(sqrWave_wFiltered, 'b-');





