function signalFilt = lowpassfilt(signal, Fs, lowf)
% low pass filter in the frequency field

l_s = length(signal);
signalFFT = fft(signal);

% filter in frequency field
hz=linspace(0, Fs, l_s);
filterlow = (1-1./(1+exp(-hz+lowf))); 

lp_signalFFT = signalFFT.*filterlow;
signalFilt = real(ifft(lp_signalFFT));

return