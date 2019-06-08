function y = filt50Hz(data1D, Fs)

%this function was written to convert the noisy original s(t) EEG data into
%less noisy x(t) data, to remove the known noise at 50Hz

f0 = 50;                %50 Hz noise from AC source in Germany
fn = Fs/2;              %Nyquist rate
freqRatio = f0/fn;      %normalize
notchWidth = 0.1;       %small width

%compute zeros
notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%and poles
notchPoles = (1-notchWidth) * notchZeros;

%figure()
%zplane(notchZeros.', notchPoles.');

b = poly( notchZeros ); %# Get moving average filter coefficients
a = poly( notchPoles ); %# Get autoregressive filter coefficients

% figure()
% freqz(b,a,Fs)

%#filter signal x
y = filter(b,a,data1D);


end