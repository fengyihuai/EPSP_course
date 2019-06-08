% iwthresh.m test
clc, clear, close all;
% Generate signal and set threshold. 
y = linspace(-1,1,100); 
thr = 0.4;

% Perform hard thresholding. 
ythard = iwthresh(y,'h',thr);

% Perform soft thresholding. 
ytsoft = iwthresh(y,'s',thr);

% Using some plotting commands,
% the following figure is generated.

figure, subplot(131), plot(y), grid minor, title('Original Signal');
subplot(132), plot(ythard), ylim([-1 1]), grid minor, title('Hard Threshold');
subplot(133), plot(ytsoft), ylim([-1 1]), grid minor, title('Soft Threshold');