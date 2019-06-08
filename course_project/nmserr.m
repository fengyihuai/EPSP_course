function  NMSE = nmserr(z1, z2)
% reference: Akhtar, M.T., W. Mitsuhashi and C.J. James, Employing spatially 
% constrained ICA and wavelet denoising, for automatic removal of artifacts 
% from multichannel EEG data. Signal Processing, 2012. 92(2): p. 401-416.

if length(z1) ~= length(z2)
    error('input z1, z2 must have the same length!')
end
n = length(z1);
sum_ratio = sum((z1-z2).^2)/sum(z1.^2);
%  mean:E{} denotes mathematical expectation
NMSE = 20 * log10(mean(sum_ratio));

end