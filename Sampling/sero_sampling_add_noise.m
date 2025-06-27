function S = sero_sampling_add_noise(S, SNR)
% function S = sero_sampling_add_noise(S, SNR)
% ADD NOISE
% This application has the following interpretation.
% Given that the true signal is 1, the high resolution image has a signal
% to noise ration equal to SNR. The low resolution image has an snr that is
% higher by a factor of the aspect ratio.
S = sqrt((S + randn(size(S))./SNR).^2+(randn(size(S))./SNR).^2);

