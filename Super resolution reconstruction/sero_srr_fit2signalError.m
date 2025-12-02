function [Ez, D] = sero_srr_fit2signalError(m, TR, B, W, S)
% function Ez = sero_srr_fit2signalError(m, TR, B, W, S)
%
% m are the estimated high resolution parameters nz x 4, where nz is the
% number of high resolution positions along z
% TR is the TR matrix, n x nz, where n is the number of shots
% B is the b-value vector, n x 1
% W is the slice profile matrix, n x nz
% S is the measured signal, n x 1

if nargin < 1 % example
    snr = 100;
    [ts, n, z, a] = sero_srr_default_image_pars();

    % create object
    T = sero_srr_generate_structure(z, 1);

    % Create SERO sampling matrix
    [~, TR,  B,  W ] = sero_srr_sampling_v3(ts, n, a, z);
    TR(isnan(TR)) = 0;

    % Get signal
    S   = sero_srr_fit2data(T, TR,  B,  1, W);

    % Add noise to signal
    S   = sero_srr_add_noise(S,  snr);

    % Estimate parameters
    m  = sero_srr_data2fit_1d (S,  TR,  B,  W,  1, [], [], []);

    m(45) = m(45)*2;
    %     m(70) = m(70)/20;
    %     m(140)= m(140)*40;

    % Run function
    tic
    [Ez, D] = sero_srr_fit2signalError(m, TR, B, W, S);
    toc
    % get position
    p  = sero_srr_tr2position(TR);

    % Plot
    clf
    subplot(2,1,1)
    plot(p, D, 'k.')
    hold on
    plot(Ez, '.k-')
    legend('Bias vs pos', 'RMSE')
    title('Error vs position')

    subplot(2,1,2)
    hist(D, 50)
    title(['Signal error distribution (\mu = ' num2str(mean(D)) ' \sigma = ' num2str(std(D)) ')'])
    return
end

% Estimated signal
E  = sero_srr_fit2data(m, TR, B, 1, W);

% Signal difference
D  = E-S;

% Average position of signals
p  = sero_srr_tr2position(TR);

% Position based weight
pp = repmat(p, 1, size(TR,2)) - repmat(1:size(TR,2), size(TR,1), 1);
s  = 2;
ww = exp( -pp.^2 / (2*s^2) );

% Weighted squared difference along z-direction
Ez = sqrt((D.^2)' * (ww ./ sum(ww,1)));






