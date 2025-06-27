clear, clc

% snr at signal = 1
snr = 100;

rng(42)

% Get defailt imaging parameters
[timeStep, nSamples, z, ar] = sero_sampling_default_image_pars();

% create line phantom object
T = sero_sampling_generate_structure(z, 1);

% Create SERO sampling matrix
[~, TR,  B,  W ] = sero_sampling_sero_gen1(timeStep, nSamples, ar, z);

TR(isnan(TR)) = 0;


% Get signal
S   = sero_sampling_fit2data(T, TR,  B,  1, W);

% Add noise to signal 
N   = sero_sampling_add_noise(S,  snr);

% do fit
m  = sero_recon_data2fit_1d_reg_v2(N,  TR,  B,  W,  1, 1e-2, 1);


%% Plot result

pnam = {'s0', 'D', 'T1', 'V'};

clf
for i = 1:4
    subplot(2,2,i)
    plot(1:z, T(:,i), 'b.'); hold on
    plot(1:z, mean(m(:,i,:),3), 'k')
    % plot(1:z, mean(mr(:,i,:),3), 'r')
    title(pnam{i})
    xlim([0 z+1])
    legend('Ground Truth', 'SERO')
end

fh = gcf;
fh.Theme = 'light';

% 1/m(1,5)