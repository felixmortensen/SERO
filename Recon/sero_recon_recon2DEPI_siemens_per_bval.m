function volumes = sero_recon_recon2DEPI_siemens_per_bval(fn_data, fn_seq, xps_fn, rec_delay)
% volumes = pulseq_recon2DEPI_siemens_per_bval(fn_data, fn_seq, xps_fn, rec_delay)
% Reads Siemens raw data and performs ghost correction per unique b-value.
% Based on pulseq_recon2DEPI_siemens.
% All processing downstream (reshaping, FFT, etc.) remains unchanged.

% Set default arguments
if nargin < 2 || isempty(fn_seq)
    fn_seq = pulseq_fn_data2seq(fn_data);
end
if nargin < 3 || isempty(xps_fn)
    error('An xps file containing b-values must be provided!');
end
if nargin < 4
    rec_delay = 'auto';
end

%% Read data and sequence
% Use a special function that converts Siemens raw data to manageable
% format. Please see: https://github.com/CIC-methods/FID-A/blob/master/inputOutput/mapVBVD/README.md

twix_obj = mapVBVD(fn_data);
if iscell(twix_obj)
    data_matrix = single(twix_obj{end}.image.unsorted());
else
    data_matrix = single(twix_obj.image.unsorted());
end
clear twix_obj

seq_obj = mr.Sequence();
seq_obj.read(fn_seq, 'detectRFuse');
ktr_adc_global = seq_obj.calculateKspacePP();

nADC = size(data_matrix, 1);         %number of frequency samples
[ktr_adc_temp, t_adc] = seq_obj.calculateKspacePP();
dwell_time = (t_adc(nADC) - t_adc(1)) / (nADC - 1);
clear t_adc

bvals = load(xps_fn).xps.b;           
nAcq = size(data_matrix, 3);

% Replicate bvals so that each image's b-value repeats for each phase encoding step
repFactor = nAcq / numel(bvals);       
k = ones(1, repFactor);
bRep = k' * bvals';    
bRep = squeeze(bRep(:));
unique_b = unique(bRep, 'stable'); %preserver original acquitisiton order
ktr_adc_all = cell(length(unique_b), 1);  %to store ghost-corrected trajectories

% Ghost correction per unique b-value.
switch rec_delay
    case 0
        fprintf('Skipping ghost correction (no delay estimation)!\n');
        for i = 1:length(unique_b)
            ktr_adc_all{i} = ktr_adc_global; 
        end

    case 'auto'
        fprintf('Performing ghost correction per unique b-value!\n');
        for i = 1:length(unique_b)
            idx = find(bRep == unique_b(i)); %idx over total acquisitions
            data_subset = data_matrix(:, :, idx);
            
            % Separate odd and even echoes over the whole block.
            oddIdx  = 1:2:size(data_subset, 3);
            evenIdx = 2:2:size(data_subset, 3);

            if isempty(evenIdx) || isempty(oddIdx)
                warning('Not enough echoes for b = %g', unique_b(i));
                ktr_adc_all{i} = ktr_adc_global;  %fallback to global trajectory
                continue;
            end
            
            % Estimate ghost delay from the entire set of acquisitions for this b-value.
            A = ifftshift(ifft(ifftshift(data_subset(end:-1:1, :, evenIdx), 1)), 1);
            B = ifftshift(ifft(ifftshift(data_subset(:, :, oddIdx), 1)), 1);
            
            tmp = A .* conj(B);
            tmp = tmp(2:end, :, :) .* conj(tmp(1:end-1, :, :));
            phase = angle(sum(tmp(:)));
            
            rec_delay_est = double(phase / (4 * pi) * nADC * dwell_time);
            fprintf('b = %g, delay = %g s\n', unique_b(i), rec_delay_est);
            
            % Compute the ghost-corrected trajectory for this b-value
            ktr_adc_subset = seq_obj.calculateKspacePP('trajectory_delay', rec_delay_est);
            ktr_adc_all{i} = ktr_adc_subset;
        end

    otherwise
        fprintf('Applying manual delay of %g s!\n', rec_delay);
        for i = 1:length(unique_b)
            ktr_adc_all{i} = seq_obj.calculateKspacePP('trajectory_delay', rec_delay);
        end
end

ghostCorrectedData = data_matrix;

%% Reconstruction: One Volume Per b-value

% Preallocate a cell array to hold each reconstructed volume.
volumes = cell(length(unique_b), 1);

% Construct volume for each b-value
for i = 1:length(unique_b)
    idx = find(bRep == unique_b(i));  %these are the acquisitions for b-value i
    data_subset = ghostCorrectedData(:, :, idx);
    
    nPE = nAcq/length(bvals); %phase-encoding lines per image
    nImages = numel(idx) / nPE;

    if mod(numel(idx), nPE) ~= 0
        warning('Acquisition count for b = %g is not an integer multiple of %d.', unique_b(i), nPE); %ensure number of acquisitions (current b‑value) is divisible by nPE
        continue;
    end
    
    ktr_adc_subset = ktr_adc_all{i};  %get the ghost-corrected trajectory for this b-value
    
    %% Calculate parameters from the corrected trajectory for this group.
    k_last = ktr_adc_subset(:, end);
    k_2last = ktr_adc_subset(:, end - nADC);
    delta_ky = k_last(2) - k_2last(2);
    Ny_post = round(abs(k_last(2) / delta_ky));

    if k_last(2) > 0
         Ny_pre = round(abs(min(ktr_adc_subset(2, :)) / delta_ky));
    else
         Ny_pre = round(abs(max(ktr_adc_subset(2, :)) / delta_ky));
    end

    Nx = 2 * max([Ny_post, Ny_pre]);
    Ny = Nx;
    Ny_sampled = Ny_pre + Ny_post + 1;
    nkxx = Nx;  %number of samples along readout.
    
    % Create a uniform kx grid using this group's parameters.
    kxmin = min(ktr_adc_subset(1, :));
    kxmax = max(ktr_adc_subset(1, :));
    kxmax1 = kxmax / (Nx/2 - 1) * (Nx/2);  %adjust for FFT center
    kmaxabs = max(kxmax1, -kxmin);
    kxx = ((-Nx/2):(Nx/2-1)) / (Nx/2) * kmaxabs;
    

    % Reshape the data into [nADC x nCoils x nPE x nImages]
    data_subset_reshaped = reshape(data_subset, [nADC, size(data_subset, 2), nPE, nImages]);
        
    ktr_adc2 = reshape(ktr_adc_subset, [size(ktr_adc_subset,1), nADC, size(ktr_adc_subset,2)/nADC]); %reshape the ghost-corrected trajectory for the readout.
    
    interpData = zeros(nkxx, size(data_subset,2), nPE, nImages, 'like', data_subset_reshaped);

    for img = 1:nImages
        for pe = 1:nPE
            overall_idx = (img - 1) * nPE + pe; %index corresponding to the current echo line
            for c = 1:size(data_subset,2)
                interpData(:, c, pe, img) = interp1( ktr_adc2(1, :, overall_idx), data_subset_reshaped(:, c, pe, img), kxx, 'spline', 0);
            end
        end
    end
    
    if size(interpData,1) > nkxx
        interpData((nkxx+1):end,:,:,:) = []; %remove excess rows in readout dimension
    end

    %% Reshape for multiple slices or repetitions
    imgVolume = ifftshift(ifft(ifftshift(interpData, 1)), 1);

    if Ny_sampled < Ny
        imgVolume = padarray(imgVolume, [0 0 (Ny - Ny_sampled) 0], 'pre'); %zero padding for phase encodign lines
    end

    imgVolume = ifftshift(ifft(ifftshift(imgVolume, 3), [], 3), 3);
    
    volumes{i} = imgVolume;
end
