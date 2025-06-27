function E = sero_seq_trb2pos_v2(tr, b, dz, spacing)
% function E =  sero_seq_trb2pos_with_gap(tr, b, dz, gap_factor)
% tr = position-wise repetition times in seconds
% b  = shot-wise b-values in s/m2 with arbitrary scale
% dz = step distance in meters
% spacing = factor to determine the gap between positions (1 for full gap)

if nargin < 1
    [~, tr, b] = srp_sampling_conventional(0.15,1000,1,50);
    dz = 1.5e-3;
    spacing = 0;
end

if nargin < 4
    spacing = 0;
end

% Number of unique positions and shots
n_pos = size(tr,2);
n_shot = size(tr,1);

% Position list
pos_l = (1:n_pos) * dz * (1 + spacing) - dz/2;
pos_l = pos_l - min(pos_l) - range(pos_l)/2;

% Slice center
pos = (tr > 0) .* pos_l;
pos(pos == 0) = nan;
mpos = nanmean(pos, 2);
mthc = sum(tr>0, 2)*dz;

% Relative diffusion gradient amplitudes
gradientAmp = ones(n_shot, 3) .* sqrt(b./max(b)) / sqrt(3);

% Relative spoiler gradient amplitudes
spoilerAmp = ones(n_shot, 3);

% Output sampling scheme
E = [mpos mthc gradientAmp spoilerAmp];
end