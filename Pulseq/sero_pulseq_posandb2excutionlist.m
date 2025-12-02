function exel = sero_pulseq_posandb2excutionlist(posi, sthi, bval, dz)
% function exel = sero_pulseq_posandb2excutionlist(posi, sthi, bval, dz)
% posi = position indices
% sthi = slice thickness in number of voxel indices
% bval = shot-wise b-values in s/m2 with arbitrary scale
% dz   = step distance in meters

% Number of unique positions and shots
n_shot = size(posi,1);

% Position list
pos = posi(:) * dz;
pos = pos - min(pos)-range(pos)/2; % center on middle of interval

% Thickness list
slth = sthi(:) * dz;

% Relative diffusion gradient amplitudes
gradientAmp = ones(n_shot, 3) .* sqrt(bval./max(bval)) / sqrt(3);

% Relative spoiler gradient amplitudes
spoilerAmp = ones(n_shot, 3);

% Output sampling scheme
exel = [pos slth gradientAmp spoilerAmp];
