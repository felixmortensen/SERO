function [M, z, sliceThick] = sero_pulseq_sliceProfile(rf_wf, rf_t, flipAngle, gSlice, sliceThick, do_rewind)
% function [M, z, sliceThick] = sero_pulseq_sliceProfile(rf_wf, rf_t, flipAngle, gSlice, sliceThick, do_rewind)
%
% By Filip Sz but very inspired by https://larsonlab.github.io/MRI-education-resources/
%
% rf_wf       RF pulse in [mT] formated as 1xn
% rf_t        RF time vector in [ms] formated as 1xn
% flipAng     RF flip angle in [degrees]
% gSlice      Gradient amplitude in [mT/m]
% sliceThick  Slice thickness in [m]

if size(rf_wf,1)>1 || size(rf_t,1)>1
    error('Check input size!')
end

if isempty(sliceThick)
    sliceThick = 12e-3;
    do_recalc = 1;
else
    do_recalc = 0;
end

if nargin < 6
    do_rewind = 0;
end

gamma = 42.58; % kHz/mT
dt    = rf_t(2)-rf_t(1); % ms
rf_wf = (flipAngle * pi / 180) * rf_wf / sum(rf_wf) / (2*pi * gamma * dt);

% find center of mass time of RF
rf_ct = rf_wf*rf_t'/sum(rf_wf);
ff    = (max(rf_t)-rf_ct)/dt;

N  = numel(rf_wf);
BW = sliceThick * gamma * gSlice; % kHz
df = linspace(-2*BW, 2*BW, 200); % kHz
z  = df / gSlice / gamma; % m

M  = repmat([0 0 1]', [1 length(df)]);

for n = 1:N
    for f = 1:length(df)
        M(:, f) = bloch_rotate( M(:,f), dt, [real(rf_wf(n))  imag(rf_wf(n))  df(f)/gamma]);
    end
end

if do_rewind
    for f = 1:length(df)
        M(:, f) = bloch_rotate( M(:,f), dt, [0  0  -df(f)/gamma*ff]);
    end
end


if do_recalc
    Mxy = sqrt(M(1,:).^2 + M(2,:).^2);
    [pk,loc,sliceThick_new] = findpeaks(abs(Mxy),z,'MinPeakProminence',.4);
    
    if sliceThick_new>sliceThick
        error('Estimated ST is larger than assumed. Check!')
    end
    sliceThick = sliceThick_new;
end
