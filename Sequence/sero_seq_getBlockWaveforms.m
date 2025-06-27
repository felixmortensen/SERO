function [gwf, rf, dt, ti, rwf] = sero_seq_getBlockWaveforms(seq, blockInd)
% function [gwf, rf, dt, ti, rwf] = sero_seq_getBlockWaveforms(seq, blockInd)

if nargin < 2 || isempty(blockInd) || all(blockInd == 0)
    blockInd = [1 numel(seq.blockEvents)];
end

if isscalar(blockInd)
    blockInd = [1 1]*blockInd;
end

dt      = seq.gradRasterTime;
[wave_data, tfp_excitation, tfp_refocusing] = sero_seq_getWaveDataAndTimes(seq, blockInd);

ttot = 0;
for i = 1:numel(wave_data)
    ttot = max([ttot max(wave_data{i}(1,:))]);
end

ti   = linspace(0, ttot, round(ttot/seq.gradRasterTime))';


for i = 1:numel(wave_data)
    t = wave_data{i}(1,:);
    g = wave_data{i}(2,:);

    if isempty(t) | isempty(g)
        wf(:,i) = zeros(size(ti));
        continue
    end

    if i < 4
        g = mr.convert(g, 'Hz/m', 'mT/m');
        g = g/1e3; % to T/m
    end

    [ut, uind] = unique(t);

    if numel(ut)<numel(t)
        warning(['Time points are not unique in wave_data{' num2str(i) '}!']);
        t = t(uind);
        g = g(uind);
    end

    try
        wf(:,i) = interp1(t, g, ti);
    catch me
        error(me.message)
        wf(:,i) = zeros(size(ti));
    end

end

wf(isnan(wf)) = 0;

gwf = real(wf(:,1:3));
rwf = wf(:,4:5);

rf  = nan(size(ti));

if ~isempty(tfp_excitation) && ~isempty(tfp_refocusing)
    rf(:) = 0;

    rfTimes = sort([tfp_excitation(1,:) tfp_refocusing(1,:)]);

    for i = 1:numel(rfTimes)
        rf(ti>rfTimes(i)) = (-1)^(i+1);
    end

    rf(ti<tfp_excitation(1,1)) = nan;
end