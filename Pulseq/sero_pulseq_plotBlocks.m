function h = sero_pulseq_plotBlocks(seq, blockInd)
% function h = sero_pulseq_plotBlocks(seq, blockInd)

if nargin < 2 || isempty(blockInd)
    blockInd = [1 numel(seq.blockEvents)];
end

[gwf, rf, dt, ti, rwf] = sero_pulseq_getBlockWaveforms(seq, blockInd);

col = pulseq_colXYZR();

tStart = sum(seq.blockDurations(1:(blockInd(1))))-seq.blockDurations(blockInd(1));


% PLOT GWF
hold on
for i = 1:size(gwf,2)
    h(i) = plot(ti+tStart, gwf(:,i), '-', 'color', col(i,:));
end

yl = get(gca, 'YLim');


% PLOT DEPHASING DIR
hdd = plot(ti+tStart, rf*yl(2)/2, 'k--');


% PLOT RFWF
y  = real(rwf);
y  = y/max(abs(y(:)))*yl(2);
hrf1 = plot(ti+tStart, y(:,1), '-', 'color', [col(4,:) 0.1]);
hrf2 = plot(ti+tStart, y(:,2), '-', 'color', [col(4,:) 1]);


% ADD LABELS
ylabel('g(t) [mT/m] and rf(t) [a.u.]')
xlabel('Time [s]')

h   = cat(1, h(:), hdd, hrf1, hrf2);

indEx = sero_pulseq_getRFBlockInd(seq, blockInd);

for i = 1:numel(indEx)
    rf90 = seq.getBlock(indEx(i));
    trf = sum(seq.blockDurations(1:indEx(i)));
    text(trf, yl(2)*0.55, ['\Deltaz = ' num2str(rf90.rf.freqOffset/rf90.gz.waveform(2)*1e3, '%0.1f') ' mm'], 'hori', 'center')
end



