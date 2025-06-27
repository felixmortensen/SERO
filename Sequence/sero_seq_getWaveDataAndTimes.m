function [wave_data, tfp_excitation, tfp_refocusing] = sero_seq_getWaveDataAndTimes(seq, blockInd)
% function [wave_data, tfp_excitation, tfp_refocusing] = sero_seq_getWaveDataAndTimes(seq, blockInd)

[wave_data, tfp_excitation, tfp_refocusing] = seq.waveforms_and_times(true, blockInd);

[indEx, indRef, indSat] = sero_seq_getRFBlockInd(seq, blockInd);
indAll = sort([indEx(:); indRef(:); indSat(:)]);

startEx  = find(indAll==indEx(1),1);
startRef = find(indAll==indRef(1),1);

if isempty(indSat)
    stepSize  = 2;
else
    stepSize  = 3;
end


if isempty(tfp_refocusing)
    tmp            = tfp_excitation;
    tfp_excitation = tmp(:, startEx:stepSize:end);
    tfp_refocusing = tmp(:,startRef:stepSize:end);
else
    % error('situation not recognized')
end


tr = wave_data{4}(1,:);

rwf = [];
for i = 1:numel(indAll)
    tmp = seq.getBlock(indAll(i));
    rwf = [rwf; 0; tmp.rf.signal(:)];
end

nc = min([numel(tr) numel(rwf)]);

wave_data{5}(1,:) = tr(1:nc);
wave_data{5}(2,:) = rwf(1:nc);

end