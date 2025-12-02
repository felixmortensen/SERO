function [gwf, rf, dt, ti, rwf] = sero_pulseq_getShotWaveforms(seq, n)
% function [gwf, rf, dt, ti, rwf] = sero_pulseq_getShotWaveforms(seq, n)

if ~isscalar(n)
    error('This function is intended for a single shot at a time!')
end

blockInd = sero_pulseq_shotToBlockInd(seq, n);

[gwf, rf, dt, ti, rwf] = sero_pulseq_getBlockWaveforms(seq, blockInd);