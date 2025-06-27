function B = sero_seq_calculateBtensor(seq)
% function B = pulseq_calculateBtensor(seq)
%
% Returns b-tensor in 1x6 "Voight-like" format
% This function may not be general!

if ischar(seq) || isstring(seq)
    tmp = mr.Sequence;
    tmp.read(seq);
    seq = tmp;
end

ind   = pulseq_getExcitationBlockInd(seq);
nShot = numel(ind);
B     = zeros(nShot, 6);

for i = 1:nShot
    [gwf, rf, dt] = pulseq_getShotWaveforms(seq, i);

    tmp = fwf_gwf_to_btens(gwf, rf, dt, 2*pi()*seq.sys.gamma);

    % Here we save the b-tensor in a voight-like format that is compatible
    % with the md-dmri framework.
    B(i,:) = tmp([1 5 9 2 3 6]) .* [1 1 1 sqrt(2) sqrt(2) sqrt(2)];
end

