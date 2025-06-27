function h = sero_seq_plot_pns(seq, blockInd, hw)
% function h = sero_seq_plot_pns(seq, blockInd, hw)
%
% This function requires the SAFE nerve stimulation framework:
% https://github.com/filip-szczepankiewicz/safe_pns_prediction

if nargin < 2 || isempty(blockInd)
    blockInd = [1 numel(seq.blockEvents)];
end

if nargin < 3
    hw = safe_example_hw;
end

[gwf, rf, dt] = sero_seq_getBlockWaveforms(seq, blockInd);

pns = safe_gwf_to_pns(gwf, rf, dt, hw, 1);
h = safe_plot(pns, dt);
