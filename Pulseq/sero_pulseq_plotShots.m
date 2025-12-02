function h = sero_pulseq_plotShots(seq, n)
% function h = sero_pulseq_plotShots(seq, n)
%
% seq is the pulseq sequence object
% n is the 1x2 shot interval vector [start end].

if nargin < 2
    n = 1;
end

if numel(n) == 1 %#ok<ISCL>
    n = [n n];
end

indBlock = sero_pulseq_shotToBlockInd(seq, n);

h = sero_pulseq_plotBlocks(seq, indBlock);