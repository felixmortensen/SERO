function [blockInd, indEx, indRef, indSat] = sero_pulseq_shotToBlockInd(seq, n)
% function [blockInd, indEx, indRef, indSat] = sero_pulseq_shotToBlockInd(seq, n)

if numel(n) == 1 %#ok<ISCL>
    n = [n n];
end

if all(n==0)
    n = [1 numel(seq.blockEvents)];
end

[indEx, indRef, indSat] = sero_pulseq_getRFBlockInd(seq);

nShot = numel(indEx);

indStart = min([indSat; indEx],[],1);

if n(2) < nShot
    ind2 = indStart(n(2)+1)-1;
else
    ind2 = numel(seq.blockEvents);
end

blockInd = [indStart(n(1)) ind2];