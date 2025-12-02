function ind = sero_pulseq_getExcitationBlockInd(seq)
% function ind = sero_pulseq_getExcitationBlockInd(seq)
% 
% Function returns a list of indices that correspond to excitation events.

nBlock = numel(seq.blockEvents);

ind = [];

for i = 1:nBlock
    if all(seq.blockEvents{i}([2 5 7]))
        ind(end+1) = i;
    end
end
