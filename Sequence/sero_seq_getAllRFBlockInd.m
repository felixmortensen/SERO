function ind = sero_seq_getAllRFBlockInd(seq)
% function ind = sero_seq_getAllRFBlockInd(seq)
% 
% Function returns a list of indices that correspond to RF events.

nBlock = numel(seq.blockEvents);

ind = [];

for i = 1:nBlock
    if all(seq.blockEvents{i}([2 5]))
        ind(end+1) = i;
    end
end
