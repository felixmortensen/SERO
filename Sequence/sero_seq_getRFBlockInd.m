function [indEx, indRef, indSat] = sero_seq_getRFBlockInd(seq, blockInd)
% function [indEx, indRef, indSat] = sero_seq_getRFBlockInd(seq, blockInd)

if nargin < 2
    blockInd = [];
end

indAll = sero_seq_getAllRFBlockInd(seq);
indEx  = sero_seq_getExcitationBlockInd(seq);

nShot  = numel(indEx);
indRef = zeros(size(indEx));

for i = 1:nShot
    % Assume that the pulse after the excitation is a refocusing pulse
    indRef(i) = indAll(find(indAll>indEx(i), 1));
end

% Check if saturation is used
if numel(indAll)==(numel(indEx)*3)
    indSat = zeros(size(indEx));

    for i = 1:nShot
        % Assume that the pulse before the excitation is a saturation pulse
        indSat(i) = indAll(find(indAll<indEx(i), 1, 'last'));
    end

else
    indSat = [];
end

% Mask
if ~isempty(blockInd)
    indEx (indEx <blockInd(1) | indEx >blockInd(2)) = [];
    indRef(indRef<blockInd(1) | indRef>blockInd(2)) = [];
    indSat(indSat<blockInd(1) | indSat>blockInd(2)) = [];
end