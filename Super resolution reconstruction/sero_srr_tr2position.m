function p = sero_srr_tr2position(TR, w)
% function p = sero_srr_tr2position(TR, w)

if nargin < 2
    w = sum(TR>0, 2);
end

ind = 1:size(TR,2);

p = (TR>0) * ind' ./ w;

