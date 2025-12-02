function [wght, w] = sero_srr_weight_f(TR, edg)
% function [wght, w] = sero_srr_weight_f(TR, edg)

if nargin < 2
    w = max(sum(TR>0,2));
    edg = 1 * w;
end

p = sero_srr_tr2position(TR);

mid = (max(p)+min(p))/2;
dfe = mid - abs(p-mid);
wght = dfe/edg;
wght(wght>1) = 1;
