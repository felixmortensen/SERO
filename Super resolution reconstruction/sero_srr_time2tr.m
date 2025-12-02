function tr = sero_srr_time2tr(t)
% function tr = sero_srr_time2tr(t)

tr = nan(size(t));

for i = 1:size(t,2)
    ind = find(~isnan(t(:,i)));
    tr(ind, i) = diff([-inf; t(ind,i)]);
end