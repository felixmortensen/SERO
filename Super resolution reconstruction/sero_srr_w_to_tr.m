function [tr, t] = sero_srr_w_to_tr(w, t_shot)
% function [tr, t] = sero_srr_w_to_tr(w, t_shot)

t = w .* (1:size(w,1))' * t_shot;

tr = nan(size(t));

for i = 1:size(t,2)
    
    ind = find(t(:,i));
    
    ds = diff([-inf; ind]); % FIXME!!!
    
    tr(ind, i) = ds * t_shot;
end

tr(isnan(tr)) = 0;