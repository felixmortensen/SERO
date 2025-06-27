function [t, tr, b, w, pos] = sero_sampling_sero_gen1(ts, n, a, z, in_seed)
%Generate first generation SERO sampling scheme with predefined positions,
%TR, b, w

if nargin < 1
    [ts, n, z, a] = sero_sampling_default_image_pars();
    [t, tr, b, w] = sero_sampling_sero_gen1(ts, n, a, z);
    
    clf
    sero_sampling_plot_sampling_stats_v2(tr, b, w)
    
    return
end

keep_seed = rng;

if nargin < 5 || isempty(in_seed)
    rng(42);
end




bv = [1.2 0.8 0.4 0.1];
% bv = [2 1.4 0.7 0.1];

b  = repmat(bv', ceil(n/numel(bv)), 1);
b  = b(1:n);
b  = b(randperm(numel(b)));

t = nan(n,z);

n_history = 10;

xold = ones(1,n_history)*inf;

x0 = repmat(1:(z-a+1), ceil(n/(z-a+1)),1);

x0 = x0(:);

x0l = x0(randperm(numel(x0)));

pos = zeros(n, 1);

for i = 1:n
    
    x0 = x0l(i);
    
    for j = 1:n_history
        if any( abs(x0-xold)<(a/2) )
            tmp =  x0l(i:end);
            x0l(i:end) = tmp(randperm(numel(tmp)));
            x0 = x0l(i);
        end
    end
    
    ind = x0:(x0+a-1);

    pos(i) = (mean(ind)-1/2)-z/2;
    
    t(i, ind) = i*ts;
    
    xold(1)=[];
    xold(end+1) = x0;
end


tr = nan(size(t));

for i = 1:size(t,2)
    
    ind = find(~isnan(t(:,i)));
    
    ds = diff([-inf; ind]); % FIXME!!!
    
    tr(ind, i) = ds * ts;
end

w = double(tr>0);

rng(keep_seed)

