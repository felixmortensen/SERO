function [t, tr, b, w] = sero_sampling_sampling_conventional(ts, n, a, z, do_interleaved, start_pos)

if nargin < 1
    [ts, n, z] = sero_sampling_default_image_pars();
    [t, tr, b, w] = sero_sampling_sampling_conventional(ts, n, 1, z, 0);
    a=1;
    clf
    sero_sampling_plot_sampling_stats_v2(tr, b, w)
    return
end

if nargin < 5
    do_interleaved = 0;
end

if nargin < 6
    start_pos = 1;
end

if start_pos > a
    error('makes no sense')
end

% bv = [1 .5 0];
bv = [2 1.5 1 0.5 0];

b = repmat(bv, z, ceil(n/z));
b = b(1:n);
b = b(:);

t = nan(n,z);

pos = start_pos:a:(z-a+1);

x0l = repmat(pos, ceil(n/numel(pos)),1);

if do_interleaved
    x0l = x0l(:, [1:2:end 2:2:end]); % interleaved
end

x0l = x0l';
x0l = x0l(:);

for i = 1:n
    
    x0 = x0l(i);
    
    ind = x0:(x0+a-1);
    
    t(i, ind) = i*ts;
end


tr = nan(size(t));

for i = 1:size(t,2)
    
    ind = find(~isnan(t(:,i)));
    
    ds = diff([-10; ind]); % FIXME!!!
    
    tr(ind, i) = ds * ts;
end

w = double(tr>0);

% Fix steady state
tr(~isnan(tr)) = ts*z/a;