function [t, tr, b, w] = sero_sampling_sampling_conventional_slice_shift(ts, n, a, z, do_interleaved)

if nargin < 1
    [ts, n, z, a] = sero_sampling_default_image_pars();
    [t, tr, b, w] = sero_sampling_sampling_conventional_slice_shift(ts, n, a, z, 1);
    
    clf
    sero_sampling_plot_sampling_stats_v2(tr, b, w)
    
    return
end

if nargin < 5
    do_interleaved = 0;
end

t  = [];
tr = [];
b  = [];

for i = 1:a
    
    start_pos = i;
    [ct, ctr, cb] = sero_sampling_sampling_conventional(ts, round(n/a), a, z, do_interleaved, start_pos);
    
    t  = [t; ct];
    tr = [tr; ctr];
    b  = [b; cb];
    
end

w = double(tr>0);