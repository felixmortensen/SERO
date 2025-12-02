function [t, tr, b, w] = sero_srr_sampling_conventional_slice_shift(ts, n, a, z, do_interleaved)

if nargin < 1
    [ts, n, z, a] = sero_srr_default_image_pars();
    [t, tr, b, w] = sero_srr_sampling_conventional_slice_shift(ts, n, a, z, 1);
    
    clf
    sero_srr_plot_sampling_stats_v2(tr, b, w)
    
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
    [ct, ctr, cb] = sero_srr_sampling_conventional(ts, round(n/a), a, z, do_interleaved, start_pos);
    
    t  = [t; ct];
    tr = [tr; ctr];
    b  = [b; cb];
    
end

w = double(tr>0);