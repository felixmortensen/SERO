function [t, tr, b, w] = sero_sampling_sampling_direct(ts, n, z, do_interleaved)

if nargin < 1
    [ts, n, z] = sero_sampling_default_image_pars();
    [t, tr, b, w] = sero_sampling_sampling_direct(ts, n, z, 0);
    
    clf
    sero_sampling_plot_sampling_stats(tr, b, w)
    
    return
end

if nargin < 4
    do_interleaved = 0;
end

a = 1;

[t, tr, b, w] = sero_sampling_sampling_conventional(ts, n, a, z, do_interleaved);