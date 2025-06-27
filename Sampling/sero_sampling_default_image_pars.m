function [ts, n, z, a] = sero_sampling_default_image_pars()

ttot = 300;
ts = 0.15;
n = round(ttot/ts);
z = 50;
a = 4;