function fn_par = sero_recon_fn_nii2xps(fn_nii)
% function fn_par = sero_recon_fn_nii2xps(fn_nii)

[a,b] = msf_fileparts(fn_nii);
fn_par = [a filesep b '_xps.mat'];