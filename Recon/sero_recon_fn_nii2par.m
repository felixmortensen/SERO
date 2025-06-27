function fn_par = sero_recon_fn_nii2par(fn_nii)
% function fn_seq = sero_recon_fn_data2nii(fn_data)

[a,b] = msf_fileparts(fn_nii);
fn_par = [a filesep b '_par.nii.gz'];
