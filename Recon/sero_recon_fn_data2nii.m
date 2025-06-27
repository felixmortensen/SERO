function fn_seq = sero_recon_fn_data2nii(fn_data)
% function fn_seq = sero_recon_fn_data2nii(fn_data)

[a,b] = fileparts(fn_data);
fn_seq = [a filesep b '.nii.gz'];
