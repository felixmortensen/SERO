function fn_seq = sero_srr_fn_data2nii_aco(fn_data)
% function fn_seq = sero_srr_fn_data2nii_aco(fn_data)

[a,b] = fileparts(fn_data);
fn_seq = [a filesep b '_aco.nii.gz'];
