clear, clc

%raw data to stack of 2D images

fn_dat = 'location of raw data file';
fn_seq = 'SERO_example_sequence.seq';
fn_xps = 'location of xps file';
fn_nii = '2D_example.nii';

sero_pulseq_data2nii_adaptiveCombine(fn_dat, fn_seq, fn_nii, fn_xps)

%%

%2D stack to SRR

new_nii = 'SRR_example.nii';
mask_fn = 'location of mask file';

regstr = 1e-2;
do_t1 = 1;
mode = 1;

sero_srr_nii2fit(fn_nii, fn_xps, mask_fn, regstr, do_t1, mode, new_nii)