function fn_xps = sero_pulseq_fn_seq2xps(fn_seq)
% function fn_xps = sero_pulseq_fn_seq2xps(fn_seq)

[a,b,c] = fileparts(fn_seq);

fn_xps = fullfile(a, [b '_xps.mat']);


