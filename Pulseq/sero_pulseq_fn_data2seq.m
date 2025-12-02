function fn_seq = sero_pulseq_fn_data2seq(fn_data)
% function fn_seq = sero_pulseq_fn_data2seq(fn_data)

[a,b] = fileparts(fn_data);
fn_seq = [a filesep b '.seq'];
