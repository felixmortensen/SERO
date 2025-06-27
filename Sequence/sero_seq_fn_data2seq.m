function fn_seq = sero_seq_fn_data2seq(fn_data)
% function fn_seq = pulseq_fn_data2seq(fn_data)

[a,b] = fileparts(fn_data);
fn_seq = [a filesep b '.seq'];
