function this_path = sero_seq_goToFolder()
% function this_path = pulseq_goToFolder()

this_path = fileparts(mfilename('fullpath'));
cd(this_path)
