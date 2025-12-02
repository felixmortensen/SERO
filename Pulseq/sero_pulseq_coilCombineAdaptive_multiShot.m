function data_out =  sero_pulseq_coilCombineAdaptive_multiShot(data)
% function data_out =  sero_pulseq_coilCombineAdaptive_multiShot(data)
% Assume data is in format [xvox, coil, yvox, shot].
% Needs to be in format [coil, yvox, xvox, shot] for converter.

siz      = size(data);
data_out = zeros(siz([3 1 4]));

for i = 1:siz(4) % loop over shots
    tmp = permute(data(:,:,:,i), [2 1 3]);
    data_out(:,:,i) = sero_pulseq_coilCombineAdaptive(tmp);
    disp(['Adaptive combine of shot ' num2str(i) '/' num2str(siz(4))])
end


