function xps = sero_seq_SERO_seq2xps(seq)
% function xps = sero_seq_SERO_seq2xps(seq)
%
% This function creates an experimental parameter structure from the
% sequence object or sequence file path. Note that the funciton is mainly
% adapted fro the SERO project and may therefore lack general function.

if isstring(seq) || ischar(seq)
    fn_seq = seq;
    seq = mr.Sequence();
    seq.read(fn_seq);
end


% Get rf events and set indices based on spin echo design
[indEx, indRef, indSat] = sero_seq_getRFBlockInd(seq);
nShot = numel(indEx);


% Gradient encoding
gwf = cell(nShot,1);
rf  = cell(nShot,1);
dt  = cell(nShot,1);

for i = 1:nShot
    [gwf{i}, rf{i}, dt{i}] = sero_seq_getShotWaveforms(seq, i);
end

xps = fwf_xps_from_gwfl(gwf, rf, dt, 2*pi*seq.sys.gamma);


% get position encoding
[~, ~, ~, ~, ~, ~, slicepos, t_slicepos] = seq.calculateKspacePP();


% Decide on behavior which depends on if seq is structure direct from code
% or loaded from file. They behave differently for some reason.
if size(slicepos,2) == nShot
    i_start = 1;
    i_step  = 1;

elseif isempty(indSat)
    i_start = 1;
    i_step  = 2;

else
    i_start = 2;
    i_step  = 3;

end

slicepos   = slicepos(3, i_start:i_step:end);
t_slicepos = t_slicepos(i_start:i_step:end);

[upos_abs, ~, pos_ind] = uniquetol(slicepos, 0.001);

upos = round((upos_abs-min(upos_abs))*10000)/10000;
pos  = upos(pos_ind);
dz   = median(diff(upos));

[slt, W] = sero_seq_seq2sliceThickness(seq, dz);

ar   = round(slt/dz);
npos = round(max((pos(:)+slt(:))/dz));

t    = nan(xps.n, npos);
w    = zeros(xps.n, npos);

for i = 1:xps.n
    tmp = upos(pos_ind(i));
    ind = (1:ar(i))+round(tmp/dz);
    t(i, ind) = t_slicepos(i);

    wInd = floor(numel(W{i})/2-ar(i)/2)+1;
    wInd = wInd:(wInd+ar(i)-1);
    w(i, ind) = W{i}(wInd);
end

tr = sero_recon_time2tr(t);


% Export
xps.t  = t;
xps.tr = tr;
xps.w  = w;
xps.seqHash = seq.signatureValue;