function fn_nii = sero_pulseq_data2nii_adaptiveCombine(fn_data, fn_seq, fn_nii, xps_fn, do_norm)
% function fn_nii = sero_pulseq_data2nii_adaptiveCombine(fn_data, fn_seq, fn_nii, do_norm)

if nargin < 2
    fn_seq = [];
end

if nargin < 3 || isempty(fn_nii)
    fn_nii = sero_pulseq_fn_data2nii_aco(fn_data);
end

if nargin < 4 || isempty(xps_fn)
    error('Provide xps file containing b-values!')
end

if nargin < 5
    do_norm = 1;
end

img_mc = sero_pulseq_recon2DEPI_siemens_per_bval(fn_data, fn_seq, xps_fn, 'auto');

xps = load(xps_fn);
bvals = xps.xps.b;
nImagesTotal = numel(bvals);

unique_b = unique(bvals, 'stable'); %preserve original acquisition order

% Preallocate the final 4D array.
[nx, nCoils, nPE, ~] = size(img_mc{1});
fullVolume = zeros(nx, nCoils, nPE, nImagesTotal, 'like', img_mc{1}); %ensure same data type as volumes

for i = 1:numel(unique_b)
   
    im_idx = find(bvals == unique_b(i));

    if numel(im_idx) ~= size(img_mc{i}, 4)
        error('Mismatch: Found %d images with b = %g but pulseqstructed volume has %d images.', ... %dimensionality check
             numel(im_idx), unique_b(i), size(img_mc{i}, 4));
    end

    fullVolume(:,:,:,im_idx) = img_mc{i};
end

img_mc = fullVolume; %dimensions [nx, nCoils, nPE, nImagesTotal]

if do_norm
    img_mc = img_mc / mean(abs(img_mc(:)));
end

img_co = sero_pulseq_coilCombineAdaptive_multiShot(img_mc);

mdm_nii_write(single(img_co), fn_nii)