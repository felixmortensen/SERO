function fn_par = sero_srr_nii2fit_voxel(nii_fn, xps_fn, mask_fn, reg_str, do_t1, mode, new_nii, opt)
%function fn_par = sero_srr_nii2fit_voxel(nii_fn, xps_fn, mask_fn, reg_str, do_t1, mode, new_nii, opt)


if isempty(xps_fn) || nargin < 2
    xps_fn = mdm_fn_nii2xps(nii_fn);
end

if isempty(mask_fn) || nargin < 3
    mask_fn = [];
end

if isempty(reg_str) || nargin < 4
    reg_str=0;
end

if isempty(do_t1) || nargin < 5
    do_t1 = 1;
end

if nargin < 8
    opt = srp_data2fit_opt();
end

%Load and unpack experimental pars
xps = sero_srr_xps_load(xps_fn);

TR = xps.tr;
B  = xps.b;
W  = xps.w;
xps.lambda = reg_str;

TR(isnan(TR)) = 0;
TR(isinf(TR)) = 10;

%Load signal data
I   = mdm_nii_read(nii_fn);
I = abs(I);
siz = size(I);


% Mask
% Load mask and create weight vectors
if ~isempty(mask_fn)
    M    = mdm_nii_read(mask_fn);
    %WGHT = srp_mask2wght(M, TR, 2); %This doens't work
    WGHT = ones(size(TR));
else
    M    = ones([siz(1) siz(2) size(TR,2)]);
    WGHT = ones(size(TR));
end

% Prelocate parameter output
P = zeros([siz(1) siz(2) size(TR,2) 4]);

% Start fit
n_col = sum(M,3);
n     = sum(n_col(:)>0);
c = 1;

for i = 1:siz(1)
    for j = 1:siz(1)
        for k = 1:50
            % mask = squeeze(M(i,j,k));
            % 
            % if ~any(mask)
            %     continue
            % end

            m  = TR(:,k) > 0;

            Sobs = double(squeeze(I(i,j,m)));

            %slice-profile weights, TRs, b-values for those shots
            tr_k = TR(m,k);         
            b_k  = B(m);
            w_k  = W(m,k);

            x0 = sero_srr_setGuess(Sobs, tr_k, b_k, w_k, do_t1);

            T = sero_srr_data2fit_1d_reg_v2(Sobs, tr_k, b_k, w_k, do_t1, reg_str, 0, x0); %mode = 0

            P(i,j,k,1:4) = single(T(1,1:4));

            disp([num2str(c) '/' num2str(n)])

            c = c+1;
        end
    end
end

mdm_nii_write(flip(P,3), new_nii);
disp(['File ' new_nii ' saved!'])
end
