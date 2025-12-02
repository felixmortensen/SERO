function fn_par = sero_srr_nii2fit_average(nii_fn, xps_fn, mask_fn, do_t1, opt)
% function fn_par = sero_srr_nii2fit(nii_fn, xps_fn, mask_fn, do_t1, opt)

if isempty(xps_fn) || nargin < 2
    xps_fn = mdm_fn_nii2xps(nii_fn);
end

if isempty(mask_fn) || nargin < 3
    mask_fn = [];
end

if nargin < 5
    opt = sero_srr_data2fit_opt();
end

% Load and unpack experimental pars
xps = sero_srr_xps_load(xps_fn);

TR = xps.tr;
B  = xps.b*1e-9;
W  = xps.w;

TR(isnan(TR)) = 0;


% Load signal data
I   = mdm_nii_read(nii_fn);
siz = size(I);


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
    for j = 1:siz(2)

        % m = squeeze(M(i,j,:));
        % 
        % if ~any(m)
        %     continue
        % end

        % wght = squeeze(WGHT(i,j,:));

        S = double(abs(squeeze(I(i,j,:))));

        tic
        T = sero_srr_setGuess(S, TR, B, W, do_t1, 1/25);
        toc

        for k = 1:size(P,4)
            P(i,j,:,k) = T(:,k);
        end

        disp([num2str(c) '/' num2str(n)])

        c = c+1;

    end
end

fn_par = sero_srr_fn_nii2par_avg(fn_nii);
mdm_nii_write(flip(single(P),3), fn_par);
disp(['File ' new_nii ' saved!'])



