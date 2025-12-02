function [SNRmean, STD, E] = sero_srr_snr_check_invivo(fn_nii, fn_fit, fn_xps, fn_msk, do_t1)
% function SNR = sero_srr_snr_check_invivo(fn_nii, fn_fit, fn_xps, fn_msk, do_t1)


%read files
I   = mdm_nii_read(fn_nii); I = abs(I);
xps = load(fn_xps);
P   = mdm_nii_read(fn_fit);
S0  = P(:,:,:,1);
M   = mdm_nii_read(fn_msk);

tr = xps.tr;
tr(isnan(tr)) = 0;

P = flip(P,3);

%forward model
E = srp_fit2data_3D( ...
    P(:,:,:,1), P(:,:,:,2), P(:,:,:,3), P(:,:,:,4), ...
    xps.tr, xps.b, xps.tr>0, do_t1, []);


M = P(:,:,:,2) > 0.5 & P(:,:,:,2) < 1.1;

STD = zeros(size(P,[1 2 3]));
ES = STD;

for i = 1:size(S0,1)
    for j = 1:size(S0,2)
        for k = 1:size(S0,3)
            % if ~M(i,j,k)
            %     continue
            % end

            m    = tr(:,k) > 0;
            Sest = squeeze(E(i,j,m>0));
            Sobs = squeeze(I(i,j,m>0));
            if any(isnan(Sest))
                continue
            end
            STD(i,j,k) = std(Sobs - Sest);
            ES(i,j,k) = mean(Sest(:));

            SNRmean(i,j,k) = sqrt(mean((Sest/STD(i,j,k)).^2));
        end
    end
    i
end

mdm_nii_write(flip(SNRmean,3), 'SNRmean.nii');
end
