function xps = sero_recon_xps_load(fn_xps)
% function xps = sero_recon_xps_load(fn_xps)

load(fn_xps)

if exist('xps', 'var')
    %do nothing

else
    xps.tr = tr;
    xps.b = b;

    try
        xps.w = w;
    catch
        
    end
end


if ~isfield(xps, 'w')
    xps.w = xps.tr>0;
    warning('Filling in weight matrix from TR>0')
end