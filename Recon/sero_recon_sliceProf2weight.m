function [weight, bins] = sero_recon_sliceProf2weight(slp, z, dz, sliceThick, do_plot)
% function [weight, bins] = sero_recon_sliceProf2weight(slp, z, dz, sliceThick, do_plot)

ar = round(sliceThick/dz);

z_mm  = z*1e3;
dz_mm = dz*1e3;

if iseven(ar)
    bins   = (ceil(min(z)/dz):1:ceil(max(z)/dz))*dz_mm;
    weight = zeros(numel(bins)-1, 1);
else
    bins   = ((ceil(min(z)/dz)+1/2):1:ceil(max(z)/dz))*dz_mm;
    weight = zeros(numel(bins)-1, 1);
end


for i = 1:(numel(bins)-1)

    ind = z_mm>=bins(i) & z_mm<=bins(i+1);
    weight(i) = abs(mean(slp(ind)));

end


if do_plot
    stairs(bins(1:end-1), weight)
    ylabel('|M_{xy}|')
    xlabel('Position [mm]')
    title('Discrete slice profile')
end