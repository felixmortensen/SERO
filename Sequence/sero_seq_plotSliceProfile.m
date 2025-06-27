function sliceThick = sero_seq_plotSliceProfile(rfShape, rfTime, gAmp, flipAngle, sliceThick, dz, do_rewind)
% function pulseq_plotSliceProfile(rfShape, rfTime, gAmp, flipAngle, sliceThick, dz, do_rewind)

[M, z, sliceThick] = pulseq_sliceProfile(rfShape, rfTime*1000, flipAngle, mr.convert(gAmp, 'Hz/m', 'mT/m'), sliceThick, do_rewind);

Mc  = M(1,:) + M(2,:)*1i;
Mxy = sqrt(M(1,:).^2 + M(2,:).^2);
Mz  = M(3,:);

Mopt = zeros(size(z));
Mopt(abs(z)<=(sliceThick/2)) = 1;

[weight, bins] = srp_sliceProf2weight(Mc, z, dz, sliceThick, 0);

%% PLOT
subplot(2,2,1)
plot(rfTime*1000, rfShape)
title('RF-pulse')
ylabel('Amplitude [a.u.]')
xlabel('Time [ms]')


subplot(2,2,2)
plot(z*1e3, Mxy, z*1e3, Mz, z*1e3, Mopt, 'k:')
title(['Slice profile (ST= ' num2str(sliceThick*1e3) ' mm)'])
ylabel('Magnetization [1]')
xlabel('Position [mm]')
xlim([-1 1]*sliceThick*1e3)

subplot(2,2,3)
plot(z*1e3, abs(Mc), z*1e3, angle(Mc))
title('Slice profile')
ylabel('Magnetization [1]')
xlabel('Position [mm]')


subplot(2,2,4)
xb = bins(1:end-1);
yb = ones(size(xb)) * abs(mean(Mc))/sliceThick*range(z);
stairs(xb, weight); hold on;
plot(xb, yb)
title('Discrete slice profile')
ylabel('Magnetization [1]')
xlabel('Position [mm]')

for i = 1:numel(bins(1:end-1))
    x = mean([bins(i) bins(i+1)]);
    if abs(x)>sliceThick*1e3
        continue
    end
    text(x, weight(i), num2str(weight(i)*100, '%0.0f'), 'hori', 'center', 'vert', 'bottom')
end

ylim([-0.2 max(weight)+0.2])
xlim([-1 1]*sliceThick*1e3)

%
%
% plot(1:n_el, wsl, 'o--', 'MarkerFaceColor','r'); hold on
% plot(2:(n_el-1), wsl(2:(end-1)), 'o-', 'MarkerFaceColor','w'); hold off
% ylim([0 1.1])
% xlim([.5 n_el+.5])
% title('Signal weight around center')
% ylabel('Signal fraction [1]')
% xlabel('Index')