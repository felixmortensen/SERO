function val = sero_pulseq_roundToRaster(val, dt)
% function val = pulseq_roundToRaster(val, dt)

val = round(val/dt)*dt;