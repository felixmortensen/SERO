function val = sero_seq_roundToRaster(val, dt)
% function val = pulseq_roundToRaster(val, dt)

val = round(val/dt)*dt;