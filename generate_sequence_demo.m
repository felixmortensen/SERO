clear,clc

fn_seq = 'SERO_example_sequence.seq';

% Basic settings
dz = 1.75e-3;
nPos = 50;
ar = 4; 

% Hardware limits
lims = sero_seq_sysLims_prisma();

% Sequence options
opt = sero_seq_SERO_opt();
opt.bMax = 1.4e9;
opt.nShot = 1000;
opt.TE = 70e-3;
opt.Nx = 120;
opt.Ny = 120;
opt.FOV = opt.Nx * dz;

% Create sampling scheme
[posi, b, sthi] = sero_sampling_sero_gen2(opt.nShot, ar, nPos);

exel = sero_seq_posandb2excutionlist(posi, sthi, b, dz);

% Create and write the seq structure
seq = sero_seq_SERO_create_v4(exel, opt, lims);
write(seq, fn_seq);

% Calculate xps
xps = sero_seq_SERO_seq2xps(seq);
mdm_xps_save(xps, sero_seq_fn_seq2xps(fn_seq));


%% PLOT
figure(1); clf
subplot(2,2,1)
sero_seq_plotShots(seq, 1)

subplot(2,2,3)
hw = safe_hw_prisma_xr_sh05;
blInd = sero_seq_shotToBlockInd(seq, [1 5]);
sero_seq_plot_pns(seq, blInd, hw)

subplot(1,2,2)
imagesc(xps.tr)


figure(2); clf
fn_asc = "";
ascData = mr.Siemens.readasc(fn_asc);
sero_seq_plotMechSpectrum(seq, [1 100], 'shot', [], [], ascData);
