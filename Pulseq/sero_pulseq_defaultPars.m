function par = sero_pulseq_defaultPars()
% function par = sero_pulseq_defaultPars()

par.tShot = 0.15;          % step time [s]
par.nShot = 1000;       % number of shots/samples
par.nPos = 50;          % number of positions/slices
par.ar = 4;             % aspect ratio
par.dz = 1.5e-3;         % step length in z [m]
par.fov = 180e-3;       % in-plane FOV [m]
par.Nx = 120;           % Image matrix in x [1]
par.Ny = par.Nx;    % Image matrix in y [1]
par.slth = 6e-3;           % slice width [m]
par.sliceGap = 0;          %slice gap factor

par.maxbValue = 2000;    % b-value [s/mm^2]
par.TE = 80e-3;         % Echo time [s]
par.spoilFac = 20;      % Spoiler strength compared to edge of k-space
par.pe_enable = 1;      % a flag to disable phase encoding
par.ro_os = 1;          % frequency oversampling factor (in contrast to the product sequence we don't really need it)
par.roTime = 6.0e-4;    % Time spent reading each line (regulates bandwidth) [s]
par.PFFactor = 0.5;     % Fraction of lines to read before TE. 1: full sampling 0: start with ky=0
par.do_fatsat = 0;      % Boolean if do fat saturation

par.rf_dur    = 3e-3; % [s]
par.rf_apod   = 0.5;    % decay to zero [0-1]
par.rf_tbwp   = 4;      % edge sharpness