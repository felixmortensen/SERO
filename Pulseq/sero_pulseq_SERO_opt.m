function opt = sero_pulseq_SERO_opt(opt)
% function opt = sero_pulseq_SERO_opt(opt)
%
% 1. Call this functions without input arguments to get default structure -
% e.g. opt = sero_pulseq_SERO_opt()
% 2. Change the options as you like
% e.g. opt.TE = 100e-3;
% 3. Call this function again with the options as input to updated derived
% parameters and perform a check
% e.g. opt = pulseq_SERO_opt(opt)

% Defined pars
if nargin < 1
    opt.sliceGap = 0;
    opt.tShot = 0.15;
    opt.FOV = 0.220;
    opt.Nx = 110;
    opt.Ny = 110;
    opt.bMax = 1.4e9;
    opt.TE = 80e-3;
    opt.SpoilerFactor = 30;
    opt.freqOverSamp = 1;
    opt.ReadoutTime = 0.6e-3;
    opt.PFFactor = 0.5;
    opt.do_fatsat = 1;

    opt.rf_dur = 6e-3;
    opt.rf_apod = .3;
    opt.rf_tbwp = 10;
end

% Derived pars
% none so far

% Check consistency
% This is better solved with an object, but good enough for now as it
% mainly protects agains adding new fields that are not used.
% Add more checks for consistency
if numel(fieldnames(opt)) ~= 15
    error('opt is inconsistent!')
end

% For example
if opt.bMax < 0
    error('b-value must be 0 or positive number!')
end