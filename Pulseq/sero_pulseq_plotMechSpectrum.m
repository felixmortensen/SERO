function sero_pulseq_plotMechSpectrum(seq, n, mode, fInterval, tRes, ascData)
% function sero_pulseq_plotMechSpectrum(seq, n, mode, fInterval, tRes, ascData)
%
% Function to plot the mechanical spectrum from seq based on the Matlab
% function "pspectrum".
% n is 1x2 index vector that defines start and end index.
% mode defines if n refers to shots or block number.
% finterval is the frequency interval to include in the analysis
% tRes is the temporal resolution to include in the analysis.
% ascData is the hardware specification specific to Siemens that can be
% compiled by calling ascData=mr.Siemens.readasc(fn_ASC). This structure
% contains frequency intervals that define where the system has mechanical
% resonances, i.e., should be avoided.

if nargin < 2 || isempty(n)
    n = 0;
end

if nargin < 3 || isempty(mode)
    mode = 'block';
end

if nargin < 4 || isempty(fInterval)
    fInterval = [0 5000];
end

if nargin < 5 || isempty(tRes)
    tRes = .1;
end

if nargin < 6
    ascData = [];
end


switch mode
    case 'block'
        blockInd = n;
    case 'shot'
        blockInd = sero_pulseq_shotToBlockInd(seq, n);
end


dt  = seq.gradRasterTime;
gwf = sero_pulseq_getBlockWaveforms(seq, blockInd);
swf = diff([gwf; 0 0 0],1,1)/dt;


% Gradient waveform power spectrum
[gps, gf] = pspectrum(gwf, 1/dt, 'power', 'FrequencyLimits', fInterval);
for i = 1:3
    [gps2d(:,:,i), gf2d, gt2d] = pspectrum(gwf(:,i), 1/dt, 'spectrogram', 'FrequencyLimits', fInterval, 'TimeResolution', tRes);
end

% d/dt of gwf (slew rate waveform) power spectrum
[sps, sf] = pspectrum(swf, 1/dt, 'power', 'FrequencyLimits', fInterval);
for i = 1:3
    [sps2d(:,:,i), sf2d, st2d] = pspectrum(swf(:,i), 1/dt, 'spectrogram', 'FrequencyLimits', fInterval, 'TimeResolution', tRes);
end


% plot
subplot(2,2,1)
plot(gf/1e3, gps/max(gps(:)))
xlabel('Frequency [kHz]')
ylabel('Power [a.u.]')
title('Average power of g(t)')

ylim([0 1.1])
xlim(fInterval/1e3)
plot_dangerzone(ascData)

subplot(2,2,2)
imagesc(gt2d, gf2d/1e3,  log(sum(gps2d,3)))
xlabel('Time [s]')
ylabel('Frequency [kHz]')
title('Log(Power(g(t))) vs time')

subplot(2,2,3)
plot(sf/1e3, sps/max(sps(:)))
xlabel('Frequency [kHz]')
ylabel('Power [a.u.]')
title('Average power of dg(t)/dt')

ylim([0 1.1])
xlim(fInterval/1e3)
plot_dangerzone(ascData)

subplot(2,2,4)
imagesc(st2d, sf2d/1e3,  log(sum(sps2d,3)))
xlabel('Time [s]')
ylabel('Frequency [kHz]')
title('Log(Power(dg(t)/dt)) vs time')
end


function plot_dangerzone(ascData)

if isempty(ascData)
    return
end

s      = ascData.asGPAParameters(1).sGCParameters;
mid    = s.aflAcousticResonanceFrequency/1e3;
bwz    = s.aflAcousticResonanceBandwidth/1e3;
nZones = numel(mid);

hold on

amp = get(gca, 'YLim');

for i = 1:nZones

    if all([mid(i) bwz(i)]==0)
        continue
    end

    x = mid(i) + [-1 1 1 -1]*bwz(i)/2;
    y = [[1 1]*amp(1) [1 1]*amp(2)];

    patch(x, y, [.8 .2 .2], 'facealpha', 0.3, 'linestyle', 'none');
    plot([1 1]*mid(i), amp, 'k--')
end

end