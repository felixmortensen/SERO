function sero_seq_writeSERO_v2(sampling, fileName, varargin)
% Return a diffusion weighted pulse sequence based on SERO and writeEpiDiffusionRS
% SERO assumes that the user suppies an event schedule that contains slice
% positioning, diffusion encoding vector, and a spoiler gradient vector.
%modifcation etc

%Requires first generation SERO sampling scheme (sero_recon_sampling_v3)

% Check if the output filename is provided, and that it is entered
% correctly
if nargin < 2 || ~ischar(fileName) || ~endsWith(fileName, '.seq')
    error('Invalid output filename. It must be a string ending with ".seq".');
else
    o_name = fileName;
end

%% Set basic imaging and diffusion parameters by unpacking params

% Get the parameters
params = pulseq_InitialiseParameters_v2(varargin{:});

% Unpack parameters for position matrix
ts        = params.ts;            % step time [s]
nShot     = params.nShot;       % number of shots/samples
nPos      = params.nPos;        % number of positions/slices
ar        = params.ar;
dz        = params.dz;            %step length [m]

%choose sampling scheme
if strcmp(sampling, 'sero')

    [~, tr, b] = srp_generateSERO(ts, nShot, ar, nPos);
    E = pulseq_trb2pos(tr, b, dz);

elseif strcmp(sampling, 'direct')

    [~, tr, b] = srp_sampling_conventional(ts, nShot, 1, nPos, 1);
    E = pulseq_trb2pos(tr, b, dz);

elseif strcmp(sampling, 'slice')

    [~, tr, b] = srp_sampling_conventional_slice_shift(ts, nShot, ar, nPos, 1);
    E = pulseq_trb2pos(tr, b, dz);

else
    error('Invalid sampling scheme. Choose either sero, direct or slice!')
end

%unpack remaining imaging parameters
fov       = params.FOV;          % in-plane FOV [m]
Nx        = params.Nx;             % Image matrix in x [1]
Ny        = Nx;              % Image matrix in y [1]
slth      = params.SliceWidth;   %slice width [m]
fov_z     = range(E(:,1))+slth;            % FOV in slice direction [m], difference between last and first slice + one half slicewidth at the top and one half at the bottom
hirez     = params.HighResolutionGrid;            % Intended high resolution grid
t_shot    = params.t_shot;

maxbValue = params.Max_b_value;             % b-value [s/mm^2]
TE        = params.Echo_time;           % Echo time [s]
spoilFac  = params.SpoilerFactor;             % Spoiler strength compared to edge of k-space

pe_enable = params.PhaseEncoding_enable;               % a flag to disable phase encoding
ro_os     = params.Frequency_oversampling_factor;               % frequency oversampling factor (in contrast to the product sequence we don't really need it)
roTime    = params.ReadoutTime;          % Time spent reading each line (regulates bandwidth) [s]
PFFactor  = params.PFFactor;               % Fraction of lines to read before TE. 1: full sampling 0: start with ky=0

rf_faFactor = params.rf_faFactor;
do_fatsat = params.do_fatsat;
do_rfSpoil = params.do_rfSpoil;

%% Tweak the slice profile
% The slice profile in the default settings is very crappy. We can improve
% this at the cost of some echo time (lonmger RF pulses).

switch 1
    case 1 % Original
        rf_dur = 3e-3; % [s]
        apod   = 0.5;    % decay to zero [0-1]
        tbwp   = 4;      % edge sharpness
    case 2 % better
        rf_dur = 6e-3;
        apod   = .35;
        tbwp   = 8;
    case 3 % better yet
        rf_dur = 8e-3;
        apod   = .35;
        tbwp   = 15;
end


%% Set system limits
lims = pulseq_sysLims_prisma();

%% Create sequence object and populate it
seq = mr.Sequence(lims.img); % Create a new sequence object

% Set rf durations
tRFex  = rf_dur;
tRFref = rf_dur;
tRFsat = 12e-3;

% Create fat-sat pulse
sat_ppm = -3.45;
sat_freq = sat_ppm*1e-6*lims.img.B0*lims.img.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims.img,'Duration',tRFsat,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
gz_fs = mr.makeTrapezoid('z',lims.img,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm


% Create 90 degree slice selection pulse and gradient -- including rewinder!
[rf90 , gz90, gzReph] = mr.makeSincPulse(pi/2,'system',lims.img,'Duration',tRFex,...
    'SliceThickness',slth,'PhaseOffset',pi/2,'apodization',apod,'timeBwProduct',tbwp);

[~, gzr_t, gzr_a] = mr.makeExtendedTrapezoidArea('z',gz90.amplitude,0, gzReph.area+0.5*gz90.amplitude*gz90.fallTime,lims.img);
gz90n=mr.makeExtendedTrapezoid('z','system',lims.img,'times',[0 gz90.riseTime gz90.riseTime+gz90.flatTime+gzr_t]+gz90.delay, ...
    'amplitudes', [0 gz90.amplitude gzr_a]);

% Create 180 degree slice refocusing pulse and gradients
[rf180, gz180] = mr.makeSincPulse(pi*rf_faFactor,'system',lims.img,'Duration',tRFref,...
    'SliceThickness',slth,'apodization',apod,'timeBwProduct',tbwp,'use','refocusing');


% define the output trigger to play out with every slice excitatuion
trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

% Define other gradients and ADC events
deltak = 1/fov;
kWidth = Nx*deltak;

% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(deltak/lims.img.maxSlew)/10e-6/2)*10e-6*2;

% we round-up the duration to 2x the gradient raster time
% the split code below fails if this really makes a trpezoid instead of a triangle...
gy = mr.makeTrapezoid('y',lims.img,'Area',-deltak,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area
extra_area = blip_dur/2*blip_dur/2*lims.img.maxSlew; % check unit!;
gx = mr.makeTrapezoid('x',lims.img,'Area',kWidth+extra_area,'duration',roTime+blip_dur);
actual_area = gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude = gx.amplitude/actual_area*kWidth;
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea = gx.amplitude*gx.flatTime;

% Prepare spoilers
spoilArea = spoilFac * gx.area;
gSpl = mr.makeTrapezoid('x', 'Area', spoilArea, 'system', lims.spl);

% calculate ADC
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively.
adcDwellNyquist=deltak/gx.amplitude/ro_os;
% round-down dwell time to 100 ns
adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
adcSamples=floor(roTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% FOV positioning requires alignment to grad. raster... -> TODO

% split the blip into two halves and produce a combined synthetic gradient
gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims.img);
[gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims.img);

% pe_enable support
gy_blipup.waveform=gy_blipup.waveform*pe_enable;
gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

% phase encoding and partial Fourier
Ny_pre=round(PFFactor*Ny/2-1); % PE steps prior to ky=0, excluding the central line
Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
Ny_meas=Ny_pre+Ny_post;

% Pre-phasing gradients
gxPre = mr.makeTrapezoid('x',lims.img,'Area',-gx.area/2);
gyPre = mr.makeTrapezoid('y',lims.img,'Area',Ny_pre*deltak);
[gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);

% EPI revinders (not sure if these help or hurt)
gxPost = gxPre; gxPost.amplitude = -gxPre.amplitude;
gyPost = gyPre; gyPost.amplitude = -(Ny_post-1) * gy_blipdown.area / (gyPost.riseTime+gyPost.flatArea) * 2; % FSz: I dont quite understand where the 2 comes from

% relax the PE prepahser to reduce stimulation
gyPre = mr.makeTrapezoid('y',lims.img,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
gyPre.amplitude=gyPre.amplitude*pe_enable;

% Calculate delay times
durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx);
rfCenterInclDelay=rf90.delay + mr.calcRfCenter(rf90);
rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
delayTE1   =ceil((TE/2 - mr.calcDuration(rf90,gz90n) + rfCenterInclDelay - rf180centerInclDelay)/lims.img.gradRasterTime)*lims.img.gradRasterTime;
delayTE2tmp=ceil((TE/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter)/lims.img.gradRasterTime)*lims.img.gradRasterTime;
assert(delayTE1>=0);
gxPre.delay=0;
gyPre.delay=0;
delayTE2=delayTE2tmp-mr.calcDuration(gxPre,gyPre);
[gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
assert(delayTE2>=0);

% diffusion weithting calculation
% delayTE2 is our window for small_delta
% delayTE1+delayTE2-delayTE2 is our big delta
% we anticipate that we will use the maximum gradient amplitude, so we need
% to shorten delayTE2 by gmax/max_sr to accommodate the ramp down
small_delta=delayTE2-ceil(lims.diff.maxGrad/lims.diff.maxSlew/lims.diff.gradRasterTime)*lims.diff.gradRasterTime;
big_delta=delayTE1+mr.calcDuration(rf180,gz180);
% we define bFactCalc function below to eventually calculate time-optimal
% gradients. for now we just abuse it with g=1 to give us the coefficient

g=sqrt(maxbValue*1e6/pulseq_calculateBval(1,small_delta,big_delta));
gr=ceil(g/lims.diff.maxSlew/lims.diff.gradRasterTime)*lims.diff.gradRasterTime;
gDiff=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims.diff);
assert(mr.calcDuration(gDiff)<=delayTE1);
assert(mr.calcDuration(gDiff)<=delayTE2);

% Make 3-axis gradient objects for convenience
chn = {'x', 'y', 'z'};
for i = 1:3
    gDiff_s.(chn{i}) = gDiff;
    gDiff_s.(chn{i}).channel = chn{i};

    gSpl_s.(chn{i}) = gSpl;
    gSpl_s.(chn{i}).channel = chn{i};
end


%% Compose sequence event blocks
for e_ind = 1:size(E,1)

    % Save start time for the shot
    t_start = sum(seq.blockDurations);

    % Fat saturation
    if do_fatsat
        seq.addBlock(rf_fs,gz_fs);
    end

    % Set RF spoiling phase offset
    if do_rfSpoil
        rf_ph = mod(117*(e_ind^2 + e_ind + 2), 360)*pi/180;
    else
        rf_ph = 0;
    end

    % Excitation
    slicePos = E(e_ind,1);
    rf90.freqOffset=gz90.amplitude*slicePos;
    rf90.phaseOffset=pi/2-2*pi*rf90.freqOffset*mr.calcRfCenter(rf90) + rf_ph;
    seq.addBlock(rf90,gz90n,trig);

    % Diffusion prep
    for i = 1:3
        gDiff_s.(chn{i}).amplitude = gDiff.amplitude * E(e_ind, i+1);
    end

    % Diffusion before refocus
    seq.addBlock(mr.makeDelay(delayTE1), gDiff_s.x, gDiff_s.y, gDiff_s.z);

    % Refocus
    rf180.freqOffset=gz180.amplitude*slicePos;
    rf180.phaseOffset=-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180) + rf_ph;
    seq.addBlock(rf180,gz180)

    % Diffusion after refocus
    seq.addBlock(mr.makeDelay(delayTE2), gDiff_s.x, gDiff_s.y, gDiff_s.z);

    % EPI
    adc.phaseOffset = -2*pi*rf180.freqOffset*mr.calcRfCenter(rf180) + rf_ph;
    seq.addBlock(gxPre,gyPre);
    for i=1:Ny_meas
        if i==1
            seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
        elseif i==Ny_meas
            seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
        else
            seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
        end
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end

    % Rewinder after EPI
    %     seq.addBlock(gxPost,gyPost);

    % Spoiling
    for i = 1:3
        gSpl_s.(chn{i}).amplitude  = gSpl.amplitude  * E(e_ind, i+4);
    end
    seq.addBlock(mr.makeDelay(1e-3), gSpl_s.x, gSpl_s.y, gSpl_s.z);

    % Add fill time
    t_end = sum(seq.blockDurations);
    t_fill = pulseq_roundToRaster(t_shot-(t_end-t_start), lims.img.gradRasterTime);
    seq.addBlock(mr.makeDelay(t_fill));
end

%%
if nargin < 1
    seq.plot();
    return; %break the code there if the statement is fullfilled
end
%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov fov_z]);
seq.setDefinition('Name', 'epi-diff');

%Extract the file parts from o_name
[path, baseFilename, extension] = fileparts(o_name);

%Construct the new filename with the desired suffix
newFilename = fullfile(path, [baseFilename '_TR_b.mat']);

%Save TR and b in new .mat file with the filename ending on _TR_b.mat
%instead of .seq
save(newFilename, 'tr', 'b');

write(seq, o_name);

%Specify the name of the zip file
zipFilename = fullfile(path, [baseFilename '.zip']);

%Create a zip file containing both the MATLAB data file and pulse sequence file
zip(zipFilename, {newFilename, o_name});
disp(['New sequence file available in ' zipFilename])

% Delete the original files from directory and keep them only in zip file
delete(newFilename, o_name);

end

