function seq = sero_seq_SERO_create_v4(E, opt, lims)
% function seq = sero_seq_SERO_create_v4(E, opt, lims)
%
% E is the second generation SERO sampling scheme
% opt contains all pulse sequence settings
% lims are the hardware limitations

if nargin < 2
    opt = sero_seq_SERO_opt();
end

if nargin < 3
    lims = sero_seq_sysLims_prisma();
end

%% Unpack parameters from options structure
tShot     = opt.tShot;           % step time/t_shot [s]
fov       = opt.FOV;             % in-plane FOV [m]
Nx        = opt.Nx;              % Image matrix in x [1]
Ny        = opt.Ny;              % Image matrix in y [1]

maxbValue = opt.bMax*1e-6;       % b-value [s/mm^2]
TE        = opt.TE;              % Echo time [s]
spoilFac  = opt.SpoilerFactor;   % Spoiler strength compared to edge of k-space

ro_os     = opt.freqOverSamp;    % frequency oversampling factor (in contrast to the product sequence we don't really need it)
roTime    = opt.ReadoutTime;     % Time spent reading each line (regulates bandwidth) [s]
PFFactor  = opt.PFFactor;        % Fraction of lines to read before TE. 1: full sampling 0: start with ky=0
do_fatsat = opt.do_fatsat;       % Use fat saturation or not

rf_dur    = opt.rf_dur;
rf_apod   = opt.rf_apod;
rf_tbwp   = opt.rf_tbwp;

% Set rf durations
tRFex   = rf_dur;
tRFref  = rf_dur;
tRFsat  = 12e-3;
sat_ppm = -3.45;


% Set derived paramters
nShot    = size(E,1);
slicePos = E(:,1);
slThic   = E(:,2);
fov_z    = range(slicePos+slThic);  % FOV in slice direction [m], difference between last and first slice + one half slicewidth at the top and one half at the bottom
[uSlThic, ~, ind_st]  = unique(slThic);

if numel(uSlThic)>1
    error('this is not supported yet');
end

% Vestigial parameters
do_phaseEnc = true;               % a flag to disable phase encoding
rf_ph       = 0;


%% Create sequence object and populate it
seq = mr.Sequence(lims.img); % Create a new sequence object

% Create fat-sat pulse
sat_freq = sat_ppm*1e-6*lims.img.B0*lims.img.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims.spl,'Duration',tRFsat,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
gz_fs = mr.makeTrapezoid('z',lims.spl,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create 90 degree slice selection pulse and gradient (including rewinder!)
% and create 180 degree slice refocusing pulse and gradients
for i = 1:numel(uSlThic)

    slth                            = uSlThic(i);
    
    [rf90(i), gz90(i), gzReph(i)]   = mr.makeSincPulse(pi/2,'system',lims.img,'Duration',tRFex,...
                                                        'SliceThickness',slth,'PhaseOffset',pi/2,...
                                                        'apodization',rf_apod,'timeBwProduct',rf_tbwp);

    [~, gzr_t, gzr_a]               = mr.makeExtendedTrapezoidArea('z',gz90(i).amplitude,0, ...
                                                                    gzReph(i).area+0.5*gz90(i).amplitude*gz90(i).fallTime,...
                                                                    lims.img);
    
    gz90n(i)                        = mr.makeExtendedTrapezoid('z','system',lims.img,'times',...
                                                                [0 gz90(i).riseTime gz90(i).riseTime+gz90(i).flatTime+gzr_t]+gz90(i).delay, ...
                                                                'amplitudes', [0 gz90(i).amplitude gzr_a]);

    [rf180(i), gz180(i)]            = mr.makeSincPulse(pi,'system',lims.img,'Duration',tRFref,...
                                                        'SliceThickness',slth,'apodization',rf_apod,...
                                                        'timeBwProduct',rf_tbwp,'use','refocusing');
end


% define the output trigger to play out with every slice excitatuion
trig = mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'


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
gxAmpReset = gx.amplitude;
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


% split the blip into two halves and produce a combined synthetic gradient
gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims.img);
[gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims.img);


% pe_enable support
gy_blipup.waveform=gy_blipup.waveform*do_phaseEnc;
gy_blipdown.waveform=gy_blipdown.waveform*do_phaseEnc;
gy_blipdownup.waveform=gy_blipdownup.waveform*do_phaseEnc;


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
gyPre.amplitude=gyPre.amplitude*do_phaseEnc;


% Calculate delay times
durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx);
rfCenterInclDelay = rf90(1).delay + mr.calcRfCenter(rf90(1));
rf180centerInclDelay = rf180(1).delay + mr.calcRfCenter(rf180(1));
delayTE1   = ceil((TE/2 - mr.calcDuration(rf90(1),gz90n) + rfCenterInclDelay - rf180centerInclDelay)/lims.img.gradRasterTime)*lims.img.gradRasterTime;
delayTE2tmp = ceil((TE/2 - mr.calcDuration(rf180(1),gz180(1)) + rf180centerInclDelay - durationToCenter)/lims.img.gradRasterTime)*lims.img.gradRasterTime;
assert(delayTE1>=0);
gxPre.delay = 0;
gyPre.delay = 0;
delayTE2 = delayTE2tmp-mr.calcDuration(gxPre,gyPre);
[gxPre,gyPre] = mr.align('right',gxPre,'left',gyPre);
assert(delayTE2>=0);

% diffusion weithting calculation
% delayTE2 is our window for small_delta
% delayTE1+delayTE2-delayTE2 is our big delta
% we anticipate that we will use the maximum gradient amplitude, so we need
% to shorten delayTE2 by gmax/max_sr to accommodate the ramp down
small_delta = delayTE2-ceil(lims.diff.maxGrad/lims.diff.maxSlew/lims.diff.gradRasterTime)*lims.diff.gradRasterTime;
big_delta = delayTE1+mr.calcDuration(rf180(1),gz180(1));
% we define bFactCalc function below to eventually calculate time-optimal
% gradients. for now we just abuse it with g=1 to give us the coefficient

g = sqrt(maxbValue*1e6/sero_seq_calculateBval(1,small_delta,big_delta));
gr = ceil(g/lims.diff.maxSlew/lims.diff.gradRasterTime)*lims.diff.gradRasterTime;
gDiff = mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims.diff);
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
for e_ind = 1:nShot

    % select current objects
    this_ind   = ind_st(e_ind);
    this_rf90  = rf90(this_ind);
    this_gz90  = gz90(this_ind);
    this_gz90n = gz90n(this_ind);
    this_rf180 = rf180(this_ind);
    this_gz180 = gz180(this_ind);

    % Save start time for the shot
    t_start = sum(seq.blockDurations);

    % Fat saturation
    if do_fatsat
        seq.addBlock(rf_fs,gz_fs);
    end

    % Excitation
    this_rf90.freqOffset  = this_gz90.amplitude*slicePos(e_ind,1);
    this_rf90.phaseOffset = pi/2 - 2*pi*this_rf90.freqOffset * mr.calcRfCenter(this_rf90) + rf_ph;
    seq.addBlock(this_rf90, this_gz90n, trig);

    % Diffusion prep
    for i = 1:3
        gDiff_s.(chn{i}).amplitude = gDiff.amplitude * E(e_ind, i+2);
    end

    % Diffusion before refocus
    seq.addBlock(mr.makeDelay(delayTE1), gDiff_s.x, gDiff_s.y, gDiff_s.z);

    % Refocus
    this_rf180.freqOffset  = this_gz180.amplitude*slicePos(e_ind,1);
    this_rf180.phaseOffset = -2*pi*rf180(1).freqOffset*mr.calcRfCenter(this_rf180) + rf_ph;
    seq.addBlock(this_rf180, this_gz180);

    % Diffusion after refocus
    seq.addBlock(mr.makeDelay(delayTE2), gDiff_s.x, gDiff_s.y, gDiff_s.z);
    seq.addBlock(mr.makeDelay(1e-3));

    % EPI
    adc.phaseOffset = -2*pi*rf180(1).freqOffset*mr.calcRfCenter(rf180(1)) + rf_ph;
    seq.addBlock(gxPre,gyPre);

    gx.amplitude = gxAmpReset;

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

    % Spoiling
    for i = 1:3
        gSpl_s.(chn{i}).amplitude  = gSpl.amplitude  * E(e_ind, i+5);
    end
    seq.addBlock(mr.makeDelay(1e-3), gSpl_s.x, gSpl_s.y, gSpl_s.z);

    % Add fill time
    t_end = sum(seq.blockDurations);
    t_fill = sero_seq_roundToRaster(tShot-(t_end-t_start), lims.img.gradRasterTime);
    seq.addBlock(mr.makeDelay(t_fill));
end


% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov fov_z]);
seq.setDefinition('Name', 'epi-diff');



