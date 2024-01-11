%This script is for internal use only. It may be removed in the future.

% Helper script creates data for help examples for
% phased.RangeDopplerResponse
%
% Radar system is based on featured example:
%    "Waveform Design to Improve Performance of an Existing Radar System"
%
% Creates 2 data sets, one using rectangular waveform and the other using a
% FMCW waveform.
%
% The designed system uses coherent integration instead of noncoherent
% integration (which was used in the featured example). This is because
% Doppler processing is a form of coherent integration.

%   Copyright 2016 The MathWorks, Inc.

rng default;

%% Design specifications
fc = 10e9;           % Center frequency
prop_speed = physconst('LightSpeed');   % Propagation speed
pd = 0.9;            % Probability of detection
pfa = 1e-6;          % Probability of false alarm
max_range = 5000;    % Maximum unambiguous range
min_rcs = 1;         % Required target radar cross section
range_res = 50;      % Required range resolution

%% Design systems and generate radar data
[rx,mf,fs,noisepower] = ...
    getRectdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res);
rectdata = struct('fs',fs,'propspeed',prop_speed,'fc',fc, ...
    'rxdata',rx,'mfcoeffs',mf,'noisepower',noisepower);

[rx,xref,slope,fs,noisepower] = ...
    getFMCWdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res);
fmcwdata = struct('fs',fs,'propspeed',prop_speed,'fc',fc, ...
    'rxdata',rx,'refsig',xref,'sweepslope',slope,'noisepower',noisepower);

save('RangeDopplerResponseExampleData.mat','rectdata','fmcwdata');

%% LFM Example
clear; figure(1); clf;

% Example:
%   Calculate range-Doppler response of a pulsed rectangular radar signal
%   using matched filter approach. The signal includes three target
%   returns. Two of them are around 2000 m away and the third is around
%   3500 m away. In addition, two of them are stationary relative to the
%   radar while the other one is moving away from radar at around 100 m/s.
    
% Load example data
load('RangeDopplerResponseExampleData','rectdata');
fs = rectdata.fs;
propspeed = rectdata.propspeed;
fc = rectdata.fc;
rxdata = rectdata.rxdata;
mfcoeffs = rectdata.mfcoeffs;
noisepower = rectdata.noisepower;

% Create range-Doppler response for matched filter processing. Interpolate
% to 1024 Doppler bins
rngdopresp = phased.RangeDopplerResponse(...
    'DopplerFFTLengthSource','Property',...
    'DopplerFFTLength',1024,...
    'DopplerOutput','Speed',...
    'OperatingFrequency',fc,...
    'SampleRate',fs,...
    'PropagationSpeed',propspeed);

% Perform range-Doppler processing
[resp,rng_grid,dop_grid] = rngdopresp(rxdata,mfcoeffs);

% Calculate noise power after range-Doppler processing
mfgain = mfcoeffs'*mfcoeffs;
noisepower_proc = mfgain*noisepower;

imagesc(dop_grid,rng_grid,pow2db(abs(resp).^2/noisepower_proc));
xlabel('Speed (m/s)'); ylabel('Range (m)'); title('Range Doppler Map');
set(get(colorbar,'YLabel'),'String','SNR (dB)'); caxis([0 60]);

%% FMCW Example
clear; figure(2); clf;

% Example:
%   Plot the range-Doppler response of an FMCW signal. The
%   signal is not dechirped. The signal contains the return
%   from one target which is around 2200 m away and has a
%   normalized Doppler frequency of about -0.36 relative to the
%   radar.

% Load example data
load('RangeDopplerResponseExampleData','fmcwdata');
fs = fmcwdata.fs;
propspeed = fmcwdata.propspeed;
rxdata = fmcwdata.rxdata;
refsig = fmcwdata.refsig;
sweepslope = fmcwdata.sweepslope;

% Create range-Doppler response for dechirp and FFT processing
rngdopresp = phased.RangeDopplerResponse(...
    'RangeMethod','FFT',...
    'SweepSlope',sweepslope,...
    'DechirpInput',true,...
    'SampleRate',fs,...
    'PropagationSpeed',propspeed);

% Plot range-Doppler response of processed data
plotResponse(rngdopresp,rxdata,refsig, ...
    'Unit','db','NormalizeDoppler',true);

%% Helpers
function [rx_data,mfcoeffs,fs,noisepower] = getRectdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res)

%% Waveform and waveform processing
pulse_bw = prop_speed/(2*range_res);    % Pulse bandwidth
pulse_width = 1/pulse_bw;               % Pulse width
prf = prop_speed/(2*max_range);         % Pulse repetition frequency
fs = 2*pulse_bw;                        % Sampling rate
waveform = phased.RectangularWaveform(...
    'PulseWidth',pulse_width,...
    'PRF',prf,...
    'SampleRate',fs);

% Matched filter coefficients
mfcoeffs = getMatchedFilter(waveform);

% Coherent pulse integration
num_pulse_int = 128;
int_gain_dB = pow2db(num_pulse_int); % Coherent integration gain

%% Receiver
receiver = phased.ReceiverPreamp(...
    'Gain',20,...
    'NoiseFigure',0,...
    'SampleRate',fs,...
    'EnableInputPort',true,...
    'SeedSource','Property',...
    'Seed',2007);

%% Transmitter
tx_gain = 20;

snr_min = albersheim(pd, pfa); % One pulse
snr_min = snr_min-int_gain_dB; % Improvement after coherent pulse integration

lambda = prop_speed/fc;
peak_power = radareqpow(lambda,max_range,snr_min,pulse_width,...
    'RCS',min_rcs,'Gain',tx_gain);

transmitter = phased.Transmitter(...
    'Gain',tx_gain,...
    'PeakPower',peak_power,...
    'InUseOutputPort',true);

%% Targets
tgtpos = [[2000;0;0],[2000;0;0],[3500;0;0]];
tgtvel = [[100;0;0],[0;0;0],[0;0;0]];
tgtrcs = [1.6 2.2 1.05];

%% Generate data
rx_data = simulateSystem(fc,waveform,transmitter,receiver,num_pulse_int, ...
    tgtpos,tgtvel,tgtrcs);

noise_bw = fs/2;
noisepower = noisepow(noise_bw,...
    receiver.NoiseFigure,receiver.ReferenceTemperature);
end

function [rx_data,xref,sweepslope,fs,noisepower] = getFMCWdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res)

%% Waveform and waveform processing
sweep_bw = range2bw(range_res,prop_speed);
sweep_time = 5.5*range2time(max_range,prop_speed);
sweepslope = sweep_bw/sweep_time;
fs = sweep_bw;
waveform = phased.FMCWWaveform(...
    'SweepBandwidth',sweep_bw,...
    'SweepTime',sweep_time,...
    'SampleRate',fs);

% Reference signal for dechirp processing
xref = waveform();

% Coherent pulse integration
num_pulse_int = 128;
int_gain_dB = pow2db(num_pulse_int); % Coherent integration gain

%% Receiver
receiver = phased.ReceiverPreamp(...
    'Gain',20,...
    'NoiseFigure',0,...
    'SampleRate',fs,...
    'EnableInputPort',false,...
    'SeedSource','Property',...
    'Seed',2007);

%% Transmitter
tx_gain = 20;

snr_min = albersheim(pd, pfa); % One pulse
snr_min = snr_min-int_gain_dB; % Improvement after coherent pulse integration

lambda = prop_speed/fc;
peak_power = radareqpow(lambda,max_range,snr_min,sweep_time,...
    'RCS',min_rcs,'Gain',tx_gain);

transmitter = phased.Transmitter(...
    'Gain',tx_gain,...
    'PeakPower',peak_power,...
    'InUseOutputPort',true);

%% Targets
prf = 1/sweep_time;
fdnorm = -0.36;
v = -lambda/2*fdnorm*prf;
tgtpos = [2200;0;0];
tgtvel = [v;0;0];
tgtrcs = 1.6;

%% Generate data
rx_data = simulateSystem(fc,waveform,transmitter,receiver,num_pulse_int,tgtpos,tgtvel,tgtrcs);

noise_bw = fs/2;
noisepower = noisepow(noise_bw,...
    receiver.NoiseFigure,receiver.ReferenceTemperature);
end

function rxdata = simulateSystem(fc,waveform,transmitter,receiver,num_pulse_int,tgtpos,tgtvel,tgtrcs)

fs = waveform.SampleRate;

%% Transmit and receive antennas
antenna = phased.IsotropicAntennaElement(...
    'FrequencyRange',[5e9 15e9]);

radiator = phased.Radiator(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

collector = phased.Collector(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

%% Radar platform
sensormotion = phased.Platform(...
    'InitialPosition',[0; 0; 0],...
    'Velocity',[0; 0; 0]);

%% Targets
tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel);
target = phased.RadarTarget('MeanRCS',tgtrcs,'OperatingFrequency',fc);

%% Environment
channel = phased.FreeSpace(...
    'SampleRate',fs,...
    'TwoWayPropagation',true,...
    'OperatingFrequency',fc);

% Update sensor and target positions
if isa(waveform,'phased.RectangularWaveform')
    time_step = 1/waveform.PRF;
elseif isa(waveform,'phased.FMCWWaveform')
    time_step = waveform.SweepTime;
end

%% Simulate system
% Build data cube
pulse = waveform();
num_fast_time = length(pulse);
rxdata = zeros(num_fast_time,num_pulse_int);
for m = 1:num_pulse_int
    [sensorpos,sensorvel] = sensormotion(time_step);
    [tgtpos,tgtvel] = tgtmotion(time_step);
    
    % Calculate the target angles as seen by the sensor
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);
    
    % Simulate propagation of pulse in direction of targets
    pulse = waveform();
    [txsig,txstatus] = transmitter(pulse);
    txsig = radiator(txsig,tgtang);
    txsig = channel(txsig,sensorpos,tgtpos,sensorvel,tgtvel);
    
    % Reflect pulse off of targets
    tgtsig = target(txsig);
    
    % Receive target returns at sensor
    rxsig = collector(tgtsig,tgtang);
    if receiver.EnableInputPort
        rxpulse = receiver(rxsig,~(txstatus>0));
    else
        rxpulse = receiver(rxsig);
    end
    rxdata(:,m) = rxpulse;
end
end
