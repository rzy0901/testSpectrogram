%This script is for internal use only. It may be removed in the future.

% Helper script creates data for help examples for phased.RangeResponse
%
% Radar system is based on featured example:
%    "Waveform Design to Improve Performance of an Existing Radar System"
%
% Creates 2 data sets, one using LFM waveform and the other using a FMCW
% waveform.
%
% Only one pulse (sweep) of data is created for each waveform to avoid
% needing pulse integration in help example.

%   Copyright 2016 The MathWorks, Inc.

rng default;

%% Design specifications
fc = 10e9;           % Center frequency
prop_speed = physconst('LightSpeed');   % Propagation speed
pd = 0.9;            % Probability of detection
pfa = 1e-6;          % Probability of false alarm
max_range = 8000;    % Maximum unambiguous range
min_rcs = 1;         % Required target radar cross section
range_res = 50;      % Required range resolution

%% Design systems and generate radar data
[rx,mf,fs,noisepower] = ...
    getLFMdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res);
lfmdata = struct('fs',fs,'propspeed',prop_speed, ...
    'rxdata',rx,'mfcoeffs',mf,'noisepower',noisepower);

[rx,xref,slope,fs,noisepower] = ...
    getFMCWdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res);
fmcwdata = struct('fs',fs,'propspeed',prop_speed, ...
    'rxdata',rx,'refsig',xref,'sweepslope',slope,'noisepower',noisepower);

save('RangeResponseExampleData.mat','lfmdata','fmcwdata');

%% LFM Example
clear; figure(1); clf;

% Example:
%   Calculate range response of a pulsed LFM radar signal. The signal
%   includes three target returns. Two are around 2000 m away and the
%   third is around 3500 m away.
    
% Load example data
load('RangeResponseExampleData','lfmdata');
fs = lfmdata.fs;
propspeed = lfmdata.propspeed;
rxdata = lfmdata.rxdata;
mfcoeffs = lfmdata.mfcoeffs;
noisepower = lfmdata.noisepower;

% Create range response for matched filter processing
rngresp = phased.RangeResponse(...
    'SampleRate',fs,...
    'PropagationSpeed',propspeed);

% Perform range processing
[resp,rng_grid] = rngresp(rxdata,mfcoeffs);

% Calculate noise power after range processing
mfgain = mfcoeffs'*mfcoeffs;
noisepower_proc = mfgain*noisepower;

subplot(2,1,1);
plot(rng_grid,pow2db(abs(rxdata).^2./noisepower));
ylabel('SNR (dB)'); title('Before Range Processing'); xlim(rng_grid([1 end]));

subplot(2,1,2);
plot(rng_grid,pow2db(abs(resp).^2./noisepower_proc));
ylabel('SNR (dB)'); title('After Range Processing'); xlim(rng_grid([1 end]));
xlabel('Range (m)');

%% FMCW Example
clear; figure(2); clf;

% Example:
%   Plot the range response of an FMCW signal. The signal is not dechirped.
%   The signal contains the return from one target which is around 2200 m
%   away.

% Load example data
load('RangeResponseExampleData','fmcwdata');
fs = fmcwdata.fs;
propspeed = fmcwdata.propspeed;
rxdata = fmcwdata.rxdata;
refsig = fmcwdata.refsig;
sweepslope = fmcwdata.sweepslope;

% Create range response for dechirp and FFT processing
rngresp = phased.RangeResponse(...
    'RangeMethod','FFT',...
    'SweepSlope',sweepslope,...
    'DechirpInput',true,...
    'SampleRate',fs,...
    'PropagationSpeed',propspeed);

% Plot range response of processed data
plotResponse(rngresp,rxdata,refsig,'Unit','db');

%% Helpers
function [rx_data,mfcoeffs,fs,noisepower] = getLFMdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res)

%% Waveform and waveform processing
pulse_bw = prop_speed/(2*range_res);    % Pulse bandwidth
pulse_width = 20/pulse_bw;              % Pulse width
prf = prop_speed/(2*max_range);         % Pulse repetition frequency
fs = 2*pulse_bw;                        % Sampling rate
waveform = phased.LinearFMWaveform(...
    'SweepBandwidth',pulse_bw,...
    'PulseWidth',pulse_width,...
    'PRF',prf,...
    'SampleRate',fs);

% Matched filter coefficients
mfcoeffs = getMatchedFilter(waveform);

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
snr_min = albersheim(pd, pfa);
lambda = prop_speed/fc;
peak_power = radareqpow(lambda,max_range,snr_min,pulse_width,...
    'RCS',min_rcs,'Gain',tx_gain);

transmitter = phased.Transmitter(...
    'Gain',tx_gain,...
    'PeakPower',peak_power,...
    'InUseOutputPort',true);

%% Targets
tgtpos = [[2000;0;0],[2100;0;0],[3500;0;0]];
tgtvel = [[0;0;0],[0;0;0],[0;0;0]];
tgtrcs = [1.6 2.2 1.05];

%% Generate data
rx_data = simulateSystem(fc,waveform,transmitter,receiver,tgtpos,tgtvel,tgtrcs);

noise_bw = fs/2;
noisepower = noisepow(noise_bw,...
    receiver.NoiseFigure,receiver.ReferenceTemperature);
end

function [rx_data,xref,sweepslope,fs,noisepower] = getFMCWdata(fc,prop_speed,pd,pfa,max_range,min_rcs,range_res)

%% Waveform and waveform processing
sweep_bw = prop_speed/(2*range_res);
sweep_time = 5.5*range2time(max_range,prop_speed);
fs = sweep_bw;
sweepslope = sweep_bw/sweep_time;
waveform = phased.FMCWWaveform(...
    'SweepBandwidth',sweep_bw,...
    'SweepTime',sweep_time,...
    'SampleRate',fs);

% Reference signal for dechirp processing
xref = waveform();

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
snr_min = albersheim(pd, pfa);
lambda = prop_speed/fc;
peak_power = radareqpow(lambda,max_range,snr_min,sweep_time,...
    'RCS',min_rcs,'Gain',tx_gain);

transmitter = phased.Transmitter(...
    'Gain',tx_gain,...
    'PeakPower',peak_power,...
    'InUseOutputPort',true);

%% Targets
tgtpos = [2200;0;0];
tgtvel = [0;0;0];
tgtrcs = 1.6;

rx_data = simulateSystem(fc,waveform,transmitter,receiver,tgtpos,tgtvel,tgtrcs);

%% Generate data
noise_bw = fs/2;
noisepower = noisepow(noise_bw,...
    receiver.NoiseFigure,receiver.ReferenceTemperature);
end

function rxpulse = simulateSystem(fc,waveform,transmitter,receiver,tgtpos,tgtvel,tgtrcs)

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

%% Simulate system
% Update sensor and target positions
if isa(waveform,'phased.LinearFMWaveform')
    time_step = 1/waveform.PRF;
elseif isa(waveform,'phased.FMCWWaveform')
    time_step = waveform.SweepTime;
end
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
end
