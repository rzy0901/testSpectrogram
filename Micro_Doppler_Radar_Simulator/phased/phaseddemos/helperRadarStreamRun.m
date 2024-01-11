function helperRadarStreamRun
% This function helperRadarStreamRun is only in support of
% RadarStreamExample. It may be removed in a future release.

%   Copyright 2013-2014 The MathWorks, Inc.

%#codegen

coder.extrinsic('helperRadarStreamDisplay');

fs = 6e6;
bw = 3e6;
c = 3e8;
fc = 10e9;
prf = 18750;

waveform = phased.LinearFMWaveform(...
    'SampleRate',fs,'PulseWidth',20/bw,...
    'PRF',prf,'SweepBandwidth',bw);

transmitter = phased.Transmitter(...
    'PeakPower',2240,'InUseOutputPort',true);

radiator = phased.Radiator(...
    'Sensor',phased.IsotropicAntennaElement,...
    'PropagationSpeed',c,'OperatingFrequency',fc);

collector = phased.Collector(...
    'Sensor',phased.IsotropicAntennaElement,...
    'PropagationSpeed',c,'OperatingFrequency',fc);

receiver = phased.ReceiverPreamp(...
    'SampleRate',fs,'EnableInputPort',true);

sensormotion = phased.Platform;

tgtpos = [[2000.66;0;0],[6532.63;0;0],[6845.04;0;0]];
tgtvel = [[0;0;0],[300;0;0],[-300;0;0]];
tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel);

tgtrcs = [2.2 1.1 1.05];
target = phased.RadarTarget(...
    'MeanRCS',tgtrcs,'PropagationSpeed',c,'OperatingFrequency',fc);

channel = phased.FreeSpace(...
    'PropagationSpeed',c,'OperatingFrequency',fc,...
    'TwoWayPropagation',true,'SampleRate',fs);

match_sig = getMatchedFilter(waveform);
matchedfilter =  phased.MatchedFilter(...
    'Coefficients',match_sig);

fast_time_grid = 0:1/fs:1/prf-1/fs;
range_gates = c*fast_time_grid/2; 
max_range = 8e3;
lambda = c/fc;

rng_loss = 2*fspl(range_gates,lambda); % factor 2 for round trip
ref_loss = 2*fspl(max_range,lambda);   % reference loss for maximum range

tvg = phased.TimeVaryingGain(...
    'RangeLoss',rng_loss,...
    'ReferenceLoss',ref_loss);

num_pulse_int = 10;
threshold = 4e-11;

num_pulse_samples = numel(fast_time_grid);
rx_pulses = complex(zeros(num_pulse_samples,num_pulse_int)); % pre-allocate 
mf_pulses = complex(zeros(num_pulse_samples,num_pulse_int)); % pre-allocate 
detect_pulse = zeros(num_pulse_samples,1);

for m = 1:100*num_pulse_int
    
    % Update sensor and target positions
    [sensorpos,sensorvel] = sensormotion(1/prf);
    [tgtpos,tgtvel] = tgtmotion(1/prf);

    % Calculate the target angles as seen by the sensor
    [~,tgtang] = rangeangle(tgtpos,sensorpos);
    
    % Simulate propagation of pulse in direction of targets
    pulse = waveform();
    [pulse,txstatus] = transmitter(pulse);
    txsig = radiator(pulse,tgtang);
    txsig = channel(txsig,sensorpos,tgtpos,sensorvel,tgtvel);
    
    % Reflect pulse off of targets
    tgtsig = target(txsig);
    
    % Receive target returns at sensor
    rxsig = collector(tgtsig,tgtang);
    nn = mod(m-1,num_pulse_int)+1;
    rx_pulses(:,nn) = receiver(rxsig,~(txstatus>0));
    
    % Detection processing
    mf_pulses(:,nn) = matchedfilter(rx_pulses(:,nn));
    mf_pulses(:,nn) = tvg(mf_pulses(:,nn));
    
    % Perform pulse integration every 'num_pulse_int' pulses
    if nn == num_pulse_int                        % detection 
        detect_pulse = pulsint(mf_pulses,'noncoherent');
    end
    
    helperRadarStreamDisplay(pulse,abs(rx_pulses(:,nn)),...
        abs(mf_pulses(:,nn)),detect_pulse,...
        sqrt(threshold)*ones(num_pulse_samples,1))
end
end
