%This script is for internal use only. It may be removed in the future.

% Helper script creates data for help examples for
% phased.RangeEstimator and phased.DopplerEstimator
%
% Radar system is based on featured example:
%    "Automotive Adaptive Cruise Control Using FMCW and MFSK Technology"

%   Copyright 2016 The MathWorks, Inc.

rng default

%% Setup scenario
% Radar parameters
%   These numbers are taken from: helperslexFMCWMultiTargetsDOAParam.m
fs = 150e6;
propspeed = 3e8;
fc = 77e9;
Tpulse = 5.5*400/propspeed;

peakpower = 0.00316227766016838;
txgain = 36.0042142909402;

Nf = 4.5;
rxgain = 42.0042142909402;

radarpos = [0;0;0];
radarvel = [0;0;0];

% Target parameters
Ntgts = 3;
tgtpos = zeros(Ntgts);
tgtpos(1,:) = [500 520 750];

tgtvel = zeros(3,Ntgts);
tgtvel(1,:) = [-60 20 40];

tgtrcs = db2pow(30)*[1 1 1];

%% Transmitter
bw = fs/2;
pulsewidth = 1/bw;
pulsesamps = pulsewidth*fs;
pulsewidth = ceil(pulsesamps)/fs;
prf = 1/Tpulse;

waveform = phased.RectangularWaveform(...
    'SampleRate',fs,...
    'PulseWidth',pulsewidth,...
    'PRF',prf);

matchingcoeff = getMatchedFilter(waveform);

numPulse = 10;
numRng = fs/prf;
numDop = 32;

rangedopresp = phased.RangeDopplerResponse(...    
    'SampleRate',fs,...
    'PropagationSpeed',propspeed,...
    'DopplerFFTLengthSource','Property',...
    'DopplerFFTLength',numDop,...
    'DopplerOutput','Speed',...
    'OperatingFrequency',fc);

sig = waveform();
bwrms = rmsbw(sig,fs);
rngrms = bw2range(bwrms,propspeed);

transmitter = phased.Transmitter(...
    'PeakPower',peakpower,...
    'Gain',txgain,...
    'InUseOutputPort',true);

txarray = phased.IsotropicAntennaElement;

radiator = phased.Radiator(...
    'Sensor',txarray,...
    'PropagationSpeed',propspeed,...
    'OperatingFrequency',fc);

sensormotion = phased.Platform(radarpos,radarvel);

%% Channel
channel = phased.FreeSpace(...
    'SampleRate',fs,...    
    'PropagationSpeed',propspeed,...
    'OperatingFrequency',fc,...
    'TwoWayPropagation',true);

%% Receiver
rxarray = clone(txarray);

collector = phased.Collector(...
    'Sensor',rxarray,...
    'PropagationSpeed',propspeed,...
    'OperatingFrequency',fc);

receiver = phased.ReceiverPreamp(...
    'SampleRate',fs,...
    'Gain',rxgain,...
    'NoiseFigure',Nf);

%% Targets
tgtmotion = phased.Platform(tgtpos,tgtvel);

target = phased.RadarTarget(...
    'PropagationSpeed',propspeed,...
    'OperatingFrequency',fc,...
    'MeanRCS',tgtrcs);

%% Signal processing
rangeestimator = phased.RangeEstimator(...
    'ClusterInputPort',false,...
    'VarianceOutputPort',true,...
    'NoisePowerSource','Input port',...
    'RMSResolution',rngrms);

dopestimator = phased.DopplerEstimator(...
    'ClusterInputPort',false,...
    'VarianceOutputPort',true,...
    'NoisePowerSource','Input port',...
    'NumPulses',numPulse);
    

%% Build data cube
dt = 1/prf;
cube = zeros(numRng,numPulse);
for n = 1:numPulse
    
    [sensorpos,sensorvel] = sensormotion(dt);
    [tgtpos,tgtvel] = tgtmotion(dt);
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);
    
    % Simulate propagation of pulse in direction of targets
    sig = waveform();
    [txsig,txstatus] = transmitter(sig);
    txsig = radiator(txsig,tgtang);
    txsig = channel(txsig,sensorpos,tgtpos,sensorvel,tgtvel);
    
    % Reflect pulse off of targets
    tgtsig = target(txsig);
    
    % Receive target returns at sensor
    rxcol = collector(tgtsig,tgtang);
    rxsig = receiver(rxcol);
    cube(:,n) = rxsig;
end

%% Range and Doppler processing
[resp,rnggrid,spdgrid] = rangedopresp(cube,matchingcoeff);
mfgain = matchingcoeff'*matchingcoeff;
dopgain = numPulse;
noisebw = fs;
noisepower = noisepow(noisebw,receiver.NoiseFigure,receiver.ReferenceTemperature);
noisepowerprc = mfgain*dopgain*noisepower;

%% Substitute for CFAR detector
detidx = NaN(Ntgts,2);

tgtrng = rangeangle(tgtpos,radarpos);
tgtspd = radialspeed(tgtpos,tgtvel,radarpos,radarvel);
% tgtdop = 2*speed2dop(tgtspd,propspeed/fc);
for m = 1:numel(tgtrng)
    [~,iMin] = min(abs(rnggrid-tgtrng(m)));
    detidx(m,1) = iMin;
    [~,iMin] = min(abs(spdgrid-tgtspd(m)));
    detidx(m,2) = iMin;
end
detidx = detidx';
noise = noisepowerprc*ones(1,size(detidx,2));

%% Compute estimates
[rngest,rngvar] = rangeestimator(resp,rnggrid,detidx,noise);
[spdest,spdvar] = dopestimator(resp,spdgrid,detidx,noise);

%% Save off variables for help examples
save('RangeDopplerEstimatorData.mat', ...
    'rngrms','numPulse','resp','rnggrid','spdgrid','detidx','noise');

%% Helpers
function bwrms = rmsbw(x,Fs)
%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, 2nd ed.,
%       McGraw-Hill Professional Engineering, 2014

[Pxx,W] = periodogram(x,[],[],Fs,'centered');
dW = diff(W(1:2));

Es = sum(Pxx*dW);
bwrms = sqrt(sum(W.^2.*Pxx*dW)./Es);
end
