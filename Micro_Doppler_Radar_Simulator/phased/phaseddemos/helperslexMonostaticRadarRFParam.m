function helperslexMonostaticRadarRFParam
% This function is only in support of slexMonostaticRadarExample. 
% It may be removed in a future release.

%   Copyright 2014 The MathWorks, Inc.

    [propSpeed, fc, pulseBw, prf, fs, txGain, peakPower, ...
     matchingCoeff, metersPerSample, rangeOffset, rangeLoss, ...
     referenceLoss, target1Rcs, target1Pos, target1Vel] = calcParams();
    
    % Environment
    paramRadarRF.propSpeed = propSpeed;
    paramRadarRF.fc = fc;
    % Waveform parameters
    paramRadarRF.pulseBw = pulseBw;
    paramRadarRF.prf = prf;
    paramRadarRF.fs = fs;
    % Transmitter parameters
    paramRadarRF.txGain = txGain;
    paramRadarRF.peakPower =  peakPower;
    % Matched filter parameters
    paramRadarRF.matchingCoeff = matchingCoeff;
    % Time varying gain parameters 
    paramRadarRF.metersPerSample = metersPerSample;
    paramRadarRF.rangeOffset = rangeOffset;
    paramRadarRF.rangeLoss = rangeLoss;
    paramRadarRF.referenceLoss = referenceLoss;
    % Radar parameters
    paramRadarRF.target1Rcs = target1Rcs;
    paramRadarRF.target1Pos = target1Pos;
    paramRadarRF.target1Vel = target1Vel;

    assignin('base','paramRadarRF',paramRadarRF);

end

function [propSpeed, fc, pulseBw, prf, fs, txGain, peakPower, ...
          matchingCoeff, metersPerSample, rangeOffset, rangeLoss, ...
          referenceLoss, target1Rcs, target1Pos, target1Vel]  = calcParams()  
    % Environment
    propSpeed = physconst('LightSpeed');   % Propagation speed
    fc = 10e9;           % Operating frequency
    lambda = propSpeed/fc;


    % Constraints
    maxRange = 5000;    % Maximum unambiguous range
    rangeRes = 50;      % Required range resolution
    pd = 0.9;            % Probability of detection
    pfa = 1e-6;          % Probability of false alarm
    tgtRcs = 1;         % Required target radar cross section
    numPulseInt = 10;  % Integrate 10 pulses at a time


    % Waveform parameters
    pulseBw = propSpeed/(2*rangeRes);    % Pulse bandwidth
    pulseWidth = 1/pulseBw;               % Pulse width
    prf = propSpeed/(2*maxRange);         % Pulse repetition frequency
    fs = 2*pulseBw;    

    % Transmitter parameters
    snrMin = albersheim(pd, pfa, numPulseInt);
    txGain = 20;
    peakPower =  ...
        radareqpow(lambda,maxRange,snrMin,pulseWidth,...
                   'RCS',tgtRcs,'Gain',txGain);

    % Matched filter parameters
    hwav = phased.RectangularWaveform(...
        'PulseWidth',1/pulseBw,...
        'PRF',prf,...
        'SampleRate',fs);
    matchingCoeff = getMatchedFilter(hwav);

    % Delay introduced due to filter
    matchingDelay = size(matchingCoeff,1)-1;

    % Time varying gain parameters 
    fastTimeGrid = unigrid(0,1/fs,1/prf,'[)');
    rangeGates = propSpeed*fastTimeGrid/2; 
    metersPerSample = rangeGates(2);
    rangeOffset = -rangeGates(2)*matchingDelay;
    rangeLoss = 2*fspl(rangeGates,lambda);
    referenceLoss = 2*fspl(maxRange,lambda);

    %Radar parameters
    target1Rcs = 0.6;
    target1Pos = [1988.66;0;0];
    target1Vel = [ 0; 0; 0 ];

end
