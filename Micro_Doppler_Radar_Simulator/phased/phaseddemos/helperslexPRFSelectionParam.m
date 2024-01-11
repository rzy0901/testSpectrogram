function helperslexPRFSelectionParam
% This function is only in support of slexMonostaticRadarExample. 
% It may be removed in a future release.

%   Copyright 2014 The MathWorks, Inc.

    [propSpeed, fc, pulseBw, prf, fs, txGain, peakPower, ...
     matchingCoeff, metersPerSample, rangeOffset, rangeLoss, ...
     referenceLoss, target1Rcs, target1Pos, target1Vel] = calcParams();
    
    % Environment
    paramPRFSel.propSpeed = propSpeed;
    paramPRFSel.fc = fc;
    % Waveform parameters
    paramPRFSel.pulseBw = pulseBw;
    paramPRFSel.prf = prf;
    paramPRFSel.fs = fs;
    % Transmitter parameters
    paramPRFSel.txGain = txGain;
    paramPRFSel.peakPower =  peakPower;
    % Matched filter parameters
    paramPRFSel.matchingCoeff = matchingCoeff;
    % Time varying gain parameters 
    paramPRFSel.metersPerSample = metersPerSample;
    paramPRFSel.rangeOffset = rangeOffset;
    paramPRFSel.rangeLoss = rangeLoss;
    paramPRFSel.referenceLoss = referenceLoss;
    % Radar parameters
    paramPRFSel.target1Rcs = target1Rcs;
    paramPRFSel.target1Pos = target1Pos;
    paramPRFSel.target1Vel = target1Vel;

    assignin('base','paramPRFSel',paramPRFSel);

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
    fs = 2*pulseBw;   
    
    prf_search = propSpeed/(2*maxRange);  % Pulse repetition frequency
    lockRange = 2000;                   % Lock in target
    prf_lock = propSpeed/(2*lockRange);
    prf = [prf_search prf_lock];

    % Transmitter parameters
    snrMin = albersheim(pd, pfa, numPulseInt);
    txGain = 20;
    peakPower =  ...
        radareqpow(lambda,maxRange,snrMin,pulseWidth,...
                   'RCS',tgtRcs,'Gain',txGain);

    % Matched filter parameters
    hwav = phased.RectangularWaveform(...
        'PulseWidth',1/pulseBw,...
        'PRF',prf_search,...
        'SampleRate',fs);
    matchingCoeff = getMatchedFilter(hwav);

    % Delay introduced due to filter
    matchingDelay = size(matchingCoeff,1)-1;

    % Time varying gain parameters 
    fastTimeGrid = unigrid(0,1/fs,1/prf_search,'[)');
    rangeGates = propSpeed*fastTimeGrid/2; 
    metersPerSample = rangeGates(2);
    rangeOffset = -rangeGates(2)*matchingDelay;
    rangeLoss = 2*fspl(rangeGates,lambda);
    referenceLoss = 2*fspl(maxRange,lambda);

    %Radar parameters
    target1Rcs = [1 0.5];
    target1Pos = [2005 3500;0 0;0 0];
    target1Vel = [ -300 -200; 0 0; 0 0 ];

end
