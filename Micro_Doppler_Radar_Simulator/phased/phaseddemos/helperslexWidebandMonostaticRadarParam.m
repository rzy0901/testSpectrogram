function helperslexWidebandMonostaticRadarParam
% This function is only in support of slexWidebandMonostaticRadarExample. 
% It may be removed in a future release.

%   Copyright 2016 The MathWorks, Inc.

% Set random number generator to a known state so that results are
% reproducible
rng default;

[propSpeed, fc, groundCoefficient, numSubbands, fs, txGain, ...
    peakPower, pulseBw, pulseWidth, prf, numPulseInt, numDopBins, ...
    stretchRef, stretchSpan, stretchSlope, decimateFactor, padLength, ...
    stretchMetersPerSamples, stretchRngLim, stretchFs, ...
    numCFARGuard, numCFARTrain, cfarThreshold, cfarRngIdx, cfarSpdIdx, ...
    cfarRngLim, cfarSpeedLim, sensorPos, sensorVel, ...
    targetPos, targetVel, rcsPattern, rcsFreq, rcsAz, rcsEl] = calcParams();

    % Environment
    paramWidebandRadar.propSpeed = propSpeed;
    paramWidebandRadar.fc = fc;
    paramWidebandRadar.groundCoefficient = groundCoefficient;
    paramWidebandRadar.numSubbands = numSubbands;
    % Transmitter parameters
    paramWidebandRadar.txGain = txGain;
    paramWidebandRadar.peakPower =  peakPower;
    % Waveform parameters
    paramWidebandRadar.pulseBw = pulseBw;
    paramWidebandRadar.prf = prf;
    paramWidebandRadar.pulseWidth = pulseWidth;
    paramWidebandRadar.fs = fs;
    % Pulse integration parameters
    paramWidebandRadar.numPulseInt = numPulseInt;
    paramWidebandRadar.numDopBins = numDopBins;
    % Stretch processor parameters
    paramWidebandRadar.stretchRef = stretchRef;
    paramWidebandRadar.stretchSpan = stretchSpan;
    paramWidebandRadar.stretchSlope = stretchSlope;
    paramWidebandRadar.decimateFactor = decimateFactor;
    paramWidebandRadar.padLength = padLength;
    paramWidebandRadar.stretchMetersPerSamples = stretchMetersPerSamples;
    paramWidebandRadar.stretchRngLim = stretchRngLim;
    paramWidebandRadar.stretchFs = stretchFs;
    % CFAR parameters
    paramWidebandRadar.numCFARGuard = numCFARGuard;
    paramWidebandRadar.numCFARTrain = numCFARTrain;
    paramWidebandRadar.cfarThreshold = cfarThreshold;
    paramWidebandRadar.cfarRngIdx = cfarRngIdx;
    paramWidebandRadar.cfarSpdIdx = cfarSpdIdx;
    % Sensor parameters
    paramWidebandRadar.sensorPos = sensorPos;
    paramWidebandRadar.sensorVel = sensorVel;
    % Target parameters
    paramWidebandRadar.targetPos = targetPos;
    paramWidebandRadar.targetVel = targetVel;
    paramWidebandRadar.rcsPattern = rcsPattern;
    paramWidebandRadar.rcsFreq = rcsFreq;
    paramWidebandRadar.rcsAz = rcsAz;
    paramWidebandRadar.rcsEl = rcsEl;
    % Scope parameters
    paramWidebandRadar.rangeLim = cfarRngLim;
    paramWidebandRadar.radialSpeedLim = cfarSpeedLim;

    assignin('base','paramWidebandRadar',paramWidebandRadar);

end

function [propSpeed, fc, groundCoefficient, numSubbands, fs, txGain, ...
    peakPower, pulseBw, pulseWidth, prf, numPulseInt, numDopBins, ...
    stretchRef, stretchSpan, stretchSlope, decimateFactor, padLength, ...
    stretchMetersPerSamples, stretchRngLim, stretchFs, ...
    numCFARGuard, numCFARTrain, cfarThreshold, cfarRngIdx, cfarSpdIdx, ...
    cfarRngLim, cfarSpeedLim, sensorPos, sensorVel, ...
    targetPos, targetVel, rcsPattern, rcsFreq, rcsAz, rcsEl] = calcParams()
    % Environment
    propSpeed = physconst('LightSpeed');    % Propagation speed
    groundCoefficient = -0.7;
    numSubbands = 64;
    
    % Radar design parameters
    maxRng = 5000;      % Required max detection range
    pd = 0.9;           % Probability of detection
    pfa = 1e-6;         % Probability of false alarm
    minTgtRcs = 1;      % Required target radar cross section
    numPulseInt = 10;   % Integrate 10 pulses at a time
    fc = 10e9;          % Operating frequency
    pulseBw = 0.10*fc;  % Wideband system with bandwidth set to 10% of
                        % center frequency
    pulseDuty = 0.2;    % Pulse duty cycle, fraction of PRI length
    
    lambda = propSpeed/fc;
    fs = 2*pulseBw;
    rangeRes = propSpeed/(2*pulseBw);

    % Sensor parameters
    sensorPos = [0;0;100];
    sensorVel = [0;0;0];
    
    % Waveform parameters
    pulseBw = propSpeed/(2*rangeRes);       % Pulse bandwidth
    prf = propSpeed/(2*maxRng);             % Pulse repetition frequency
    prf = fs/floor(fs/prf);                 % Must have integer number of
    priLength = fs/prf;                     % samples in a PRI
    pulseWidth = pulseDuty/prf;

    % Transmitter parameters
    snrMin = albersheim(pd, pfa, numPulseInt);
    txGain = 20;
    peakPower =  ...
        radareqpow(lambda,maxRng,snrMin,pulseWidth,...
                   'RCS',minTgtRcs,'Gain',txGain);

    % Target parameters
    radialSpeed = 100;
    targetPos = [3000;0;100];
    targetVel = [sensorVel(1)+radialSpeed;0;0];
    
    % Wideband target RCS model
    targetScatPos = [-0.5 -0.5 0.5 0.5;0.5 -0.5 0.5 -0.5;0 0 0 0];
    
    numRCSFreq = 4;
    azgrid = -180:20:180;
    elgrid = 0;
    rcsFreq = ((0:numRCSFreq-1)/numRCSFreq-0.5)*fs+fc;
    [cylRCS,rcsAz,rcsEl] = helperCylinderRCSPattern(propSpeed,rcsFreq,azgrid,elgrid);
    
    nf = numel(rcsFreq);
    naz = numel(rcsAz);
    nel = numel(rcsEl);
    rcsPattern = zeros(nel,naz,nf);
    for k = 1:nf
        for m = 1:nel
            pos = targetScatPos*rcsFreq(k)/fc;
            sv = steervec(pos,[rcsAz;rcsEl(m)*ones(1,naz)]);
            % sv is squared due to round trip in a monostatic scenario
            rcsPattern(m,:,k) = abs(sqrt(cylRCS(m,:,k)).*sum(sv.^2)).^2;
        end
    end
    
    rcsPattern = rcsPattern/max(rcsPattern(:)); % Normalize RCS to 1sqm
    rcsPattern = minTgtRcs*rcsPattern; % Scale to desired peak RCS
    
    % Stretch processor parameters
    stretchWindow = 200/2*[-1 1]+3000;

    stretchSlope = pulseBw/pulseWidth;
    stretchRef = mean(stretchWindow);
    stretchSpan = diff(stretchWindow);
    stretchSpanMax = propSpeed*(fs/2)/(2*stretchSlope);

    decimateFactor = floor(stretchSpanMax/stretchSpan);
    stretchPRILength = ceil(priLength/decimateFactor);
    padLength = stretchPRILength*decimateFactor;
    stretchFs = fs/decimateFactor;

    % Range-Doppler parameters
    numDopBins = 32;
    xcpi = zeros(stretchFs/prf,numDopBins);
    rngdopresp = phased.RangeDopplerResponse(...
        'RangeMethod','FFT',...
        'PropagationSpeed',propSpeed,...
        'SampleRate',stretchFs,...
        'SweepSlope',stretchSlope,...
        'DopplerOutput','Speed',...
        'OperatingFrequency',fc);
    [~,stretchRng,radialSpeedBins] = rngdopresp(xcpi);
    stretchRng = stretchRng+stretchRef;
    stretchMetersPerSamples = diff(stretchRng(1:2));
    stretchRngLim = [min(stretchRng) max(stretchRng)];
    
    % CFAR parameters
    nGuardRng = 4;
    nTrainRng = 16;
    nCUTRng = 1+nGuardRng+nTrainRng;
    iCUTRng = interp1(stretchRng,1:stretchPRILength,...
        stretchRef+stretchSpan/16*[-1 1],'nearest'); % Zoomed in range
    iCUTRng = iCUTRng(1):iCUTRng(2);
    iCUTRng = iCUTRng(iCUTRng>=nCUTRng & iCUTRng<=stretchPRILength-nCUTRng+1);
    cfarRngLim = [stretchRng(iCUTRng(1)) stretchRng(iCUTRng(end))];
    cfarRngLim = [min(cfarRngLim) max(cfarRngLim)];
    
    nGuardDop = 3;
    nTrainDop = 3;
    nCUTDop = 1+nGuardDop+nTrainDop;
    iCUTDop = 1:numDopBins;
    iCUTDop = iCUTDop(iCUTDop>=nCUTDop & iCUTDop<=numDopBins-nCUTDop+1);
    cfarSpeedLim = [radialSpeedBins(iCUTDop(1)) radialSpeedBins(iCUTDop(end))];
    cfarSpeedLim = [min(cfarSpeedLim) max(cfarSpeedLim)];

    cfarRngIdx = iCUTRng([1 end]);
    cfarSpdIdx = iCUTDop([1 end]);
    
    numCFARGuard = [nGuardRng nGuardDop];
    numCFARTrain = [nTrainRng nTrainDop];

    cfarThreshold = db2pow(snrMin);
    
end
