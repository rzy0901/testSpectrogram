function helperslexBeamformerInterferenceParam
% This function is only in support of slexBeamscanMVDRDOAExample. 
% It may be removed in a future release.

%   Copyright 2014-2016 The MathWorks, Inc.

% Environment
prop_speed = physconst('LightSpeed');   % Propagation speed
fc = 100e6;           % Operating frequency
lambda = prop_speed/fc; % Wavelength
paramBeamformerI.propSpeed = prop_speed;
paramBeamformerI.fc = fc;

% Antenna
paramBeamformerI.Antenna = phased.ULA('NumElements',10,'ElementSpacing',0.5*lambda);

% Pulse
fs  = 1000; %1khz
paramBeamformerI.fs = fs;
prf = 1/.3;
paramBeamformerI.prf = prf;
paramBeamformerI.samplesPerFrame = fs/prf;

% LCMV Constraint Matrix
steeringvec = phased.SteeringVector('SensorArray',paramBeamformerI.Antenna);
paramBeamformerI.cMatrix = steeringvec(fc,[43 45 47]);

assignin('base','paramBeamformerI',paramBeamformerI);
