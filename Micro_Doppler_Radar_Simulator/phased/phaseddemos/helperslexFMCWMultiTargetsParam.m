function helperslexFMCWMultiTargetsParam
% This function helperslexFMCWMultiTargetsParam is only in support of
% slexFMCWMultiTargetsExample. It may be removed in a future release.

%   Copyright 2014 The MathWorks, Inc.

paramFMCWMT.Fs = 150e6;
paramFMCWMT.bw = 150e6;
paramFMCWMT.T = 5.5*400/3e8;
paramFMCWMT.ppow = 0.00316227766016838;
paramFMCWMT.TxGain = 36.0042142909402;
paramFMCWMT.RadarVel = [ 100*1000/3600; 0; 0];
paramFMCWMT.RadarPos = [0;0;0];
paramFMCWMT.Fc = 77e9;
paramFMCWMT.RCS = 50;
paramFMCWMT.CarVel = [ 60*1000/3600; 0; 0];
paramFMCWMT.CarPos = [50;0;0];
paramFMCWMT.C = 3e8;
paramFMCWMT.NF = 4.5;
paramFMCWMT.RxGain = 42.0042142909402;
paramFMCWMT.slope = paramFMCWMT.bw/paramFMCWMT.T;
paramFMCWMT.lambda = paramFMCWMT.C/paramFMCWMT.Fc;
paramFMCWMT.TruckRCS = 1000;
paramFMCWMT.TruckVel = [ 130*1000/3600; 0; 0];
paramFMCWMT.TruckPos = [150;0;0];
paramFMCWMT.NRangeFFT = 2048;
paramFMCWMT.NDopplerFFT = 256;
paramFMCWMT.NBuffer = 64;

assignin('base','paramFMCWMT',paramFMCWMT)