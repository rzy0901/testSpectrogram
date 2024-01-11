function helperslexFMCWParam
% This function helperslexFMCWParam is only in support of
% slexFMCWExample. It may be removed in a future release.

%   Copyright 2014 The MathWorks, Inc.

paramFMCW.Fs = 150e6;
paramFMCW.T = 5.5*400/3e8;
paramFMCW.ppow = 0.00316227766016838;
paramFMCW.TxGain = 36.0042142909402;
paramFMCW.RadarVel = [ 100*1000/3600; 0; 0];
paramFMCW.RadarPos = [0;0;0];
paramFMCW.Fc = 77e9;
paramFMCW.RCS = 100;
paramFMCW.CarVel = [ 96*1000/3600; 0; 0];
paramFMCW.CarPos = [43;0;0];
paramFMCW.C = 3e8;
paramFMCW.NF = 4.5;
paramFMCW.RxGain = 42.0042142909402;

assignin('base','paramFMCW',paramFMCW)