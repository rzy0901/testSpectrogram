function helperslexFMCWRFParam
% This function helperslexFMCWParam is only in support of
% slexFMCWExample. It may be removed in a future release.

%   Copyright 2014 The MathWorks, Inc.

% In an ACC setup, the maximum range the radar needs to monitor is around
% 200 m and the system needs to be able to distinguish two targets that are
% 1 meter apart. From these requirements, one can compute the waveform
% parameters.
paramFMCWRF.Fc = 77e9;
paramFMCWRF.C = 3e8;
paramFMCWRF.lambda = paramFMCWRF.C/paramFMCWRF.Fc;

% The sweep bandwidth can be determined according to the range resolution (1m).
paramFMCWRF.bw = range2bw(1,paramFMCWRF.C);

% The sweep time can be computed based on the time needed for the signal to
% travel the unambiguous maximum range. In general, for an FMCW radar 
% system, the sweep time should be at least 5 to 6 times the round trip 
% time. This example uses a factor of 5.5.
% The range of the radar is 200m, for a total round trip of 400m.
% The sweep slope is calculated using both sweep bandwidth and sweep time.
paramFMCWRF.T = 5.5*range2time(200,paramFMCWRF.C); 
paramFMCWRF.slope = paramFMCWRF.bw/paramFMCWRF.T;

paramFMCWRF.Fs = 150e6;

% Car and Radar position and speed
paramFMCWRF.CarVel = [ 120*1000/3600; 0; 0];
paramFMCWRF.CarPos = [50;0;0];
paramFMCWRF.RadarVel = [ 90*1000/3600; 0; 0];
paramFMCWRF.RadarPos = [0;0;0];

% car_state = 1;
% if car_state == 1 % car is very far away with very different speed
%     paramFMCWRF.CarPos = [200;0;0];
%     paramFMCWRF.CarVel = [ 0*1000/3600; 0; 0];
%     paramFMCWRF.RadarVel = [ 230*1000/3600; 0; 0];
% end
% if car_state == 2 % car is close by at the same speed
%     paramFMCWRF.CarPos = [5;0;0];
%     paramFMCWRF.CarVel = [ 50*1000/3600; 0; 0];
%     paramFMCWRF.RadarVel = [ 50*1000/3600; 0; 0];
% end
paramFMCWRF.RCS = db2pow(min(10*log10(norm(paramFMCWRF.CarPos))+5,20));

paramFMCWRF.lna_gain = 42;         % in dB
paramFMCWRF.pa_gain = 30;          % in dB
paramFMCWRF.da_gain = 0;           % in dB
paramFMCWRF.mix_loss = 0;          % in dB
paramFMCWRF.vga_gain = 0;    

% Linear noiseless RF frontend
IncludeRFimperfections = 1;
paramFMCWRF.HarmonicOrder = 1;
paramFMCWRF.pa_ip3 = inf;          % in dBm
paramFMCWRF.pa_1dBc = inf;         % in dBm
paramFMCWRF.pa_psat = inf;         % in dBm
paramFMCWRF.lna_ip3 = inf;         % in dBm
paramFMCWRF.Risolation = 1e30;     % in Ohm
paramFMCWRF.mix_ip2 = inf;         % in dBm
paramFMCWRF.NoiseTemp = 1e-12;     % in K
paramFMCWRF.lna_nf = 0;            % in dB
paramFMCWRF.pa_nf = 0;             % in dB
paramFMCWRF.da_nf = 0;             % in dB
paramFMCWRF.mix_nf = 0;            % in dB
paramFMCWRF.vga_nf = 0;            % in dB

if IncludeRFimperfections == 1
    paramFMCWRF.HarmonicOrder = 3;
    paramFMCWRF.pa_ip3 = 23;       % in dB
    paramFMCWRF.pa_1dBc = 12.2;    % in dB
    paramFMCWRF.pa_psat = 14.9;    % in dB
    paramFMCWRF.lna_ip3 = -18;     % in dB
    paramFMCWRF.Risolation = 150e4; % in Ohm 
    paramFMCWRF.mix_ip2 = 65.3;    % in dB
    %paramFMCWRF.mix_ip2 = 48.3;    % in dB
    paramFMCWRF.NoiseTemp = 290;   % in K
    paramFMCWRF.lna_nf = 3.0;      % in dB
    paramFMCWRF.pa_nf = 8;         % in dB
    paramFMCWRF.da_nf = 4.0;       % in dB
    paramFMCWRF.mix_nf = 5;        % in dB
    paramFMCWRF.vga_nf = 6;        % in dB
end

assignin('base','paramFMCWRF',paramFMCWRF)