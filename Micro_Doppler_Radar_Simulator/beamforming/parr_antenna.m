% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     parr_antenna.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       September 2015 created
%
%  *************************************************************************************
%    Description:
% 
%    function returns amplitudes weighted by gain coefficients calculated
%    for phased antenna array model
%
%    [amg] = parr_antenna(am, az, el)
%
%    Inputs:
%
%       1. am   - array of input amplitudes
%       2. az   - azimuth   angles array [deg]
%       3. el   - elevation angles array [deg]
%       4. paa  - structure with PAA parameters
%
%    Outputs:
%
%       1. amg  - array of output amplitudes weighted by antenna gain coefficients
%
%  *************************************************************************************/
function [ amg ] = parr_antenna(am, az, el, W, paa)

Nray = length(am);
U = zeros(1, Nray);

az_rad = az * pi / 180;
el_rad = el * pi / 180;
% Calculate antenna gain value for given directions of rays
for i = 1 : Nray    
    U(i) = pattern(el_rad(i), az_rad(i), W, paa);
end

% Total power in whole space 
power = 0;
for el = 0 : 0.1 : pi/2
    power = power + simp(0, pi*2, el, W, 1e-1, paa);
end
power = power * 0.1;

% Calculate antenna gain
g = sqrt(4 * pi * U / power);

amg = am .* g.';
