% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     ant_gain.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       November 2015 created
%
%  *************************************************************************************
%    Description:
%
%    function returns amplitudes weighted by antenna gain for target phased antenna array space position
%
%    [amg] = ant_gain(ant_type,hpbw,am,az,el,az_max,el_max)
%
%    Inputs:
%
%       1. am       - amplitudes array 
%       2. az       - TX/RX azimuths array
%       3. el       - TX/RX elevations array
%       4. az_max   - azimuth angle of a ray with maximum power
%       5. el_max   - elevation angle of a ray with maximum power
%       6. paa_struct - keeps parameters of PAA
%
%    Outputs:
%
%       1. amg - output amplitudes weighted by antenna gain coefficients
%
%  *************************************************************************************/
function [amg] = ant_gain(am, az, el, az_max, el_max, paa_struct)

% Unpack parameters
Nx = paa_struct.Nx;
Ny = paa_struct.Ny;
dx = paa_struct.dx;
dy = paa_struct.dy;
lyam = paa_struct.lyam;
        
% Calculation of weight vector for ray w/ maximum power
k = 2 * pi / lyam;
a0 = k * sin(el_max * pi / 180) * cos(az_max * pi / 180) * dx;
a1 = k * sin(el_max * pi / 180) * sin(az_max * pi / 180) * dy;

% Calculate weight coefficients for PAA
W = zeros(1, Nx * Ny);
for nx = 0 : Nx - 1
   for ny = 0 : Ny - 1
       W(ny + 1 + Ny * nx) = exp(-1j * ( a0 * nx + a1 * ny )) / sqrt(Nx*Ny);
   end
end

% Apply antenna pattern
amg = parr_antenna(am, az, el, W, paa_struct);
end