% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     pattern.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       September 2015 created
%
%  *************************************************************************************
%    Description:
% 
%    Calculates value of Antenna Pattern for given direction
%
%    out = pattern(el, az, W)
%
%    Inputs:
%
%       1. el   - elevation angles array [deg]
%       2. az   - azimuth   angles array [deg]
%       3. W    - weight coefficients
%       4. paa  - structure with PAA parameters
%
%    Outputs:
%
%       1. out  - value of antenna pattern in specified direction
%
%  *************************************************************************************/
function out = pattern(el, az, W, paa)

Nx = paa.Nx;
Ny = paa.Ny;
dx = paa.dx;
dy = paa.dy;
lyam = paa.lyam;

k = 2 * pi / lyam;

% Calculate value of antenna pattern in specified direction
out = 0;
for nx = 0 : Nx - 1
    for ny = 0 : Ny - 1
        phase = exp(1j * ( k * sin(el) * cos(az) * dx * nx + k * sin(el) * sin(az) * dy * ny ));
        out = W(ny + 1 + Ny * nx) * phase + out;
    end
end
out = abs(out)^2;

end

