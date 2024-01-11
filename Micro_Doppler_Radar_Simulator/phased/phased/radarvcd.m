function [vcp, vcpangles] = radarvcd(freq, rfs, anht, varargin)
%radarvcd Vertical coverage diagram for radar
%   [VCP, VCPANGLES] = radarvcd(FREQ, RFS, ANHT) calculates the vertical
%   coverage pattern VCP at elevation angles VCPANGLES for a radar system
%   operating at a carrier frequency FREQ (in Hz), free space range RFS and
%   antenna height of ANHT. By default, the unit of RFS is km and the unit
%   of ANHT is m. FREQ must be less than or equal to 10e9.
%   
%   [...] = radarvcd(..., 'RangeUnit', RNGUNIT) specifies the unit of range
%   RFS as one of 'nmi' | 'mi' | 'km' | 'ft' | 'm'. The default value is
%   'km'.
% 
%   [...] = radarvcd(..., 'HeightUnit', HTUNIT) specifies the unit of
%   height ANHT as one of 'nmi' | 'mi' | 'km' | 'ft' | 'm'. The default
%   value is 'm'.
% 
%   [...] = radarvcd(..., 'Polarization', pol) specifies the polarization
%   of the transmitted wave. The parameter, pol, can be one of 'H' | 'V',
%   where the default value is 'H'. 'H' indicates horizontal polarization
%   and 'V' indicates vertical polarization.
% 
%   [...] = radarvcd(..., 'SurfaceDielectric', EPSC) specifies the complex
%   dielectric constant of the reflecting surface. The default value of
%   EPSC depends on the value of FREQ. It uses a sea water dielectric
%   model in [1].
% 
%   [...] = radarvcd(..., 'SurfaceRoughness', HTVAR) specifies the
%   roughness (height variation) of the reflecting surface. The default
%   value of HTVAR is 0, indicating a smooth surface. The roughness is
%   modeled as a sinusoid and HTVAR is the crest-to-trough height. The unit
%   of height is specified by HTUNIT.
% 
%   [...] = radarvcd(..., 'AntennaPattern', PAT, 'PatternAngles', PATEL)
%   specifies the antenna elevation pattern and corresponding elevation
%   angles (in degrees). Both PAT and PATANG must be column vectors and
%   they must have same size. PATANG must be between -90 and 90 degrees. In
%   general, to properly compute the coverage, the pattern should be
%   specified from -90 to 90 degrees.
% 
%   [...] = radarvcd(..., 'TiltAngle', TILTANG) specifies the tilt angle,
%   in degrees, of the antenna with respect to the surface.
% 
%   [...] = radarvcd(..., 'MaxElevation', MAXEL) specifies the maximum
%   elevation angle, in degrees, for which the vertical coverage pattern
%   will be calculated. The default value of MAXEL is 60.
%
%   radarvcd(...) plots the radar coverage diagram.
%
%   % Example:
%   %   Plot the radar vertical coverage pattern assuming the antenna has a
%   %   sinc pattern. The frequency is 100 MHz, the antenna height is 20
%   %   feet, and the range is 100 nautical miles. Assume the surface is
%   %   smooth, the antenna is not tilted and in horizontal polarization.
%
%   pat_angles = linspace(-90,90,361)';
%   pat_u = 1.39157/sind(90/2)*sind(pat_angles);
%   pat = sinc(pat_u/pi);
% 
%   freq = 100e6;     % MHz
%   ant_height = 20;  % ft
%   rng_fs = 100;     % Nautical mile
%   tilt_ang = 0;  % degrees
%   radarvcd(freq,rng_fs,ant_height,...
%       'RangeUnit','nmi','HeightUnit','ft','AntennaPattern',pat,...
%       'PatternAngles',pat_angles,'TiltAngle',tilt_ang);
%
%   See also phased, blakechart

%   Copyright 2012-2018 The MathWorks, Inc.

%   References
%     [1] Blake, L.V., Machine Plotting of Radar Vertical-Plane Coverage
%     Diagrams, Naval Research Laboratory, 1970 (NRL Report 7098)

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

coder.extrinsic('unitsratio');

phased.internal.narginchk(3, 21, nargin);

c = physconst('LightSpeed');
lambda = c/freq;
[rngunit, htunit, epsc, pol, anpat, patangles, ...
    tilt, htsd, vcpelmax] = parseInput(freq, rfs, anht, varargin{:});

% Make sure rfs is row vector
if iscolumn(rfs)    
    rfsC = reshape(rfs, 1, length(rfs));
else
    rfsC = rfs;
end
rfsC = sort(rfsC);  % Lo to hi rfs

% Convert anht and hstd to meter
anht = anht*coder.internal.const(unitsratio('m', htunit));
htsd = htsd*coder.internal.const(unitsratio('m', htunit));

thdel = 4.92/coder.internal.const(unitsratio('ft','m'))/...
    (freq*1e-6 * anht);   % Code snippet from NRL Report 7098, p38
thetad = (0:thdel:phased.internal.deg2rad(vcpelmax)); % Direct ray elevation angles
patangles = phased.internal.deg2rad(patangles);
tilt = phased.internal.deg2rad(tilt);

ae = effearthradius;
phi = sqrt((tan(thetad)/3).^2 + 2*anht/(3*ae)) - tan(thetad)/3;  % Eq25

% Compute vertical coverage pattern
pdCE = pathdiff_CE(thetad, anht, phi); % Curved earth, Eq 24
pdFE = pathdiff_FE(thetad, anht);  % Flat earth, Eq 23
ang_idx = find(pdCE>=lambda/4,1);
if ~isempty(ang_idx)
    start_ang_idx = ang_idx(1); %for codegen
else
    start_ang_idx = 1;
end
idxH = abs(pdCE - pdFE) < 0.01*lambda; % Higher elevation
idxL = ~idxH;                          % Lower elevation
idxHLimit = idxH(start_ang_idx:end);
pd = [pdCE(idxL) pdFE(idxH)];
pdLimit = pd(start_ang_idx:end);
thetadLimit = thetad(start_ang_idx:end);
psigrz = thetadLimit + phi(start_ang_idx:end);  % Eq 26
thetar = -psigrz; % Angle of reflected wave, Fig 1
thetar(idxHLimit) = -thetadLimit(idxHLimit);

[rhoref, phiref] = reflectioncoeff(epsc, psigrz, pol); % Eq 29, 30
r = roughness(lambda, htsd/(2*sqrt(2)), psigrz);  % Eq 20
D = divfactor(anht, thetadLimit); % Eq 27


D(idxHLimit) = 1;
alpha = (2 * pi * pdLimit / lambda) + phiref;  % Eq 9

% Propagation factor Eq 8 
theta1 = rem((thetadLimit - tilt), pi/2);
theta2 = rem((thetar - tilt), pi/2);
apatd = antfactor(anpat, patangles, theta1);
apatr = antfactor(anpat, patangles, theta2);
grc = (r .* D .* rhoref .* apatr./apatd);   % Eq 7
pf = abs(apatd).*(sqrt(1 + grc.^2 + 2*grc.*cos(alpha))); % Eq 8

% Vertical coverage pattern, same unit as range
% vcp is a matrix with each column is a separate pattern
pattern = pf.' * rfsC;  
patternang = phased.internal.rad2deg(thetadLimit'); 
% Choose rmax, hmax
if nargout % No output args
    vcp = pattern;
    vcpangles = patternang;
else    
    if ~isempty(coder.target)
        coder.internal.assert(false,'phased:rocsnr:invalidCodegenOutput','radarvcd');
    end
    % EXPO = 0.143859/unitsratio(htunit, 'km');  % per htunit
    blakechart(pattern, patternang,'RangeUnit', rngunit,...
        'HeightUnit', htunit, ...
        'ScalePower', 1);
end
end
 
% ----------------------------------
%           HELPER FUNCTIONS       |
% ----------------------------------
function [rngunitOut, htunit, epsc, pol, anpat, patangles, tilt, htsd, vcpelmax] = parseInput(freq, rfs, ha, varargin)

eml_assert_no_varsize(1:nargin, freq, rfs, ha, varargin{:});
funName = 'radarvcd';

defaultRangeUnit = 'km';
defaultHeightUnit = 'm';
defaultSurfaceDielectric = seacmplxdielectric(freq);
defaultPolarization = 'H';
[defaultAntennaPattern, defaultPatternAngles] = sincpattern();
defaultTiltAngle = 0;
defaultSurfaceRoughness = 0;
defaultMaxElevation = 60;

if isempty(coder.target)
    par = inputParser;
    par.addParameter('RangeUnit', defaultRangeUnit);
    par.addParameter('HeightUnit', defaultHeightUnit);
    par.addParameter('Polarization', defaultPolarization);
    par.addParameter('SurfaceDielectric', defaultSurfaceDielectric);
    par.addParameter('SurfaceRoughness', defaultSurfaceRoughness);
    par.addParameter('AntennaPattern', defaultAntennaPattern); % Vector
    par.addParameter('PatternAngles', defaultPatternAngles);  % Vector
    par.addParameter('TiltAngle', defaultTiltAngle);
    par.addParameter('MaxElevation', defaultMaxElevation); 

    par.parse(varargin{:});
    rngunit = par.Results.RangeUnit;
    htunit = par.Results.HeightUnit;
    pol = par.Results.Polarization;
    epsc = par.Results.SurfaceDielectric;
    htsd = par.Results.SurfaceRoughness;
    anpat = par.Results.AntennaPattern;
    patangles = par.Results.PatternAngles;
    tilt = par.Results.TiltAngle;
    vcpelmax = par.Results.MaxElevation;
else
     parms = struct('RangeUnit', uint32(0), ...
                    'HeightUnit', uint32(0), ...
                    'Polarization', uint32(0), ...
                    'SurfaceDielectric', uint32(0), ...
                    'SurfaceRoughness', uint32(0), ...
                    'AntennaPattern', uint32(0), ...
                    'PatternAngles', uint32(0), ...
                    'TiltAngle', uint32(0), ...
                    'MaxElevation', uint32(0));

     
     pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
     rngunit = eml_get_parameter_value(pstruct.RangeUnit,defaultRangeUnit,varargin{:});
     htunit = eml_get_parameter_value(pstruct.HeightUnit,defaultHeightUnit,varargin{:});
     pol = eml_get_parameter_value(pstruct.Polarization,defaultPolarization,varargin{:});
     epsc = eml_get_parameter_value(pstruct.SurfaceDielectric,defaultSurfaceDielectric,varargin{:});
     htsd = eml_get_parameter_value(pstruct.SurfaceRoughness,defaultSurfaceRoughness,varargin{:});
     anpat = eml_get_parameter_value(pstruct.AntennaPattern,defaultAntennaPattern,varargin{:});
     patangles = eml_get_parameter_value(pstruct.PatternAngles,defaultPatternAngles,varargin{:});
     tilt = eml_get_parameter_value(pstruct.TiltAngle,defaultTiltAngle,varargin{:});
     vcpelmax = eml_get_parameter_value(pstruct.MaxElevation,defaultMaxElevation,varargin{:});

end
    

% Validate optional parameters
validateattributes(freq,{'double'},{'positive','nonzero','finite','scalar','real','<=',10e9}, ...
    funName, 'Frequency', 1);

validateattributes(rfs,{'double'},{'positive','nonzero','finite', 'vector','real'}, ...
    funName, 'Free Space Range', 2);

validateattributes(ha,{'double'},{'positive','nonzero','scalar','finite','real'}, ...
    funName, 'Antenna Height', 3);

rngunit = validatestring(rngunit,{'km','nmi','mi', 'ft','m'}, funName,'RangeUnit');
if strcmp(rngunit,'nmi')
    rngunitOut = 'nm';
else
    rngunitOut = rngunit;
end

htunit = validatestring(htunit,{'km','nmi','mi', 'ft','m'}, funName,'HeightUnit');
if strcmp(htunit,'nmi')
    htunit = 'nm';
end

validateattributes(epsc,{'double'},{'nonempty','finite', 'scalar'}, ...
    funName, 'SurfaceDielectric');

pol = validatestring(pol,{'H', 'V'}, funName,'Polarization');

validateattributes(anpat,{'double'},{'nonempty','vector','finite','real'}, ... 
    funName, 'AntennaPattern');

validateattributes(patangles,{'double'},{'vector', '<=', 90, '>=', -90, ...
    'nrows', size(anpat, 1), 'ncols', size(anpat, 2)}, funName, 'PatternAngles');

validateattributes(tilt,{'double'},{'nonempty',  'finite', 'scalar', 'real', ...
    '>', -90,  '<', 90}, funName,'TiltAngle');

validateattributes(htsd,{'double'},{'nonempty', 'finite', 'scalar', 'real', '>=', 0}, ...
    funName, 'SurfaceRoughness');

validateattributes(vcpelmax,{'double'},{'finite', 'scalar', 'real', '<=', 90, '>', 0}, ...
    funName, 'MaxElevation');
end

function [rhoref, phiref] = reflectioncoeff(epsc, psigrz, pol)
% seareflectioncoeff  Returns reflection coefficient of a surface with
% dielectric constant of epsc 
% LAMBDA Wavelength
% EPSC Complex dielectric constant
% POL Polarization
% PSIGRZ Grazing Angle
% RHOREF Magnitude of complex reflection coefficient
% PHIREF Phase of complex reflection coefficient
if (pol ==  'H')
    % Horizontal polarization
    refc = (sin(psigrz) - sqrt(epsc - cos(psigrz).^2))./ ...
           (sin(psigrz) + sqrt(epsc - cos(psigrz).^2));
elseif (pol ==  'V')
    % Vertical polarization
    refc = (epsc*sin(psigrz) - sqrt(epsc - cos(psigrz).^2))./ ...
           (epsc*sin(psigrz) + sqrt(epsc - cos(psigrz).^2));
end
rhoref = abs(refc);
phiref = -angle(refc);
end

function epsc = seacmplxdielectric(freq)
% cmplxdielectric Returns complex dielectric constant of sea
% FREQ Carrier frequency

% Ref: NRL Report 7098 Eq 32.

freqMhz = freq/1e6;   % Convert to MHZ, used formula is for MHz
if (freqMhz <= 1500)
    epsr = 80;
    sigma = 4.3;
elseif (freqMhz > 1500) && (freqMhz <= 3000)
    epsr = 80 - 0.00733*(freqMhz - 1500);
    sigma = 4.30 + 0.00148*(freqMhz - 1500);
else   % (freqMhz > 3000) && (freqMhz <= 10000)
    epsr = 69 - 0.00243*(freqMhz - 3000);
    sigma = 6.52 + 0.001314*(freqMhz - 3000);
end
lambda = physconst('LightSpeed')/freq;
epsc = epsr - 1i*60*lambda*sigma;
end

function r = roughness(lambda, H, psigrz)
% roughness Returns roughness factor of surface assuming Gaussian height
% distribution. 
% FREQ Carrier frequency
% HSTD Std. dev of height in meters
% PSIGRZ Grazing angle

% Ref NRL Report 7098 eq 20
r = exp(-2*((2*pi*H*sin(psigrz))/lambda).^2);
end

function D = divfactor(ha, thetad)
% divfactor Returns divergence factor 
% HA Antenna height
% THETAD Elevation angle of direct ray path

% Ref: NRL Report 7098, Eq 27, 28
ae = effearthradius;
zeta = sqrt(ae/(2*ha)) * tan(thetad);
D = sqrt((1 + (2 * zeta ./ sqrt(zeta.^2 + 3)))/3);
end

function antresp = antfactor(anpat, patangles, theta)
% antresp Returns interpolated antenna response at theta based on
% provided response samples anpat and patangles
antresp = interp1(patangles, anpat, theta);
end

function delFE = pathdiff_FE(thetad, ha)
% Path difference for Flat Earth assumption
delFE = 2*ha*sin(thetad);
end

function delCE = pathdiff_CE(thetad, ha, phi)
% Path difference for Curved Earth assumption
ae = effearthradius;
delCE = sqrt(ha^2 + ae * (ae + ha) * phi.^2) .* sin(thetad + phi).^2 * 2;
end

function y = sinclocal(x)
i=find(x==0);                                                              
x(i)= 1;      % From LS: don't need this is /0 warning is off                           
y = sin(pi*x)./(pi*x);                                                     
y(i) = 1;   
end

function [apat, ath] = sincpattern()
% Default antenna pattern is sinc pattern
% Eq 38
ath = linspace(-pi/2, pi/2, 361);
HPBW = phased.internal.deg2rad(10);
k = 1.39157/sin(HPBW/2);
u = k*sin(ath);
apat = sinclocal(u/pi);
ath = phased.internal.rad2deg(ath);   % Return in degrees
end


% [EOF]
