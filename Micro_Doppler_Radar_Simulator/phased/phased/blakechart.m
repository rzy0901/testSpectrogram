function varargout = blakechart(vcp, vcpangles, varargin)
%blakechart Blake chart for radar
%   blakechart(VCP, VCPANGLES) plots vertical coverage pattern, VCP, for a
%   radar system on a Blake chart which is a range-height-angle chart. This
%   function uses the CRPL exponential model for refractive index of the
%   atmosphere. VCP can be a matrix whose columns represent individual
%   vertical coverage patterns. VCPANGLES is a column vector whose number
%   of rows is the same as number of rows of VCP. Each entry in VCPANGLES
%   specifies the elevation angle (in degrees) at which the vertical
%   coverage pattern is measured.
%
%   blakechart(VCP, VCPANGLES, RMAX, HMAX) specifies the range limit RMAX
%   and height limit HMAX for the Blake chart, respectively. The unit of
%   RMAX and HMAX are determined by the values of RangeUnit and HeightUnit
%   parameters. The default unit for both of them are 'km'.
% 
%   blakechart(...,'RangeUnit', RUNIT, 'HeightUnit', HUNIT) specifies
%   the range and height units in RUNIT and HUNIT, as one of 'km' | 'm' |
%   'mi' | 'nmi' | 'ft', where the default values for both of them are
%   'km'.
% 
%   blakechart(..., 'ScalePower', SPOW) specifies the range and height axis
%   scale power as a scalar between 0 and 1. The default value for SPOW is
%   1/4.
% 
%   blakechart(..., 'SurfaceRefractivity', NS, 'RefractionExponent', REXP)
%   specifies the surface refractivity, NS (in 1e-6/km), and the exponent
%   factor, REXP (in 1/km), for the atmospheric refraction model. The
%   default value of NS is 313 and the default value of REXP is 0.143859.
%
%   The atmospheric refraction model is given by
%
%   N(h) = 1 + NS*exp(-REXP*h)
%
%   where h is the height in kilometers and N is the index of refraction.
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
%   tilt_ang = 0;     % degrees
%   [vcd, vcd_ang] = radarvcd(freq,rng_fs,ant_height,...
%       'RangeUnit','nmi','HeightUnit','ft','AntennaPattern',pat,...
%       'PatternAngles',pat_angles,'TiltAngle',tilt_ang);
%
%   blakechart(vcd, vcd_ang, 200, 600e3, 'RangeUnit','nmi',...
%       'HeightUnit','ft')
%   
%   See also phased, radarvcd.

%   Copyright 2012-2018 The MathWorks, Inc.

%   References:
%
%   [1] Radio Ray (Radar) Range-Height-Angle Charts, L. V. Blake, NRL,
%       1968 
%   [2] Bean, B.R. & et.al. , CRPL exponential reference
%       atmosphere, US National Bureau of Standards, 1959

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

CURVECOLOR = [0.75 0.75 0.75];

% Check number of arguments is between 2 and 14 and validate
phased.internal.narginchk(2, 14, nargin);
% Check for either 2 or 4 required arguments
if (nargin == 3) % Error case
    error('phased:blakechart:badArgCombination', ...
        getString(message('phased:blakechart:badArgCombination')));
end

if ((nargin  >= 4) && isnumeric(varargin{1}) && isnumeric(varargin{2}))
    rmax = varargin{1};
    hmax = varargin{2};
    validatevcp(vcp,vcpangles);
    validateboundarymax(rmax,hmax);
    [spow, Ns, expo, rngunit, htunit] = parseInput(varargin{3:end});
else 
    % rmax, hmax not specified, so estimate!
    validatevcp(vcp, vcpangles);
    [spow, Ns, expo, rngunit, htunit] = parseInput(varargin{:});
    [rmax, hmax] = estrmaxhmax(vcp, rngunit, htunit);
end   

% Convert to SI units. All computation done in SI units
Ns = Ns*1e-6;  % To standard unit From N-unit/million, surface index
rmaxkm = rmax*unitsratio('km', rngunit);
hmaxkm = hmax*unitsratio('km', htunit);
vcpkm = vcp*unitsratio('km', rngunit);
vcpangles = phased.internal.deg2rad(vcpangles);
cond = rmaxkm >= hmaxkm;
if ~cond
    coder.internal.assert(cond, 'phased:blakechart:RngHtErr','RMAX','HMAX');
end

% Compute ellipticity factor 
E = (rmaxkm/hmaxkm)^spow;
% Define el-angle range (true angle)
theta_0 = get_elangle();           % True angle
theta_1 = true2chart(theta_0, E);  % Angle on chart

% Computation begins
[htX, htY, hts] = getheightXY();
[bdx, bdy, intrth0] = getboundingxy(htX(end,:), htY(end,:));
[rngX, rngY, rngs] = getrangeXY();
[angX, angY] = getangleXY();

% Begin plots on [0 1 0 1] axis
newplot();
set(gca, 'Color', 'none');
axis([0 1 0 1]);
hold on;

h_patch = patch(bdx, bdy, 'w', 'EdgeColor', 'k', 'EdgeAlpha', 0.5,'HitTest','off');
hasbehavior(h_patch,'legend',false);

% Plot const height loci
if (slantrange(hmaxkm, 0, Ns, expo) >= rmaxkm)
    for row = 1:size(htY,1) % Each row -> each const ht curve
        % For each height curve
        %idxc = find(htX(row,:).^2 + htY(row, :).^2 > R2, 1, 'last');
        idxc = find(htX(row,:) > cos(theta_0), 1, 'last');
        % Interpolation to make the plot smooth at edge
        htX(row, idxc) = cos(theta_0(idxc));  
        htY(row, idxc) =  interp1(htX(row,:), htY(row,:), ...
            cos(theta_0(idxc)), 'pchip');  
        htY(row, 1:idxc-1) = NaN;
    end
end

h_heightcurve = plot(htX(1:end-1, :)', htY(1:end-1,:)', 'Color', CURVECOLOR);
for m_line = 1:numel(h_heightcurve)
    hasbehavior(h_heightcurve(m_line),'legend',false);
end

% Plot const range curves
h_constrange = plot(rngX, rngY, 'Color', CURVECOLOR, 'HitTest', 'off', 'HandleVis', 'off','HitTest','off');
for m_line = 1:numel(h_constrange)
    hasbehavior(h_constrange(m_line),'legend',false);
end

% Plot angle pokes
h_angpoke = line(angX, angY, 'Color', CURVECOLOR, 'HitTest', 'off', 'HandleVis', 'off','HitTest','off');
for m_line = 1:numel(h_angpoke)
    hasbehavior(h_angpoke(m_line),'legend',false);
end
h_angpoke_end = line(angX(:, end), angY(:, end), 'Color', 'k','HitTest','off');
for m_line = 1:numel(h_angpoke_end)
    hasbehavior(h_angpoke_end(m_line),'legend',false);
end

% ---------------------------------------------------------------
%     Draw pattern factor if any specified%
[~, c] = size(vcpkm);
for idx = c:-1:1
    vcpangc = true2chart(vcpangles, E);
    [vcpX(:,idx), vcpY(:,idx)] = rngth2xy(vcpkm(:, idx), vcpangc, rmaxkm, spow, E);
end
h_vcp = plot(vcpX, vcpY);
hold off;
if nargout
    varargout{1} = h_vcp;
end

%----------- Ticks and Labels -----------------------------------
title(getString(message('phased:blakechart:Title')));

% Angle labels
tx = angX(2, :);
ty = angY(2, :);
[angth, angr] = cart2pol(tx, ty);
angr = angr + 0.02;   % Position the label slightly outside the boundary
[angx, angy] = pol2cart(angth, angr);

for idx = 1:length(angth)
    anglbl = phased.internal.rad2deg(chart2true(angth(idx), E));
    if (anglbl < 1 && anglbl ~= 0)
        text(angx(idx), angy(idx), num2str(anglbl, '%.1f'));
    else
        text(angx(idx), angy(idx), num2str(anglbl, '%.0f'));
    end
    
end
% text(1.02, 0, '0');    
% text(0, 1.02, '90');    

% Range and height axis labels
xnorm = (rngs/rmaxkm).^spow;
ynorm = (hts/hmaxkm).^spow;
if E <= 5
    tidx = [1 15 19:2:length(xnorm) length(xnorm)];
    tidy = [1 20 25 29:2:length(ynorm) length(ynorm)];
elseif  (E > 5 && E < 50)
    tidx = [1 10:2:18 19:2:length(xnorm) length(xnorm)];
    tidy = [1 10:2:18 19:2:length(ynorm) length(ynorm)];
else
    tidx = [1:2:length(xnorm) length(xnorm)];
    tidy = [1:2:length(ynorm) length(ynorm)];
end
tidx = unique(tidx);
tidy = unique(tidy);

XTicks = xnorm(tidx);
XTickLabel = rngs(tidx)*unitsratio(rngunit, 'km');

YTicks = ynorm(tidy);
YTickLabel = hts(tidy)*unitsratio(htunit, 'km');

set(gca, 'XTick', XTicks);
set(gca, 'XTickLabel', XTickLabel);
if strcmp(rngunit,'nm')
    rngunit = 'nmi';
end
xlabel(getString(message('phased:blakechart:XLabel', rngunit)));

set(gca, 'YTick', YTicks);
set(gca, 'YTickLabel', YTickLabel);
if strcmp(htunit,'nm')
    htunit = 'nmi';
end
ylabel(getString(message('phased:blakechart:YLabel', htunit)));


% Nested helper functions
% --------------------------
    function [htX, htY, hts] = getheightXY()               
        hts = get_heightsamples(hmax, E)*unitsratio('km',htunit);
        htY = zeros(length(hts), length(theta_0));
        htX = zeros(length(hts), length(theta_0));
        % Compute constant height loci, call slantrange()
        for m = 1:length(hts)
            % Calculate range values for constant height
            srange = zeros(size(theta_0));
            for k = 1:length(theta_0)
                srange(k) = slantrange(hts(m), theta_0(k), Ns, expo);                
            end    
            [htX(m, :), htY(m, :)] = rngth2xy(srange, theta_1, rmaxkm, spow, E);
     
        end
    end

    function [bdx, bdy, intrth0] = getboundingxy(bdhtx, bdhty)
        % Solves for the xy coordinate of the boundary of the chart
        if (slantrange(hmaxkm, 0, Ns, expo) >= rmaxkm)
            % When range corresponding to hmax exceeds rmax
            % Solve for point where loci for const hmax intersects with const rmax
            if (hmaxkm == rmaxkm) % fzero fails to solve in this case
                intrth0 = pi/2;
            else  % Use fzero to solve
                % make sure answer is between 0 and pi/2
                thsol = fzero(@(th)(slantrange(hmaxkm, th, Ns, expo) - rmaxkm), 0);
                intrth0 = rem(abs(thsol), pi/2);
            end
            intrth1 = atan(E*tan(intrth0));    % Angle on chart
            [~, yintr] = rngth2xy(rmaxkm, intrth1, rmaxkm, spow, E);
            bdx = [cos(linspace(0, intrth0)) bdhtx(bdhty > yintr)  0]; 
            bdy = [E*sin(linspace(0, intrth0)) bdhty(bdhty > yintr) 0];
        elseif (slantrange(hmaxkm, 0, Ns, expo) < rmaxkm)
            % Solve for angle where the const range curve meets y = 1 line
            intrth0 = abs(fzero(@(th)(E*sin(th) - 1), 0));  % Returns angle on chart
            crngth = [theta_0(theta_0 < intrth0) intrth0];
            bdx = [cos(crngth) 0 0 ];
            bdy = [E*sin(crngth) 1 0];       
        end
    end

    function [angX, angY] = getangleXY()
        % Returns X,Y coor for angle spokes
        % Set index for angle to be marked
        % theta_0 is 181 samples
        if E <= 5 
            idxa = [1 101:10:171 181];
        elseif (E > 5 && E <= 25)
            idxa = [1:10:101 111:10:151 181];
        elseif (E > 25 && E <= 50)  
            idxa = [1:10:101 111 121 136 151 171 181];
        elseif (E > 50 && E<= 75)
            idxa = [1:10:101 111 121 136 151 181];
        else
            idxa = [1:11 21:10:101 111 181];
        end
        
        if (slantrange(hmaxkm, 0, Ns, expo) >= rmaxkm)
            %theta_0s = theta_0(ds+1:ds:end-ds);    % Sampled theta_0
            theta_0s = theta_0(idxa);
            tha0 = theta_0s(theta_0s <= intrth0); % part with const rng boundary
            st = find(theta_0s > intrth0, 1, 'first');
            thb0 = theta_0(idxa(st:end));
            thsz = length(tha0) + length(thb0);
            
            angX = [zeros(1, thsz) ; cos(tha0) htX(end, idxa(st:end))];
            angY = [zeros(1, thsz) ; E*sin(tha0) htY(end, idxa(st:end))];
        elseif (slantrange(hmaxkm, 0, Ns, expo) < rmaxkm)
            theta_0s = theta_0(idxa);    % Sampled theta_0
            tha0 = theta_0s(theta_0s <= intrth0); % part with const rng boundary
            thb0 = theta_0s(theta_0s > intrth0); % part with y = 1 boundary
            thb1 = atan(E*tan(thb0));  % Convert to angles on chart for this portion
            thsz = length(tha0) + length(thb1);
            
            angX = [zeros(1, thsz); cos(tha0) 1./tan(thb1)];
            angY = [zeros(1, thsz); E*sin(tha0) ones(1, length(thb1))];
        end
    end

    function [rngX, rngY, rngs] = getrangeXY()
        % Range related calcs
        rngs = get_rangesamples(rmax, E)*unitsratio('km', rngunit);
        rngX = cos(theta_0')*((rngs/rmaxkm).^spow);
        rngY =  E*sin(theta_0')*((rngs/rmaxkm).^spow);
        if (slantrange(hmaxkm, 0, Ns, expo) >= rmaxkm)
            % Const range
            for col = 1:length(rngs)
                chtbdY = interp1(htX(end, :), htY(end, :), rngX(:, col)); % Boundary
                temp = rngY(:, col);
                temp(temp > chtbdY) = chtbdY(temp > chtbdY);
                rngY(:, col) = temp;
            end            
        elseif (slantrange(hmaxkm, 0, Ns, expo) < rmaxkm)
            % Const Range
            rngY(rngY > 1) = 1;
        end
    end
end % blakechart

% Validation Helper functions
% -----------------------------------------+
function [spow, Ns, expo, rngunit, htunit] = parseInput(varargin)

funName = 'blakechart';

% Define default values for optional inputs.
defaultRangeUnit = 'km';
defaultHeightUnit = 'km';
defaultScalePower = 1/4;
defaultRefractivity = 313;
defaultExponent = 0.143859;  % exponent for height in 1/km

par = inputParser;    % Parser object
par.addParameter('RangeUnit', defaultRangeUnit);
par.addParameter('HeightUnit', defaultHeightUnit);
par.addParameter('ScalePower',defaultScalePower);
par.addParameter('SurfaceRefractivity', defaultRefractivity);
par.addParameter('RefractionExponent', defaultExponent);

par.parse(varargin{:});
rngunit = par.Results.RangeUnit;
htunit = par.Results.HeightUnit;
spow = par.Results.ScalePower;
Ns = par.Results.SurfaceRefractivity;
expo = par.Results.RefractionExponent;

% Validate optional parameters
[rngunit, htunit] = validateoptions(rngunit, htunit, funName,spow, Ns, expo);
end

function validateboundarymax(rmax,hmax)
    validateattributes(rmax,{'double'},{'positive','nonzero','scalar','finite', 'real'}, ...
        'blakechart', 'RMAX');    
    validateattributes(hmax,{'double'},{'positive','nonzero','scalar','finite', 'real'}, ...
        'blakechart', 'HMAX');
end

function validatevcp(vcp, vcpangles)
    % vcp and vcpangles are validated in the main function
    validateattributes(vcp,{'double'},{'finite', 'real'}, ...
        'blakechart', 'VCP');
    validateattributes(vcpangles,{'double'},{'finite', 'real', ...
        'column', '<=', 90, 'ncols', 1, 'nrows', size(vcp, 1)}, ...
        'blakechart', 'VCPANGLES');
    % allow negative angle, could start from a small negative angle
end

function [rngunit, htunit] = validateoptions(rngunit, htunit, funName,spow, Ns, expo)

    validateattributes(spow,{'double'},{'positive', 'nonzero', 'scalar','finite', 'real', ...
        '>', 0, '<=', 1}, funName,'ScalePower');

    validateattributes(Ns,{'double'},{'positive', 'nonzero','scalar','finite', 'real'},...
        funName,'Refractivity');

    validateattributes(expo,{'double'},{'positive','nonzero', 'scalar', 'finite', 'real'},...
        funName,'Exponent');

    rngunit = validatestring(rngunit,{'km','nmi','mi', 'ft', 'm'}, funName,'RangeUnit');
    if strcmp(rngunit,'nmi')
        rngunit = 'nm';
    end

    htunit = validatestring(htunit,{'km','nmi','mi', 'ft', 'm'}, funName,'HeightUnit');
    if strcmp(htunit,'nmi')
        htunit = 'nm';
    end
end

% Calculation helper functions
% -----------------------------------------+
function srange = slantrange(height, theta, Ns, d)
%RANGEHEIGHTANGLE Coordinate on range-height-angle chart
%   SRANGE = RANGEHEIGHTANGLE(HEIGHT, THETA) Returns slant range (SRANGE)
%   for given HEIGHT and ANGLE on a range-height-angle chart
%   HEIGHT Height from earth surface in km
%
%   Ref: Radio Ray (Radar) Range-Height-Angle Charts, L. V. Blake, NRL,
%   1968 Eq 4, p 11

% All calculation are done in SI units
earthrad = effearthradius(0)*unitsratio('km', 'm');  % in km
etaz = atmosrefractiveindex(0, Ns, d);
srange = integral(@integrand, 0, height);

    function y = integrand(h)
        % Integrand function
        nuh = atmosrefractiveindex(h, Ns, d);
        %nuh = atmosrefractiveindex(h);
        y = nuh./ ...
            sqrt(1 - (etaz * cos(theta)./((nuh .* (1 + h/earthrad)))).^2);
    end
end

function refindex = atmosrefractiveindex(height, Ns, expconst)
%ATMOSREFRACTIVEINDEX Exponential Atmospheric Refractive Index
%   ATMOSREFRACTIVEINDEX Returns atmospheric refractive index at HEIGHT for
%   given NS, DECAYCONST and HEIGHT. This function uses CRPL exponential
%   model to compute the refractive index.
%
%   Ref: CRPL exponential reference atmosphere, Bean, B.R. and Thayer, GD,
%   US Department of the Commerce, National Bureau of Standards, 1959

%   Eq 3, p 6, NRL Report 6650

%   NS Surface refractivity of earth, scalar
%   EXPCONST Exponential constant for the refractive index model
%   HEIGHT Height from earth surface

% CRPL exponential atmospheric model   
refindex = 1 + Ns * exp(-expconst * height);
end

function [xcor, ycor] = rngth2xy(rng, th1, rmaxkm, pow, E)
% Maps range and el_theta to x, y coor
% th1 is angle on chart
parma = (rng/rmaxkm).^pow;
parmb = E * parma;
% Compute x, y coordinates
xcor = parma./sqrt(1 + (tan(th1)/E).^2);
ycor = parmb.*sqrt(1 - (xcor./parma).^2);
end

function theta_ch = true2chart(theta_tr, E)
% Converts true angles to angles on chart
% theta_ch = theta_tr if E = 1

% Eq 13, p 12, NRL Report 6650

theta_ch = atan(E*tan(theta_tr));
end

function theta_tr = chart2true(theta_ch, E)
% Converts true angles to angles on chart
% theta_ch = theta_tr if E = 1

% Eq 13, p 12, NRL Report 6650

theta_tr = atan(tan(theta_ch)/E);
end

% Sample generation helper functions
% -----------------------------------------+
function elangle = get_elangle()
% Returns properly sample elevation angles in deg for RHAchart generation
% Based on code snippet from NRL report 7098
% Needs optimization
idx = 1;
for p = 1:2
    if p == 1
        st = 1;
        ed = 100;
    elseif p == 2
        st = 11;
        ed = 91;
    end
    for k = st:ed
        n = (k - 1)*10^(p - 1);
        elev = n/10;
        elangle(idx) = phased.internal.deg2rad(elev); %#ok<AGROW>
        idx = idx + 1;
    end
end
end

function htsamp = get_heightsamples(hmax, E)
% Returns height sample points in range [0 hmax] for constant height curves
% Sample points for height scale
expe = fix(log10(hmax));
if E <= 5
    exps = expe - 3; % 3 decade scale
elseif  (E > 5 && E < 50)
    exps = expe - 2;
else
    exps = expe - 1;
end

ht0 = 10^exps;   % Starting point
ht = zeros(1,9*(expe-exps)+1);
ht(1) = ht0;
for k = exps:(expe-1)
    ht((k-exps)*9+(1:9)+1) = linspace(2*10^k, 10^(k+1), 9);
end
htx = ht(end) + 10^(k+1);
htsamp = [ht htx:0.5*10^(k+1):hmax];
if htsamp(end) < hmax
    htsamp = [htsamp hmax];
end
end

function rngsamp = get_rangesamples(rmax, E)
% Returns range sample points in range [0 rmax] for constant height curves
% Sample points for height scale
expe = fix(log10(rmax));

if E <= 5
    exps = expe - 3; % 3 decade scale
elseif  (E > 5 && E < 50)
    exps = expe - 2; % 2 decade
else
    exps = expe - 1; % 1 decade
end

rng0 = 10^exps;   % Starting point
rng = zeros(1,9*(expe-exps)+1);
rng(1) = rng0;
for k = exps:(expe-1)
    rng((k-exps)*9+(1:9)+1) = linspace(2*10^k, 10^(k+1), 9);
end
rngt = rng(end) + 0.5*10^(k+1);
rngsamp = [rng rngt:10^(k+1):rmax];
if rngsamp(end) ~= rmax
    % Make sure rmax is included in range samples
    rngsamp = [rngsamp rmax];
end
end

function [rmax, hmax] = estrmaxhmax(vcp, rngunit, htunit)
% Estimate rmax and hmax limits for blakechart based on values of vcp 
% rmax is estimated using the maximum range covered of vcp
% hmax is estimated using the mean range coverage of vcp
% This function needs to be optimized to better estimate hmax.
 rlim = max(vcp(:, end)); 
 hlim = mean(vcp(:,end))*unitsratio(htunit, rngunit);  
 rmax = nextfactor(rlim);
 hmax = nextfactor(hlim);
end

function nf = nextfactor(x)
% Returns the next 10 multiple after x
% Eg nextfactor(13) returns 20, nextfactor(123) returns 200
ex = fix(log10(x));
fac = 10^ex;  
re = rem(x, fac);
y = (x - re)/fac; 
nf = (y + 1) * fac;
end

% [EOF]
