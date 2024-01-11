function g = taylortaperc(pos,D,varargin)
%TAYLORTAPERC   Taylor nbar taper for circular aperture
%   W = TAYLORTAPERC(POS,D) returns the value of Taylor nbar taper at
%   positions given by POS in a circular aperture with a diameter of D. POS
%   can be either a 2-row matrix or a a 3-row matrix. If POS is a 2-row
%   matrix, each column in the matrix specifies the position of the
%   corresponding element, in [x;y] form, of the array in the array plane.
%   If POS is a 3-row matrix, its columns represent the positions of
%   sensors, in [x;y;z] format, in the array. Note that POS should use the
%   same unit as D. W is a column vector containing the resulting taper.
%
%   W = TAYLORTAPERC(POS,D,NBAR) specifies the number of nearly
%   constant-level sidelobes adjacent to the mainlobe, NBAR, as a scalar.
%   NBAR must be an integer greater than or equal to one. The default value
%   of NBAR is 4.
%
%   W = TAYLORTAPERC(POS,D,NBAR,SLL) specifies the maximum sidelobe level,
%   SLL (in dB), relative to the mainlobe peak as a negative value. The
%   default value of SLL is -30.
%
%   W = TAYLORTAPERC(POS,D,NBAR,SLL,CPOS) specifies the center of array in
%   CPOS. The format of CPOS is the same as used in POS. The default center
%   of the array is assumed to be at the center of all specified positions.
%
%   NBAR should satisfy NBAR >= 2*A^2+0.5, where A is equal to
%   acosh(10^(-SLL/20))/pi, otherwise the sidelobe level specified is not
%   guaranteed.
% 
%   % Example:
%   %   Apply taylor taper to a circular aperture array. Assume nbar is 2,
%   %   and sidelobe level is -25 dB. Plot the array power pattern at 300 
%   %   MHz.
%
%   % create circular aperture
%   radius = 5; dist = 0.5;
%   pos = getElementPosition(phased.URA(2*radius/dist,dist));
%
%   % remove all elements outside the circle
%   pos(:,sum(pos.^2) > radius^2) = [];
%   ant = phased.ConformalArray('ElementPosition',pos);
%
%   % apply taper and plot pattern
%   ant.Taper = taylortaperc(pos,2*radius,2,-25,[0;0;0]);
%   pattern(ant,3e8,-1:0.01:1,0,'CoordinateSystem','uv','Type','powerdb')
%
%   See also phased, TAYLORWIN.

%   Copyright 2016 The MathWorks, Inc.

%   Reference
%   [1] Taylor, T., Design of Circular Aperture for Narrow Beamwidth and 
%       Low Sidelobes. IRE Trans. on Antennas and Propagation, Jan, 1960
%   [2] Hansen, R. C., Tables of Taylor Distributions for Circular Aperture
%       Antennas. IRE Trans. on Antenna and Propagation, Jan, 1960
%   [3] Hansen, R. C., Array Pattern control and Synthesis, Proceedings of
%       the IEEE, January, 1992

%#codegen

narginchk(2,5);
if size(pos,1) == 2
    validateattributes(pos,{'double'},{'2d','finite','nonnan','nonempty',...
        'real','nrows',2},'taylortaperc','pos');
else
    sigdatatypes.validate3DCartCoord(pos,'taylortaperc','pos');
end

sigdatatypes.validateDistance(D,'taylortaperc','D',{'scalar','positive','finite'});

if nargin < 3 || isempty(varargin{1})
    nbar = 4;
else
    nbar = sigdatatypes.validateIndex(varargin{1},'taylortaperc','NBAR',{'scalar','<=',10});
end

if nargin < 4 || isempty(varargin{2})
    sll = -30;
else
    sll = varargin{2};
    validateattributes(sll,{'double'},{'scalar','real','<',0},...
        'taylortaperc','SLL');
end

% Center position
if size(pos,1) == 2
    if nargin < 5 || isempty(varargin{3})
        posc = mean(pos,2);
    else
        posc = varargin{3};
        validateattributes(posc,{'double'},{'real','nonnan','nonempty',...
            'finite','column','nrows',2},'taylornbar','POSC');
    end

    p = hypot(pos(1,:)-posc(1,:),pos(2,:)-posc(2,:))/D;
else
    if nargin < 5 || isempty(varargin{3})
        posc = mean(pos,2);
    else
        posc = varargin{3};
        sigdatatypes.validate3DCartCoord(posc,'taylornbar','POSC',...
            {'column'});
    end

    p = hypot(hypot(pos(1,:)-posc(1,:),pos(2,:)-posc(2,:)),pos(3,:)-posc(3,:))/D;
end

cond = any(p>0.5);
if cond
    coder.internal.errorIf(cond,'phased:system:array:ElementPositionOutOfBound');
end

n = 1:nbar-1;
A = acosh(10^(-sll/20))/pi;
sigma = besselj1zeros(nbar)/sqrt(A^2+(nbar-1/2)^2);

zn = sigma*sqrt(A^2+(n-1/2).^2);
un = besselj1zeros(n);
F_num = prod(1-bsxfun(@rdivide,un.^2,zn(:).^2),1);
temp = 1 - bsxfun(@rdivide,un.^2,un(:).^2);
temp(1:nbar:end) = 1;
F_dec = prod(temp,1);
F = [1 -besselj(0,pi*un).*F_num./F_dec];

xi = 2*p(:).';  % position normalized by radius
un = besselj1zeros([0 n]);
g = 2/(pi^2)*sum(bsxfun(@rdivide,...
    bsxfun(@times,F(:),besselj(0,pi*un(:)*xi)),...
    besselj(0,pi*un(:)).^2));
g = g(:);


function u = besselj1zeros(n)
%besselj1zeros
%   U = besselj1zeros(N) return Nth zeros for Bessel function of the first
%   kind, J_n(x), where n=1.

%   [1] Taylor, T., Design of Circular Aperture for Narrow Beamwidth and 
%       Low Sidelobes. IRE Trans. on Antennas and Propagation, Jan, 1960

uparam = [0 1.2196699 2.2331306 3.2383154 4.2410628 5.2427643 6.2439216 ...
    7.2447598 8.2453948 9.2458927 10.2462933];

u = uparam(n+1);
