function varargout = polsignature(sMat,type,epsilon,tau)
%polsignature   Polarization signature
%   RESP = polsignature(RCSMAT) returns the normalized synthesized radar
%   cross section polarization signature, RESP (in square meters), of the
%   scattering cross section, RCSMAT. RCSMAT must be a 2x2 matrix in the
%   form of [rcs_hh rcs_hv;rcs_vh rcs_vv] where rcs_hv specifies the radar
%   cross section when the transmit antenna is vertically polarized and the
%   receiving antenna is horizontally polarized.
%
%   The signature, RESP, is represented as a function of the transmit
%   antenna polarization. The transmitting antenna polarization is
%   characterized by the ellipticity angle and the tilt angle. RESP is a
%   matrix whose elements correspond to the signature at given ellipticity
%   angle and tilt angle pair.
%
%   RESP = polsignature(RSMAT,TYPE) specifies the polarization signature
%   type as one of 'c' | 'x', where 'c' indicates co-polarized signature
%   and 'x' indicates cross-polarized signature. The default value of this
%   parameter is 'c'.
%
%   RESP = polsignature(RSMAT,TYPE,EPSILON) specifies the transmit antenna
%   polarization's ellipticity angle (in degrees) as a length-Q vector.
%   EPSILON must be between -45 and 45 degrees. The default value of
%   EPSILON is from -45 degrees to 45 degrees, with a step of 1 degree.
%   
%   RESP = polsignature(RSMAT,TYPE,EPSILON,TAU) specifies the transmit
%   antenna polarization's tilt angle (in degrees) as a length-P vector.
%   TAU must be between -90 and 90 degrees. The default value of TAU is
%   from -45 degrees to 45 degrees, with a step of 1 degree.
%   
%   The signature, RESP, is represented as a function of the transmitting
%   antenna polarization. The transmitting antenna polarization is
%   characterized by the ellipticity angle, EPSILON and the tilt angle,
%   TAU. RESP is a PxQ matrix whose elements correspond to the signature at
%   given ellipticity angle and tilt angle pair.
%
%   polsignature(...) plots the polarization signature.
%
%   % Example:
%   %   Plot the co-polarization signature of the a dihedral, whose RCS 
%   %   matrix is given by [0 1;1 0].
%
%   polsignature([0 1;1 0],'c')
%
%   See also phased, polellip, polloss, stokes.

%   Copyright 2012-2016 The MathWorks, Inc.

%   References
%   [1] Fawwaz Ulaby & Charles Elachi, Radar Polarimetry for Geoscience
%       Applications, Artech House, 1990
%   [2] Jong-sen Lee & Eric Pottier, Polarimetric Radar Imaging: From
%       Basics to Applications, CRC Press, 2009

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(1,4,nargin);

if nargin < 2
    type = 'c';
end
    
type = validatestring(type,{'c','x'},'polarizationResponse','Type');
% 'c': col-pol, 'x': cross-pol

validateattributes(sMat,{'double'},{'finite',...
    'size',[2 2]},'polsignature','RCSMAT');
% sMat: 2x2 scattering matrix, BSA

if nargin < 3
    epsilon = -45:45;
end
sigdatatypes.validateAngle(epsilon,'polsignature','EPSILON',...
    {'vector','>=',-45,'<=',45});

if nargin < 4
    tau = -90:90;
end
sigdatatypes.validateAngle(tau,'polsignature','TAU',...
    {'vector','>=',-90,'<=',90});

eml_assert_no_varsize(1:nargin,sMat,type,epsilon,tau);

[egrid,tgrid] = meshgrid(epsilon,tau);

Ei = ellipseAnglesToField(tgrid(:).',egrid(:).');
Es = sMat*Ei;

NumPts = numel(tau)*numel(epsilon);
pResp = complex(zeros(NumPts,1));

if strncmp(type,'c',1)
    % co-pol receiving field, same polarization as incident
    Er = Ei;
    respType = 'Co-Pol';
else
    % cross-pol receiving field, orthogonal polarization of incident
    % in terms of polarization ratio, Pr = -1/conj(Pi)
    Er = conj(flipud(Ei));
    Er(1,:) = -Er(1,:);
    respType = 'Cross-Pol';
end

for m = 1:NumPts
    pResp(m) = Er(:,m).'*Es(:,m);
end
pResp = abs(pResp);  % return magnitude, converted to power if needed in plot
pResp = pResp/max(pResp);
pResp = reshape(pResp,numel(tau),[]);

if nargout
    varargout{1} = pResp;
    varargout{2} = epsilon;
    varargout{3} = tau;
else
    cond = isempty(coder.target);
    if ~cond
        coder.internal.assert(cond,'phased:phased:invalidCodegenOutput','polsignature');
    end
    p = phased.internal.PolarizationResponse('OrientationAngle',tau,...
        'EllipticityAngle',epsilon,'Pattern',pResp,'ResponseType',respType);
    plot(p,'Units','Power');
end

function eField = ellipseAnglesToField(tau, epsilon)
%ELLIPSEANGLESTOFIELD Convert ellipse angles to field vector

% eField vector in the form of [Ex; Ey]
% tau: tilt angle [-pi/2, pi/2]  in degrees
% epsilon: ellipticity [0, pi/4] in degrees

% tau, epsilon must be row vectors

% if tau > 45
%     tau = tau - pi/2;
% elseif tau < -45
%     tau = tau + pi/2;
% end

tau = phased.internal.deg2rad(tau);
epsilon = phased.internal.deg2rad(epsilon);

% return normalized field, absolute phase is omitted (Mott07, eq 1.30)
eField = [cos(tau).*cos(epsilon)-1i.*sin(tau).*sin(epsilon);...
    sin(tau).*cos(epsilon)+1i.*cos(tau).*sin(epsilon)];


% [EOF]
