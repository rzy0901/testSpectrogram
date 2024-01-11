function rho = polloss(fv_t,fv_r,pos_r,ax_r,pos_t,ax_t)
%polloss    Polarization loss
%   RHO = polloss(FV_T,FV_R) returns the polarization loss, RHO (in dB),
%   between the polarized field of a transmitter, FV_T, and the polarized
%   field of a receiver, FV_R. Both FV_T and FV_R are length-2 column
%   vectors in the form of [H;V] where H is the horizontal polarization and
%   V is the vertical polarization. FV_T is measured in the transmitter's
%   local coordinate system, at the direction toward the receiver. FV_R is
%   measured in the receiver's local coordinate system, at the direction
%   toward the transmitter.
%
%   RHO = polloss(FV_T,FV_R,POS_R) specifies the position of the receiver
%   (in meters), POS_R, as a length-3 column vector. The position is
%   specified in the [x;y;z] form in the global coordinate system.
%
%   RHO = polloss(FV_T,FV_R,POS_R,AXES_R) specifies the orthonormal axes,
%   AXES_R, that defines the receiver's local coordinate system as a 3x3
%   matrix. Each column defines the x, y, and z axes respectively.
%
%   RHO = polloss(FV_T,FV_R,POS_R,AXES_R,POS_T) specifies the position of
%   the transmitter (in meters), POS_T, as a length-3 column vector. The
%   position is specified in the [x;y;z] form in the global coordinate
%   system.
%
%   RHO = polloss(FV_T,FV_R,POS_R,AXES_R,POS_T,AXES_T) specifies the
%   orthonormal axes, AXES_T, that defines the transmitter's local
%   coordinate system as a 3x3 matrix. Each column defines the x, y, and z
%   axes respectively.
%
%   % Example:
%   %   Determine the polarization loss between a pair of transmitting and
%   %   receiving antenna. Assuming the transmitting antenna is a vertical 
%   %   dipole at origin. The receiving antenna is also a vertical dipole, 
%   %   located 100 meters away along the x-axis. However, the receiving
%   %   dipole is rotated 45 degrees around x-axis.
%
%   fv_t = [0;1]; %vertical polarization
%   fv_r = [0;1]; %vertical polarization
%   pos_r = [100;0;0];
%   ax_r = rotx(45)*azelaxes(0,0); %dipole rotated around x axis
%   rho = polloss(fv_t,fv_r,pos_r,ax_r)
%
%   See also phased, polellip, stokes.

%   Copyright 2012 The MathWorks, Inc.

%   References
%   [1] Harold Mott, Polarization in Antennas and Radar, John Wiley & Sons,
%   1986

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,6,nargin);
if nargin < 6
    ax_t = eye(3);
end

if nargin < 5
    pos_t = [0;0;0];
end

if nargin < 4
    ax_r = eye(3);
end

if nargin < 3
    pos_r = [0;0;0];
end
eml_assert_no_varsize(1:nargin,fv_t,fv_r,pos_r,ax_r,pos_t,ax_t);
validateattributes(fv_t,{'double'},{'finite','2d','nrows',2},...
    'polloss','FV_T');

validateattributes(fv_r,{'double'},{'finite','2d','nrows',2},...
    'polloss','FV_R');

sigdatatypes.validate3DCartCoord(pos_r,'polloss','POS_R',{'ncols',1});
sigdatatypes.validate3DCartCoord(ax_r,'polloss','AXES_R',{'ncols',3});

sigdatatypes.validate3DCartCoord(pos_t,'polloss','POS_T',{'ncols',1});
sigdatatypes.validate3DCartCoord(ax_t,'polloss','AXES_T',{'ncols',3});

% convert field in tx to global coordinates

rcoord_t = global2localcoord(pos_r,'rs',pos_t,ax_t);
fv_t_local = sph2cartvec([fv_t/norm(fv_t);0],...    % no radial components
    rcoord_t(1),rcoord_t(2)); 
fv_t_global = phased.internal.local2globalvec(fv_t_local,ax_t);

% convert filed in rx to global coordinates

tcoord_r = global2localcoord(pos_t,'rs',pos_r,ax_r);
fv_r_local = sph2cartvec([fv_r/norm(fv_r);0],...
    tcoord_r(1),tcoord_r(2));
fv_r_global = phased.internal.local2globalvec(fv_r_local,ax_r);

rho = abs(fv_t_global.'*fv_r_global).^2;   % same coord, no conjugate
rho = -pow2db(rho);







% [EOF]
