function w = cbfweights(pos_in,ang_in,N)
%cbfweights Conventional beamformer weights for a sensor array 
%   W = cbfweights(POS,ANG) returns the narrowband conventional beamformer
%   weights of the sensor array defined in POS, for the directions
%   specified in ANG (in degrees). W is an NxM matrix, where N is the
%   number of elements in the sensor array and M is the number of
%   beamforming directions. Each column of W contains the weights of the
%   array for the corresponding direction specified in ANGLE. W is computed
%   by dividing the corresponding steering vector by N.
%
%   POS represents the locations of elements in the sensor array, specified
%   in the unit of signal wavelength. All elements in the sensor array are
%   assumed to be isotropic. POS can be either a 1xN vector, a 2xN matrix,
%   or a 3xN matrix, where N is the number of elements in the sensor array.
%   If POS is a 1xN vector, then it represents y-coordinates of elements in
%   a linear array that is along y axis. If POS is a 2xN matrix, then the
%   array lies in the yz plane. In this case, each column of POS represents
%   the [y;z] coordinates of the corresponding element. If POS is a 3xN
%   matrix, then the array has arbitrary shape. In this case, each column
%   of POS represents the [x;y;z] coordinates of the corresponding element.
%
%   ANG represents the beamforming directions. ANG can be either a 1xM
%   vector or a 2xM matrix, where M is the number of desired beamforming
%   directions. If ANG is a 2xM matrix, each column specifies the direction
%   in the space in [azimuth; elevation] form (in degrees). The azimuth
%   angle must be between -180 and 180 degrees and the elevation angle must
%   be between -90 and 90 degrees. The azimuth angle is defined in the xy
%   plane; it is the angle measured from the x axis, which is also the
%   array normal direction, toward the y axis. The elevation angle is
%   defined as the angle from the xy plane toward the z axis. If ANG is a
%   1xM vector, then it contains the azimuth angles of all the directions,
%   and the corresponding elevation angles are assumed to be 0.
%
%   W = cbfweights(POS,ANG,N) returns the weights obtained with N-bit phase
%   shifters. Specifying N as 0 indicates no quantization
%   effect in phase shifters.
%
%   % Example:
%   %   Calculate the beamforming weights of a 10-element half-wavelength
%   %   spacing ULA in the direction of 30 and 45 degrees azimuth. 
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   w = cbfweights(elementPos,[30 45]);
%
%   See also phased, steervec, mvdrweights, lcmvweights,
%   phased.PhaseShiftBeamformer

%   Copyright 2012-2015 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, pp. 53, Wiley, 2002

%#codegen

phased.internal.narginchk(2,3,nargin);
if nargin < 3
    N = 0;
end

[pos,ang,N_ele] = parseInput(pos_in,ang_in,N);

w = phased.internal.steeringvec(pos,1,1,ang,N)./N_ele;


%---------------------------------
function [pos,ang,N_ele] = parseInput(pos_in,ang_in,N)

eml_assert_no_varsize(1:nargin,pos_in,ang_in,N);
if size(pos_in,1) == 1
    pos = [zeros(1,size(pos_in,2));pos_in;zeros(1,size(pos_in,2))];
elseif size(pos_in,1) == 2;
    pos = [zeros(1,size(pos_in,2));pos_in];
else
    pos = pos_in;
end

sigdatatypes.validate3DCartCoord(pos,'cbfweights','POS');

N_ele = size(pos,2);

if size(ang_in,1) == 1
    ang = [ang_in;zeros(1,size(ang_in,2))];
else
    ang = ang_in;
end
sigdatatypes.validateAzElAngle(ang,'cbfweights','ANG');

validateattributes(N,{'double'},...
    {'scalar','integer','nonnegative'},'cbfweights','N');



% [EOF]
