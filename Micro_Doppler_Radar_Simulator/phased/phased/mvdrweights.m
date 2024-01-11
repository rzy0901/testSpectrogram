function w = mvdrweights(pos_in,ang_in,cov,N)
%mvdrweights   MVDR beamformer weights for a sensor array
%   W = mvdrweights(POS,ANG,COV) returns the narrowband minimum variance
%   distortionless response (MVDR) beamformer weights, W, of the sensor
%   array defined in POS, for the directions specified in ANG (in degrees)
%   and the covariance matrix, COV. W is an NxM matrix, where N is the
%   number of elements in the sensor array and M is the number of
%   beamforming directions. Each column of W contains the weights of the
%   array for the corresponding direction specified in ANG.
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
%   COV is an NxN matrix representing the covariance matrix across the N
%   sensor elements for all incoming signals.
%
%   W = mvdrweights(POS,ANG,COV,N) returns the weights obtained with N-bit
%   phase shifters. Specifying N as 0 indicates no quantization effect in
%   phase shifters.
%
%   % Example:
%   %   Calculate the MVDR beamforming weights of a 10-element 
%   %   half-wavelength spacing ULA in the direction of 30 and 45 degrees
%   %   azimuth. Assume there are two signals coming from azimuth 60 and
%   %   -60 degrees, respectively. The noise is white across all elements
%   %   and the SNR is 10 dB.
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   Sn  = sensorcov(elementPos,[-60 60],db2pow(-10));
%   w = mvdrweights(elementPos,[30 45],Sn);
%
%   plotangl = -90:90;
%   vv = steervec(elementPos,plotangl);
%   plot(plotangl,mag2db(abs(w'*vv)))
%   grid on
%   xlabel('Azimuth Angle (degrees)');
%   ylabel('Normalized Power (dB)');
%   legend('30 deg','45 deg');
%
%   See also phased, steervec, phased.MVDRBeamformer, cbfweights,
%   lcmvweights

%   Copyright 2012 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, pp. 442, Wiley, 2002


%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,4,nargin);
if nargin < 4
    N = 0;
end

[pos,ang] = parseInput(pos_in,ang_in,cov,N);

Constraints =  phased.internal.steeringvec(pos,1,1,ang,N);

convinvAS = cov\Constraints;

w = bsxfun(@rdivide,convinvAS,sum(conj(Constraints).*convinvAS));

%---------------------------------
function [pos,ang,N_elem,N_angl] = parseInput(pos_in,ang_in,scov_in,N)
eml_assert_no_varsize(1:nargin, pos_in,ang_in,scov_in,N);
if size(pos_in,1) == 1
    pos = [zeros(1,size(pos_in,2));pos_in;zeros(1,size(pos_in,2))];
elseif size(pos_in,1) == 2;
    pos = [zeros(1,size(pos_in,2));pos_in];
else
    pos = pos_in;
end

sigdatatypes.validate3DCartCoord(pos,'mvdrweights','POS');
N_elem = size(pos,2);

if size(ang_in,1) == 1
    ang = [ang_in;zeros(1,size(ang_in,2))];
else
    ang = ang_in;
end
sigdatatypes.validateAzElAngle(ang,'mvdrweights','ANG');
N_angl = size(ang,2);


validateattributes(scov_in,{'double'},{'finite','size',[N_elem N_elem]},...
        'mvdrweights','COV');
tol = 10*eps(max(abs(diag(scov_in))));   % based on stats cholcov
cond = any(any(abs(scov_in - scov_in') > tol));
if cond
    coder.internal.errorIf(cond,...
                       'phased:phased:notHermitian','COV');
end

% Check for positive semi definite
[~,sED] = eig((scov_in+scov_in')/2);  % ensure Hermitian
sED = diag(sED);
tol = eps(max(abs(sED))); % based on stats cholcov
sED(abs(sED)<=tol)=0;
cond = any(sED<0);
if cond
    coder.internal.errorIf(cond,...
                       'phased:phased:notPositiveSemiDefinite','SCOV');
end

validateattributes(N,{'double'},...
    {'scalar','integer','nonnegative'},'mvdrweights','N');


% [EOF]
