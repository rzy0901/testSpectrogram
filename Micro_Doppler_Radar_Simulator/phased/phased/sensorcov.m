function Rx = sensorcov(pos_in,ang_in,ncov_in,scov_in)
%sensorcov Covariance matrix of received signal
%   R = sensorcov(POS, ANG) returns the received covariance matrix, R, of
%   narrowband plane wave signals at a sensor array. POS represents the
%   element positions of the sensor array. ANG represents the incoming
%   directions of each plane wave signal.
%
%   POS represents the locations of elements in the sensor array, specified
%   in the unit of signal wavelength. All elements in the sensor array are
%   assumed to be isotropic. POS can be either a 1xN vector, a 2xN matrix,
%   or a 3xN matrix, where N is the number of elements in the sensor array.
%   If POS is a 1xN vector, then it represents y-coordinates of elements in
%   a linear array that is along the y axis. If POS is a 2xN matrix, then
%   the array lies in the yz plane. In this case, each column of POS
%   represents the [y;z] coordinates of the corresponding element. If POS
%   is a 3xN matrix, then the array has arbitrary shape. In this case, each
%   column of POS represents the [x;y;z] coordinates of the corresponding
%   element.
%
%   ANG represents the directions of the incoming signals. ANG can be
%   either a 1xM vector or a 2xM matrix, where M is the number of incoming
%   signals. If ANG is a 2xM matrix, each column specifies the direction in
%   the space in [azimuth; elevation] form (in degrees). The azimuth angle
%   must be between -180 and 180 degrees and the elevation angle must be
%   between -90 and 90 degrees. The azimuth angle is defined in the xy
%   plane; it is the angle measured from the x axis, which is also the
%   array normal direction, toward the y axis. The elevation angle is
%   defined as the angle from the xy plane toward the z axis. If ANG is a
%   1xM vector, then it contains the azimuth angles of all the directions,
%   and the corresponding elevation angles are assumed to be 0.
%
%   R is an NxN matrix covariance matrix of the received signal at the
%   sensor array.
%
%   R = sensorcov(POS,ANG,NCOV) specifies the noise characteristics in NCOV
%   as either a nonnegative scalar, a 1xN vector, or an NxN positive
%   definite matrix. If NCOV is a nonnegative scalar, it represents the
%   noise power (in watts) of the white noise across all sensor elements.
%   If NCOV is a 1xN vector, it represents the noise power (in watts) of
%   each receiving sensor. If NCOV is an NxN matrix, then it represents the
%   covariance matrix for the noise across all sensor elements. The default
%   value of NCOV is 0.
%
%   R = sensorcov(POS,ANG,NCOV,SCOV) also specifies the incoming signal
%   characteristics in SCOV as either a positive scalar, a 1xM vector, or
%   an MxM positive semidefinite matrix. If SCOV is a scalar, it represents
%   the power (in watts) of the incoming signals. All incoming signals are
%   assumed to be uncorrelated and share the same power level. If SCOV is a
%   1xM matrix, then it represents the power of individual incoming
%   signals. However, all incoming signals are still uncorrelated to each
%   other. If SCOV is an MxM matrix, then it represents the covariance
%   matrix for all incoming signals. The default value of SCOV is 1.
%
%
%   % Example:
%   %   Calculate the covariance matrix of a received signal at an 
%   %   10-element half-wavelength spacing ULA. Assume there are two 
%   %   signals coming from azimuth 30 and 60 degrees, respectively. The 
%   %   noise is white across all elements and the SNR is 10 dB. Using
%   %   the covariance matrix, estimate the direction of arrival of the
%   %   received signals.
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%
%   Rx = sensorcov(elementPos,[30 60],db2pow(-10));
%   doa = rootmusicdoa(Rx,2)
%
%   See also phased, phased.SteeringVector, sensorsig

%   Copyright 2012-2016 The MathWorks, Inc.

%   References
%   [1] Harry Van Trees, Optimum Array Processing, Wiley-Interscience, 2002

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(2,4,nargin);

if nargin < 3
    ncov_in = 0;
end

if nargin < 4
    scov_in = 1;
end

[pos,ang,ncov,scov] = parseInput(pos_in,ang_in,ncov_in,scov_in);

sv = phased.internal.steeringvec(pos,1,1,ang);
Rx = sv*scov*sv' + ncov;
Rx = (Rx+Rx')/2;  % ensure Hermitian

%---------------------------------
function [pos,ang,ncov,scov] = parseInput(pos_in,ang_in,ncov_in,scov_in)
eml_assert_no_varsize(1:nargin,pos_in,ang_in,ncov_in,scov_in);
if size(pos_in,1) == 1
    pos = [zeros(1,size(pos_in,2));pos_in;zeros(1,size(pos_in,2))];
elseif size(pos_in,1) == 2;
    pos = [zeros(1,size(pos_in,2));pos_in];
else
    pos = pos_in;
end
sigdatatypes.validate3DCartCoord(pos,'sensorcov','POS');

N_elem = size(pos,2);

if size(ang_in,1) == 1
    ang = [ang_in;zeros(1,size(ang_in,2))];
else
    ang = ang_in;
end
sigdatatypes.validateAzElAngle(ang,'sensorcov','ANG');

N_ang = size(ang,2);

if isscalar(ncov_in)
    validateattributes(ncov_in,{'double'},{'scalar','real',...
        'finite','nonnegative'},'sensorcov','NCOV');
    ncov = ncov_in*eye(N_elem);
elseif isrow(ncov_in)
    validateattributes(ncov_in,{'double'},{'real',...
        'finite','positive','size',[1 N_elem]},'sensorcov','NCOV');
    ncov = diag(ncov_in);
else
    validateattributes(ncov_in,{'double'},{'finite','size',[N_elem N_elem]},...
        'sensorcov','NCOV');
    tol = 10*eps(max(abs(diag(ncov_in))));  % based on stats cholcov
    cond = any(any(abs(ncov_in - ncov_in') > tol));
    if cond
        coder.internal.errorIf(cond,...
             'phased:phased:notHermitian','NCOV');
    end
    % Check for positive definite
    [~,nFlagPD] = chol(ncov_in);
    cond = logical(nFlagPD);
    if cond
        coder.internal.errorIf(cond,...
             'phased:phased:notPositiveDefinite','NCOV');
    end
    ncov = ncov_in;
end



if isscalar(scov_in)
    validateattributes(scov_in,{'double'},{'scalar','real',...
        'finite','positive'},'sensorcov','SCOV');
    scov = scov_in*eye(N_ang);
elseif isrow(scov_in)
    validateattributes(scov_in,{'double'},{'real',...
        'finite','positive','size',[1 N_ang]},'sensorcov','SCOV');
    scov = diag(scov_in);
else
    validateattributes(scov_in,{'double'},{'finite','size',[N_ang N_ang]},...
        'sensorcov','SCOV');
    tol = 10*eps(max(abs(diag(scov_in))));   % based on stats cholcov
    cond = any(any(abs(scov_in - scov_in') > tol));
    if cond
        coder.internal.errorIf(cond,...
             'phased:phased:notHermitian','SCOV');
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

    scov = scov_in;
end



% [EOF]
