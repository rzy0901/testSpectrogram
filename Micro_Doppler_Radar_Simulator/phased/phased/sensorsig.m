function [xOut, Rx_t, Rx] = sensorsig(pos_in,Nsnapshots,ang_in,varargin)
%sensorsig Received signal at sensor array
%   X = sensorsig(POS,NS,ANG) simulates the received narrowband plane wave
%   signals at a sensor array. POS represents the element positions of the
%   sensor array. NS is a positive integer scalar indicating the number of
%   snapshots of the simulated signal. ANG represents the incoming
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
%   X is an NSxN matrix that contains the received signal at the sensor
%   array. Each column of X represents the received signal at the
%   corresponding element.
%
%   The input signals are assumed to be constant-modulus signals with
%   random phases.
%
%   X = sensorsig(POS,NS,ANG,NCOV) specifies the noise characteristics in
%   NCOV as either a nonnegative scalar, a 1xN vector, or an NxN positive
%   definite matrix. If NCOV is a nonnegative scalar, it represents the
%   noise power (in watts) of the white noise across all sensor elements.
%   If NCOV is a 1xN vector, it represents the noise power (in watts) of
%   each receiving sensor. If NCOV is an NxN matrix, then it represents the
%   covariance matrix for the noise across all sensor elements. The default
%   value of NCOV is 0.
%
%   X = sensorsig(POS,NS,ANG,NCOV,SCOV) also specifies the incoming signal
%   characteristics in SCOV as either a positive scalar, a 1xM vector, or
%   an MxM positive semidefinite matrix. If SCOV is a scalar, it represents
%   the power (in watts) of the incoming signals. All incoming signals are
%   assumed to be uncorrelated and share the same power level. If SCOV is a
%   1xM matrix, then it represents the power of individual incoming
%   signals. However, all incoming signals are still uncorrelated to each
%   other. If SCOV is an MxM matrix, then it represents the covariance
%   matrix for all incoming signals. The default value of SCOV is 1.
%
%   X = sensorsig(...,'Taper',TAPER) specifies the array taper as either a
%   scalar or an Nx1 vector, where N is the number of elements in the
%   array. If TAPER is a scalar, all elements in the array use the same
%   taper. If TAPER is a vector, its entries specify the taper applied to
%   the corresponding elements in the array.
%
%   [X,RT] = sensorsig(...) returns the theoretical covariance matrix of
%   the received signal in an NxN matrix, RT.
%
%   [X,RT,R] = sensorsig(..) also returns the sample covariance matrix of
%   the received signal in an NxN matrix, R. R is derived from X.
%
%   % Example:
%   %   Simulate the received signal at an 8-element half-wavelength
%   %   spacing ULA. Assume there are two signals coming from azimuth 30
%   %   and 60 degrees, respectively. The noise is white across all 
%   %   elements and the SNR is 15 dB. Simulate 100 snapshots and use the
%   %   data to estimate the direction of arrivals. Compare to the 
%   %   estimated direction of arrival when the 7th element is attenuated 
%   %   by 20 dB. 
%
%   N = 8;     % Elements in array
%   d = 0.5;   % sensor spacing half wavelength
%   k = 100;   % number of snapshots
%   elementPos = (0:N-1)*d;
%
%   % simulate signal
%   [x,~,RS] = sensorsig(elementPos,k,[55 60],db2pow(-15));
%   
%   % Use root MUSIC to get the direction of arrival estimation
%   rootmusicdoa(RS,2)
%   
%   x(:,7) = x(:,7)*db2mag(-20); % attenuate the 7th element
%   RS = (x.'*conj(x))/k; % Calculate the sample covariance matrix
%   rootmusicdoa(RS,2)
%
%   See also phased, phased.SteeringVector, sensorcov

%   Copyright 2012-2018 The MathWorks, Inc.

%   References
%   [1] Harry Van Trees, Optimum Array Processing, Wiley-Interscience, 2002

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,7,nargin);

[pos,ang,ncov,scov,taper,nonoiseflag,N_ele] = parseInput(...
    pos_in,Nsnapshots,ang_in,varargin{:});

if isscalar(scov)
    sU = sqrt(scov);
else
    %chol needs var size inputs for codegen
    scov_chol = scov;
    coder.varsize('scov_chol',size(scov));
    [sUv, sFlagPD] = chol(scov_chol);
    if sFlagPD
        [sEV,sED] = eig((scov+scov')/2);  % ensure Hermitian
        sED = diag(sED);
        tol = eps(max(abs(sED))); % based on stats cholcov
        sED(abs(sED)<=tol)=0;
        cond = any(sED<0);
        if cond
            coder.internal.errorIf(cond,...
                'phased:phased:notPositiveSemiDefinite','SCOV');
        end
        posidx = (sED>0);
        sUt = sEV(:,posidx)*diag(sqrt(sED(posidx)));
        sU = sUt';
    else
        N_ang = size(ang,2);
        sUf = sUv(1:N_ang,1:N_ang);
        sU = sUf;
    end
end
x_in = exp(1i*2*pi*rand(Nsnapshots,size(sU,1)))*sU;
sv = phased.internal.steeringvec(pos,1,1,ang);
sv = bsxfun(@times,sv,taper);

x = x_in*sv.';

if ~nonoiseflag
    ncov_chol = ncov;
    coder.varsize('ncov_chol',size(ncov));

    [nU,nFlagPD] = chol(ncov_chol);

    cond = logical(nFlagPD);
    if cond
        coder.internal.errorIf(cond,...
             'phased:phased:notPositiveDefinite','NCOV');
    end

    noise = phased.internal.cwgn(1,Nsnapshots,N_ele)*nU;
    xOut = x + noise;
else
    xOut = x;
end

if nargout > 1
    Rx_t = sv*scov*sv' + ncov;
    Rx_t = (Rx_t+Rx_t')/2;  % ensure Hermitian
end

if nargout > 2
    if Nsnapshots < N_ele
        coder.internal.warning('phased:phased:notEnoughSnapshots','R','NS',sprintf('%d',N_ele));
    end
    Rx = xOut.'*conj(xOut)/Nsnapshots;
    Rx = (Rx+Rx')/2;  % ensure Hermitian
end

%---------------------------------
function [pos,ang,ncov,scov,taper,nonoiseflag,N_ele] = ...
    parseInput(pos_in,Nsnapshots,ang_in,varargin)

if coder.target('MATLAB')
    p = inputParser;
    addOptional(p,'ncov',0);
    addOptional(p,'scov',1);

    addParameter(p,'taper',1);

    parse(p,varargin{:});

    ncov_in = p.Results.ncov;
    scov_in = p.Results.scov;
    taper_in = p.Results.taper;

else
    if nargin > 3
        if ischar(varargin{1})
            ncov_in = 0;
            scov_in = 1;
            taper_in = eml_parse_param(varargin{:});
        else
            if isempty(varargin{1})
                ncov_in = 0;
            else
                ncov_in = varargin{1};
            end
            if nargin > 4
                if ischar(varargin{2})
                    scov_in = 1;
                    taper_in = eml_parse_param(varargin{2:end});
                else
                    if isempty(varargin{2})
                        scov_in = 1;
                    else
                        scov_in = varargin{2};
                    end
                    taper_in = eml_parse_param(varargin{3:end});
                end
            else
                scov_in = 1;
                taper_in = 1;
            end
        end
    else
        ncov_in = 0;
        scov_in = 1;
        taper_in = 1;
    end
   
end

eml_assert_no_varsize(6,pos_in,Nsnapshots,ang_in,ncov_in,scov_in,taper_in); 

if size(pos_in,1) == 1
    pos = [zeros(1,size(pos_in,2));pos_in;zeros(1,size(pos_in,2))];
elseif size(pos_in,1) == 2
    pos = [zeros(1,size(pos_in,2));pos_in];
else
    pos = pos_in;
end
sigdatatypes.validate3DCartCoord(pos,'sensorsig','POS');

N_ele = size(pos,2);

sigdatatypes.validateIndex(Nsnapshots,'sensorsig','NS',{'scalar'});

if size(ang_in,1) == 1
    ang = [ang_in;zeros(1,size(ang_in,2))];
else
    ang = ang_in;
end
sigdatatypes.validateAzElAngle(ang,'sensorsig','ANG');

N_ang = size(ang,2);

nonoiseflag = false;
if isscalar(ncov_in)
    validateattributes(ncov_in,{'double'},{'scalar','real',...
        'finite','nonnegative'},'sensorsig','NCOV');
    ncov = ncov_in*eye(N_ele);
    if ~ncov_in
        nonoiseflag = true;
    end
elseif isrow(ncov_in)
    validateattributes(ncov_in,{'double'},{'real',...
        'finite','positive','size',[1 N_ele]},'sensorsig','NCOV');
    ncov = diag(ncov_in);
else
    validateattributes(ncov_in,{'double'},{'finite','size',[N_ele N_ele]},...
        'sensorsig','NCOV');
    tol = 10*eps(max(abs(diag(ncov_in))));  % based on stats cholcov
    cond = any(any(abs(ncov_in - ncov_in') > tol));
    if cond
        coder.internal.errorIf(cond,...
             'phased:phased:notHermitian','NCOV');
    end
    ncov = ncov_in;
end

if isscalar(scov_in)
    validateattributes(scov_in,{'double'},{'scalar','real',...
        'finite','positive'},'sensorsig','SCOV');
    scov = scov_in*eye(N_ang);
elseif isrow(scov_in)
    validateattributes(scov_in,{'double'},{'real',...
        'finite','positive','size',[1 N_ang]},'sensorsig','SCOV');
    scov = diag(scov_in);
else
    validateattributes(scov_in,{'double'},{'finite','size',[N_ang N_ang]},...
        'sensorsig','SCOV');
    tol = 10*eps(max(abs(diag(scov_in))));   % based on stats cholcov
    cond = any(any(abs(scov_in - scov_in') > tol));
    if cond
        coder.internal.errorIf(cond,...
             'phased:phased:notHermitian','SCOV');
    end
    scov = scov_in;
end

if isscalar(taper_in)
    taper = taper_in*ones(N_ele,1);
else
    taper = taper_in;
end
validateattributes(taper,{'double'},{'finite','nonempty','size',[N_ele 1]},...
    'sensorsig','TAPER');

function taper_in = eml_parse_param(varargin)

params = struct('Taper',uint32(0));
pstruct = eml_parse_parameter_inputs(params,[],varargin{:});
taper_in = eml_get_parameter_value(pstruct.Taper,1,varargin{:});

% [EOF]
