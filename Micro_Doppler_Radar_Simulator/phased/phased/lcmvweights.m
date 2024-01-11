function w = lcmvweights(constraints,response,cov)
%lcmvweights   LCMV beamformer weights for a sensor array
%   W = lcmvweights(CONSTR,RESP,COV) returns the narrowband linear
%   constraint minimum variance (LCMV) beamformer weights, W, of a the
%   sensor array, performed on the signal covariance matrix, COV, with
%   desired responses, RESP, corresponding to the constraint matrix,
%   CONSTR.
%
%   CONSTR represents the constraint matrix used for the LCMV beamformer.
%   CONSTR is as an NxK matrix. Each column of the matrix is a constraint
%   and N is the number of elements in the sensor array. K must be less
%   than or equal to N.
%
%   RESP represents the desired response used for LCMV beamformer. RESP is
%   a length-K column vector where K is the number of constraints. Each
%   element in the vector defines the desired response of the constraint
%   specified in the corresponding column of CONSTR.
%
%   COV is an NxN matrix representing the covariance matrix across the N
%   sensor elements for all incoming signals.
%
%   % Example:
%   %   Calculate the LCMV beamforming weights of a 10-element 
%   %   half-wavelength spacing ULA with a response of 1 in the direction 
%   %   of -20 degrees azimuth and 0 in the direction of 0 and 20 degrees
%   %   azimuth. Assume there are two signals coming from azimuth 60 and 
%   %   -60 degrees, respectively. The noise is white across all elements
%   %   and the SNR is 10 dB.
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   sv = steervec(elementPos,[-20 0 20]);
%   resp = [1 0 0]';
%   Sn  = sensorcov(elementPos,[-60 60],db2pow(-10));
%
%   w = lcmvweights(sv,resp,Sn);
%   plotangl = -90:90;
%   vv = steervec(elementPos,plotangl);
%   plot(plotangl,mag2db(abs(w'*vv)))
%   grid on
%   xlabel('Azimuth Angle (degrees)');
%   ylabel('Normalized Power (dB)');
%
%   See also phased, steervec, cbfweights, mvdrweights,
%   phased.LCMVBeamformer

%   Copyright 2012 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, pp. 527, Wiley, 2002

%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,3,nargin);

parseInput(constraints,response,cov);

% the covariance matrix is defined as E{x.'*conj(x)}
AS = constraints;
covinvAS = cov\AS;
ASHcovinvAS = AS'*covinvAS;

% Ensure Hermitian by averaging ASHcovinvAS
w  = (covinvAS/((ASHcovinvAS+ASHcovinvAS')/2))*response;     %LCMV weights


%---------------------------------
function  parseInput(constraints_in,response_in,scov_in)
eml_assert_no_varsize(1:nargin, constraints_in,response_in,scov_in);
validateattributes( constraints_in, { 'double' }, { '2d', 'finite', 'nonempty' }, ...
                    'lcmvweights', 'CONSTR');
cond =  any(all(constraints_in==0));
if cond
    coder.internal.errorIf(cond, ...
         'phased:beamformer:SMI:expectedNonZero');
end

numConstr = size(constraints_in,2);
numElem = size(constraints_in,1);
cond = numConstr > numElem;
if cond
    coder.internal.errorIf(cond, ...
         'phased:LCMVBeamformer:TooManyConstraints');
end

validateattributes( response_in, { 'double' }, { 'column', 'finite', 'nonempty', ...
                   'size',[numConstr 1]}, 'lcmvweights', 'RESP');

validateattributes(scov_in,{'double'},{'finite','size',[numElem numElem]},...
        'lcmvweights','COV');
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


% [EOF]
