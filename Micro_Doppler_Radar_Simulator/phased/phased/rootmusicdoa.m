function ang = rootmusicdoa(R, Nsig, varargin)
%rootmusicdoa   Root MUSIC direction of arrival (DOA)
%   ANG = rootmusicdoa(R,NSIG) computes using Root MUSIC the vector of
%   estimated direction of arrivals, ANG, in degrees of NSIG signals whose
%   covariance matrix estimate is given by the positive definite matrix R.
%   R is an NxN matrix representing the covariance matrix of the input to a
%   uniform linear array (ULA) with half wavelength spacing between
%   elements.  Exact conjugate-symmetry of R is ensured by forming (R+R')/2
%   inside the function.
%
%   NSIG is a positive integer scalar specifying the number of
%   signals to estimate.  NSIG should be less than N.
%
%   ANG = rootmusicdoa(R,NSIG,'ElementSpacing',DIST) returns the estimated
%   direction of arrivals for a ULA with element spacing, DIST, specified
%   as a positive real scalar in terms of signal wavelength.
%   
%   % Examples:
%
%   % Example 1:
%   %   Calculate the direction of arrival of 4 uncorrelated signals  
%   %   impinging on a 10-element half-wavelength spacing ULA.  Assume the 
%   %   signals coming from azimuth 0,-25,45, and 60 degrees,
%   %   respectively. The noise is white across all elements and the SNR 
%   %   is 5 dB.
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   angles = [0 -25 45 60];
%   Nsig = 4;
%   % Received covariance matrix
%   R = sensorcov(elementPos,angles,db2pow(-5));
%   doa = rootmusicdoa(R,Nsig)
%
%   % Example 2:
%   %   Calculate the direction of arrival of a square wave impinging on a
%   %   10-element half-wavelength spacing ULA at 30 degrees azimuth.
%   %   The noise is white across all elements and the SNR is -5 dB.
%
%   K = 1000; %number of snapshots
%   t=0:K-1;
%   x = square(2*pi*.01*t);
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   sv = steervec(elementPos,30);
%   xc = sv*x;  % collected signal
%   % Add noise
%   sigmaN = db2pow(-5);
%   noise = (randn(N,K)+ 1i*randn(N,K))*(sigmaN/2)^.5;
%   xcn = xc+noise;
%   R = (xcn*xcn')/K;
%   rootmusicdoa(R,1)
%
%   See also phased, espritdoa, spsmooth, aictest, mdltest

%   Copyright 2012-2018 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, pp. 1159, Wiley, 2002

%#ok<*EMCA>
%#codegen

phased.internal.narginchk(2,4,nargin);

[elSpacing, eigenvects,diagEigenVals]= ...
    parseInput(R,Nsig,varargin{:});

N = size(R,1);
%Sort eigenvectors
[~,indx1] = sort(diagEigenVals,'descend');
 eigenvects = eigenvects(:,indx1);
% Separate the signal and noise eigenvectors
noise_eigenvects = eigenvects(:,Nsig+1:end);

% Form a polynomial D
% D consists of a sum of polynomials given by the product of the noise
% subspace eigenvectors and the reversed and conjugated version.
D = complex(zeros(2*N-1,1));
for i = 1:N-Nsig
    D = D + conv(noise_eigenvects(:,i),conj(flipud(noise_eigenvects(:,i))));
end

% Take the angle of the NSIG roots of D closest to the unit circle.
roots_D = roots(D);
roots_D1 = roots_D(abs(roots_D) < 1);
[~,indx] = sort(abs(abs(roots_D1)-1));
sorted_roots = roots_D1(indx);
doas = angle(sorted_roots(1:Nsig));

%Convert estimated angle in sin-space to degrees. This method is valid for
%ULA only.

u = doas/(2*pi*elSpacing);
% check whether all elements of u are within [-1,1]
idx = find(abs(u)<=1);
if  length(idx) <Nsig && isempty(coder.target)
        warning(message('phased:phased:internal:AbstractULASubspaceDOA:InvalidPsi',Nsig));
end
if isempty(idx)
    ang = zeros(1,0);
else
    ang = asind(u(idx));
    % convert to row vector
    ang = ang(:).';
end

%--------------------------------------
function [elSpacing,eigenvects,diagEigenVals] = ...
    parseInput(R,numSignals,varargin)


eml_assert_no_varsize(1:nargin, R, numSignals, varargin{:});
validateattributes(R,{'double'},{'finite','square'},...
        'rootmusicdoa','R');
tol = 10*eps(max(abs(diag(R))));   % based on stats cholcov
cond = any(any(abs(R - R') > tol));
if cond
    coder.internal.errorIf(cond,...
         'phased:phased:notHermitian','R');
end

% Check for positive semi definite
[eigenvects,sEDArg] = eig((R+R')/2);  % ensure Hermitian
sED = diag(real(sEDArg));
diagEigenVals =sED;
tol = eps(max(abs(sED))); % based on stats cholcov
sED(abs(sED)<=tol)=0;
cond = any(sED<0);
if cond
    coder.internal.errorIf(cond,...
         'phased:phased:notPositiveSemiDefinite','R');
end
                   

M = size(R,1);

validateattributes(numSignals,{'double'},...
    {'scalar','positive','finite','integer','<',M},'rootmusicdoa','NSIG');

defaultElSpacing = 0.5;

if isempty(coder.target)
     p = inputParser;
     p.addParameter('ElementSpacing',defaultElSpacing);
     p.parse(varargin{:});
     elSpacing = p.Results.ElementSpacing;
 else
     parms = struct('ElementSpacing',uint32(0));
     pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
     elSpacing = ...
         eml_get_parameter_value(pstruct.ElementSpacing,defaultElSpacing,varargin{:});
end

validateattributes(elSpacing,{'double'},...
    {'scalar','positive','finite'},'rootmusicdoa','DIST');


% [EOF]

