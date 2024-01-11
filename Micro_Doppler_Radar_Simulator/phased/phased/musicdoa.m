function [ang,spec,specAng] = musicdoa(R, Nsig, varargin)
%musicdoa   MUSIC direction of arrival (DOA)
%   ANG = musicdoa(R,NSIG) computes, using MUSIC, the vector of estimated
%   directions of arrival, ANG, of NSIG signals whose covariance matrix
%   estimate is given by the positive definite matrix R. Angle units are in
%   degrees. R is an NxN matrix representing the covariance of the input to
%   a uniform linear array (ULA) with half-wavelength spacing between
%   elements. Exact conjugate symmetry of R is ensured by forming (R+R')/2
%   inside the function.
%
%   NSIG is a positive integer scalar specifying the number of
%   signals to estimate. NSIG should be less than N.
%
%   [ANG,SPEC,SPECANG] = musicdoa(R,NSIG) returns the spatial spectrum,
%   SPEC, and corresponding broadside scan angles, SPECANG, of the NSIG
%   signals whose covariance matrix estimate is given by R. Detected
%   sources appear as peaks in the spatial spectrum.
%
%   [...] = musicdoa(R,NSIG,...,'ScanAngles',SCANANG) searches a vector of
%   broadside scan angles, SCANANG, and returns the NSIG angles
%   corresponding to the highest peaks in the MUSIC spatial spectrum. The
%   angles in SCANANG must lie in the range [-90,90]. By default, SCANANG
%   is the vector -90:90.
%
%   [...] = musicdoa(R,NSIG,...,'ElementSpacing',DIST) returns the
%   estimated spatial spectrum and directions of arrival for a ULA with
%   element spacing, DIST, specified as a positive real scalar in terms of
%   signal wavelength.
%
%   % Examples:
%
%   % Example 1:
%   %   Calculate the directions of arrival of 4 uncorrelated signals  
%   %   impinging on a 10-element ULA with half-wavelength spacing. Assume 
%   %   the signals are coming from azimuth 0, -25, 45, and 60 degrees,
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
%   doa = musicdoa(R,Nsig)
%
%   % Example 2:
%   %   Repeat example 1 and plot the spatial spectrum of the signals. 
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   angles = [0 -25 45 60];
%   Nsig = 4;
%   % Received covariance matrix
%   R = sensorcov(elementPos,angles,db2pow(-5));
%   [doa,spec,specang] = musicdoa(R,Nsig);
%   figure;
%   plot(specang,spec)
%   xlabel('Broadside Angles (degrees)')
%   ylabel('Magnitude Spectrum')
%   title('MUSIC Spatial Spectrum')
%
%   % Example 3:
%   %   Calculate the direction of arrival of a square wave impinging on a
%   %   10-element ULA with half-wavelength spacing at 30 degrees azimuth.
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
%   musicdoa(R,1)
%
%   See also phased, rootmusicdoa, espritdoa, spsmooth

%   Copyright 2016 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, pp. 1159, Wiley, 2002

%#ok<*EMCA>
%#codegen

phased.internal.narginchk(2,6,nargin);

[elSpacing, scanAng, eigenvects, diagEigenVals]= ...
    parseInput(R,Nsig,varargin{:});

N = size(R,1);

% Sort eigenvectors
[~,indx1] = sort(diagEigenVals,'descend');
 eigenvects = eigenvects(:,indx1);
 
% Separate the signal and noise eigenvectors
noise_eigenvects = eigenvects(:,Nsig+1:end);

% Compute steering vectors. Broadside angles are equal to azimuth angles
% for zero elevation.
elementPos = (0:N-1)*elSpacing;
sv = steervec(elementPos,scanAng);
 
% Calculate the spatial spectrum. Add a small positive constant to prevent
% division by zero.
D = sum(abs((sv'*noise_eigenvects)).^2,2)+eps(1); % 9.44 in [1]
spec = sqrt(1./D).';
specAng = scanAng;

% Find DOA
[~,locs] = findpeaks(spec,'SortStr','descend');
D = min(Nsig,length(locs));
assert(D <= Nsig);
if D>0
    ang = scanAng(locs(1:D));
else
    ang = zeros(1,0);
end

%--------------------------------------
function [elSpacing,scanAng,eigenvects,diagEigenVals] = ...
    parseInput(R,numSignals,varargin)

eml_assert_no_varsize(1:nargin, R, numSignals, varargin{:});
validateattributes(R,{'double'},{'finite','square'},...
        'musicdoa','R');
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
    {'scalar','positive','finite','integer','<',M},'musicdoa','NSIG');

defaultElSpacing = 0.5;
defaultScanAng = -90:90;

if isempty(coder.target)
     p = inputParser;
     p.addParameter('ElementSpacing',defaultElSpacing);
     p.addParameter('ScanAngles',defaultScanAng);
     p.parse(varargin{:});
     elSpacing = p.Results.ElementSpacing;
     scanAng = p.Results.ScanAngles;
else
    parms = struct( ...
        'ElementSpacing',uint32(0), ...
        'ScanAngles',uint32(0));
    poptions = struct( ...
        'PartialMatching','unique', ...
        'StructExpand',false);
    pstruct = coder.internal.parseParameterInputs(parms,poptions,varargin{:});
    elSpacing = coder.internal.getParameterValue(...
      pstruct.ElementSpacing,defaultElSpacing,varargin{:});
    scanAng = coder.internal.getParameterValue(...
      pstruct.ScanAngles,defaultScanAng,varargin{:});
end

validateattributes(elSpacing,{'double'},...
  {'scalar','positive','finite'},'musicdoa','DIST');

validateattributes(scanAng,{'double'},...
  {'vector','finite','<=',90,'>=',-90},'musicdoa','SCANANG');

scanAng = scanAng(:)';

% [EOF]