function ang = espritdoa(R, Nsig, varargin)
%espritdoa   ESPRIT direction of arrival (DOA)
%   ANG = espritdoa(R,NSIG) computes using TLS ESPRIT the vector of
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
%   ANG = espritdoa(R,NSIG,'ElementSpacing',ELSP) returns the estimated
%   direction of arrivals for a ULA with element spacing, ELSP, specified
%   as a positive real scalar in terms of signal wavelength. The default
%   value of ELSP is 0.5.
%   
%   ANG = espritdoa(R,NSIG,'RowWeighting',W) returns the estimated
%   direction of arrivals for a ULA with the row weighing factor, W, for
%   signal subspace eigenvectors as a positive integer scalar. The default
%   is 1. This property controls the weights applied to the selection
%   matrices. In most cases the higher value the better. However it can
%   never be greater than (Ns-1)/2 where Ns is the number of elements of
%   the subarray. The number of subarray element, Ns, is equal to the
%   number of array elements, N, minus one.
%
%   % Example:
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
%   doa = espritdoa(R,Nsig)
%
%   See also phased, rootmusicdoa, spsmooth, aictest, mdltest
%
%   Copyright 2012-2018 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, pp. 1176, Wiley, 2002

%#ok<*EMCA>
%#codegen

phased.internal.narginchk(2,8,nargin);

[elSpacing, eigenvects,diagEigenVals,saSpacing,rweight,M]= ...
    parseInput(R,Nsig,varargin{:});

%Sort eigenvectors
[~,indx] = sort(diagEigenVals,'descend');
eigenvects = eigenvects(:,indx);

D = Nsig;
% check eigenvalues against source dimension. ESPRIT does not work
% when there are less then D non-zero eigenvalues.
D_act = sum(diagEigenVals>eps(max(abs(diagEigenVals))));
if D_act < D
    coder.internal.errorIf(D_act < D ,...
        'phased:phased:ESPRITEstimator:NotEnoughRank', D, D_act);
end

% Row weighting
Ns = M-saSpacing; %number of elements in a subarray
ms = rweight;
w = min(ms,Ns-ms+1);                             % Eq 9.133 in [1]
weights = diag(sqrt([1:w-1 w*ones(1,Ns-2*(w-1)) w-1:-1:1])); % Eq 9.132 in [1]
O = zeros(Ns,saSpacing);

% Selection Matrices
Js1 = [weights O]; % Eq 9.134 in [1]
Js2 = [O weights]; % Eq 9.135 in [1]

% Selecting subarray signal subspaces
Us1 = Js1*eigenvects(:,1:D);
Us2 = Js2*eigenvects(:,1:D);
% TLS-ESPRIT
C = [Us1';Us2']*[Us1 Us2];    % Eq. (9.123) in [1]
[U,~,~] = svd(C);             % C is 2*D x 2*D
V12 = U(1:D,D+1:2*D);         % D x D
V22 = U(D+1:2*D,D+1:2*D);     % D x D
psi = -V12/V22;               % Eq. (9.122) in [1]
psieig = eig(psi);
%   Extract angle information estimated from two subarrays based on the
%   distance of the phase center between subarrays.
doas = 1/saSpacing*angle(psieig);

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
function [elSpacing,eigenvects,diagEigenVals,saSpacing,rweight,numElem] = ...
    parseInput(R,numSignals,varargin)

eml_assert_no_varsize(1:nargin, R, numSignals, varargin{:});
validateattributes(R,{'double'},{'finite','square'},...
        'espritdoa','R');
tol = 10*eps(max(abs(diag(R))));   % based on stats cholcov
cond = any(any(abs(R - R') > tol));
if cond
    coder.internal.errorIf(cond,...
         'phased:phased:notHermitian','R');
end

% Check for positive semi definite
[eigenvects,sED] = eig((R+R')/2);  % ensure Hermitian
sED = diag(sED);
diagEigenVals = sED;
tol = eps(max(abs(sED))); % based on stats cholcov
sED(abs(sED)<=tol)=0;
cond = any(sED<0);
if cond
    coder.internal.errorIf(cond,...
          'phased:phased:notPositiveSemiDefinite','R');
end
                   

numElem = size(R,1);

defaultElSpacing = 0.5;
defaultSaSpacing = 1;
defaultW = 1;
if isempty(coder.target)
     p = inputParser;
     p.addParameter('ElementSpacing',defaultElSpacing);
     p.addParameter('SubarrayOffset',defaultSaSpacing);
     p.addParameter('RowWeighting',defaultW);
     p.parse(varargin{:});
     elSpacing = p.Results.ElementSpacing;
     saSpacing = p.Results.SubarrayOffset;
     rweight =  p.Results.RowWeighting;
 else
     parms = struct('ElementSpacing',uint32(0),...
         'SubarrayOffset',uint32(0),...
         'RowWeighting',uint32(0));
     pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
     elSpacing = ...
         eml_get_parameter_value(pstruct.ElementSpacing,defaultElSpacing,varargin{:});
     saSpacing = ...
         eml_get_parameter_value(pstruct.SubarrayOffset,defaultSaSpacing,varargin{:});
     rweight = ...
         eml_get_parameter_value(pstruct.RowWeighting,defaultW,varargin{:});
end

  validateattributes(elSpacing,{'double'},...
    {'scalar','positive','finite'},'espritdoa','ELSP');
%Nelement of subarray = Number of elements of array - subarray spacing
numElemSubArray = numElem-saSpacing;
%eq 9.133
validateattributes(rweight,{'double'},...
    {'scalar','positive','finite','integer','<=',ceil(numElemSubArray/2)},'espritdoa','W'); 

validateattributes(saSpacing,{'double'},...
    {'scalar','positive','finite','integer','<',numElem-1},'espritdoa','DS');

%p. 1171: number of elements in each subarray >= numSignals +1

validateattributes(numSignals,{'double'},...
    {'scalar','positive','finite','integer','<',numElemSubArray},'espritdoa','NSIG');
% [EOF]

