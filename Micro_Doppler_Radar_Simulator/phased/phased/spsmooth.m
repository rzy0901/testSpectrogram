function RSM = spsmooth(R, L, fb)
%spsmooth   Spatial smoothing of a covariance matrix
%   RSM = spsmooth(R,L) computes the average covariance matrix,
%   RSM, of the NxN covariance matrix, R,over L maximum overlapped
%   subarrays. L is a positive integer less than N. RSM is an
%   (N-L+1)x(N-L+1) matrix.
%   
%   RSM = spsmooth(R,L,'fb') performs forward-backward averaging on
%   R.
%   
%   % Example:
%   %   Calculate the direction of arrival of 2 correlated signals  
%   %   impinging on a 10-element half-wavelength spacing ULA.  Assume the 
%   %   signals coming from azimuth 0 and -25 degrees, respectively and 
%   %   that the signal coming from -25 degrees is attenuated by 1/5.
%   %   The noise is white across all elements and the SNR is 5 dB.
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   angles = [0 -25];
%   ac = [1 1/5]; %Attenuation coefficient
%   %Signal cov matrix
%   scov = ac'*ac;
%   Nsig = 2;
%   % Received covariance matrix
%   R = sensorcov(elementPos,angles,db2pow(-5),scov);
%   % Perform spatial smoothing to be able to detect
%   % L-1 coherent signals. We choose L = 3.
%   L = 3;
%   RSM = spsmooth(R,L,'fb');
%   doasm = rootmusicdoa(RSM,Nsig);
%   %DOA without smoothing
%   doa =  rootmusicdoa(R,Nsig);
%
%   See also phased, espritdoa, rootmusicdoa, aictest, mdltest

%   Copyright 2012 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, pp. 719, Wiley, 2002
%
%   [2] D.A. Linebarger, R.D DeGroat and E.M. Dowling. Efficient
%       direction-finding methods employing forward/backward averaging.
%       IEEE Trans. Signal Process. vol.SP-42, pp.2136-2145, August 1994

%#ok<*EMCA>
%#codegen

phased.internal.narginchk(2,3,nargin);

eml_assert_no_varsize(1:2, R, L);
validateattributes(R,{'double'},{'finite','square'},...
        'spsmooth','R');     
M = size(R,1);
validateattributes(L,{'double'},...
    {'scalar','positive','finite','integer','<',M},'spsmooth','L');
fbFlag = false;
if  nargin == 3
    validatestring(fb,{'fb'},'spsmooth','',3);
    fbFlag = true;
end
%Spatial smoothing
%L is number of subarray, M is size of subarray
M = size(R,1)-L+1;
if isreal(R)
    RSM = zeros(M);
else
    RSM = complex(zeros(M));
end
for k = 1:L
    RSM = RSM+ R(k:k+M-1,k:k+M-1);
end
RSM = RSM/L;
if fbFlag
    %Spatial smoothing forward backward
    if isreal(R)
        RSMB = zeros(M);
    else
        RSMB = complex(zeros(M));
    end
    for k = 1:L
        RSMB = RSMB+rot90(R(k:k+M-1,k:k+M-1).',2);
    end
    RSMB = RSMB/L;
    RSM = (RSM+RSMB)/2;
end

% [EOF]

