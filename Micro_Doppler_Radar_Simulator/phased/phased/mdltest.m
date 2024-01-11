function NSIG = mdltest(x, fb)
%mdltest   Minimum description length test
%   NSIG = mdltest(X) estimates the number of signals, NSIG, from the
%   received signal, X, using the Minimum Description Length test. X 
%   is a complex KxN matrix, where K is the number of snapshots and N
%   is the number of sensors at the receiver.
%   
%   NSIG = mdltest(X,'fb') estimates the number of signals by performing
%   forward-backward averaging on the input.
%
%   % Example:
%   %   Estimate the number of signals from a 300 snapshots collected    
%   %   input of 2 signals impinging on a 10-element half-wavelength 
%   %   spacing ULA. Assume the signals coming from azimuth 0 and -25 
%   %   degrees, respectively.The noise is white across all elements  
%   %   and the SNR is 5 dB.
%
%   N = 10;     % Elements in array
%   d = 0.5;    % sensor spacing half wavelength
%   elementPos = (0:N-1)*d;
%   angles = [0 -25];
%   % Received covariance matrix
%   k = 300;
%   x = sensorsig(elementPos,k,angles,db2pow(-5));
%   Nsig = mdltest(x)
%
%   See also phased, espritdoa, rootmusicdoa, spsmooth, aictest

%   Copyright 2012 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCA>
%#codegen

phased.internal.narginchk(1,2,nargin);
eml_assert_no_varsize(1,x);
validateattributes(x,{'double'},{'finite','nonempty'},...
        'mdltest','x');     

fbFlag = false;
if  nargin == 2
    validatestring(fb,{'fb'},'mdltest','',2);
    fbFlag = true;
end

K = size(x,1);
R = 1/K*(x'*x);
if fbFlag
  R = spsmooth(R,1,'fb');
end
[~, eigenvalsOut] = eig(R);
eigenvals = real(eigenvalsOut);
eigenvals = sort(diag(eigenvals),'descend');
eigenvals(eigenvals<0) = 0;
        
NSIG = phased.internal.mdltest(eigenvals,K,fbFlag);


% [EOF]

