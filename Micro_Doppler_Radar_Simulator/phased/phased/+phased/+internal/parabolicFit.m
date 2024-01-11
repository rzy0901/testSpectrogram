function deltak = parabolicFit(x,k0)
%This function is for internal use only. It may be removed in the future.

%PARABOLICFIT   Perform parabolic fit for peak detection 
%   DELTAK = phased.internal.parabolicFit(X,K0) performs the parabolic fit
%   to do the peak detection in X around index K0. X can be a matrix whose
%   columns are independent channels and K0 specifies the initial peak
%   locations in X. K0 can be a matrix so each channel may have multiple
%   peaks. The number of columns in K0 must match the number of columns in
%   X. 
%
%   The function performs a 3-point parabolic fit around each entry in K0
%   and then use the fitted parabolic function to find the true peak
%   location. The differences between the true peak locations and initial
%   peak locations, K0, are returned in DELTAK. This process can often be
%   used to get a better peak location and reduce straddle loss.
%
%   % Examples:
%
%   % Example 1:
%   %   Perform a parabolic fit to find the peak location in a spectrum.
%
%   fs = 100; N = 128; t = (0:N-1).'/fs; NFFT = 128; 
%   f = 30; x = exp(1i*2*pi*f*t)+(randn(N,1)+1i*randn(N,1));  
%   xfabs = abs(fft(x.*hann(N),NFFT)); freq = (0:NFFT-1)*fs/NFFT;
%   [~,k0] = max(xfabs);
%   dk = phased.internal.parabolicFit(xfabs,k0);
%   f_est = freq(k0)+dk*fs/N
%
%   % Example 2:
%   %   Fit a parabolic function
%   
%   t = (0:5).';   x = -(t-2.8).^2;
%   [~,k0] = max(x);
%   dk = phased.internal.parabolicFit(x,k0);
%   pkloc = t(k0)+dk*mean(diff(t))
%
%   See also phased.

%   Copyright 2015 The MathWorks, Inc.

% References: 
% [1] Mark Richards, Fundamentals of Radar Signal Processing, pp266 
%     eq 5.103

%#ok<*EMCA>
%#codegen

% x: columns are channels
% k0: can be a matrix, column # should match x

idx = bsxfun(@plus,k0,(0:size(x,2)-1)*size(x,1));
deltak = 0.5*(x(idx-1)-x(idx+1))./...
    (x(idx+1)+x(idx-1)-2*x(idx));  
deltak(abs(deltak)>=1) = 0;  % the delta should never be more than 1 sample
