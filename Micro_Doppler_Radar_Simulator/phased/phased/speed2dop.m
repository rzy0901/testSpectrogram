function dp = speed2dop(sp,lambda)
%speed2dop Convert speed to Doppler shift
%   DOP = speed2dop(V,LAMBDA) converts the speed V (in m/s) to the
%   corresponding Doppler frequency shift DOP (in Hz). LAMBDA is a scalar
%   representing the signal wavelength (in meters). The function assumes
%   the Doppler shift is associated with one way wave propagation.
%
%   % Example:
%   %   Calculate the Doppler shift corresponding to the speed of 300 m/s.
%   %   The signal wavelength is 3 m.
%
%   sp = 300; wavelen = 3;
%   dp = speed2dop(sp,wavelen)
%
%   See also phased, dop2speed, dopsteeringvec.

%   Copyright 2010 The MathWorks, Inc.

%   Reference 
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%       McGraw-Hill, 2001 
%   [2] Theodore Rappaport, Wireless Communications Principles & Practices, 
%       Prentice Hall, 1996

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,2,nargin);

validateattributes(sp,{'double'},{'real','finite'},'speed2dop','SPEED');
validateattributes(lambda,{'double'},{'scalar','positive','finite'},'speed2dop','LAMBDA');

dp = sp/lambda;

% [EOF]
