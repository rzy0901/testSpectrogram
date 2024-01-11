function sp = dop2speed(dp,lambda)
%dop2speed Convert Doppler shift to speed
%   V = dop2speed(DOP,LAMBDA) converts the Doppler frequency shift DOP (in
%   Hz) to the corresponding speed V (in m/s). LAMBDA is a scalar
%   representing the signal wavelength (in meters). The function assumes
%   the Doppler shift is associated with one way wave propagation.
%
%   This function supports single and double precision for input data and
%   arguments. If the input data, DOP, is single precision, the output
%   data is single precision. If the input data, DOP, is double precision,
%   the output data is double precision. The precision of the output is
%   independent of the precision of the arguments.
%
%   % Example:
%   %   Calculate the speed corresponding to a 300Hz Doppler shift. The
%   %   signal wavelength is 3 m.
%
%   dp = 300; wavelen = 3;
%   sp = dop2speed(dp,wavelen)
%
%   See also phased, speed2dop, dopsteeringvec.

%   Copyright 2010 The MathWorks, Inc.

%   Reference 
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%       McGraw-Hill, 2001 
%   [2] Theodore Rappaport, Wireless Communications Principles & Practices, 
%       Prentice Hall, 1996

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(2,2,nargin);
eml_assert_no_varsize(2,dp,lambda);

validateattributes(dp,{'double','single'},{'real','finite'},'dop2speed','DOP');
validateattributes(lambda,{'double','single'},{'scalar','positive','finite'},'dop2speed','LAMBDA');

lambda = cast(lambda,class(dp));
% doppler frequency is defined as dp = sp/lambda
sp = dp*lambda;

% [EOF]
