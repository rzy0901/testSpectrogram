function g = aperture2gain(Ae,lambda)
%aperture2gain  Convert effective aperture to gain
%   G = aperture2gain(A,LAMBDA) converts the effective aperture A (in
%   square meters) of an antenna to the corresponding gain G (in dB) when
%   the antenna is used to capture an electromagnetic wave with wavelength
%   LAMBDA (in meters). A can be a vector but LAMBDA must be a scalar. G
%   has the same dimensionality as A with each entry in G representing the
%   gain for the corresponding effective aperture in A.
%
%   % Example:
%   %   An antenna has an effective aperture of 3 square meters. Find its
%   %   gain when used to capture an electromagnetic wave with a wavelength
%   %   of 10 cm.
%
%   g = aperture2gain(3,0.1)
%
%   See also phased, gain2aperture.

%   Copyright 2010 The MathWorks, Inc.

%   Reference
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 2001

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,2,nargin);
eml_assert_no_varsize(2,Ae,lambda);
sigdatatypes.validateDistance(lambda,'aperture2gain','LAMBDA',...
    {'scalar','positive'});
sigdatatypes.validateArea(Ae,'aperture2gain','A',...
    {'vector'});
g = 4*pi*Ae./(lambda^2);
g = pow2db(g);

end
