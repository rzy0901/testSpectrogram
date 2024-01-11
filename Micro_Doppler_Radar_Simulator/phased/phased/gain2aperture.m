function Ae = gain2aperture(g, lambda)
%gain2aperture  Convert gain to effective aperture
%   A = gain2aperture(G,LAMBDA) converts the gain G (in dB) of an antenna
%   to the corresponding effective aperture A (in square meters) when the
%   antenna is used to capture an electromagnetic wave with wavelength
%   LAMBDA (in meters). G can be a vector but LAMBDA must be a scalar. A
%   has the same dimensionality as G with each entry in A representing the
%   effective aperture for the corresponding gain in G.
%
%   % Example:
%   %   An antenna has a gain of 3 dB. Calculate its effective aperture 
%   %   when used to capture an electromagnetic wave with a wavelength of 
%   %   10 cm.
%
%   a = gain2aperture(3,0.1)
%
%   See also phased, aperture2gain.

%   Copyright 2010 The MathWorks, Inc.

%   Reference
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 2001


%#codegen 
%#ok<*EMCA

phased.internal.narginchk(2,2,nargin);
eml_assert_no_varsize(2, g, lambda);
validateattributes(g,{'numeric'},{'real','nonempty','vector'},...
    'gain2aperture','G');
sigdatatypes.validateDistance(lambda,'gain2aperture','LAMBDA',{'scalar','positive'});

Ae = db2pow(g)*lambda^2/(4*pi);

% [EOF]
