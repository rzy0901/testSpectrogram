function stemp = systemp(nf,reftemp)
%systemp Receiver system noise temperature
%   STEMP = systemp(NF) calculates the effective system noise temperature,
%   STEMP (in kelvin), based on the noise figure, NF (in dB). The reference
%   temperature is assumed to be 290 kelvin.
%
%   STEMP = systemp(NF,REFTEMP) specifies the reference temperature (in
%   kelvin) as a nonnegative scalar, REFTEMP, where the default is 290.
%
%   % Example:
%   %   Calculate the system noise temperature of a receiver with a 300 K
%   %   reference temperature and a 5 dB noise figure.
%
%   stemp = systemp(5,300)
%
%   See also phased, phased.ReceiverPreamp, noisepow.

%   Copyright 2010-2012 The MathWorks, Inc.

%   Reference 
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%       McGraw-Hill, 2001 

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(1,2,nargin);
if nargin < 2
    reftemp = 290;
end

sigdatatypes.validateTemperature(reftemp,'systemp','REFTEMP',{'scalar'});
validateattributes(nf,{'double'},...
    {'nonempty','nonnan','scalar','nonnegative'},'systemp','NF');

stemp = reftemp * db2pow(nf);


% [EOF]
