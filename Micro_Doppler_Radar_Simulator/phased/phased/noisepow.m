function npow = noisepow(nbw,nf,reftemp)
%noisepow Noise power at the receiver
%   NPOWER = noisepow(NBW,NF,REFTEMP) calculates the noise power NPOWER (in
%   Watts) based on the noise bandwidth, NBW (in Hz), the noise figure, NF
%   (in dB), and the reference temperature, REFTEMP (in kelvin).
%
%   % Example:
%   %   Calculate the noise power of a receiver whose noise bandwidth is 10
%   %   kHz, noise figure is 1 dB, and reference temperature is 300 K.
%
%   npower = noisepow(10e3,1,300)
%
%   See also phased, phased.ReceiverPreamp, systemp.

%   Copyright 2010 The MathWorks, Inc.

%   Reference 
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%       McGraw-Hill, 2001 

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(1,3,nargin);

if nargin < 3
    reftemp = 290;
end
if nargin < 2
    nf = 0;
end
eml_assert_no_varsize(1:nargin,nbw,nf,reftemp);
sigdatatypes.validateFrequency(nbw,'noisepow','NBW',{'scalar'});
sigdatatypes.validateTemperature(reftemp,'noisepow','REFTEMP',{'scalar'});
validateattributes(nf,{'double'},...
    {'nonempty','nonnan','scalar','nonnegative'},'noisepow','NF');

B = physconst('Boltzmann');
npow =  B * systemp(nf,reftemp) * nbw;


% [EOF]
