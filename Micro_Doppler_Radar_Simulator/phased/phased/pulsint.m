function sigout = pulsint(sigin, method)
%pulsint  Pulse integration
%   Y = pulsint(X) performs the noncoherent (video) integration on the
%   received pulses X. Each column of X is considered as one pulse and
%   pulsint noncoherently integrates all available pulses. Y is the
%   resulting integrated pulse.
%   
%   Y = pulsint(X,METHOD) performs the pulse integration using the
%   specified method METHOD as one of 'coherent' | 'noncoherent'. The
%   default value of METHOD is 'noncoherent'.
%
%   If METHOD is 'coherent', the output Y is given by
%   
%   Y = X_1 + X_2 + ... + X_n.
%
%   If METHOD is 'noncoherent', the output Y is given by
%
%   Y = sqrt(|X_1|^2 + ... + |X_n|^2).
%
%   % Example:
%   %   Noncoherently integrate 10 pulses.
%
%   x = repmat(sin(2*pi*(0:99)'/100),1,10)+0.1*randn(100,10);   
%   y = pulsint(x);
%      
%   subplot(211), plot(abs(x(:,1)));
%   ylabel('Magnitude');
%   title('First Pulse');
%   subplot(212), plot(abs(y));
%   ylabel('Magnitude');
%   title('Integrated Pulse');
%
%   See also phased, phased.MatchedFilter.

%   Copyright 2007-2011 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(1,2,nargin);

if nargin < 2
    method = 'noncoherent';
end

method = validatestring(method,{'coherent','noncoherent'},...
    'pulsint','METHOD');

% Convert input signal to columns
if isvector(sigin)
    sigin = sigin(:);
end

if method(1) == 'n' 
% Noncoherent pulse integration
    sigout = sum(abs(sigin).^2,2);
    sigout = sqrt(sigout);
else
    sigout = sum(sigin,2);
end

% [EOF]
