function sv = steeringvec(pos,freq,c,ang,N)
%This function is for internal use only. It may be removed in the future.

%STEERINGVEC Steering vector for a sensor array
%   SV = phased.internal.steeringvec(POS,FREQ,C,ANG) returns the steering
%   vector of a sensor array at given frequencies and angles. 
%
%   POS is a 3-row matrix representing the positions of the elements.
%   Each column of POS is in the form of [x;y;z] (in meters). C is a scalar
%   representing the propagation speed (in m/s). ANG is a 2-row matrix
%   representing the incident directions. Each column of ANG is in the form
%   of [azimuth;elevation] form (in degrees). Azimuth angles must be within
%   [-180 180] and elevation angles must be within [-90 90]. FREQ is a
%   length-L row vector with each entry represents a frequency (in Hz).
%
%   TAU has a dimension of NxMxL where N is the number of columns in POS
%   and M is the number of columns in ANG.
%   
%   SV = phased.internal.steeringvec(POS,FREQ,C,ANG,N) returns the steering
%   vector with N-bit phase shifters. 
%
%   Example:
%   %   Calculate the steering vector of a 4-element ULA in the direction 
%   %   of 30 degrees azimuth. Assume the frequency is 300 MHz.
%
%   ha = phased.ULA(4);
%   sv = phased.internal.steeringvec(getElementPosition(ha),3e8,3e8,[30;0])

%   Copyright 2015 The MathWorks, Inc.

%#codegen

if nargin<5
    N = 0;
end
tau = phased.internal.elemdelay(pos,c,ang);
if N == 0
    if isscalar(freq)
        sv = exp(-1i*2*pi*freq*tau);
    else
        freqtau = reshape(tau(:)*freq,...
            size(pos,2),size(ang,2),[]);
        sv = exp(-1i*2*pi*freqtau);
    end
else
    if isscalar(freq)
        freqtau = phased.internal.quantizePhase(-freq*tau,N,1);
        sv = exp(1i*2*pi*freqtau);
    else
        freqtau = reshape(tau(:)*freq,...
            size(pos,2),size(ang,2),[]);
        freqtau = phased.internal.quantizePhase(-freqtau,N,1);
        sv = exp(1i*2*pi*freqtau);
    end
end

% [EOF]
