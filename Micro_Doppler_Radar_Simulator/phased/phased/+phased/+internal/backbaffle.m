function g = backbaffle(g,ang)
%This function is for internal use only. It may be removed in the future.

%BACKBAFFLE Back baffle the response
%   G = phased.internal.backbaffle(R,ANG) back baffle the response R. R is
%   a length-M vector the original response before the back baffling. ANG
%   is a 2xM matrix whose columns are the angles where the corresponding
%   responses in R are measured. Each column in ANG is in the form of
%   [azimuth; elevation] (in degrees). G is a length-M vector containing
%   the response after the back baffling. After back baffling, responses
%   corresponding to the angles that are beyond +/- 90 degrees in azimuth
%   are set to 0.
%
%   Note that the responses should be represented in linear scale.
%
%   % Example:
%   %   Back baffle a uniform response.
%   resp = [1 1 1 1 1];
%   ang = [-91 -90 0 90 91;0 0 0 0 0];
%   resp = phased.internal.backbaffle(resp,ang);

%   Copyright 2010-2016 The MathWorks, Inc.

% set all azimuth beyond +/- 90 to 0
% note elevation +/- 90 should not be set to 0 even if azimuth is beyond
% +/- 90 because all elevation 90 is the same point.
% Need tolerance to deal with numerical round-offs.

%#codegen

idx = phased.internal.isAngleAtBack(ang);
g(idx) = 0;


% [EOF]
