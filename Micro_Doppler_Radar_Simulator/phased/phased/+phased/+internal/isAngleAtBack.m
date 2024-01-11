function idx = isAngleAtBack(ang)
%This function is for internal use only. It may be removed in the future.

%ISANGLEATBACK Test if an angle is in back lobe
%   IND = phased.internal.isAngleAtBack(ANG) returns whether angles
%   specified in ANG that are at the backside. ANG is a 2xM matrix whose
%   columns are the angles in the form of [azimuth; elevation] (in
%   degrees). IND is a length-M logical vector indicating whether the
%   corresponding angle is at the back side. For angles that are beyond +/-
%   90 degrees in azimuth, the corresponding entries in IND are set to
%   true.
%
%   % Example:
%   %   Check if given angles are at the back side.
%   ang = [-91 -90 0 90 91;0 0 0 0 0];
%   ind = phased.internal.isAngleAtBack(ang);

%   Copyright 2016 The MathWorks, Inc.

%#codegen

% set all azimuth beyond +/- 90 to 0
% note elevation +/- 90 should not be set to 0 even if azimuth is beyond
% +/- 90 because all elevation 90 is the same point.
% Need tolerance to deal with numerical round-offs.

idx = ((ang(1,:)<(-90-sqrt(eps))) | (ang(1,:)>(90+sqrt(eps)))) ...
    & (abs(abs(ang(2,:))-90)>sqrt(eps));

end