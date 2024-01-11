function az_el = incident2azel(inc_ang,refaxes)
%This function is for internal use only. It may be removed in the future.

%INCIDENT2AZEL Covert incident angle to az/el format
%   AZEL = phased.internal.incident2azel(INC_ANG,REFAXES) converts the
%   incident angles specified in INC_ANG to the az/el angles in the
%   coordinate system defined in REFAXES.
%
%   INC_ANG is a 2-row matrix representing the input angles specified in
%   the original coordinates. Each column is an angle in the form of
%   [azimuth; elevation] (in degrees), where azimuth angles are within
%   [-180 180] and elevation angles are within [-90 90]. REFAXES is a 3x3
%   matrix with each column representing one of three orthonormal axes for
%   the new 3-D coordinate system. These axes are represented in the
%   original coordinate system. AZEL has the same dimensions as INC_ANG
%   with each column representing the corresponding angle in the new
%   coordinate system in the form of [azimuth; elevation].
%
%   Example:
%   %   Convert the direction of [20;30] to a new coordinate system where
%   %   the axes are rotated through 20 degrees azimuth and 30 degrees
%   %   elevation.
%
%   ang_azel = phased.internal.incident2azel([20;30],...
%           phased.internal.rotazel(eye(3),[20;30]))

%   Copyright 2011 The MathWorks, Inc.

%#codegen

%This function is for internal use only. It may be removed in the future.
LclCoord = phased.internal.global2localcoord(...
    [inc_ang; ones(1,size(inc_ang,2))],'ss',[0; 0; 0],...
    refaxes);
% obtain the angle
az_el = LclCoord(1:2,:);


% [EOF]
