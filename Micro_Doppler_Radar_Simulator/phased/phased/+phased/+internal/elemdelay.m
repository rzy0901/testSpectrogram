function tau = elemdelay(pos,c,ang)
%This function is for internal use only. It may be removed in the future.

%ELEMDELAY Delay among sensor elements in a sensor array
%   TAU = phased.internal.elemdelay(POS,C,ANG) returns the delay among
%   sensor elements in a sensor array for a given direction specified in
%   ANG. 
%   
%   POS is a 3-row matrix representing the positions of the elements.
%   Each column of POS is in the form of [x;y;z] (in meters). C is a scalar
%   representing the propagation speed (in m/s). ANG is a 2-row matrix
%   representing the incident directions. Each column of ANG is in the form
%   of [azimuth;elevation] form (in degrees). Azimuth angles must be within
%   [-180 180] and elevation angles must be within [-90 90]. TAU is an NxM
%   matrix where N is the number of columns in POS and M is the number of
%   columns in ANG.
%
%   Example:
%   %   Calculate the delay for a 4-element ULA in the direction of 30
%   %   degrees azimuth.
%
%   ha = phased.ULA(4);
%   tau = phased.internal.elemdelay(getElementPosition(ha),3e8,[30;0])

%   Copyright 2011 The MathWorks, Inc.

%#codegen

azang = ang(1,:);
elang = ang(2,:);

% angles defined in local coordinate system
incidentdir = [-cosd(elang).*cosd(azang);...
    -cosd(elang).*sind(azang);...
    -sind(elang)];
tau = pos.'*incidentdir/c;


% [EOF]
