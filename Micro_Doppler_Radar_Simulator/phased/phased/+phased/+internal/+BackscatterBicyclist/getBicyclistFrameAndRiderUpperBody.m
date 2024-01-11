function [frameAndRiderUpperBody,mounts] = getBicyclistFrameAndRiderUpperBody(spacing)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

% Setup bicycle frame
[frameScatterers,rearWheelMountPos,pedalMountPos] = phased.internal.BackscatterBicyclist.getBicyclistFrame(spacing);

% Define cyclist body 
[riderUpperBodyScatterers,hipMountPos] = phased.internal.BackscatterBicyclist.getBicyclistRiderUpperBody(frameScatterers,spacing);

% Assign stationary bicycle scatterers
frameAndRiderUpperBody = [frameScatterers...
 riderUpperBodyScatterers];

% Assign mount outputs
mounts.FrontWheel = [0; 0; 0]; 
mounts.RearWheel = rearWheelMountPos;
mounts.Pedals = pedalMountPos;
mounts.Hip = hipMountPos; 

end