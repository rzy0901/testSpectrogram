function A = global2localcoordjac(tgtpos, sensorpos, laxes)
%This function is for internal use only. It may be removed in the future.

% global2localcoord Calculates the Jacobian of the global2localcoord 
%   A = global2localcoordjac(tgtpos, sensorpos, laxes) when
%   tgtpos and sensorpos are in global rectangular coordinates and the
%   Jacobian is of the correponding spherical coordinates. laxes is the
%   matrix that defines the orinetation of the sensor.

%   Copyright 2016 The MathWorks, Inc.

% internal function, no error checking is performed

%#codegen 

relpos = tgtpos - sensorpos;
relposlocal = laxes' * relpos;
xrel = relposlocal(1);
yrel = relposlocal(2);
zrel = relposlocal(3);
xysq = xrel^2 + yrel^2;
xyzsq = xrel^2 + yrel^2 + zrel^2;
A = zeros(3,3);

if xyzsq == 0 % The target and the sensor are colocated
    A = laxes * [zeros(2,3); ones(1,3)];
elseif xysq == 0 % The target has the same (x,y) but not the same z
    x = -1 / zrel * 180 / pi;
    A = laxes * [zeros(1,3); [x x 0]; [0 0 1]];
else % The normal case
    % Since the local coordinates depend on the global coordinates, need to
    % use the chain rule to account for the partial derivatives of local
    % coordinates with respect to global coordinates.
    A(1,:) = laxes * [-yrel;        xrel;       0]/xysq;
    A(2,:) = laxes * [-xrel*zrel;  -yrel*zrel;  xysq]/sqrt(xysq)/xyzsq;
    A(3,:) = laxes * [xrel;         yrel;       zrel]/sqrt(xyzsq);
    A(1:2, :) = A(1:2, :) * 180 / pi; %convert to degrees
end