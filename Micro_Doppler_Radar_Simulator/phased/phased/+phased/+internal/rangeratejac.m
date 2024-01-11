function jacobianRangeRateRow = rangeratejac(state,sensorpos, sensorvel, planar)
%This function is for internal use only. It may be removed in the future.

% rangeratejac Calculate the Jacobian of the range rate with respect to the
% state, in local coordinates defined by sensorpos and sensor vel. planar
% is a flag = 1 if the state is of a planar motion and 0 if 3D motion. It's
% only needed if state has 6 states.
% This function is for use only with CV, CA and CT motion models.
% state = target state
% sensorpos = [xr;yr;zr] - position of the sensor
% sensorvel = [vxr;vyr;vzr] - velocity of the sensor
% planar = 1 if state has 6 states and it's a CA model

%   Copyright 2016 The MathWorks, Inc.

% internal function, no error checking is performed

%#codegen 

% First, process the state correctly
switch numel(state)
    case 2 %state = [x;vx]
        tgtpos = [state(1); 0; 0];
        tgtvel = [state(2); 0; 0];
        elementsToTake = 1:2; %which elements of the final vector to return
    case 3 %state = [x;vx;ax]
        tgtpos = [state(1); 0; 0];
        tgtvel = [state(2); 0; 0];
        elementsToTake = 1:3; %which elements of the final vector to return
    case 4 %state = [x;vx;y;vy]
        tgtpos = [state(1); state(3); 0];
        tgtvel = [state(2); state(4); 0];
        elementsToTake = [1:2,4:5]; %which elements of the final vector to return
    case 5 %state = [x;vx;y;vy;omega]
        tgtpos = [state(1); state(3); 0];
        tgtvel = [state(2); state(4); 0];
        elementsToTake = [1:2,4:5]; %which elements of the final vector to return
    case 6
        if planar %state = [x;vx;ax;y;vy;ay]
            tgtpos = [state(1); state(4); 0];
            tgtvel = [state(2); state(5); 0];
            elementsToTake = 1:6; %which elements of the final vector to return
        else %state = [x;vx;y;vy;z;vz]
            tgtpos = [state(1); state(3); state(5)];
            tgtvel = [state(2); state(4); state(6)];
            elementsToTake = [1:2,4:5,7:8]; %which elements of the final vector to return
        end
    case 9 %state = [x;vx;ax;y;vy;ay;z;vz;az]
        tgtpos = [state(1); state(4); state(7)];
        tgtvel = [state(2); state(5); state(8)];
        elementsToTake = 1:9; %which elements of the final vector to return
end

% Calculate the derivatives. Note that if r == 0, the Jacobian of rdot
% should be calculated as [Delta_xdot, Delta_ydot, Delta_zdot] / rdot,
% and if rdot == 0 as well, it is simply 1 for each relative velocity
% compnonent
relpos = tgtpos - sensorpos;
relvel = tgtvel - sensorvel;
rnorm = norm(relpos);
relvelnorm = norm(relvel);
if rnorm > 0 % Regular calculation
    rdot2x = (relvel(1) * rnorm - (relpos'*relvel) * relpos(1) / rnorm)...
        / rnorm^2;
    rdot2xdot = relpos(1) / rnorm;
    rdot2y = (relvel(2) * rnorm - (relpos'*relvel) * relpos(2) / rnorm)...
        / rnorm^2;
    rdot2ydot = relpos(2) / rnorm;
    rdot2z = (relvel(3) * rnorm - (relpos'*relvel) * relpos(3) / rnorm)...
        / rnorm^2;
    rdot2zdot = relpos(3) / rnorm;
    jacobianRangeRateRow = [rdot2x, rdot2xdot, 0, rdot2y, rdot2ydot, 0, rdot2z, rdot2zdot, 0];
elseif relvelnorm % Take the components of the relative velocity normalized by its norm
    jacobianRangeRateRow = [0, relvel(1)/relvelnorm, 0, 0, relvel(2)/relvelnorm, 0, 0, relvel(3)/relvelnorm, 0];
else % 1 for each component of the relative velocity
    jacobianRangeRateRow = [0, 1, 0, 0, 1, 0, 0, 1, 0];
end
jacobianRangeRateRow = jacobianRangeRateRow(elementsToTake);
end