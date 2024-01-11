function initialPosition = getWireSpokedWheel(wheelDiameter,numWheelSpokes,spacing)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

%% Initial Position
% Local origin is center of wheel 
r = wheelDiameter/2;
nRadii = floor(r/spacing)+1; 
radius = r:-spacing:(r-(nRadii-1)*spacing); 
nPoints = [round(2*pi*radius(1)/spacing)*2 numWheelSpokes*ones(1,numel(radius)-1)];
initialPosition = zeros(3,sum(nPoints)); 

% Define indices 
startIndices = zeros(size(nPoints)); 
endIndices = zeros(size(nPoints)); 
for m = 1:numel(radius) 
    if m == 1
        startIndices(m) = 1;
    else
        startIndices(m) = endIndices(m-1)+1;
    end
    endIndices(m) = (startIndices(m)-1)+nPoints(m); 
end

% Define inner circles 
for m = 1:numel(radius)
    res = 360/nPoints(m); 
    theta = 0:res:(360-res);
    initialPosition(1,startIndices(m):endIndices(m)) = radius(m).*cosd(theta); % x-coordinate
    initialPosition(3,startIndices(m):endIndices(m)) = radius(m).*sind(theta); % z-coordinate
end

end


