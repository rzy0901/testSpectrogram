function [frameAndRiderUpperBody, pedals, legs, ...
    frontWheel, rearWheel, mounts] = getBicyclistParts(wheelDiameter,numWheelSpokes,spacing,...
    upperLegVec,lowerLegVec)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

%% Initialize parts in local coordinates
[frameAndRiderUpperBody,mounts] = phased.internal.BackscatterBicyclist.getBicyclistFrameAndRiderUpperBody(spacing);
[pedals,mounts] = phased.internal.BackscatterBicyclist.getBicyclistPedals(spacing,mounts);

%% Add wheels
wheelScatterers = phased.internal.BackscatterBicyclist.getWireSpokedWheel(wheelDiameter,numWheelSpokes,spacing);
frontWheel = bsxfun(@plus, wheelScatterers, mounts.FrontWheel);
rearWheel = bsxfun(@plus, wheelScatterers, mounts.RearWheel);

%% Create legs
[angUpperL,angLowerL,angUpperR,angLowerR] = ...
    phased.internal.BackscatterBicyclist.getBicyclistLegAngles(...
    spacing,mounts,upperLegVec,lowerLegVec); % Get angle of legs
[leftLeg, rightLeg] = phased.internal.BackscatterBicyclist.getBicyclistLegs(...
    upperLegVec,lowerLegVec,angUpperL,angLowerL,angUpperR,angLowerR); % Calculate position
leftLeg = squeeze(leftLeg); 
rightLeg = squeeze(rightLeg); 

% Shift legs into position about hips 
nPtsUpper = length(upperLegVec);
leftLeg(:,1:nPtsUpper) = bsxfun(@plus, leftLeg(:,1:nPtsUpper), mounts.Hip(:,1));
leftLeg(:,(nPtsUpper+1):end) = bsxfun(@plus, leftLeg(:,(nPtsUpper+1):end), leftLeg(:,nPtsUpper));
rightLeg(:,1:nPtsUpper) = bsxfun(@plus, rightLeg(:,1:nPtsUpper), mounts.Hip(:,2));
rightLeg(:,(nPtsUpper+1):end) = bsxfun(@plus, rightLeg(:,(nPtsUpper+1):end), rightLeg(:,nPtsUpper));
legs = [leftLeg rightLeg];


%% Center 
transvec = [0.634852524625840; 0; wheelDiameter/2]; % Center the bike

% Bicycle frame + rider upper body 
frameAndRiderUpperBody = bsxfun(@plus, frameAndRiderUpperBody, transvec);

% Pedals 
pedals = bsxfun(@plus, pedals, transvec);

% Front Wheel
frontWheel = bsxfun(@plus, frontWheel, transvec); 

% Rear wheel
rearWheel = bsxfun(@plus, rearWheel, transvec);

% Legs
legs = bsxfun(@plus, legs, transvec);

% Translate mounts
mounts.FrontWheel = bsxfun(@plus, mounts.FrontWheel, transvec); 
mounts.RearWheel = bsxfun(@plus, mounts.RearWheel, transvec); 
mounts.Pedals = bsxfun(@plus, mounts.Pedals, transvec); 
mounts.Feet = bsxfun(@plus, mounts.Feet, transvec);
mounts.Hip = bsxfun(@plus, mounts.Hip, transvec);

end
