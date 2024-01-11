function [initialPosition,mounts] = getBicyclistPedals(spacing,mounts) 
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

%% Initialize
pedalBarLength = 0.2; 
nPtsPedalBar = ceil(pedalBarLength/spacing)+1; 
pedalLength = 0.125; 
nPtsPedal = ceil(pedalLength/spacing)+1; % Number of points for each pedal

nPtsTotal = nPtsPedalBar + 2*nPtsPedal;
initialPosition = zeros(3,nPtsTotal); 

%% Build 
% Pedal bar
posPedalBar = zeros(3,nPtsPedalBar);
posPedalBar(1,:) = linspace(-pedalBarLength/2,pedalBarLength/2,nPtsPedalBar);
initialPosition(:,1:nPtsPedalBar) = posPedalBar; 
idxLast = nPtsPedalBar; 

% Pedals
posPedals = zeros(3,nPtsPedal*2); 
posPedals(1,:) = [posPedalBar(1,1).*ones(1,nPtsPedal) posPedalBar(1,end).*ones(1,nPtsPedal)]; 
posPedals(2,:) = [(1:nPtsPedal).*pedalLength./(nPtsPedal) -(1:nPtsPedal).*pedalLength./(nPtsPedal)]; 
initialPosition(:,idxLast+1:idxLast+2*nPtsPedal) = posPedals;

% Pedal mount position for lower leg 
mounts.Feet =  [ [posPedals(1,nPtsPedal); posPedals(2,nPtsPedal); posPedals(3,nPtsPedal)] ...
    [posPedals(1,2*nPtsPedal); posPedals(2,2*nPtsPedal); posPedals(3,2*nPtsPedal)] ]; 
mounts.IdxFeet = [idxLast+nPtsPedal,idxLast+2*nPtsPedal];

%% Mount
% Move pedals to pedal mounting position
initialPosition = bsxfun(@plus,initialPosition,mounts.Pedals);

% Move feet to the position of the pedal
mounts.Feet = bsxfun(@plus,mounts.Feet,mounts.Pedals);

