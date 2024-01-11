function [initialPosition,rearWheelMountPos,pedalMountPos] = getBicyclistFrame(spacing) 
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

%% Initialize 
headTubeForkLength = 0.564;
nPtsHeadTubeFork = ceil(headTubeForkLength/spacing)+1; 
topTubeLength = 0.704;
nPtsTopTube = ceil(topTubeLength/spacing)+1; 
handleBarConnectorLength = 0.125;
nPtsConnector = ceil(handleBarConnectorLength/spacing)+1; 
handleBarLength = 0.25;
nPtsHandleBar = ceil(handleBarLength/spacing)+1; 
seatTubeLength = 0.566; 
nPtsSeatTube = ceil(seatTubeLength/spacing)+1; 
seatLength = 0.2;
nPtsSeat = ceil(seatLength/spacing)+1; 
if mod(nPtsSeat,2)~=0
    nPtsSeat = nPtsSeat+1;
end
seatStaysLength = 0.563; 
nPtsSeatStays = ceil(seatStaysLength/spacing)+1; 
downTubeLength = 0.719;
nPtsDownTube = ceil(downTubeLength/spacing)+1; 

nPtsTotal = nPtsHeadTubeFork + nPtsTopTube + nPtsConnector + nPtsHandleBar + nPtsSeatTube + nPtsSeat + nPtsSeatStays + nPtsDownTube; 
initialPosition = zeros(3,nPtsTotal);


%% Bike Frame
% Local origin is the part of the front fork that attaches to the center of
% the front wheel 

% Head tube + forks 
headTubeForkAng = deg2rad(30.259); % Angle from z-axis measured counter-clockwise 
posHeadTubeFork = zeros(3,nPtsHeadTubeFork);
for iP = 2:nPtsHeadTubeFork
    posHeadTubeFork(1,iP) = posHeadTubeFork(1,iP-1)+headTubeForkLength/(nPtsHeadTubeFork-1)*sin(headTubeForkAng); 
    posHeadTubeFork(3,iP) = posHeadTubeFork(3,iP-1)+headTubeForkLength/(nPtsHeadTubeFork-1)*cos(headTubeForkAng); 
end
posHeadTubeFork(1,:) = -posHeadTubeFork(1,:);
idxLast = nPtsHeadTubeFork; 
initialPosition(:,1:idxLast) = posHeadTubeFork;

% Top tube
posTopTube = zeros(3,nPtsTopTube);
for iP = 1:nPtsTopTube
    posTopTube(1,iP) = posHeadTubeFork(1,end)-topTubeLength/nPtsTopTube*(iP);
    posTopTube(3,iP) = posHeadTubeFork(3,end);
end
initialPosition(:,idxLast+1:idxLast+nPtsTopTube) = posTopTube; 
idxLast = idxLast+nPtsTopTube; 

% Handle bar connector
handleBarConnectorAng = deg2rad(25); 
posHandleBarConnector = zeros(3,nPtsConnector); 
for iP = 1:nPtsConnector
posHandleBarConnector(1,iP) = posHeadTubeFork(1,end)+...
    handleBarConnectorLength/nPtsConnector*(iP).*cos(handleBarConnectorAng);
posHandleBarConnector(3,iP) = posHeadTubeFork(3,end)+...
    handleBarConnectorLength/nPtsConnector*(iP).*sin(handleBarConnectorAng);
end
initialPosition(:,idxLast+1:idxLast+nPtsConnector) = posHandleBarConnector; 
idxLast = idxLast+nPtsConnector; 

% Handle bar 
posHandleBar = zeros(3,nPtsHandleBar); 
posHandleBar(1,1:end) = posHandleBarConnector(1,end); 
posHandleBar(2,1:end) = linspace(-handleBarLength/2,handleBarLength/2,nPtsHandleBar);
posHandleBar(3,1:end) = posHandleBarConnector(3,end); 
initialPosition(:,idxLast+1:idxLast+nPtsHandleBar) = posHandleBar; 
idxLast = idxLast+nPtsHandleBar; 

% Seat tube 
seatTubeAng = deg2rad(20); % Angle from negative z axis measured counter-clockwise
posSeatTube = zeros(3,nPtsSeatTube);
for iP = 1:nPtsSeatTube
    posSeatTube(1,iP) = posTopTube(1,end)+seatTubeLength/nPtsSeatTube*(iP)*sin(seatTubeAng);
    posSeatTube(3,iP) = posTopTube(3,end)-seatTubeLength/nPtsSeatTube*(iP)*cos(seatTubeAng);
end
initialPosition(:,idxLast+1:idxLast+nPtsSeatTube) = posSeatTube; 
idxLast = idxLast+nPtsSeatTube; 
pedalMountPos = posSeatTube(:,end); 

% Seat 
posSeat = zeros(3,nPtsSeat); 
posSeat(1,1:nPtsSeat-1) = posTopTube(1,end)+linspace(-seatLength/2,seatLength/2,nPtsSeat-1);
posSeat(1,nPtsSeat) = posTopTube(1,end); % Connector of seat, x dim
posSeat(3,1:nPtsSeat-1) = posTopTube(3,end)+0.0612;
posSeat(3,nPtsSeat) = posTopTube(3,end)+0.0306; % Connector of seat, z dim
initialPosition(:,idxLast+1:idxLast+nPtsSeat) = posSeat; 
idxLast = idxLast+nPtsSeat; 

% Seat stays 
seatStaysAng = deg2rad(30); % Angle from negative z axis measured clockwise
posSeatStays = zeros(3,nPtsSeatStays);
for iP = 1:nPtsSeatStays
    posSeatStays(1,iP) = posTopTube(1,end)-seatStaysLength/nPtsSeatStays*(iP)*sin(seatStaysAng);
    posSeatStays(3,iP) = posTopTube(3,end)-seatStaysLength/nPtsSeatStays*(iP)*cos(seatStaysAng);
end
initialPosition(:,idxLast+1:idxLast+nPtsSeatStays) = posSeatStays; 
rearWheelMountPos = posSeatStays(:,end); 
idxLast = idxLast+nPtsSeatStays; 

% Down tube
downTubeAng = deg2rad(40); % Angle from x-axis measured counter-clockwise
posDownTube = zeros(3,nPtsDownTube);
for iP = 1:nPtsDownTube
    posDownTube(1,iP) = posSeatTube(1,end)+downTubeLength/(nPtsDownTube+1)*(iP)*cos(downTubeAng);
    posDownTube(3,iP) = posSeatTube(3,end)+downTubeLength/(nPtsDownTube+1)*(iP)*sin(downTubeAng);
end
initialPosition(:,idxLast+1:idxLast+nPtsDownTube) = posDownTube; 

end
