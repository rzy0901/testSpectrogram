function [initialPosition,hipMountPos] =  getBicyclistRiderUpperBody(frameScatterers,spacing)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

%% Initialize 
torsoLength = 0.5;
nPtsTorso = ceil(torsoLength/spacing)+1;
headLength = 0.125;
nPtsHead = ceil((2*pi*headLength/2)/spacing)+1; % Number of points per slice; 2 slices
hipsLength = 0.25; 
nPtsHips = ceil(hipsLength/spacing)+1; 
shoulderLength = 0.25; 
nPtsShoulders = ceil(shoulderLength/spacing)+1;
upperArmLength = 0.3; 
nPtsUpperArms = ceil(upperArmLength/spacing)+1; % Number of points for each upper arm; 2 upper arms
lowerArmLength = 0.3; 
nPtsLowerArms = ceil(lowerArmLength/spacing)+1; % Number of points for each lower arm; 2 lower arms

nPtsTotal = nPtsTorso + 2*nPtsHead + nPtsHips + nPtsShoulders + 2*nPtsUpperArms + 2*nPtsLowerArms;
initialPosition = zeros(3,nPtsTotal);

%% Identify location of seat 
[maxVal,~] = max(frameScatterers(3,:));
idxSeat = find(frameScatterers(3,:)==maxVal); 
idxSeat = median(idxSeat); 

%% Build 
% Bicyclist torso
torsoAng = 45; 
posTorso = zeros(3,nPtsTorso); 
posTorso(1,:) = frameScatterers(1,idxSeat)+torsoLength/(nPtsTorso-1).*(0:nPtsTorso-1).*cosd(torsoAng); 
posTorso(3,:) = frameScatterers(3,idxSeat)+torsoLength/(nPtsTorso-1).*(0:nPtsTorso-1).*sind(torsoAng); 
initialPosition(:,1:nPtsTorso) = posTorso; 
idxLast = nPtsTorso; 

% Bicyclist head
headAng = torsoAng; 
posHead = zeros(3,nPtsHead*2); % Times 2 for 2 slices 
headRadius = headLength/2;
headCtr = [posTorso(1,end)+headRadius*cosd(headAng) ... % x
    0 ... % y 
    posTorso(3,end)+headRadius*sind(headAng)]; % z 

% X-Z slice
ang = 0:360/nPtsHead:(360-360/nPtsHead); 
posHead(1,1:nPtsHead) = headCtr(1)+headRadius.*cosd(ang);
posHead(2,1:nPtsHead) = headCtr(2); 
posHead(3,1:nPtsHead) = headCtr(3)+headRadius.*sind(ang); 

% Y-Z slice
posHead(1,nPtsHead+1:2*nPtsHead) = headCtr(1); 
posHead(2,nPtsHead+1:2*nPtsHead) = headCtr(2)+headRadius.*cosd(ang); 
posHead(3,nPtsHead+1:2*nPtsHead) = headCtr(3)+headRadius.*sind(ang); 
initialPosition(:,idxLast+1:idxLast+2*nPtsHead) = posHead; 
idxLast = idxLast+2*nPtsHead; 

% Bicyclist hips 
posHips = zeros(3,nPtsHips); 
posHips(1,:) = frameScatterers(1,idxSeat); 
posHips(2,:) = linspace(-hipsLength/2,hipsLength/2,nPtsHips); 
posHips(3,:) = frameScatterers(3,idxSeat); 
initialPosition(:,idxLast+1:idxLast+nPtsHips) = posHips; 
idxLast = idxLast+nPtsHips; 

% Hip mount position for upper leg 
hipMountPos =  [ [posHips(1,end); posHips(2,end); posHips(3,end)] ...
    [posHips(1,1); posHips(2,1); posHips(3,1)] ]; 

% Bicyclist shoulders 
posShoulders = zeros(3,nPtsShoulders); 
posShoulders(1,:) = posTorso(1,end); 
posShoulders(2,:) = linspace(-shoulderLength/2,shoulderLength/2,nPtsShoulders); 
posShoulders(3,:) = posTorso(3,end); 
initialPosition(:,idxLast+1:idxLast+nPtsShoulders) = posShoulders; 
idxLast = idxLast+nPtsShoulders; 

% Bicyclist upper arms
upperArmAng = 45; 
posUpperArms = zeros(3,nPtsUpperArms*2); 
posUpperArms(1,:) = [posShoulders(1,1)+upperArmLength/(nPtsUpperArms-1).*(0:nPtsUpperArms-1).*cosd(upperArmAng)...
    posShoulders(1,1)+upperArmLength/(nPtsUpperArms-1).*(0:nPtsUpperArms-1).*cosd(upperArmAng)]; 
posUpperArms(2,:) = [posShoulders(2,1).*ones(1,nPtsUpperArms) posShoulders(2,end).*ones(1,nPtsUpperArms)];
posUpperArms(3,:) = [posShoulders(3,1)-upperArmLength/(nPtsUpperArms-1).*(0:nPtsUpperArms-1).*sind(upperArmAng)...
    posShoulders(3,1)-upperArmLength/(nPtsUpperArms-1).*(0:nPtsUpperArms-1).*sind(upperArmAng)];
initialPosition(:,idxLast+1:idxLast+nPtsUpperArms*2) = posUpperArms; 
idxLast = idxLast+nPtsUpperArms*2; 

% Bicyclist lower arms
lowerArmAng = 15; 
posLowerArms = zeros(3,nPtsLowerArms*2); 
posLowerArms(1,:) = [posUpperArms(1,end)+lowerArmLength/(nPtsLowerArms-1).*(0:nPtsLowerArms-1).*cosd(lowerArmAng)...
    posUpperArms(1,end)+lowerArmLength/(nPtsLowerArms-1).*(0:nPtsLowerArms-1).*cosd(lowerArmAng)]; 
posLowerArms(2,:) = [posUpperArms(2,1).*ones(1,nPtsLowerArms) posUpperArms(2,end).*ones(1,nPtsLowerArms)];
posLowerArms(3,:) = [posUpperArms(3,end)-lowerArmLength/(nPtsLowerArms-1).*(0:nPtsLowerArms-1).*sind(lowerArmAng)...
    posUpperArms(3,end)-lowerArmLength/(nPtsLowerArms-1).*(0:nPtsLowerArms-1).*sind(lowerArmAng)];
initialPosition(:,idxLast+1:idxLast+nPtsLowerArms*2) = posLowerArms; 

