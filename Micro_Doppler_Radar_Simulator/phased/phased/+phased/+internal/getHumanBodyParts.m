function [bodyjointspos,bodypartslen,bodypartssz,bodypartsax] = getHumanBodyParts(H)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2018-2019 The MathWorks, Inc.

%#codegen

% H: height
% pos: [ltoe toe lankle rankle lknee rknee lhip rhip lhand rhand lelbow
% relbow lshoulder rshoulder head neck lclorigin]

% set origin at hip
torsolength = 0.288*H;
hipheight = (0.245+0.246+0.039)*H; %upper+lower legs + foot height

% head
headlen = 0.13*H;
headpos = [0;0;H-hipheight];
% neck
necklen = 0.052*H;
shoulderheight = H-headlen-necklen;
% shoulder
shoulderlen = (0.259/2)*H;
temp = [0 0 0;[1 0 -1]*shoulderlen;[1 1 1]*(shoulderheight-hipheight)];
shoulderpos = temp(:,[1 3]);
neckpos = temp(:,2);
% elbow
upperarmlen = 0.188*H;
elbowheight = shoulderheight-upperarmlen;
elbowpos = [0 0;[1 -1]*shoulderlen;[1 1]*(elbowheight-hipheight)];
% hand
lowerarmlen = 0.145*H;
handheight = elbowheight-lowerarmlen;
handpos = [0 0;[1 -1]*shoulderlen;[1 1]*(handheight-hipheight)];
% hip
hiplen = 0.191/2*H;
hipheight = shoulderheight-torsolength;
temp = [0 0 0;[1 0 -1]*hiplen;[1 1 1]*(hipheight-hipheight)];
hippos = temp(:,[1 3]);
lcloriginpos = temp(:,2);
% knee
upperleglen = 0.245*H;
kneeheight = hipheight-upperleglen;
kneepos = [0 0;[1 -1]*hiplen;[1 1]*(kneeheight-hipheight)];
% ankle
lowerleglen = 0.246*H;
ankleheight = kneeheight-lowerleglen;
anklepos = [0 0;[1 -1]*hiplen;[1 1]*(ankleheight-hipheight)];
% toe
footlen = 0.143*H;
toepos = [[1 1]*footlen;[1 -1]*hiplen;[1 1]*(0-hipheight)];

bodyjointspos = [toepos anklepos kneepos hippos handpos ...
    elbowpos shoulderpos headpos neckpos lcloriginpos];

% body parts length 
bodypartslen = [footlen lowerleglen upperleglen hiplen lowerarmlen upperarmlen ...
    shoulderlen headlen necklen torsolength ankleheight kneeheight elbowheight ...
    shoulderheight hipheight];

% body parts represented by spheroid, in the format of [a,c]
% foot, lower leg, upper leg, hip, lower arm, upper arm, shoulder, head,
% torso
bodypartssz = [0.5 footlen/2; 0.06 lowerleglen/2; 0.07 upperleglen/2; ...
    0.07 hiplen/2; 0.05 lowerarmlen/2; 0.06 upperarmlen/2; ...
    0.06 shoulderlen/2; 0.1 headlen/2; 0.15 torsolength/2].';

% assume by default spheroid is along z axis, 
% most body parts are in this direction, except the foot, the hip, and the
% shoulder

% Build body parts axes matrix (performed this way for codegen; addresses g1932019) 
footax = roty(90);
lowerlegax = eye(3);
upperlegax = eye(3);
hipax = rotx(-90);
lowerarmax = eye(3);
upperarmax = eye(3);
shoulderax = rotx(-90);
headax = eye(3);
torsoax = eye(3);
bodypartsax = zeros(3,3,9); 
bodypartsax(:,:,1) = footax;
bodypartsax(:,:,2) = lowerlegax;
bodypartsax(:,:,3) = upperlegax;
bodypartsax(:,:,4) = hipax;
bodypartsax(:,:,5) = lowerarmax;
bodypartsax(:,:,6) = upperarmax;
bodypartsax(:,:,7) = shoulderax;
bodypartsax(:,:,8) = headax;
bodypartsax(:,:,9) = torsoax;