%This script is for internal use only. It may be removed in the future.

% Helper script creates data for command-line help for
% clusterDBSCAN.
%
%   The maximum unambiguous range and Doppler in this example is 20 m
%   and 30 Hz, respectively. Generate data with the following extended
%   targets:
%      - 1 unambiguous target located at (10, 15)
%      - 1 ambiguous target in Doppler located at (10, -30)
%      - 1 ambiguous target in range located at (20, 15)
%      - 1ambiguous target in range and Doppler located at (20, 30)
%      - 5 false alarms 

%   Copyright 2019 The MathWorks, Inc.

% Define targets
maxRange = 20;
minDoppler = -30; 
maxDoppler = 30;
targetsXY = [...
    maxRange/2, maxDoppler/2;...  % Unambiguous target
    maxRange/2, -maxDoppler;...   % Ambiguous target in Doppler
    maxRange, maxDoppler/2;...    % Ambiguous target in Range
    maxRange, maxDoppler];        % Ambiguous target in Range and Doppler

% Build unambiguous data for clustering
nDetsPerTarget = [20 10 13 17];
nTgts = size(targetsXY,1);
x = [zeros(sum(nDetsPerTarget),2); 20*rand(5,2)];
maxVal = 1;
minVal = -1;
rangeVals = maxVal-minVal;
idxStart = 1;
for iT = 1:nTgts
    nDetThisTgt = nDetsPerTarget(iT);
    idxEnd = idxStart+nDetThisTgt-1;
    x(idxStart:idxEnd,:) = ...
        (rangeVals*rand(nDetThisTgt,2)+minVal)+targetsXY(iT,:);
    idxStart = idxEnd+1;
end

% Wrap detections according to maximum unambiguous range and Doppler
x(:,1) = mod(x(:,1),maxRange);
x(:,2) = mod(x(:,2)+maxDoppler,2*maxDoppler)...
    - maxDoppler;

saveName = fullfile(matlabroot,'toolbox','phased','phased','dataClusterDBSCAN.mat');
save(saveName,'x','maxRange','minDoppler','maxDoppler');