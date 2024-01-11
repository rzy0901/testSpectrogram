function [leftLegScatterersPos, rightLegScatterersPos] = getBicyclistLegs(upperLegVec,lowerLegVec,angUpperL,angLowerL,angUpperR,angLowerR)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen


%% Left leg
% Left leg scatterers
leftLegScatterersPos = calcLegScattererPos(angUpperL,angLowerL,...
    upperLegVec,lowerLegVec);


%% Right leg
% Right leg scatterers
rightLegScatterersPos = calcLegScattererPos(angUpperR,angLowerR,...
    upperLegVec,lowerLegVec);
end


%% Supporting functions
function legScatterersPos = calcLegScattererPos(angUpper,angLower,upperLegVec,lowerLegVec)
% Setup
nPtsUpper = length(upperLegVec);
nPtsLower = length(lowerLegVec);
nPtsTotal = nPtsUpper + nPtsLower; 
legScatterersPos = zeros(3,nPtsTotal,length(angUpper));

% Upper leg
if coder.target('MATLAB')
    cosUp = repmat(cos(angUpper),nPtsUpper,1);
    sinUp = repmat(sin(angUpper),nPtsUpper,1);
    legScatterersPos(1,1:nPtsUpper,:) = upperLegVec.'.*cosUp;
    legScatterersPos(3,1:nPtsUpper,:) = upperLegVec.'.*sinUp;
else
    for ii = 1:length(angUpper)
        cosUp = cos(angUpper(ii));
        sinUp = sin(angUpper(ii));
        legScatterersPos(1,1:nPtsUpper,ii) = upperLegVec.*cosUp;
        legScatterersPos(3,1:nPtsUpper,ii) = upperLegVec.*sinUp;
    end
end

% Lower leg
if coder.target('MATLAB')
    cosLow = repmat(cos(angLower),nPtsLower,1);
    sinLow = repmat(sin(angLower),nPtsLower,1);
    legScatterersPos(1,nPtsUpper+1:end,:) = -lowerLegVec.'.*cosLow;
    legScatterersPos(3,nPtsUpper+1:end,:) = -lowerLegVec.'.*sinLow;
else
    for ii = 1:length(angUpper)
        cosLow = cos(angLower(ii));
        sinLow = sin(angLower(ii));
        legScatterersPos(1,nPtsUpper+1:end,ii) = -lowerLegVec.*cosLow;
        legScatterersPos(3,nPtsUpper+1:end,ii) = -lowerLegVec.*sinLow;
    end
end

end

