function [rcvAng,srcAng,pathLength,numRefl,TxPosImage,TxVelImage]=imageMethod(ChannelDepth,TxPos,RxPos,TxVel,numPaths)
%This class is for internal use only. It may be removed in the future.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

% Paths consist of a direct path and pairs of images (top and bottom
% reflection). For odd numPaths, the number of image pairs is
% (numPaths-1)/2. If an even number of paths is requested, compute
% numPaths/2 image pairs. The last image will be discarded.
if mod(numPaths,2)
  numImagePairs = (numPaths-1)/2; % num of sets of 2 images
  rcvAng = zeros(2,numPaths);
  srcAng = zeros(2,numPaths);
else
 numImagePairs = numPaths/2;
 rcvAng = zeros(2,numPaths+1);
 srcAng = zeros(2,numPaths+1);
end

% Source and receiver positions
zSrc = TxPos(3);
zRcv = RxPos(3);
zChanTop = 0;
zChanBot = -ChannelDepth;

% Transmit/Receiver distance in x-y plane
dXY=norm(TxPos(1:2)-RxPos(1:2),2);

% Compute the receiver angle and path length for each set of images.
% thetaTop and lengthTop will contain paths that are first incident on the
% top surface. Compute the z coordinate of each image.
zTop=zeros(1,numImagePairs);
zBot=zeros(1,numImagePairs);
if numPaths > 1
    zTop(1) = zSrc + 2*(zChanTop-zSrc);
    zBot(1) = zSrc + 2*(zChanBot-zSrc);
end

for imageIndx=2:numImagePairs
    % Alternate between creating images wrt top and bottom of channel.
    if(mod(imageIndx,2)~=0)
        zTop(imageIndx) = zTop(imageIndx-1) + 2*(zChanTop - zTop(imageIndx-1));
        zBot(imageIndx) = zBot(imageIndx-1) + 2*(zChanBot - zBot(imageIndx-1));
    else
        zTop(imageIndx) = zTop(imageIndx-1) + 2*(zChanBot - zTop(imageIndx-1));
        zBot(imageIndx) = zBot(imageIndx-1) + 2*(zChanTop - zBot(imageIndx-1));
    end   
end

% Assign image positions
TxPosImage = repmat(TxPos,1,2*numImagePairs+1);
TxPosImage(3,2:end) = [zTop zBot];

% Compute the elevation angle wrt the receiver of each path.
thetaTop=atan((zTop-zRcv)/dXY)*180/pi;
thetaBot=atan((zBot-zRcv)/dXY)*180/pi;

% Compute the path length.
lengthTop=sqrt((zTop-zRcv).^2 + dXY^2); 
lengthBot=sqrt((zBot-zRcv).^2 + dXY^2); 

% Create vectors of arrival (receiver) angles, and path
% lengths. Sources angles are equal in magnitude to receiver angles but have
% sign corresponding the first reflection interface (top or bottom).
[directR,directA] = rangeangle(RxPos,TxPos); 
[~,directArcv] = rangeangle(TxPos,RxPos); 
rcvAng(1,:) = directArcv(1);
srcAng(1,:) = directA(1);
rcvAng(2,:) = [directArcv(2) thetaTop thetaBot];
srcAng(2,:) = [directA(2) abs(rcvAng(2,2:numImagePairs+1)) -1*abs(rcvAng(2,numImagePairs+2:end))];
pathLength= [directR lengthTop lengthBot];
           
% Compute the number of top and bottom reflections for each path.
c_bounce = zeros(1,2*numImagePairs);
c_bounce(1:2:2*numImagePairs-1) = 1:numImagePairs;
c_bounce(2:2:2*numImagePairs) = 1:numImagePairs;
c_bounce = c_bounce(1:numImagePairs);
numTopRefl = [0 c_bounce 0 c_bounce(1:end-1)]; % Number of top bounces for each path 
numTopRefl = numTopRefl(1:2*numImagePairs+1);
numBotRefl = [0 0 c_bounce(1:end-1) c_bounce]; % Number of bottom bounces for each path 
numBotRefl = numBotRefl(1:2*numImagePairs+1);

% Velocity of images in z direction is determined by the number of bounces
TxVelImage = repmat(TxVel,1,2*numImagePairs+1);
TxVelImage(3,:) = TxVelImage(3,:).*(-1).^(numTopRefl+numBotRefl);

%Output requested number of paths
rcvAng = rcvAng(:,1:numPaths);
srcAng = srcAng(:,1:numPaths);
pathLength = pathLength(:,1:numPaths);
numRefl = [numTopRefl(1:numPaths); numBotRefl(1:numPaths)];
TxPosImage = TxPosImage(:,1:numPaths);
TxVelImage = TxVelImage(:,1:numPaths);
