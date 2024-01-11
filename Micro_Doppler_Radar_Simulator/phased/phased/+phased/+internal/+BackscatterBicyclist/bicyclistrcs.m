function rcsmat = bicyclistrcs()
%This function is for internal use only. It may be removed in the future.
%   Use mean values of bicycle sectors from reference [3] to build default rcs
%   vector for the bicyclist radar target.
%
%   The default radar cross section (RCS) of the bicyclist is defined using
%   the same method as [1], using measurements from [3]. It has been shown
%   that the RCS value exhibits a dependency on the viewing angle. Thus,
%   the individual point scatterers get an RCS value approximated by the
%   angle of view of the whole cyclist.

%   Copyright 2019 The MathWorks, Inc.

%   Reference
%   [1] Stolz, M. et al. "Multi-Target Reflection Point Model of Cyclists
%   for Automotive Radar." 2017 European Radar Conference (EURAD),
%   Nuremberg, 2017, pp. 94-97.
%
%   [2] Chen, V., D. Tahmoush, and W. J. Miceli. Radar Micro-Doppler
%   Signatures: Processing and Applications. The Institution of Engineering
%   and Technology: London, 2014.
%
%   [3] Belgiovane, D., and C. C. Chen. "Bicycles and  Human Riders
%   Backscattering at 77 GHz for Automotive Radar." 2016 10th European
%   Conference on Antennas and Propagation (EuCAP), Davos, 2016, pp. 1-5.

%#codegen

% Define RCS values for 1-180
frontRearRCS = 10^(-11.57/10); % 1-30 deg (front), 151-180 deg (rear)
sideRCS = 10^(-8.1/10); % 61-120 deg
obliqueRCS = 10^(-13.24/10); % 31-60 deg, 121-150 deg
rcs = [obliqueRCS frontRearRCS obliqueRCS sideRCS obliqueRCS frontRearRCS obliqueRCS sideRCS obliqueRCS frontRearRCS obliqueRCS];
ang = [-225 -180 -135 -90 -45 0 45 90 135 180 225]; 
angInterp = -180:180; 
rcsmat = interp1(ang,rcs,angInterp,'pchip'); 

end