function azel = dirvec2azel(dirvec)
%This function is for internal use only. It may be removed in the future.

%DIRVEC2AZEL    Convert directional vectors to azimuth/elevation angles
%   AZEL = phased.internal.dirvec2azel(DIRVEC) converts the direction
%   vectors specified in DIRVEC to the corresponding angle AZEL in the form
%   of [azimuth; elevation] (in degrees). DIRVEC is a 3-row matrix and AZEL
%   is a 2-row matrix whose columns contain the angles of
%   corresponding directional vectors in DIRVEC.
%
%   This function supports single and double precision for input data.
%   If the input data, DIRVEC, is single precision, the output data
%   is single precision. If the input data, DIRVEC, is double precision, 
%   the output data is double precision.
%
%   The directional vector is defined as [cos(el)cos(az); cos(el)sin(az);
%   sin(el)].
%
%   % Example:
%   %   Compute the angles corresponding to the direction vector of
%   %   [0.5; 0.5; 0.7071].
%
%   dirvec = phased.internal.dirvec2azel([0.5;0.5;0.7071])
%
%   See also phased, arbazel2azel, azel2dirvec

%   Copyright 2018 The MathWorks, Inc.

%#ok<*EMCA>
%#codegen

azel = zeros(2,size(dirvec,2),class(dirvec));
nzidx = dirvec(1,:)~=0;
[th,phi] = cart2sph(dirvec(1,nzidx),dirvec(2,nzidx),dirvec(3,nzidx));
azel(:,nzidx) = [th;phi];
temp = dirvec(3,~nzidx);
azel(2,~nzidx) = asin(max(-1,min(temp,1)));  % guard valid sin value
temp = dirvec(2,~nzidx)./cos(azel(2,~nzidx));
azel(1,~nzidx) = asin(max(-1,min(temp,1)));  % guard valid sin value
azel = phased.internal.rad2deg(azel);


