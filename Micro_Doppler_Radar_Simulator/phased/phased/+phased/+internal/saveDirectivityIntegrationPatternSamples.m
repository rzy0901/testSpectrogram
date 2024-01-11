function saveDirectivityIntegrationPatternSamples
% This function is for internal use only. It may be removed in the future.

%getDirectivityIntegrationPatternSamples
%   [az,el,ds] = getDirectivityIntegrationPatternSamples returns the
%   sampling points for calculating integrated field pattern. These points
%   are roughly uniformly distributed over the sphere.

load([matlabroot ('/toolbox/antenna/antenna/+em/@EmStructures/spherenew.mat')]);
c_s = bsxfun(@rdivide,Center_s,sqrt(sum(Center_s.^2)));
%c_s = Center_s;
[phi,theta] = cart2sph(c_s(1,:),c_s(2,:),c_s(3,:));
az_samp = wrapTo180(phased.internal.rad2deg(phi));
el_samp = phased.internal.rad2deg(theta);
area_samp = Area_s;
save IntegratedFieldSamplePoints az_samp el_samp area_samp
