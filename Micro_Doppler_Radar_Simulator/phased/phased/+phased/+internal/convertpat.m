function pat_newOut = convertpat(fh_conv,pat_old,ang_1_old,ang_2_old,ang_1,ang_2,toUVFlag)
% This function is for internal use only. It may be removed in the future.

%CONVERTPAT Convert pattern among different spaces
%   PAT_NEW = phased.internal.convertpat(...
%       ANGLECONVERT_FH,PAT_OLD,ANG1_OLD,ANG2_OLD,ANG1_NEW,ANG2_NEW,UVFLAG)
%   converts the PAT_OLD defined in the grid of ANG1_OLD (in degrees) and
%   ANG2_OLD (in degrees) into the grid of ANG1_NEW (in degrees) and
%   ANG2_NEW (in degrees) in a new coordinate system. The conversion of
%   angles are specified in the function handle ANGLECONVERT_FH, which
%   converts the angle in the new coordinates back to the original
%   coordinates.
%   
%   ANG1_OLD is a length-M vector, ANG2_OLD is a length-N vector, and
%   PAT_OLD is an NxM matrix. ANG1_NEW is a length-P vector, ANG2_NEW is a
%   length-Q vector, and PAT_NEW is a QxP matrix.
%
%   UVFlag is a logical value indicating whether the targeting coordinate
%   is in U/V space. If so, the resulting pattern is NaN for points outside
%   of the unit circle in U/V space.
%
%   The conversion first converts the angles in the new coordinates to the
%   old coordinate and then interpolate the values at those angles using
%   the existing grid in the old coordinates.
%
%   Example:
%   %   Convert a pattern from az/el space to u/v space.
%   
%   az = -90:90; el = -90:90; pat_azel = zeros(181,181);
%   u = -1:0.01:1; v = -1:0.01:1;
%   pat_uv = phased.internal.convertpat(@uv2azel,pat_azel,az,el,u,v,true);
%   surf(u,v,pat_uv);

%   Copyright 2011 The MathWorks, Inc.

%#codegen


[ang_1_gridArg,ang_2_gridArg] = meshgrid(ang_1,ang_2);
ang_1_grid = ang_1_gridArg(:);
ang_2_grid = ang_2_gridArg(:);

pat_new = nan(size(ang_1_grid));

if toUVFlag
    validIdx = find(hypot(ang_1_grid,ang_2_grid)<=1);
else
    validIdx = 1:numel(ang_1_grid);
end

ang_interp = fh_conv([ang_1_grid(validIdx).';ang_2_grid(validIdx).']).';
% Interpolate the pattern at the required old angle pair
pat_new(validIdx) = interp2(ang_1_old,ang_2_old,pat_old,...
    ang_interp(:,1),ang_interp(:,2),'nearest');
pat_newOut = reshape(pat_new,numel(ang_2),[]);

% [EOF]
