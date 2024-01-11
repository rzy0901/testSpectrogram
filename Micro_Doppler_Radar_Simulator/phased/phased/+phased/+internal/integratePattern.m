function intpat = integratePattern(pat,el,daz,del)
%This function is for internal use only. It may be removed in the future.

%integratePattern integrate radiation pattern
%   IPAT = phased.internal.integratePattern(PAT,EL,DeltaAZ,DeltaEl) returns
%   the integrated pattern, IPAT, based on the pattern, PAT, defined on an
%   azimuth and elevation grid. PAT is a pxqxr matrix whose pages are
%   patterns measured on p elevation angles and q azimuth angles. Elevation
%   angles are specified in a px1 vector, EL (in radians). DeltaAZ and
%   DeltaEl are azimuth and elevation angle steps (in radians) for each
%   measurement. IPAT is a 1xr vector whose entries are integrated patterns
%   corresponding to pages in PAT.

%   Copyright 2013 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

% del is px1 or scalar, in radian
% daz is 1xq or scalar, in radian
% el is px1, in radian
% pat is pxqxk
% intpat is 1xk

intpat = squeeze(sum(sum(bsxfun(@times,pat,(del.*cos(el))*daz),1),2)).';


% [EOF]
