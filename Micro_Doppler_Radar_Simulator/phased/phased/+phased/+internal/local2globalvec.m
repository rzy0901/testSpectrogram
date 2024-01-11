function gVec = local2globalvec(lclVec,lclAxes)
%This function is for internal use only. It may be removed in the future.

%local2globalvec Convert vector from local to global coordinate system
%   GVEC = phased.internal.local2globalvec(LVEC,LAXES) converts the vector
%   represented in the local coordinate system, LVEC, to its corresponding
%   representation, GVEC, in the global coordinate system. LVEC is a 3-row
%   matrix whose columns represent vectors in (x,y,z) form. GVEC has the
%   same dimension as LVEC whose columns contain the corresponding vectors,
%   in the form of (x,y,z), in the global coordinate system. LAXES is a 3x3
%   matrix whose columns specify the x, y, and z axes of the local
%   coordinate system, respectively.
%
%   % Example:
%   %   Find the corresponding representation of vector [1;2;3] in the
%   %   global coordinate system. Assuming the local coordinate system is
%   %   given by [0 0 1;1 0 0;0 1 0].
%
%   gvec = phased.internal.local2globalvec([1;2;3],[0 0 1;1 0 0;0 1 0])

%   Copyright 2012 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% must be rectangular to rectangular
% vectors are specified in 3xN matrices, each column is a vector

gVec = lclAxes*lclVec;


% [EOF]
