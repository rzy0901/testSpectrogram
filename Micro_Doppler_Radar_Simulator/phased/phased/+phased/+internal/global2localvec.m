function lclVec = global2localvec(gVec,lclAxes)
%This function is for internal use only. It may be removed in the future.

%global2localvec Convert vector from global to local coordinate system
%   LVEC = phased.internal.global2localvec(GVEC,LAXES) converts the vector
%   represented in the global coordinate system, GVEC, to its corresponding
%   representation, LVEC, in the local coordinate system. GVEC is a 3-row
%   matrix whose columns represent vectors in (x,y,z) form. LVEC has the
%   same dimension as GVEC whose columns contain the corresponding vectors,
%   in the form of (x,y,z), in the local coordinate system. LAXES is a 3x3
%   matrix whose columns specify the x, y, and z axes of the local
%   coordinate system, respectively.
%
%   % Example:
%   %   Find the corresponding representation of vector [1;2;3] in the
%   %   local coordinate system. Assuming the local coordinate system is
%   %   given by [0 0 1;1 0 0;0 1 0].
%
%   lvec = phased.internal.global2localvec([1;2;3],[0 0 1;1 0 0;0 1 0])

%   Copyright 2012 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% must be rectangular to rectangular
% vectors are specified in 3xN matrices, each column is a vector

lclVec = lclAxes'*gVec;


% [EOF]
