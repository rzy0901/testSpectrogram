function scopext(ext)
%SCOPEXT  <short description>

%   Copyright 2010-2011 The MathWorks, Inc.

% There exist two ways to add src-dcv pairs to Visual object: 
% available at http://inside.mathworks.com/wiki/Scope_Source_Visual_API
% We can add more visual plugin's but not needed for this project

% Source
h = ext.add('Visuals', 'Array Polar Response', 'phaseddemoscope.ArrayResponseVisual2DPolar', ...
    'Array Response Pattern');
h.Visible = false;

% [EOF]