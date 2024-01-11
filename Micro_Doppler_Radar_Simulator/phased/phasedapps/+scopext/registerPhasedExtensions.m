function registerPhasedExtensions(ext)
% registerPhasedExtensions Register extensions for phased array scopes.

% Copyright 2014 The MathWorks, Inc.

r = ext.add('Visuals', 'PhasedMatrix', 'phased.scopes.MatrixVisual');
r.Visible = false;
r = ext.add('Visuals', 'IntensityVisual', 'phased.scopes.IntensityVisual', ...
            getString(message('phased:scopes:IntensityDescription')), ...
            getString(message('phased:scopes:IntensityLabel')));
r.Visible = false;
r = ext.add('Visuals', 'ScenarioVisual', 'phased.scopes.ScenarioVisual');
r.Visible = false;
% [EOF]
