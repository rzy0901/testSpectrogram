function blkStruct = slblocks
%SLBLOCKS Defines the Simulink library block representation
%   for the Phased Array System Toolbox.

%   Copyright 2013 The MathWorks, Inc.

blkStruct.Name    = sprintf('Phased Array\nSystem\nToolbox');
blkStruct.OpenFcn = 'phasedlib';
blkStruct.MaskInitialization = '';

% Define the library list for the Simulink Library browser.
% Return the name of the library model and the name for it
Browser(1).Library = 'phasedlib';
Browser(1).Name    = 'Phased Array System Toolbox';
Browser(1).IsFlat  = 0; % Is this library "flat" (i.e. no subsystems)?

blkStruct.Browser = Browser;

% Define information for model updater
% blkStruct.ModelUpdaterMethods.fhSeparatedChecks = @spblksUpdateModel;
% blkStruct.ModelUpdaterMethods.fhDetermineBrokenLinks = @spblksBrokenLinksMapping;

% End of slblocks.m


