function aboutphasedtbx
%ABOUTPHASEDTBX About the Phased Array System Toolbox.
%   ABOUTPHASEDTBX Displays the version number of the Phased Array System
%   Toolbox and the copyright notice in a modal dialog box.

%   Copyright 2012 The MathWorks, Inc.

load('aboutphased.mat','respdb');
tlbx = ver('phased');
str = sprintf([tlbx.Name ' ' tlbx.Version '\n',...
    getString(message('phased:phased:phasedCopyright',...
    datestr(tlbx.Date,10)))]);
msgbox(str,tlbx.Name,'custom',respdb,jet(64),'modal');

% [EOF]
