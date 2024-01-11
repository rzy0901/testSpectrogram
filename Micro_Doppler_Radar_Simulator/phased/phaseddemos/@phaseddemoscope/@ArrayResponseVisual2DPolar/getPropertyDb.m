function propertyDb = getPropertyDb
%GETPROPERTYDB Get the propertyDb.

%   Copyright 2010-2011 The MathWorks, Inc.

propertyDb =uiscopes.AbstractLineVisual.getPropertyDb;
propertyDb.add('XLimits', 'mxArray', [-180 180]);
propertyDb.add('Polar', 'bool', true);

% [EOF]
