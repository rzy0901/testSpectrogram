function propertyChanged(this, eventData)
%PROPERTYCHANGED React to changes in the contained property objects.

%   Copyright 2010-2011 The MathWorks, Inc.

if ~ischar(eventData)
    eventData = get(eventData.AffectedObject, 'Name');
end

switch eventData
    case 'xlimits'
        set(this.Axes, 'XLim', this.getPropValue('XLimits'));
    case 'polar'
        keyboard
end

% [EOF]
