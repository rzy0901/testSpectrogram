function varargout = validate(hDlg)
%VALIDATE Returns true if this object is valid

%   Copyright 2010-2011 The MathWorks, Inc.

% Ask the super class to check the properties it adds.
[b, exception] = uiscopes.AbstractLineVisual.validate(hDlg);
%b = true;
%exception = MException.empty; % if b is true, the exception should be an empty exception object 
PhaseString = hDlg.getWidgetValue([hDlg.getSource.Register.Name 'DisplayPhase']);

% The source data should be non-negative value, check that the phase range
% is a valid variable
[~, errid, errmsg] = uiservices.evaluate(PhaseString);
if ~isempty(errid)
    b=false;
    exception = MException(errid,errmsg);
end

if nargout
    varargout = {b, exception};
elseif ~b
    throw(exception);
end

% [EOF]

