function onDataSourceChanged(this)
%ONDATASOURCECHANGED <short description>

%   Copyright 2010-2011 The MathWorks, Inc.

% If we aren't rendered, we don't need to do anything here.
source = this.Application.DataSource;
if ~ishghandle(this.Axes) || isempty(source)
    return;
end

this.NewDataListener = addNewDataListener(this.Application, makeOnNewData(this));

% -------------------------------------------------------------------------
function cb = makeOnNewData(this)

cb = @(hSource) onNewData(this, hSource);

% [EOF]