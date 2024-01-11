function propsSchema = getPropsSchema(hCfg, hDlg)
%GETPROPSCHEMA Get the propSchema.

%   Copyright 2010-2011 The MathWorks, Inc.

% Get the 3 tabs.
main = groupToTab(getMainTab(hCfg), 'Main');
axis = groupToTab(getAxisTab(hCfg, hDlg), 'Axis Properties');

% Put the tabs together.
propsSchema.Type = 'tab';
propsSchema.Tabs = {main, axis};

% -------------------------------------------------------------------------
function main = getMainTab(hCfg)

[buffer_lbl, buffer] = uiscopes.getWidgetSchema(hCfg, 'DisplayBuffer', 'edit', 1, 1);

grid    = uiscopes.getWidgetSchema(hCfg, 'Grid',    'checkbox', 2, 1);
legend  = uiscopes.getWidgetSchema(hCfg, 'Legend',  'checkbox', 3, 1);
compact = uiscopes.getWidgetSchema(hCfg, 'Compact', 'checkbox', 4, 1);

main.Type = 'group';
main.Name = uiscopes.message('ParametersLabel');
main.LayoutGrid = [5 2];
main.RowStretch = [0 0 0 0 1];
main.ColStretch = [0 1];
main.Items = {buffer_lbl, buffer, grid, legend, compact};

% -------------------------------------------------------------------------
function axis = getAxisTab(hCfg, hDlg)

displaylimits = uiscopes.getWidgetSchema(hCfg, 'AutoDisplayLimits', 'checkbox', 1, 1);
displaylimits.DialogRefresh = true;

% MinimumXLim and MaximumXLim are visible when the AutoDisplayLimits are
% set to false.
visState = ~uiservices.getWidgetValue(displaylimits, hDlg);
[minxlim_lbl, minxlim] = extmgr.getWidgetSchema(hCfg, 'MinXLim', ...
    uiscopes.message('TimeMinXLimLabel'), 'edit', 2, 1);
minxlim_lbl.Visible = visState;
minxlim.Visible     = visState;

[maxxlim_lbl, maxxlim] = extmgr.getWidgetSchema(hCfg, 'MaxXLim', ...
    uiscopes.message('TimeMaxXLimLabel'), 'edit', 3, 1);
maxxlim_lbl.Visible = visState;
maxxlim.Visible     = visState;

[ylabel_lbl,  ylabel]  = uiscopes.getWidgetSchema(hCfg, 'YLabel',  'edit', 4, 1);
[minylim_lbl, minylim] = uiscopes.getWidgetSchema(hCfg, 'MinYLim', 'edit', 5, 1);
[maxylim_lbl, maxylim] = uiscopes.getWidgetSchema(hCfg, 'MaxYLim', 'edit', 6, 1);

axis.Type = 'group';
axis.Name = uiscopes.message('ParametersLabel');
axis.LayoutGrid = [7 2];
axis.RowStretch = [zeros(1, 6) 1];
axis.ColStretch = [0 1];
axis.Items = {displaylimits, ...
    minxlim_lbl, minxlim, ...
    maxxlim_lbl, maxxlim, ...
    ylabel_lbl, ylabel, ...
    minylim_lbl, minylim, ...
    maxylim_lbl, maxylim};

% -------------------------------------------------------------------------
function tab = groupToTab(group, tabName)

tab.Name  = tabName;
tab.Items = {group};

% [EOF]
