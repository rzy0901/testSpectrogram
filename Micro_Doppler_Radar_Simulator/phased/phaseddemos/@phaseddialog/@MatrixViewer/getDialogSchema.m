function dlgStruct = getDialogSchema(this, name)
%getDialogSchema generates dialog schema for the MATRIX VIEWER block.
%   GETDIALOGSCHEMA(this,path) generates the DynamicDialog interface to the
%   MATRIX VIEWER block with associated phaseddialog.MatrixViewer object h and
%   full block path

% Copyright 1995-2010 The MathWorks, Inc.

% Here's the current assumption: this function is either being called to create
% a new dialog or to update an existing one.  If a dialog is being created, 
% then the ctor for this block must have just been called, and thus
% loadFromBlock has been called, and the block's members are up to date.  If 
% the dialog already exists, then the Mode=1 setting on the widgets has taken
% care of keeping the members up to date as well - thus, there's nothing we 
% need to do before figuring out which widgets we need.  If this turns out NOT
% to be the case, we may want some code like what's commented out below here:

%%%%%%%%%%%%%%   
% Parameters %
%%%%%%%%%%%%%%

%% no pop-up (combo-box) in first tab)

MaxWidgetInTab=7;
%% Image Property Tab
%% Widget-1:
CMapStr = dspGetLeafWidgetBase('edit','Colormap matrix:','CMapStr',this,'CMapStr');
CMapStr.Entries = set(this,'CMapStr')';
CMapStr.RowSpan = [1 1]; CMapStr.ColSpan = [1 1];
CMapStr.Visible = 1;
CMapStr.Tunable = 1;

%% Widget-2:
YMin = dspGetLeafWidgetBase('edit','Minimum input value:','YMin',this,'YMin');
YMin.Entries = set(this,'YMin')';
YMin.RowSpan = [2 2]; YMin.ColSpan = [1 1];
YMin.Visible = 1;
YMin.Tunable = 1;

%% Widget-3:
YMax = dspGetLeafWidgetBase('edit','Maximum input value:','YMax',this,'YMax');
YMax.Entries = set(this,'YMax')';
YMax.RowSpan = [3 3]; YMax.ColSpan = [1 1];
YMax.Visible = 1;
YMax.Tunable = 1;

%% Widget-4:
DataLimits = dspGetLeafWidgetBase('combobox','Data limits','DataLimits',this,'DataLimits');
DataLimits.Entries = {'Auto','User-defined'};
DataLimits.RowSpan = [4 4]; DataLimits.ColSpan = [1 1];
DataLimits.Visible = 1;
DataLimits.Tunable = 1;
DataLimits.DialogRefresh = 1;

%% Widget-5:
XData = dspGetLeafWidgetBase('edit','X-Data value:','XData',this,'XData');
XData.Entries = set(this,'XData')';
XData.RowSpan = [5 5]; XData.ColSpan = [1 1];
XData.Visible = 1;
XData.Tunable = 1;

%% Widget-6:
YData = dspGetLeafWidgetBase('edit','Y-Data value:','YData',this,'YData');
YData.Entries = set(this,'YData')';
YData.RowSpan = [6 6]; YData.ColSpan = [1 1];
YData.Visible = 1;
YData.Tunable = 1;

%% Widget-7:
AxisColorbar = dspGetLeafWidgetBase('checkbox','Display colorbar','AxisColorbar',this,'AxisColorbar');
AxisColorbar.RowSpan = [7 7]; AxisColorbar.ColSpan = [1 1];
AxisColorbar.Visible = 1;
AxisColorbar.Tunable = 1;

%%%% DYNAMICS CONTROL%%%%%%%%%%%%

if ~strcmp(this.DataLimits, 'User-defined')
    XData.Visible = 0; XData.Tunable = 0;
    YData.Visible = 0; YData.Tunable = 0;    
else
    XData.Visible = 1; XData.Tunable = 1;
    YData.Visible = 1; YData.Tunable = 1;    
end


%% Widget-Frame:
ImPropParameterPane = dspGetContainerWidgetBase('group','Parameters','ImPropParameterPane');
ImPropParameterPane.Items = {CMapStr,YMin,YMax,DataLimits,XData,YData,AxisColorbar};
ImPropParameterPane.Tag = 'ImPropParameterPane';
ImPropParameterPane.LayoutGrid = [MaxWidgetInTab 1];
ImPropParameterPane.RowStretch = [zeros(1,MaxWidgetInTab-1) 1];

%% Widget-Tab1:
ImagPropTab.Name = 'Image Properties';
ImagPropTab.Items = {ImPropParameterPane};

%% Axis Properties Tab
%% Widget-1:
AxisOrigin = dspGetLeafWidgetBase('combobox','Axis origin:','AxisOrigin',this,'AxisOrigin');
AxisOrigin.Entries = {'Upper left corner','Lower left corner'};
AxisOrigin.RowSpan = [1 1]; AxisOrigin.ColSpan = [1 1];
AxisOrigin.Visible = 1;
AxisOrigin.Tunable = 1;

%% Widget-2:
XLabel = dspGetLeafWidgetBase('edit','X-axis title:','XLabel',this,'XLabel');
XLabel.Entries = set(this,'XLabel')';
XLabel.RowSpan = [2 2]; XLabel.ColSpan = [1 1];
XLabel.Visible = 1;
XLabel.Tunable = 1;

%% Widget-3:
YLabel = dspGetLeafWidgetBase('edit','Y-axis title:','YLabel',this,'YLabel');
YLabel.Entries = set(this,'YLabel')';
YLabel.RowSpan = [3 3]; YLabel.ColSpan = [1 1];
YLabel.Visible = 1;
YLabel.Tunable = 1;

%% Widget-4:
ZLabel = dspGetLeafWidgetBase('edit','Colorbar title:','ZLabel',this,'ZLabel');
ZLabel.Entries = set(this,'ZLabel')';
ZLabel.RowSpan = [4 4]; ZLabel.ColSpan = [1 1];
ZLabel.Visible = 1;
ZLabel.Tunable = 1;

%% Widget-5:
FigPos = dspGetLeafWidgetBase('edit','Figure position, [x y width height]:','FigPos',this,'FigPos');
FigPos.Entries = set(this,'FigPos')';
FigPos.RowSpan = [5 5]; FigPos.ColSpan = [1 1];
FigPos.Visible = 1;
FigPos.Tunable = 1;

%% Widget-6:
AxisZoom = dspGetLeafWidgetBase('checkbox','Axis zoom','AxisZoom',this,'AxisZoom');
AxisZoom.RowSpan = [6 6]; AxisZoom.ColSpan = [1 1];
AxisZoom.Visible = 1;
AxisZoom.Tunable = 1;


%% Widget-Frame:
AxPropParameterPane = dspGetContainerWidgetBase('group','Parameters','AxPropParameterPane');
AxPropParameterPane.Items = {AxisOrigin,XLabel,YLabel,ZLabel,FigPos,AxisZoom};
AxPropParameterPane.Tag = 'AxPropParameterPane';
AxPropParameterPane.LayoutGrid = [MaxWidgetInTab 1];
AxPropParameterPane.RowStretch = [zeros(1,MaxWidgetInTab-1) 1];

%% Widget-Tab2:
AxisPropTab.Name = 'Axis Properties';
AxisPropTab.Items = {AxPropParameterPane};

%% Tab container
tabbedPane = dspGetContainerWidgetBase('tab','','tabPane');
tabbedPane.Tabs = {ImagPropTab,AxisPropTab};
tabbedPane.RowSpan = [2 2];
tabbedPane.ColSpan = [1 1];

dlgStruct = this.getBaseSchemaStruct(tabbedPane);
idx = findstr(this.Block.Name,sprintf('\n'));
blkName = this.Block.Name;
blkName(idx)=' ';
dlgStruct.DialogTitle = ['Sink Block Parameters: ' blkName];

% [EOF]
