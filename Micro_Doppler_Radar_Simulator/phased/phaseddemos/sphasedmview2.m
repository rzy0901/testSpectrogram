function sphasedmview2(varargin)
%SDSPMVIEW2 Matrix Viewer block Level-2 MATLAB S-function

% Copyright 1995-2014 The MathWorks, Inc.

% What's in the Figure userdata:
% ------------------------------
% Main scope figure handles:
%   fig_data.blockHandle        block handle
%   fig_data.hfig         handle to figure
%
%   fig_data.main.haxis         handle to axes
%   fig_data.main.himage        image handles
%   fig_data.main.axiszoom.on   P/V cell-array pairs to turn on zoom
%   fig_data.main.axiszoom.off     and turn it off when requested
%   fig_data.main.axiszoom.cbar    and move it to make room for colorbar
%   fig_data.main.colorbar.h    colorbar image
%   fig_data.main.colorbar.hax  colorbar axis handle
%
% Handles to menu items:
%   - appearing only in figure menu:
%       fig_data.menu.recpos     record position
%
%   - appearing in both figure and context menu:
%       fig_data.menu.top          top-level Axes and Lines in Figure
%       fig_data.menu.context      context menu
%       fig_data.menu.axiszoom     2x1, [fig;context] (checkable)
%       fig_data.menu.axiscolorbar 2x1, [fig;context]
%       fig_data.menu.autoscale
%
%
% What's in the Block userdata:
% -----------------------------
%   block_data.autoscaling  indicates autoscale computation in progress
%   block_data.hfig         handle to figure
%   block_data.haxis        handle to axes
%   block_data.himage       image handle
%   block_data.params       structure of cached block dialog parameters
%   block_data.hcolorbar    colorbar image handle
%   block_data.haxcolorbar  colorbar axis handle
%
% Block Parameters
% ----------------------------
% block.DialogPrm(1):  CMapStr: Cell array containing a string
% representation of a color map
% block.DialogPrm(2):  YMin:         Minimum input data value
% block.DialogPrm(3):  YMax:         Maximum input data value
% block.DialogPrm(4):  Data Limits: 'Auto','User-defined'
% block.DialogPrm(5):  XData: X-Data value
% block.DialogPrm(6):  YData: Y-Data value
% block.DialogPrm(7):  AxisColorbar: checkbox
% block.DialogPrm(8):  AxisOrigin: Lower left cornerUpper left corner
% block.DialogPrm(9):  XLabel: Cell array. x-axis label for main image
% block.DialogPrm(10):  YLabel: Cell array. y-axis label for main image
% block.DialogPrm(11):  ZLabel: Cell array. color bar scaling label
% block.DialogPrm(12):  FigPos:       figure position
% block.DialogPrm(13): AxisZoom:     checkbox

if nargin==1,
    mdlInitializeSizes(varargin{1});  % block
else
    %    varargin{4} => CopyFcn       = BlockCopy;
    %                   DeleteFcn     = BlockDelete;
    %                   NameChangeFcn = BlockNameChange; 
    feval(varargin{4:end});% GUI callback
end 

%% -----------------------------------------------------------
function mdlInitializeSizes(block)

% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 0;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Complexity   = 0;  % Real

% Register parameters
block.NumDialogPrms = 13; % coming from mask

% sampling mode
block.SampleTimes = [-1 0]; %Port-based sample time

% Specify if Accelerator should use TLC or call back into MATLAB file
block.SetAccelRunOnTLC(false);
block.SetSimViewingDevice(true);% no TLC required

% Specify that the block will work with context save/restore
block.SimStateCompliance = 'DefaultSimState';

% Reg methods
block.RegBlockMethod('CheckParameters',          @mdlCheckParameters);

block.RegBlockMethod('SetInputPortSamplingMode', @mdlSetInputPortFrameData);
block.RegBlockMethod('SetInputPortDimensions',   @mdlSetInputPortDimensions);
block.RegBlockMethod('SetInputPortDataType',     @mdlSetInputPortDataType);

block.RegBlockMethod('PostPropagationSetup',    @mdlPostPropSetup); %C-Mex: mdlSetWorkWidths

block.RegBlockMethod('Start',                   @mdlStart);
block.RegBlockMethod('ProcessParameters',       @mdlProcessParameters);
%block.RegBlockMethod('Update',                  @mdlUpdate); % registered in mdlPostPropSetup
block.RegBlockMethod('Terminate',               @mdlTerminate);

% see code in mdlStart

%% ---------------------------------------------------------------
function okflag = OK_TO_CHECK_PARAMS()

ss = get_param(bdroot,'simulationstatus');
okflag = strcmp(ss, 'initializing') || ...
    strcmp(ss, 'updating')     || ...
    strcmp(ss, 'running')      || ...
    strcmp(ss, 'paused');

%% ------------------------------------------------
function  mdlCheckParameters(block)

% Always call the check,
% because two things need it:
%  1 - error if it's ok_to_check_params, and
%  2 - mdlProcessParams needs to know we're okey-dokey
%
msg = CheckParams(block);
if OK_TO_CHECK_PARAMS() && ~isempty(msg)
    error(message('dsp:sdspmview2:invalidFigureParameter', msg));
end
%
% If all the vars are defined, push the changes through
% This gives us "live" parameters while the dialog is open
% and while the simulation is NOT running
%
if isempty(msg)
    mdlProcessParameters(block);
end

%% ------------------------------------------------
function mdlSetInputPortFrameData(block, idx, fd) 

block.InputPort(idx).SamplingMode = fd;

%% ------------------------------------------------
function mdlSetInputPortDataType(block, idx, dtid)

block.InputPort(idx).DatatypeID = dtid;


%% ------------------------------------------------
function mdlSetInputPortDimensions(block,idx,di)

block.InputPort(idx).Dimensions = di;

%% ------------------------------------------------
function mdlPostPropSetup(block)

% Register model update function
if (block.InputPort(1).DatatypeID < 13) %% built-in data type
     block.RegBlockMethod('Update', @mdlUpdate_builtinDT);
else %% fixed-point data type
     block.RegBlockMethod('Update', @mdlUpdate_fixedptDT);
end

%% -----------------------------------------------------------
function mdlStart(block)

% mdlStart Called before mdlUpdate in a simulation run.  
% Creates a new scope GUI if it does not exist, or restarts an 
% existing scope GUI.
%
% Updates block_data
codeGenMode = strcmp(get_param(bdroot(block.BlockHandle),'buildingrtwcode'),'on');
accelMode   = strcmp(get_param(bdroot(block.BlockHandle),'simulationmode'),'accelerator');
if (codeGenMode && ~accelMode)
    return;
end    

% To avoid contrast in color set image data to possible darkest color
blkh = block.BlockHandle;
block_data = get_param(blkh,'UserData');

if (block.InputPort(1).DatatypeID < 13) %% built-in 
     u = block.InputPort(1).Data; % save the class of input data type
else %% fixed-point
     u = block.InputPort(1).DataAsDouble; % save the class of input data type as double
end 

u(:) = block.DialogPrm(2).Data; % YMin;

% Construct new scope figure window, or bring up old one:
if isfield(block_data,'hfig') && ~isempty(block_data.hfig),
   hfig = block_data.hfig; % scope already exists
else
   hfig = [];              % scope was never run before
end

% Establish a valid scope GUI:
if ~isempty(hfig),
   % Found existing scope figure window:

   % Prepare to re-start with existing scope window:
   fig_data = restart_scope(block);

   % If things did not go well during restart, say the axis
   % was somehow deleted from the existing scope figure, etc.,
   % then hfig is left empty, and a new scope must be created.
   % Get hfig, then check if it is empty later:
   hfig = fig_data.hfig;
end

if isempty(hfig),
   % Initialize new figure window:
   % Create the scope GUI
   fig_data = create_scope(block);
end

% Get line handle:
himage = fig_data.main.himage;

% Retain the name of the figure window for use when the
% block's name changes. Name is retained in S-fcn block's
% user-data:
block_data.hfig      = fig_data.hfig;

block_data.haxis     = fig_data.main.haxis;
block_data.himage    = himage;
block_data.hcolorbar   = fig_data.main.colorbar.h;
block_data.haxcolorbar = fig_data.main.colorbar.hax;
block_data.autoscaling = []; % turn off any autoscaling, if in progress

% Set block's user data:
set_param(blkh, 'UserData', block_data);

% The following block callbacks are assumed to be set
% in the library block:
%
%   CopyFcn		      "sdspmview2([],[],[],'BlockCopy');"
%   DeleteFcn		  "sdspmview2([],[],[],'BlockDelete');"
%   NameChangeFcn     "sdspmview2([],[],[],'BlockNameChange');"

% Set menu checks according to the block_data:
SetMenuChecks(block);

% Setup scope axes:
setup_axes(block, u);  % one frame of data


%% -----------------------------------------------------------
function mdlProcessParameters(block)

blkh = block.BlockHandle;
block_data = get_param(blkh, 'UserData');

if isfield(block_data, 'hfig') && ~isempty(block_data.hfig)
  % GUI is open and a change has been made
  % Update menu checks:
  SetMenuChecks(block);

  % Handle figure position changes here:
  % Only update if the current block dialog FigPos differs from
  % the cached block_data FigPos.  The figure itself might be at
  % a different position, which we should not change UNLESS the
  % user actually made an explicit change in the mask (or, hit
  % the RecordFigPos menu item, which simply changes the mask).
  curpos = get(block_data.hfig,'Position');
  if ~isequal(curpos, block.DialogPrm(12).Data)
    set(block_data.hfig,'Position', block.DialogPrm(12).Data); % FigPos
  end
  
  % Get current image data
  u = get(block_data.himage, 'CData');
  
  setup_axes(block, u);
end
   

%% ------------------------------------------------------------
function mdlUpdate_builtinDT(block)

update_image(get_param(block.BlockHandle, 'UserData'), ...
             block.InputPort(1).Data);   %(block_data, u);

% end mdlUpdate_builtinDT

%% ------------------------------------------------------------
function mdlUpdate_fixedptDT(block)

update_image(get_param(block.BlockHandle, 'UserData'), ...
             block.InputPort(1).DataAsDouble);   %(block_data, u);

% end mdlUpdate_fixedptDT

% ---------------------------------------------------------------
function update_image(block_data, u)
% UPDATE_IMAGE Update the matrix displayed in the viewer

% u: one matrix of data
% Does not alter block_data

% If the user closed the figure window while the simulation
% was running, then hfig has been reset to empty.
%
% Allow the simulation to continue, but do not put up a new
% figure window (or error out!)
%
%
if isempty(block_data.hfig),
   return
end

% Update the image (old method):
%
% We cannot compare the new and old data, since
% Simulink has one (fixed) pointer for the data I/O,
% and "writes through" the pointer directly.  In other
% words, the HG image has its data automatically updated
% directly by Simulink.
%
% We could modify u before passing to HG in order to trigger
% a graphical update, but that would cause a deep copy of the
% entire matrix (a possibly large performance hit).  Plus,
% we couldn't add "eps" to a uint8/uint16 matrices.
% 
% xold = get(block_data.himage,'XData');
% set(block_data.himage, ...
%    'CData', u, ...
%    'XData',xold+eps, ...
%    'XData',xold);

% Update the image (new method):
set(block_data.himage, 'CData', u);

drawnow('update');

% Check if autoscaling is in progress:
if ~isempty(block_data.autoscaling),
   Autoscale(block_data.hfig);  % next frame of data
end


% ---------------------------------------------------------------
function SetMenuChecks(block)
% Called from mdlStart, mdlProcessParameters to set menu checks

blkh = block.BlockHandle;
block_data = get_param(blkh,'UserData');
fig_data   = get(block_data.hfig,'UserData');

% Update AxisZoom menu check:
%
axisZoom = block.DialogPrm(13).Data;
if axisZoom
    azoom_chkd = 'on';
    cBarVis='off';
else
    azoom_chkd = 'off';
    cBarVis='on';
end
set(fig_data.menu.axiszoom, 'Checked',azoom_chkd);

% Update Colorbar menu check:
%
showColorBar = block.DialogPrm(7).Data;
if showColorBar
    cBarCheck = 'on';
else
    cBarCheck = 'off';
end
set(fig_data.menu.axiscolorbar, ...
   'Checked', cBarCheck, ...
   'Enable',  cBarVis);

% ---------------------------------------------------------------
% Block Parameters
% ----------------------------
% block.DialogPrm(1):  CMapStr: Nx3 colormap matrix (string)
% block.DialogPrm(2):  YMin:         Minimum input data value
% block.DialogPrm(3):  YMax:         Maximum input data value
% block.DialogPrm(4):  Data Limits: 'Auto','User-defined'
% block.DialogPrm(5):  XData: X-Data value
% block.DialogPrm(6):  YData: Y-Data value
% block.DialogPrm(7):  AxisColorbar: checkbox
% block.DialogPrm(8):  AxisOrigin: Lower left cornerUpper left corner
% block.DialogPrm(9):  XLabel: x-axis label for main image
% block.DialogPrm(10):  YLabel: y-axis label for main image
% block.DialogPrm(11):  ZLabel: color bar scaling label
% block.DialogPrm(12):  FigPos:       figure position
% block.DialogPrm(13): AxisZoom:     checkbox
function msg = CheckParams(block)
% -----------------
% (skip checkboxes and popups)

msg = '';

% Check XData:
% -------------

if ~strcmp(block.DialogPrm(4).Data,'Auto')
    x = block.DialogPrm(5).Data;
    if ~isa(x,'double') || issparse(x) || ~isreal(x) || ...
            ~isequal(size(x),[1,2]) || (x(2) <= x(1)) 
        msg = 'X-Data must be a 2 element vector in the form of [xmin xmax]';
        return
    end
    x = block.DialogPrm(6).Data;
    if ~isa(x,'double') || issparse(x) || ~isreal(x) || ...
            ~isequal(size(x),[1,2]) || (x(2) <= x(1)) 
        msg = 'Y-Data must be a 2 element vector in the form of [xmin xmax]';
        return
    end
end

% Check XLabel:
% -------------
if ~ischar(block.DialogPrm(9).Data{1}),
   msg = 'X-axis label must be a string.';
   return
end

% Check YLabel:
% -------------
if ~ischar(block.DialogPrm(10).Data{1}),
   msg = 'Y-axis label must be a string.';
   return
end

% Check ZLabel:
% -------------
if ~ischar(block.DialogPrm(11).Data{1}),
   msg = 'Z-axis label must be a string.';
   return
end

% Check YMin:
% -----------
x = block.DialogPrm(2).Data;
Nx = numel(x);
if ~isa(x,'double') || issparse(x) || ~isreal(x) || (Nx ~= 1)
   msg = 'Y-minimum must be a real-valued scalar.';
   return
end
ymin = x;

% Check YMax:
% -----------
x = block.DialogPrm(3).Data;
Nx = numel(x);
if ~isa(x,'double') || issparse(x) || ~isreal(x) || (Nx ~= 1)
    msg = 'Y-maximum must be a real-valued scalar.';
    return
end
if ~isempty(x) && ~isempty(ymin) && (x <= ymin),
   msg = 'Maximum Y-axis limit must be greater than Minimum Y-axis limit.';
   return
end

% Check FigPos:
% -------------
x = block.DialogPrm(12).Data;
if ~isa(x,'double') || issparse(x) || ~isreal(x) || ...
      size(x,1)~= 1 || size(x,2)~=4,
   msg = 'Figure position must be a real-valued 1x4 vector.';
   return
end

% Check CMapStr:
% --------------
x = block.DialogPrm(1).Data{1};
if ~isa(x,'double') || (~ismatrix(x)) || (size(x,2)~=3),
   msg = 'Colormap must be an Nx3 matrix.';
   return
end
if min(min(x)) < 0 || max(max(x)) > 1.0,
   msg = 'Colormaps must contain values between 0.0 and 1.0, inclusive.';
   return
end

% ---------------------------------------------------------------
function setup_axes(block, u)
% Setup viewer x- and y-axes

% Does not alter block_data
% u = input data (one matrix)

blkh = block.BlockHandle;
block_data = get_param(blkh,'UserData');
hfig    = block_data.hfig;
hax     = block_data.haxis;
himage  = block_data.himage;
[nrows, ncols] = size(u);
haxclr  = block_data.haxcolorbar;

% Refresh display: (in case of stray markings):
% ----------------
FigRefresh(block_data.hfig);

% Setup X-axis label:
% -------------------
% Don't modify user-defined domain:
xLabel = block.DialogPrm(9).Data{1};
if ~ischar(xLabel), xLabel = 'X-Axis'; end
hxLabel = get(hax, 'XLabel');
set(hxLabel, 'String', xLabel);

% Setup Y-axis label:
% -------------------
yLabel = block.DialogPrm(10).Data{1};
if ~ischar(yLabel), yLabel='Y-Axis'; end
hyLabel = get(hax,'YLabel');
set(hyLabel, 'String', yLabel);

% Setup Colorbar label:
% ---------------------
cLabel = block.DialogPrm(11).Data{1};
if ~ischar(yLabel), cLabel='Z-Axis'; end
hyLabel = get(haxclr,'YLabel');
set(hyLabel, 'String', cLabel);

% Setup image data:
% -----------------
% NOTE: update_image() does NOT alter the xdata or ydata,
%       so it could be set up once here:

cmap = block.DialogPrm(1).Data{1};
set(hfig,'colormap', cmap);

if strcmp(block.DialogPrm(4).Data,'Auto')
    xdata = [1 ncols];
    ydata = [1 nrows];
else
    xdata =  block.DialogPrm(5).Data;
    ydata =  block.DialogPrm(6).Data;    
end
xdataWidth = (xdata(2)-xdata(1))/(ncols-1);
xLim = [xdata(1)-(xdataWidth/2) xdata(2)+(xdataWidth/2)];
ydataWidth = (ydata(2)-ydata(1))/(nrows-1);
yLim = [ydata(1)-(ydataWidth/2) ydata(2)+(ydataWidth/2)];

set(himage,'CData',u,'XData',xdata,'YData',ydata);

% Adjust axis for origin:
origin = block.DialogPrm(8).Data;
isXY=strncmp(origin,'Lower',5);
if isXY, ydir='normal';
else ydir='reverse';
end
set(hax,'YDir',ydir);

% Setup colormap scaling
% if none of ymin/ymax is empty:
% - Figure has the colormap
% - Axis has the clim/climmode
% - image has the image data
   % block_data.params.YMin
   % block_data.params.YMax
if ~isempty(block.DialogPrm(2).Data) && ...
     ~isempty(block.DialogPrm(3).Data),
  set(hax, ...
     'clim',[block.DialogPrm(2).Data block.DialogPrm(3).Data], ...
     'climmode','manual');
  set(himage,'CDataMapping','Scaled');
else
  set(hax, 'climmode','auto');
  set(himage,'CDataMapping','Direct');
end

% Update colorbar image:
% ----------------------
useColorbar = block.DialogPrm(7).Data;
axiszoom = block.DialogPrm(13).Data;

cbarANDnozoom = useColorbar & (~axiszoom);

% Colorbar (and its axis) visibility:
if cbarANDnozoom, cbarVis='on'; else cbarVis='off'; end

% Modify colorbar vertical axis limits to match
% the current scaling method in force:
N = size(cmap,1);
ylim = [1 N];
if ~isempty(block.DialogPrm(2).Data) && ...
   ~isempty(block.DialogPrm(3).Data),
  ylim = [block.DialogPrm(2).Data block.DialogPrm(3).Data];
end

set(block_data.hcolorbar, ...
   'CDataMapping','Direct', ...
   'CData',(1:N)', ...
   'YData',ylim, ...
   'vis', cbarVis);

set(haxclr, ...
   'ylim',ylim, ...
   'vis', cbarVis);


% Perform AxisZoom:
% -----------------
% Put axis into correct zoom state:
fig_data = get(hfig,'UserData');
if ~axiszoom,
   % Turn off AxisZoom:
      
   % - reset axis position
   if useColorbar,
      set(hax, fig_data.main.axiszoom.cbar{:});
   else
      set(hax, fig_data.main.axiszoom.off{:});
   end
   
else
   % Turn on AxisZoom:
   
   % - turn off top-level menus
   set(fig_data.menu.top,'vis','off');
   set(hfig,'menu','none');
   
   % - set axis position
   set(hax, fig_data.main.axiszoom.on{:});
end

set(hax, ...
    'xlim',xLim, ...
    'ylim',yLim, ...
    'zlimmode','manual');
    

% Update display with the actual line values:
update_image(block_data, u);  % one frame of data


% ---------------------------------------------------------------
function fig_data = create_scope(block)
% CREATE_SCOPE Create new scope GUI

blkh = block.BlockHandle;
% Initialize empty settings:
fig_data.main  = [];  % until we move things here
fig_data.menu  = [];

iotype = get_param(blkh,'iotype');
if strcmp(iotype,'viewer')
  fig_name = viewertitle(blkh,false);
else
  parent_path   = get_param(blkh,'parent');
  block_name_with_path = [parent_path '/' get_param(blkh,'Name')];
  fig_name = block_name_with_path;
end

hfig = figure('numbertitle', 'off', ...
   'name',         fig_name, ...
   'menubar','none', ...
   'toolbar','none', ...
   'position',     block.DialogPrm(12).Data, ...
   'nextplot',     'add', ...
   'integerhandle','off', ...
   'doublebuffer', 'off', ...
   'DeleteFcn',    @FigDelete, ...
   'HandleVisibility','callback');

% Use double-buffer when image EraseMode is set to 'none',
% AND when we're doing a forced double-update.  This would
% reduce the work to one on-screen blit instead of two...

% Axis for the image:
hax = axes('Parent',hfig, ...
   'SortMethod','childorder', ...
   'Box','on', 'ticklength',[0 0]);

% Axis for the colorbar:
haxcbar = axes('Parent',hfig, ...
   'xtick', [], ...  % turn off x ticks
   'xlim', [-0.5 1], ...
   'yaxislocation','right', ...
   'Box','on', 'ticklength',[0 0]);
   
% Set up image:
himage = image('parent',hax,'cdata',[]);

% Set up colorbar image:
hcolorbar = image('parent',haxcbar,'cdata',[],'xdata',[0 1]);

% Establish settings for all structure fields:
fig_data.blockHandle = blkh;
fig_data.hfig  = hfig;

% Store major settings:
fig_data.main.haxis   = hax;
fig_data.main.himage  = himage;

% Store settings for axis zoom:
% Cell-array contains {params, values},
% where params itself is a cell-array of Units and Position
% and values is a cell-array of corresponding values.
p = {'Units','Position'};
fig_data.main.axiszoom.off  = {p, {'Normalized',[.13 .145 .8 .8]}};
fig_data.main.axiszoom.on   = {p, {'Normalized',[ 0   0    1    1]}};
fig_data.main.axiszoom.cbar = {p, {'Normalized',[.13 .145 .645 .8]}};

% Copy colorbar data:
%
fig_data.main.colorbar.pos  = {p, {'Normalized',[.8 .145 .050 .8]}};
fig_data.main.colorbar.h    = hcolorbar;
fig_data.main.colorbar.hax  = haxcbar;

set(hax,    'Position',fig_data.main.axiszoom.cbar{2}{2});
set(haxcbar,'Position',fig_data.main.colorbar.pos{2}{2});

% Define axis menu labels:
%
if ispc
   labels = {'&Axes', '&Refresh', ...
         '&Autoscale', 'Axis &Zoom', '&Colorbar', ...
         'Save &Position'};
else
   labels = {'Axes', 'Refresh', ...
         'Autoscale', 'Axis Zoom', 'Colorbar', ...
         'Save Position'};
end
%
% Create figure AXES context menu
mAxes = uicontextmenu('parent',hfig);  % top-level Axes menu in figure
%
% submenu items:
fig_data.menu.refresh = uimenu(mAxes, 'label',labels{2}, ...
   'callback',@FigRefresh);
fig_data.menu.autoscale = uimenu(mAxes, 'label',labels{3}, ...
   'separator','on',...
   'callback', @Autoscale);
% - Create Axis Zoom item
fig_data.menu.axiszoom = uimenu(mAxes, ...
   'Label', labels{4}, ...
   'Callback', @AxisZoom);
% - Create Colorbar item
fig_data.menu.axiscolorbar = uimenu(mAxes, ...
   'Label', labels{5}, ...
   'Callback', @AxisColorbar);
% - Create Record Position item
fig_data.menu.recpos = uimenu(mAxes, 'label',labels{6}, ...
   'callback', @SaveFigPos, ...
   'separator','on');

% Store all top-level menu items in one vector
fig_data.menu.top = mAxes;

% Record figure data:
set(hfig, 'UserData', fig_data);

% Assign context menu to the axis, lines, and grid:
set([fig_data.main.haxis fig_data.main.himage], ...
   'UIContextMenu', mAxes);


% ---------------------------------------------------------------
function fig_data = restart_scope(block)
% RESTART_SCOPE Restart with existing scope window

% We want to confirm to a reasonable probability that
% the existing scope window is valid and can be restarted.

% The caller already verified that hfig is non-empty
blkh = block.BlockHandle;
block_data = get_param(blkh,'UserData');
hfig = block_data.hfig;

% We don't know if the handle points to a valid window:
if isempty(hfig) || ~ishandle(hfig),
   block_data.hfig = [];  % reset it back
   set_param(blkh,'UserData',block_data);
   fig_data = [];
   return;
end

% Something could fail during restart if the figure data was
% altered between runs ... for example, by command-line interaction.
% If errors occur, abandon the restart attempt:
try
   fig_data = get(hfig,'UserData');
   
   % In case memory (persistence) was on:
   FigRefresh(hfig);
   
   figure(hfig); % bring window forward
   
catch
   % Something failed - reset hfig to indicate error during restart:
   fig_data.hfig=[];
   block_data.hfig=[];
end

% Update data structures:
set(hfig, 'UserData',fig_data);
set_param(blkh, 'UserData',block_data);


% ---------------------------------------------------------------
function BlockNameChange %#ok - this is a block callback set in library file
% In response to the name change, we must do the following:
%
% (1) find the old figure window, only if the block had a GUI 
%     associated with it.
%     NOTE: Current block is parent of the S-function block
blkh = gcbh;
block_data = get_param(blkh, 'UserData');

% System might never have been run since loading.
% Therefore, block_data might be empty:
if isfield(block_data, 'hfig') && ~isempty(block_data.hfig),
   %isstruct(block_data),
   % (2) change name of figure window (cosmetic)
   hfig = block_data.hfig;
   
   iotype = get_param(blkh,'iotype');
   if strcmp(iotype,'viewer')
     fig_name = viewertitle(blkh,false);
   else
     parent_path   = get_param(blkh,'parent');
     block_name_with_path = [parent_path '/' get_param(blkh,'Name')];
     fig_name = block_name_with_path;
   end

   set(hfig,'name',fig_name);
   
   % (3) update figure's userdata so that the new blockname
   %     can be used if the figure gets deleted
   fig_data = get(hfig,'UserData');
   fig_data.blockHandle = blkh;
   set(hfig,'UserData',fig_data);
end



% ---------------------------------------------------------------
function CloseFigure(blkh)
% Manual (programmatic) closing of the figure window

%
block_data = get_param(blkh,'UserData');
block_data.hfig = [];
set_param(blkh, 'UserData',block_data);


% ---------------------------------------------------------------
function FigDelete(hcoNotUsed, eventStructNotUsed) %#ok
% Callback from figure window
% Called when the figure is closed or deleted

hfig = gcbf;
fig_data = get(hfig,'UserData');
if hfig ~= fig_data.hfig,
   error(message('dsp:sdspmview2:invalidHandle1'));
end

% Close the figure window
CloseFigure(fig_data.blockHandle);


% ---------------------------------------------------------------
function BlockDelete %#ok - this is a block callback set in library file
% Block is being deleted from the model

% clear out figure's close function
% delete figure manually
blkh = gcbh;
block_data = get_param(blkh,'UserData');
if isfield(block_data, 'hfig') && ~isempty(block_data.hfig),
   set(block_data.hfig, 'DeleteFcn','');
   delete(block_data.hfig);
   block_data.hfig = [];
   set_param(blkh,'UserData',block_data);
end


% ---------------------------------------------------------------
function BlockCopy %#ok - this is a block callback set in library file
% Block is being copied from the model

% clear out stored figure handle
blkh = gcbh;
block_data = get_param(blkh,'UserData');
if isstruct(block_data),
   block_data.hfig = [];
   set_param(blkh,'UserData',block_data);
end


% ---------------------------------------------------------------
function SaveFigPos(hcoNotUsed, eventStructNotUsed) %#ok
% Record the current position of the figure into the block's mask

% Get the block's name:
hfig = gcbf;
fig_data = get(hfig,'UserData');
if hfig ~= fig_data.hfig,
   error(message('dsp:sdspmview2:invalidHandle2'));
end

% Record the figure position, as a string, into the appropriate mask dialog:
FigPos = get(hfig,'Position');             % Get the fig position in pixels
blkh = fig_data.blockHandle;
set_param(blkh, 'FigPos', mat2str(FigPos)); % Record new position


% ---------------------------------------------------------------
function AxisColorbar(hcoNotUsed, eventStructNotUsed) %#ok
%function AxisColorbar_CallBack(hco, eventStruct)
% Toggle axis colorbar on and off
%
% opt is a string option and may be one of the following:
%     'toggle', 'on', 'off'
% If not passed, default is 'toggle'.
%
% hfig is the figure handle
% if missing, it is set to gcbf

hfig = gcbf;

fig_data = get(hfig, 'UserData');
blkh      = fig_data.blockHandle;
haxcbar  = fig_data.menu.axiscolorbar;

if strcmp(get(haxcbar,'Checked'),'on'),
    opt='off';
else
    opt='on';
end

% Update menu check:
set(haxcbar,'Checked',opt);

% Update block dialog setting, so param is recorded in model:
% This will indirectly update the param structure, via the
% mask dialog callbacks.
set_param(blkh, 'AxisColorbar', opt);

% ---------------------------------------------------------------
function AxisZoom(hcoNotUsed, eventStructNotUsed) %#ok
% Toggle display of zoomed-in axes
%
% opt is a string option and may be one of the following:
%     'toggle', 'on', 'off'
% If not passed, default is 'toggle'.
%
% hfig is the figure handle
% if missing, it is set to gcbf

hfig=gcbf;

fig_data = get(hfig, 'UserData');
blkh      = fig_data.blockHandle;
haxzoom  = fig_data.menu.axiszoom;

if strcmp(get(haxzoom,'Checked'),'on'),
    opt='off';
else
    opt='on';
end

% Update menu check:
set(haxzoom,'Checked',opt);

% Update block dialog setting, so param is recorded in model:
% This will indirectly update the param structure, via the
% mask dialog callbacks.
set_param(blkh, 'AxisZoom', opt);


% ---------------------------------------------------------------
function Autoscale(hco, eventStructNotUsed) %#ok
% AUTOSCALE Compute min/max y-limits for several input frames

if nargin == 2
    hfig = gcbf; % Called via callback from context menu
    fig_data = get(hfig,'UserData');
    blkh = fig_data.blockHandle;
    block_data = get_param(blkh,'UserData');

    % If an autoscale operation is currently in progress,
    % cancel it and stop:
    if ~isempty(block_data.autoscaling),
        CancelAutoscale(blkh);
        return; % EARLY RETURN
    end

    % If simulation stopped, do a one-shot autoscale and leave:
    v = get_sysparam(blkh,'simulationstatus');
    if ~strcmp(v,'running'),
        % Simulation is stopped or paused - perform simple one-shot autoscaling:
        oneshot_autoscale(hfig);
        return; % EARLY RETURN
    end


    % Begin countdown
    % This is the number of sequential frames which will be examined
    % in order to determine the min/max y-limits
    count=10;

    % Preset min and max
    ymin=+inf;
    ymax=-inf;

    % Put up an autoscale indicator
    str = ['Autoscale: ' mat2str(count)];
    htext = text('units','norm','pos',[0.5 0.5], ...
        'Color','white', ...
        'horiz','center', 'string',str);

    block_data.autoscaling = [count ymin ymax double(htext)];
    set_param(blkh, 'UserData', block_data);

else
    hfig = hco; % Called via mdlUpdate
    fig_data = get(hfig,'UserData');
    blkh = fig_data.blockHandle;
    block_data = get_param(blkh,'UserData');
    % Continue processing next frame of inputs
    % to determine autoscale limits

    count = block_data.autoscaling(1);
    ymin  = block_data.autoscaling(2);
    ymax  = block_data.autoscaling(3);
    htext = block_data.autoscaling(4);

    u = get(block_data.himage,'CData');
    if count>0,
        % Continue tracking min and max:

        count=count-1;
        ymin=double(min(ymin,min(u(:))));
        ymax=double(max(ymax,max(u(:))));

        % Update user feedback:
        set(htext,'string',['Autoscale: ' mat2str(count)]);

        block_data.autoscaling = [count ymin ymax double(htext)];
        set_param(blkh,'UserData', block_data);

    else
        % Finished computing autoscale limits

        % Remove autoscale indicator:
        delete(htext);
        htext=[];   %#ok  % reset so that terminate call deletes an empty handle

        % Turn off autoscale flag
        block_data.autoscaling = [];
        set_param(blkh, 'UserData', block_data);

        % Protect against horizontal lines:
        if (ymax==ymin),
            ymin=floor(ymin-.5);
            ymax=ceil(ymax+.5);
        end

        % Indirectly set these via the DialogApply callback:
        set_param(blkh, 'YMin', num2str(ymin));
        set_param(blkh, 'YMax', num2str(ymax));
    end
end

% ---------------------------------------------------------------
function CancelAutoscale(blkh)

% Cancel any pending autoscale operation

block_data = get_param(blkh,'UserData');

% No autoscale operation in progress:
if ~isfield(block_data,'autoscaling') || isempty(block_data.autoscaling),
   return;
end

htext = block_data.autoscaling(4);
delete(htext);
block_data.autoscaling=[];
set_param(blkh,'UserData', block_data);


% ---------------------------------------------------------------
function oneshot_autoscale(hfig)
% ONESHOT_AUTOSCALE Used when simulation is stopped
%   Cannot use multi-input autoscale, since the simulation is no longer
%   running.  Instead, we compute a one-time ymin/ymax computation, and
%   apply it to the static scope result.

fig_data = get(hfig, 'UserData');
blkh = fig_data.blockHandle;

% Get data for each line, and find min/max:
himage = fig_data.main.himage;
y = get(himage,'CData');
ymin = min(y(:));
ymax = max(y(:));

% Protect against horizontal lines:
if (ymax==ymin),
   ymin=floor(ymin-.5);
   ymax=ceil(ymax+.5);
end
% Update dialog
set_param(blkh, 'YMin', num2str(ymin));
set_param(blkh, 'YMax', num2str(ymax));

% ---------------------------------------------------------------
function mdlTerminate(block)
% TERMINATE Clean up any remaining items

blkh = block.BlockHandle;

% Cancel any pending autoscale operation:
CancelAutoscale(blkh);

% ---------------------------------------------------------------
function FigRefresh(hco, eventStructNotUsed) %#ok
% Refresh display while memory turned on

if nargin > 1
% Called via callback
    hfig=gcbf;
else % handle to figure in first (and only) argument
% Called from another local function
    hfig = hco;
end

if ~isempty(hfig),
   refresh(hfig);
end

function v = get_sysparam(sys, p)
% GET_SYSPARAM Get a system-level parameter.
%  This function is error-protected against calls on
%  subsystems or blocks, which do not have system-level
%  parameter.


% Copyright 1995-2011 The MathWorks, Inc.

if strmatch(lower(p), ...
      lower(fieldnames(get_param(sys,'objectparameters'))), 'exact'),
   % Simulink does not do param-completion, so we use 'exact' above
   
   % Object has parameter - get the corresponding value:
   v = get_param(sys,p);
else
   % No system parameter - check parent:
   parent = get_param(sys,'Parent');
   if ~isempty(parent),
      % recurse
      v = get_sysparam(get_param(sys,'parent'),p);
   else
      % param not found
      error(message('dsp:get_sysparam:invalidFcnInput', p));
   end
end
% ------------------------------------------------------------
% [EOF] sdspmview2.m

% LocalWords:  userdata hfig haxis himage axiszoom cbar hax recpos checkable
% LocalWords:  axiscolorbar autoscaling hcolorbar haxcolorbar CMap YMin YMax
% LocalWords:  XLabel YLabel ZLabel simulationstatus okey dokey buildingrtwcode
% LocalWords:  simulationmode CData DT fixedpt xold XData Nx popups Colormaps
% LocalWords:  xdata ydata YDir ymin ymax clim climmode vis zlimmode iotype
% LocalWords:  numbertitle nextplot integerhandle doublebuffer blit ticklength
% LocalWords:  xtick yaxislocation submenu blockname hco horiz ONESHOT
