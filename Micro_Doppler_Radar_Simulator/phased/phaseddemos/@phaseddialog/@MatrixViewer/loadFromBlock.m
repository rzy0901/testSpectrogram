function loadFromBlock(h)
%loadFromBlock loads parameters from Simulink block.
%   LOADFROMBLOCK(h) loads Matrix Viewer specific parameters 
%   into object h from the associated SL block.

% Copyright 1995-2003 The MathWorks, Inc.

% Image Prop Tab
h.CMapStr = h.Block.CMapStr;
h.YMin = h.Block.YMin;
h.YMax = h.Block.YMax;
h.DataLimits = h.Block.DataLimits;
h.XData = h.Block.XData;
h.YData = h.Block.YData;    
h.AxisColorbar = strcmpi(h.Block.AxisColorbar,'on');
% Axis Prop tab
h.AxisOrigin = h.Block.AxisOrigin;
h.XLabel = h.Block.XLabel;
h.YLabel = h.Block.YLabel;
h.ZLabel = h.Block.ZLabel;
h.FigPos = h.Block.FigPos;
h.AxisZoom = strcmpi(h.Block.AxisZoom,'on');

% [EOF]
