function sonarEquationCalculator
% sonarEquationCalculator Estimates target range, transmission loss,
% source level and SNR of a sonar system
%   sonarEquationCalculator launches the sonar equation calculator app. The
%   app allows the user to calculate target range, transmission loss,
%   SourceLevel, or SNR level using the sonar equations. The app supports
%   calculations for both active and passive modes.

% Copyright 2017 The MathWorks, Inc.

% Construct UI components
GUI_WIDTH = 380;
GUI_HEIGHT = 425;
EDIT_ALIGN = 'Left';
RSZ_HEIGHT_TL = 120;   % Resize for Transmission loss calculation
RSZ_HEIGHT = 150;      % Resize height for detection settings
DBLRTARROW = getString(message('phased:apps:sonareqapp:DoubleRightArrow'));
DBLLTARROW = getString(message('phased:apps:sonareqapp:DoubleLeftArrow'));

% Units factors used to convert result values to units set on app
Units.ResultRange = 1; % Convert to m by default
Units.RangeDepth = 1;  % Convert to m by default
Units.RangeFreq = 1e3; % Default units set to kHz 
Units.Depth = 1;       % Convert to m by default
Units.Range = 1;       % Convert to m by default
Units.Freq = 1e3;      % Default units set to kHz 

hGui = figure( ...
    'Name', getString(message('phased:apps:sonareqapp:SonarEqApp')), ...
    'Visible', 'off', ...
    'Toolbar', 'none', ...
    'Menubar', 'none', ...
    'NumberTitle', 'off', ...
    'DockControl', 'off', ...
    'Resize', 'off', ...
    'IntegerHandle', 'off', ...
    'Units', 'pixels', ...
    'Position', [300, 300, GUI_WIDTH, GUI_HEIGHT], ...
    'Tag', 'sonarEqApp');

% Color reference for UI controls
BKGND_COLOR = get(hGui,'DefaultUIControlBackgroundColor');
EDIT_COLOR = 'white';

% Set main gui background color to default
set(hGui, 'Color', BKGND_COLOR);

% Menu structure
hFileMenu = uimenu(hGui,'Label', getString(message('phased:apps:sonareqapp:File')), ...
    'Tag', 'FileMenu');
uimenu(hFileMenu, 'Label', getString(message('phased:apps:sonareqapp:GenReport')), ...
    'Tag', 'GenReport', 'Callback', @gen_report_Callback);
uimenu(hFileMenu, 'Label', getString(message('phased:apps:sonareqapp:GenMatlabCode')), ...
    'Tag', 'GenMCode', 'Callback', @gen_mcode_Callback);
uimenu(hFileMenu, 'Label', getString(message('phased:apps:sonareqapp:Close')), ...
    'Separator', 'on', 'Tag', 'CloseMenu', 'Callback', @close_Callback);
hHelpMenu = uimenu(hGui, 'Label', getString(message('phased:apps:sonareqapp:Help')), ...
    'Tag', 'HelpMenu');
uimenu(hHelpMenu, 'Label', getString(message('phased:apps:sonareqapp:SonarEqAppHelp')), ...
    'Tag', 'RadEqAppHelp', 'Callback', @sonareqapp_help_Callback);
uimenu(hHelpMenu, 'Label',getString(message('phased:apps:sonareqapp:PhASTHelp')) , ...
    'Tag', 'PASTHelp', 'Callback', @phased_help_Callback);
uimenu(hHelpMenu, 'Label', getString(message('phased:apps:sonareqapp:SonarEqAppAbout')), ...
    'Separator', 'on', 'Tag', 'AboutMenu', 'Callback', @about_Callback);

%------------------------------------------------------------------------%
% Sonar Equation Setting UI components
%------------------------------------------------------------------------%
hCalcpanel = uipanel(...
    'Parent', hGui, ...
    'FontSize', 10, ...
    'BorderType', 'none', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'CalcPanel');

hCalctxt = uicontrol(...
    'Parent', hCalcpanel, ...
    'Style', 'text', ...
    'String', 'Calculation:', ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'ParmCalculateLabel', ...
    'BackgroundColor', BKGND_COLOR);

calc_parms = {getString(message('phased:apps:sonareqapp:Range')), ...
    getString(message('phased:apps:sonareqapp:TransmissionLoss')),...
    getString(message('phased:apps:sonareqapp:SourceLevel')), ...
    getString(message('phased:apps:sonareqapp:SNR'))};

hCalcpopup = uicontrol(...
    'Parent', hCalcpanel, ...
    'Style', 'popupmenu', ...
    'String', calc_parms, ...
    'Enable', 'on', ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'CalculatePopup', ...
    'BackgroundColor', 'White', ...
    'Callback', @calc_popup_Callback);

% +----------% Layout manager for calculation option panel ------------------%
hLmanCalc = siglayout.gridbaglayout(hCalcpanel, 'VerticalGap',5, 'HorizontalGap', 5);
hLmanCalc.add(hCalctxt, 1, 1, 'MinimumHeight', 20, 'MinimumWidth',160, 'TopInset',4, 'LeftInset', 8);
hLmanCalc.add(hCalcpopup, 1, 2, 'MinimumWidth', 160, 'LeftInset', 10);
hLmanCalc.HorizontalWeights = [0.5 1];
hLmanCalc.clean();
hLmanCalc.update();

% Main panel
hParmpanel = uipanel(...
    'Parent', hGui, ...
    'Title', getString(message('phased:apps:sonareqapp:SonarSettings')), ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'ParamPanel');

% Dummy panel to hold common parameter.
hCommonParmpanel = uipanel(...
    'Parent', hParmpanel, ...
    'Tag', 'CommonParmPanel', ...
    'BorderType', 'none', ...
    'BorderWidth', 0, ...
    'BackgroundColor', BKGND_COLOR);

% Dummy panel for SourceLevel ans SNR ui controls
hPowSnrpanel = uipanel(...
    'Parent', hParmpanel, ...
    'Tag', 'PowSnrPanel', ...
    'BorderType', 'none', ...
    'BorderWidth', 0, ...
    'BackgroundColor', BKGND_COLOR);

% +---
% Parameter Text boxes

hNoiseLevel = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:NoiseLevel')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'NLLabel', ...
    'BackgroundColor', BKGND_COLOR);

hDItxt = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:DirectivityIndex')) ':'], ...
    'Visible','on',...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'DiLabel', ...
    'BackgroundColor', BKGND_COLOR);

hTStxt = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:TargetStrength')) ':'], ...
    'Visible','off',...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'TsLabel', ...
    'BackgroundColor', BKGND_COLOR);

hTLTx = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:TransmissionLoss')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'TransLossLabel', ...
    'BackgroundColor', BKGND_COLOR);

hRFreqTx = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:Frequency')) ':'], ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'RfreqLabel', ...
    'BackgroundColor', BKGND_COLOR);

hRDepthTx = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:Depth')) ':'], ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'RDepthLabel', ...
    'BackgroundColor', BKGND_COLOR);

hSnrtxt = uicontrol(...
    'Parent', hPowSnrpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:SNR')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'SnrLabel', ...
    'BackgroundColor', BKGND_COLOR);

hSourceLeveltxt = uicontrol(...
    'Parent', hPowSnrpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:SourceLevel')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'SourceLevelLabel', ...
    'BackgroundColor', BKGND_COLOR);

% Parameter Edit boxes

hModetxt = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:Mode')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'Tag', 'ModeLabel', ...
    'BackgroundColor', BKGND_COLOR);

Mode_parms = {getString(message('phased:apps:sonareqapp:Active')), ...
    getString(message('phased:apps:sonareqapp:Passive'))};

hModepopup = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'popupmenu', ...
    'String', Mode_parms, ...
    'Enable', 'on', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'Tag', 'ModePopup', ...
    'BackgroundColor', 'White', ...
    'Callback', @mode_popup_Callback);

hNoiseLeveledt = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'edit', ...
    'String', '73', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'NLEdit', ...
    'Callback', @noiselevel_Callback);
setappdata(hNoiseLeveledt, 'NL', get(hNoiseLeveledt, 'String'));

hDIedt = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'edit', ...
    'String','20', ...
    'Visible','on',...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'DIEdit', ...
    'Callback', @di_Callback);
setappdata(hDIedt, 'DI', get(hDIedt, 'String'));

hTSedt = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'edit', ...
    'String','25', ...
    'Visible','off',...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'TSEdit', ...
    'Callback', @ts_Callback);
setappdata(hTSedt, 'TS', get(hTSedt, 'String'));

% +----- Sonar config params ---+

hTLedt = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'edit', ...
    'String', '78', ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'TranslossEdit', ...
    'Callback', @transloss_Callback);
setappdata(hTLedt, 'transloss', get(hTLedt, 'String'));

hRFreqedt = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'edit', ...
    'String', '2', ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'RngFreqEdit', ...
    'Callback', @rng_freq_Callback);
setappdata(hRFreqedt, 'rFreq', get(hRFreqedt, 'String'));

hRDepthedt = uicontrol(...
    'Parent', hParmpanel, ...
    'Style', 'edit', ...
    'String', '10000', ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'RngDepthEdit', ...
    'Callback', @rng_depth_Callback);
setappdata(hRDepthedt, 'rDepth', get(hRDepthedt, 'String'));

hSourceLeveledt = uicontrol(...
    'Parent', hPowSnrpanel, ...
    'Style', 'edit', ...
    'String', '220', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'SourceLevelEdit', ...
    'Callback', @SourceLevel_Callback);
setappdata(hSourceLeveledt, 'peakSourceLevel', get(hSourceLeveledt, 'String'));

hSnredt = uicontrol(...
    'Parent', hPowSnrpanel, ...
    'Style', 'edit', ...
    'String', '10', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'SnrEdit', ...
    'Callback', @snr_Callback);
setappdata(hSnredt, 'snr', get(hSnredt, 'String'));

hSnrDetbtn = uicontrol(...
    'Parent', hPowSnrpanel, ...
    'Style', 'pushbutton', ...
    'String', DBLRTARROW, ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tooltip', getString(message('phased:apps:sonareqapp:SnrBtnToolTip')), ...
    'Tag', 'SnrDetPushBtn', ...
    'Callback', @snr_det_btn_Callback);

%Panel to group Sonar Transmission Loss related ui controls
hmediumPanel = uipanel(...
    'Parent', hParmpanel, ...
    'Tag', 'SonarConfigPanel', ...
    'BorderType', 'none', ...
    'BorderWidth', 0, ...
    'BackgroundColor', BKGND_COLOR);

hTLbtn = uicontrol(...
    'Parent', hmediumPanel, ...
    'Style', 'pushbutton', ...
    'String', DBLRTARROW, ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tooltip', getString(message('phased:apps:sonareqapp:TlBtnToolTip')), ...
    'Tooltip', 'Click to calculate Transmission Loss from range', ...
    'Tag', 'TLDetPushBtn', ...
    'Callback', @tl_btn_Callback);

%---------------------------------------------------%
% Units UI controls                                 %
%---------------------------------------------------%

hNoiseLevelunit = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'text', ...
    'String', ['dB//1' char(181) 'Pa'],...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'NoiselevelUnit');

hDIunit = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'text', ...
    'String', getString(message('phased:apps:sonareqapp:decibel')), ...
    'Visible','on',...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'DIUnit');

hTSunit = uicontrol(...
    'Parent', hCommonParmpanel, ...
    'Style', 'text', ...
    'String', getString(message('phased:apps:sonareqapp:decibel')), ...
    'Visible','off',...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'TSUnit');

Freq_units = {getString(message('phased:apps:sonareqapp:Hz')),...
    getString(message('phased:apps:sonareqapp:kHz')),...
    getString(message('phased:apps:sonareqapp:MHz'))};
    
hRFrequnit = uicontrol(...
    'Parent', hmediumPanel, ...
    'Style', 'popup', ...
    'String', Freq_units, ...
    'Value',2,...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'RFreqUnit',...
    'Callback', @range_freq_unit_Callback);

Range_units = {getString(message('phased:apps:sonareqapp:meter')), ...
    getString(message('phased:apps:sonareqapp:kmeter')), ...
    getString(message('phased:apps:sonareqapp:miles')), ...
    getString(message('phased:apps:sonareqapp:nautmi'))};

hRDepthunit = uicontrol(...
    'Parent', hmediumPanel, ...
    'Style', 'popup', ...
    'String', Range_units, ...
    'Value',1,...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'RDepthUnit',...
    'Callback', @range_depth_unit_Callback);

hTLunit = uicontrol(...
    'Parent', hmediumPanel, ...
    'Style', 'text', ...
    'Visible', 'off', ...
    'String', getString(message('phased:apps:sonareqapp:decibel')), ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'TransLossUnit');

hSourceLevelunit = uicontrol(...
    'Parent', hPowSnrpanel, ...
    'Style', 'text', ...
    'String', ['dB//1' char(181) 'Pa'], ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'SourceLevelUnit');

hSnrunit = uicontrol(...
    'Parent', hPowSnrpanel, ...
    'Style', 'text', ...
    'String', getString(message('phased:apps:sonareqapp:decibel')), ...
    'Visible', 'on', ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'SnrUnit');

% +------- Constants for layout ------+
EDIT_WIDTH = 85;
EDIT_HEIGHT = 25;
TEXT_HEIGHT = 25;
SNRBTN_HEIGHT = 25;
if ismac
    UNIT_WIDTH = 75;
else
    UNIT_WIDTH = 70;
end
UNIT_HEIGHT = 25;
if isunix
    % In unix, the unit popups need some top inset to align with text box.
    UNIT_INSET = 1;
elseif ispc
    UNIT_INSET = 0;
end
TOP_INSET = 5;
% +--------------------------------------------------------------+
% +-------- Layout mgr for common parameter ---------------------+
% +--------------------------------------------------------------+

hLmanCommParm = siglayout.gridbaglayout(hCommonParmpanel, ...
    'VerticalGap', 6, 'HorizontalGap', 5);
hLmanCommParm.add(hModetxt, 1, 1, 'Fill', 'Horizontal','TopInset', 5, 'LeftInset', 2);
hLmanCommParm.add(hModepopup, 1, 2:3, 'Fill','Horizontal');

hLmanCommParm.add(hNoiseLevel, 2, 1, 'Fill', 'Horizontal', ...
    'TopInset', 0, 'LeftInset', 2);
hLmanCommParm.add(hDItxt, 3, 1, 'Fill', 'Horizontal', ...
    'TopInset', 0, 'LeftInset', 2);
hLmanCommParm.add(hTStxt, 4, 1, 'Fill', 'Horizontal', ...
    'TopInset', 0, 'LeftInset', 2);
hLmanCommParm.add(hNoiseLeveledt, 2, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT);
hLmanCommParm.add(hDIedt, 3, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT);
hLmanCommParm.add(hTSedt, 4, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT);
hLmanCommParm.add(hNoiseLevelunit, 2, 3, 'MinimumWidth', UNIT_WIDTH, ...
    'TopInset', TOP_INSET);
hLmanCommParm.add(hDIunit,3, 3, 'MinimumWidth', UNIT_WIDTH, ...
    'MinimumHeight', UNIT_HEIGHT, 'TopInset', TOP_INSET);
hLmanCommParm.add(hTSunit, 4, 3, 'MinimumWidth', UNIT_WIDTH, ...
    'MinimumHeight', UNIT_HEIGHT, 'TopInset', TOP_INSET);
hLmanCommParm.HorizontalWeights = [ 1 0.009 0];

% +-------------------------------------------------------+
% Layout for SourceLevel and snr uicontrols on a dummy  panel --+
% +-------------------------------------------------------+
hLmanPowSnr = siglayout.gridbaglayout(hPowSnrpanel, 'VerticalGap', 5, ...
    'HorizontalGap', 6);
layout_SourceLevel(1);
layout_snr(3);
hLmanPowSnr.HorizontalWeights = [0 0 1 0];

%---------------------------------------------------%
% Layout manager: Parameter panel                   %
%---------------------------------------------------%
GUTTER = 10;
calc_type = get(hCalcpopup, 'Value');
hLmanParm = siglayout.gridbaglayout(hParmpanel, ...
    'VerticalGap', 5,...
    'HorizontalGap', 5);
hLmanParm.add(hCommonParmpanel,1, 1, 'MinimumHeight', 5.5*TEXT_HEIGHT, ...
    'MinimumWidth', GUI_WIDTH - GUTTER, 'TopInset', 8, 'Anchor', 'NorthWest');
hLmanParm.add(hmediumPanel, 2, 1, 'MinimumWidth', GUI_WIDTH - GUTTER, ...
    'MinimumHeight', 2.6*TEXT_HEIGHT, 'Anchor', 'NorthWest', 'TopInset', -5);
hLmanParm.add(hPowSnrpanel, 4, 1, 'MinimumHeight', 3*TEXT_HEIGHT, ...
    'MinimumWidth', GUI_WIDTH - GUTTER, 'Anchor', 'NorthEast', 'TopInset', -5);
hLmanParm.VerticalWeights = [0 0 1];
clean(hLmanParm);

% -------------- Layout Soanrconfig panel --------------------%
% Sonar configuration layout handled in the helper function
render_sonarconf(4, 3);
% Draw lines bounding the sonar config components.
% TopLine (combination of two lines with slightly different shade)
topwhitepos = 0.975;
topblackpos = topwhitepos - 0.003;
rtedgepos = 0.97;
hLineWT = annotation(hmediumPanel,'line', 'Position', [0.015 topwhitepos rtedgepos 0], ...
    'Color', [0.2 0.2 0.2]);   % Manually positioned
hLineBT = annotation(hmediumPanel,'line', 'Position', [0.015 topblackpos rtedgepos 0], ...
    'Color',[1 1 1] , 'LineWidth', 0.5);   % Manually positioned

bottomwhitepos = 0.965;
bottomblackpos = bottomwhitepos - 0.003;
% Bottom line
hLineWB = annotation(hPowSnrpanel,'line', 'Position', [0.015 bottomwhitepos rtedgepos 0], ...
    'Color',  [0 0 0]);   % Manually positioned
hLineBB = annotation(hPowSnrpanel,'line', 'Position', [0.015 bottomblackpos rtedgepos 0], ...
    'Color',[1 1 1], 'LineWidth', 0.5);   % Manually positioned

%+----------------------------------------------------------------------+%
%                         Result panel                                   %
%+----------------------------------------------------------------------+%
hResultpanel = uipanel(...
    'Parent', hGui, ...
    'FontSize', 10, ...
    'Tag', 'ResultPanel', ...
    'BorderType', 'none', ...
    'BackgroundColor', BKGND_COLOR);

hResultParmtxt = uicontrol(...
    'Parent', hResultpanel, ...
    'Style', 'text', ...
    'FontSize', 10, ...
    'Visible', 'off', ...
    'HorizontalAlignment', 'Left', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'ResultParmTxt');

hResulttxt = uicontrol(...
    'Parent', hResultpanel, ...
    'Style', 'text', ...
    'FontSize', 10, ...
    'Visible', 'off', ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'ResultTxt');

hResultTLunit = uicontrol(...
    'Parent', hResultpanel, ...
    'Style', 'text', ...
    'String', getString(message('phased:apps:sonareqapp:decibel')), ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'ResultTLUnit');

hResultSourceLevelunit = uicontrol(...
    'Parent', hResultpanel, ...
    'Style', 'text', ...
    'String', ['dB//1' char(181) 'Pa'], ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'ResultSLUnit');

hResultSnrunit = uicontrol(...
    'Parent', hResultpanel, ...
    'Style', 'text', ...
    'String', getString(message('phased:apps:sonareqapp:decibel')), ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'ResultSnrUnit');

hResultRangeunit = uicontrol(...
    'Parent', hResultpanel, ...
    'Style', 'popup', ...
    'String', Range_units, ...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'ResultRangeUnit',...
    'Callback', @result_range_unit_Callback);

% +-------------------------------+ Layout for Result panel  -------|
% +-------------------------------+
INSET = 20;
hLmanResult = siglayout.gridbaglayout(hResultpanel, 'VerticalGap', 5, ...
    'HorizontalGap', 5);
hLmanResult.add(hResultParmtxt, 1, 1, 'Fill', 'Horizontal', ...
    'LeftInset', 6);
set(hResultParmtxt, 'Visible', 'on');
hLmanResult.add(hResulttxt, 1, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'Anchor', 'Northwest', 'TopInset', TOP_INSET);
set(hResulttxt, 'Visible', 'on');
hLmanResult.add(hResultRangeunit, 1, 3);
hLmanResult.add(uicontrol('Style', 'text', 'Visible', 'off'), 1, 4, ...
    'MinimumWidth', 0.25);     % Dummy component to push unit to left.
hLmanResult.HorizontalWeights = [1 0 0 0];
% Rest of dynamic layout management is handled in render_result()

%-------------------------------------------------------------------------%
% Transmission Loss calculation ui controls
%-------------------------------------------------------------------------%
hTLpanel = uipanel( ...
    'Title', getString(message('phased:apps:sonareqapp:CalcTL')), ...
    'TitlePosition', 'lefttop', ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'BorderType', 'none', ...
    'BorderWidth', 0, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'TLPanel');

hTgtRngTx = uicontrol(...
    'Parent', hTLpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:Range')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'Tag', 'TgtRangeTxLabel', ...
    'BackgroundColor', BKGND_COLOR);

hFreqlbl = uicontrol(...
    'Parent', hTLpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:Frequency')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'freqLabel');

hDepthlbl = uicontrol(...
    'Parent', hTLpanel, ...
    'Style', 'text', ...
    'String',[getString(message('phased:apps:sonareqapp:Depth')) ':'] , ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'depthLabel');

%Edit boxes

hTgtRngTxedt = uicontrol(...
    'Parent', hTLpanel, ...
    'Style', 'edit', ...
    'String', '10000', ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'TgtRangeTxEdit', ...
    'Callback', @tgtrangetx_Callback);
setappdata(hTgtRngTxedt, 'tgtrangetx', get(hTgtRngTxedt, 'String'));

hFreqedt = uicontrol(...
    'Parent', hTLpanel, ...
    'Style', 'edit', ...
    'String', '2', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'freqEdit', ...
    'Callback', @freq_Callback);
setappdata(hFreqedt, 'freq', get(hFreqedt, 'String'));

hDepthedt = uicontrol(...
    'Parent', hTLpanel, ...
    'Style', 'edit', ...
    'String', '10000', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'depthEdit', ...
    'Callback', @depth_Callback);
setappdata(hDepthedt, 'depth', get(hDepthedt, 'String'));

%unit boxes

hTgtRngTxunit = uicontrol(...
    'Parent', hmediumPanel, ...
    'Style', 'popup', ...
    'Visible', 'off', ...
    'String', Range_units, ...
    'Value',1,...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'TgtRngTxUnit',...
    'Callback', @target_range_unit_Callback);

hFrequnit = uicontrol(...
    'Parent', hmediumPanel, ...
    'Style', 'popup', ...
    'Visible', 'off', ...
    'String', Freq_units, ...
    'Value',2,...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'HorizontalAlignment', 'Left', ...
    'Tag', 'FreqUnit',...
    'Callback', @freq_unit_Callback);

hDepthunit = uicontrol(...
    'Parent', hmediumPanel, ...
    'Style', 'popup', ...
    'Visible', 'off', ...
    'String', Range_units, ...
    'Value',1,...
    'FontSize', 10, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'DepthUnit',...
    'Callback', @depth_unit_Callback);

%-------------------------------------------------------------------------%
% Layout for Calculation of  Transmission Loss
%-------------------------------------------------------------------------%
hLmanTL = siglayout.gridbaglayout(hTLpanel, 'HorizontalGap', 5, ...
    'VerticalGap', 5);

hLmanTL.add(hTgtRngTx, 1, 1, 'Fill', 'Horizontal', 'TopInset', 15, 'LeftInset', 7);
hLmanTL.add(hFreqlbl, 2, 1, 'Fill', 'Horizontal', 'TopInset', 5, 'LeftInset', 7);
PfaLabelWidth = 183;
hLmanTL.add(hDepthlbl, 3, 1, 'MinimumWidth', PfaLabelWidth, 'TopInset', 5, ...
    'Anchor', 'West', 'LeftInset', 7);
hLmanTL.add(hTgtRngTxedt, 1, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT, ...
    'TopInset', 15,'Anchor', 'West');
hLmanTL.add(hFreqedt, 2, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT, ...
    'Anchor', 'West');
hLmanTL.add(hDepthedt, 3, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT, ...
    'Anchor', 'West');
hLmanTL.add(hTgtRngTxunit, 1, 3, 'MinimumWidth', UNIT_WIDTH, ...
    'MinimumHeight', UNIT_HEIGHT, ...
    'TopInset', 15,'Anchor', 'West');
hLmanTL.add(hFrequnit, 2, 3, 'MinimumWidth', UNIT_WIDTH, ...
    'MinimumHeight', UNIT_HEIGHT, ...
    'Anchor', 'West','TopInset', UNIT_INSET);
hLmanTL.add(hDepthunit, 3, 3, 'MinimumWidth', UNIT_WIDTH, ...
    'MinimumHeight', UNIT_HEIGHT, ...
    'Anchor', 'West','TopInset', UNIT_INSET);

%-------------------------------------------------------------------------%
% Detection settings layout ui controls
%-------------------------------------------------------------------------%
hSnrDetpanel = uipanel( ...
    'Title', getString(message('phased:apps:sonareqapp:DetSettings')), ...
    'TitlePosition', 'lefttop', ...
    'Visible', 'off', ...
    'FontSize', 10, ...
    'BorderType', 'none', ...
    'BorderWidth', 0, ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'SnrDetPanel');

hProbDetlbl = uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'text', ...
    'String', [getString(message('phased:apps:sonareqapp:ProbDetection')) ':'], ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'ProbDetLabel');

hProbFalbl = uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'text', ...
    'String',[getString(message('phased:apps:sonareqapp:ProbFalseAlarm')) ':'] , ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'ProbFaLabel');

hNumPulselbl = uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'text', ...
    'String',[getString(message('phased:apps:sonareqapp:NumPulses')) ':'] , ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'NumSampleLabel');

hSwerlingNumlbl= uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'text', ...
    'String',[getString(message('phased:apps:sonareqapp:SwerlingNum')) ':'] , ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', BKGND_COLOR, ...
    'Tag', 'SwerlingNumLabel');

% Edit inputs

hProbDetedt = uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'edit', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'ProbDetEdit', ...
    'Callback', @pd_Callback);

hProbFaedt = uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'edit', ...
    'String', '1e-3', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'ProbFaEdit', ...
    'Callback', @pfa_Callback);
setappdata(hProbFaedt, 'probfa', get(hProbFaedt, 'String'));

hNumPulseedt = uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'edit', ...
    'String', '1', ...
    'FontSize', 10, ...
    'HorizontalAlignment', EDIT_ALIGN, ...
    'BackgroundColor', EDIT_COLOR, ...
    'Tag', 'NumPulseEdit', ...
    'Callback', @numpulse_Callback);
setappdata(hNumPulseedt, 'numpulse', get(hNumPulseedt, 'String'));

hSwerlingNumpopup = uicontrol(...
    'Parent', hSnrDetpanel, ...
    'Style', 'popupmenu', ...
    'String', {'0','1','2','3','4' }, ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'right', ...
    'BackgroundColor', 'White', ...
    'Tag', 'SwerlingNumpopup', ...
    'Callback', @swerling_num_Callback);

%-------------------------------------------------------------------------%
% Layout for detection setting panel
%-------------------------------------------------------------------------%
hLmanSnrDet = siglayout.gridbaglayout(hSnrDetpanel, 'HorizontalGap', 5, ...
    'VerticalGap', 5);
hLmanSnrDet.add(hProbDetlbl, 1, 1, 'Fill', 'Horizontal', 'TopInset', 20, 'LeftInset', 7);
PfaLabelWidth = 183;
hLmanSnrDet.add(hProbFalbl, 2, 1, 'MinimumWidth', PfaLabelWidth, 'TopInset', 5, ...
    'Anchor', 'West', 'LeftInset', 7);
hLmanSnrDet.add(hNumPulselbl, 3, 1, 'Fill', 'Horizontal', 'TopInset', 5, 'LeftInset', 7);
hLmanSnrDet.add(hSwerlingNumlbl, 4, 1, 'Fill', 'Horizontal', 'TopInset', 5, 'LeftInset', 7);
hLmanSnrDet.add(hProbDetedt, 1, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT, ...
    'TopInset', 15, 'Anchor', 'West');
hLmanSnrDet.add(hProbFaedt, 2, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT, ...
    'Anchor', 'West');
hLmanSnrDet.add(hNumPulseedt, 3, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'MinimumHeight', EDIT_HEIGHT, ...
    'Anchor', 'West');
hLmanSnrDet.add(hSwerlingNumpopup, 4, 2, 'MinimumWidth', EDIT_WIDTH, ...
    'Anchor', 'Northwest');
hLmanSnrDet.HorizontalWeights = [0 1];

%-------------------------------------------------------------------------%
% Layout manager for main figure
%-------------------------------------------------------------------------%
hLmanGui = siglayout.gridbaglayout(hGui, 'VerticalGap' , 5, 'HorizontalGap', 5);
hLmanGui.add(hCalcpanel, 1, 1, 'MinimumWidth', GUI_WIDTH-GUTTER, ...
    'MinimumHeight', 30, ...
    'Anchor', 'NorthWest');
hLmanGui.add(hParmpanel, 2, 1, 'Fill', 'Vertical', ...
    'MinimumWidth', GUI_WIDTH - GUTTER, ...
    'Anchor', 'NorthWest'); %'MaximumHeight', 100, ...
hLmanGui.add(hResultpanel, 3, 1, 'MinimumWidth', GUI_WIDTH - GUTTER, ...
    'MinimumHeight', 40, ...
    'Anchor', 'Northwest');
hLmanGui.VerticalWeights = [0 1 0];
hLmanGui.clean();
hLmanGui.update();

%-------------------------------------------------------------------------%
% GUI Initialization
%-------------------------------------------------------------------------%
% Parameter initialization
calc_type = get(hCalcpopup, 'Value');
mode_type = get(hModepopup,'Value');
update_tl_Callback();
mode_popup_Callback();
render_sonarconfig_ui(calc_type,mode_type);
render_SourceLevel_snr(calc_type);
calculation_Callback();  % Initial setup
movegui(hGui,'center');
set(hGui, 'Visible', 'on', 'HandleVisibility', 'off');
if ispc
    % In windows, all the uicontrols in the GUI shift down when the GUI is
    % resized for the first time. This happens only once. As a workaround
    % for this behavior of the layout, the gui is manually resized and back
    % again.
    resize_gui('on', RSZ_HEIGHT);
    resize_gui('off', RSZ_HEIGHT);
end

%-------------------------------------------------------------------------%
% Callback functions
%-------------------------------------------------------------------------%
    function calculation_Callback(~, ~)
        % Main callback: Handles the calculation of parameters using
        % sonareq functions.
        calc_type = get(hCalcpopup, 'Value');
        mode_type = get(hModepopup,'Value');
        
        noiselevel = str2double(get(hNoiseLeveledt, 'String'));
        di = str2double(get(hDIedt,'String'));
        target_strength = str2double(get(hTSedt, 'String'));
        rng = str2double(get(hTgtRngTxedt, 'String'))* Units.Range;
        freq = str2double(get(hFreqedt, 'String'))* Units.Freq;
        depth = str2double(get(hDepthedt, 'String'))* Units.Depth;
        freq1 = str2double(get(hRFreqedt, 'String'))* Units.RangeFreq;
        depth1 = str2double(get(hRDepthedt, 'String'))* Units.RangeDepth;
        TL = str2double(get(hTLedt, 'String'));
        SL = str2double(get(hSourceLeveledt, 'String'));
        SNR = str2double(get(hSnredt, 'String'));
        
        % Calculation logic rawresult : results in default units
        if  mode_type==1  %Active
            switch calc_type
                
                case 1 % Target range
                    
                    try     % To prevent sending negative values of TL to tl2range function
                        rawresult1 = sonareqtl(SL, SNR, noiselevel,di,target_strength);
                        validateattributes(rawresult1, {'numeric'},{'positive'});
                    catch errObj
                        txt = getString(message('phased:apps:sonareqapp:RangeErrorDuetoNegTL'));
                        errordlg(txt,getString(message('phased:apps:sonareqapp:CheckInput')),'modal');
                        return;
                    end
                    rawresult = tl2range(rawresult1,freq1,depth1);
                    result = rawresult * Units.ResultRange;
                    
                case 2 % Transmission Loss
                    
                    try     % To prevent negative values of TL
                        rawresult = sonareqtl(SL, SNR, noiselevel,di,target_strength);
                        validateattributes(rawresult, {'numeric'},{'positive'});
                    catch errObj
                        txt = getString(message('phased:apps:sonareqapp:TLError'));
                        errordlg(txt,getString(message('phased:apps:sonareqapp:CheckInput')),'modal');
                        return;
                    end
                    result = rawresult;
                    
                case 3 % Peak SourceLevel
                    rawresult = sonareqsl(SNR, noiselevel,di,TL,target_strength);
                    result = rawresult;
                    
                case 4 % SNR
                    rawresult = sonareqsnr(SL, noiselevel,di, TL,target_strength);
                    result = rawresult;
            end % Switch case
            
        elseif mode_type==2   %Passive
            
            switch calc_type
                case 1  % Target range
                    try     % To prevent sending negative values of TL to tl2range function
                        rawresult1 = sonareqtl(SL, SNR, noiselevel,di);
                        validateattributes(rawresult1, {'numeric'},{'positive'});
                    catch errObj
                        txt = getString(message('phased:apps:sonareqapp:RangeErrorDuetoNegTL'));
                        errordlg(txt,getString(message('phased:apps:sonareqapp:CheckInput')),'modal');
                        return;
                    end
                    rawresult = tl2range(rawresult1,freq1,depth1);
                    result = rawresult * Units.ResultRange;
                    
                case 2 % Transmission Loss
                    
                    try     % to prevent negative values of TL
                        rawresult = sonareqtl(SL, SNR, noiselevel,di);
                        validateattributes(rawresult, {'numeric'},{'positive'});
                    catch errObj
                        txt = getString(message('phased:apps:sonareqapp:TLError'));
                        errordlg(txt,getString(message('phased:apps:sonareqapp:CheckInput')),'modal');
                        return;
                    end
                    result = rawresult;
                    
                case 3 % Peak SourceLevel
                    rawresult = sonareqsl(SNR, noiselevel,di,TL);
                    result = rawresult;
                    
                case 4 % SNR
                    rawresult = sonareqsnr(SL, noiselevel,di, TL);
                    result = rawresult;
            end % Switch case
        end
        
        % Save result value as guidata, to be used by get_parm_struct()
        guidata(hGui, rawresult);
        render_result(calc_type, result);  % Render result
    end % calculation callback

    function calc_popup_Callback(~,~)
        % Callback for main calculation popup
        calc_type = get(hCalcpopup, 'Value');
        render_sonarconfig_ui(calc_type,mode_type);
        render_SourceLevel_snr(calc_type);
        calculation_Callback();
    end

    function mode_popup_Callback(~,~)
        mode_type = get(hModepopup, 'Value');
        if mode_type == 1       % Active
            hTStxt.Visible = 'on';
            hTSedt.Visible = 'on';
            hTSunit.Visible = 'on';
        elseif mode_type == 2   % Passive
            hTStxt.Visible = 'off';
            hTSedt.Visible = 'off';
            hTSunit.Visible = 'off';
        end
        calculation_Callback();
    end

    function noiselevel_Callback(hObject, eventdata) %#ok<*INUSD>
        % Noise level callback
        noiselevel = str2double(get(hNoiseLeveledt, 'String'));
        try
            validateattributes(noiselevel, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan', 'nonempty', 'real'});
        catch errObj %#ok<*NASGU>
            holetxt = get(hNoiseLevel, 'String');
            txt=getString(message('phased:apps:sonareqapp:RealValueErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hNoiseLeveledt, 'String', getappdata(hNoiseLeveledt, 'NL'));
            return;
        end
        setappdata(hNoiseLeveledt, 'NL', get(hNoiseLeveledt, 'String'));
        calculation_Callback();
    end

    function di_Callback(hObject, eventdata)
        % Dirctivity index callback
        directivityindex = str2double(get(hDIedt, 'String'));
        di_validate_option = {'scalar', 'finite', 'nonnan', 'nonempty'};
        %err_msg_id = 'Enter proper value';
        try
            validateattributes(directivityindex, {'numeric'}, di_validate_option);
        catch errObj
            holetxt = get(hDItxt, 'String');
            txt= getString(message('phased:apps:sonareqapp:RealValueErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hDIedt, 'String', getappdata(hDIedt, 'DI'));
            return;
        end
        setappdata(hDIedt, 'DI', get(hDIedt, 'String'));
        calculation_Callback();
    end

    function ts_Callback(hObject, eventdata)
        % Target strength callback
        tgtrcs = str2double(get(hTSedt, 'String'));
        rcs_validate_option = {'scalar', 'real', 'nonnan', 'finite', 'nonempty'};
        %err_msg_id = 'Enter proper value';
        %         end
        try
            validateattributes(tgtrcs, {'numeric'}, rcs_validate_option);
        catch errObj
            holetxt = get(hTStxt, 'String');
            txt=getString(message('phased:apps:sonareqapp:RealValueErr',holetxt(1:end-1)));
            errordlg(txt,  getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hTSedt, 'String', getappdata(hTSedt, 'TS'));
            return;
        end
        setappdata(hTSedt, 'TS', get(hTSedt, 'String'));
        calculation_Callback();
    end

    function rng_freq_Callback(hObject, eventdata)
        % Callback for frequency
        rfreq = str2double(get(hRFreqedt, 'String'));
        try
            validateattributes(rfreq, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan', 'positive','nonempty', 'real'});
        catch errObj
            holetxt = get(hRFreqTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:NzPosNumErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hRFreqedt, 'String', getappdata(hRFreqedt, 'rFreq'));
            return;
        end
        % Check if frequency is less than 2MHz
        rfreq_valid = rfreq;
        cond1 = ((get(hRFrequnit, 'Value') == 1) && (rfreq_valid > 2000000));
        cond2 = ((get(hRFrequnit, 'Value') == 2) && (rfreq_valid > 2000));
        cond3 = ((get(hRFrequnit, 'Value') == 3) && (rfreq_valid > 2));
        if (cond1||cond2||cond3)
            holetxt = get(hRFreqTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinFreqErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
%             set(hRFreqedt, 'String', getappdata(hRFreqedt,'rfreq'));
            set(hRFreqedt, 'String', num2str(2));
            set(hRFrequnit, 'Value', 2);  % Reset it back to kHz
            Units.RangeFreq = 1e3;
        end
        setappdata(hRFreqedt, 'rFreq', get(hRFreqedt, 'String'));
        calculation_Callback();
    end

    function rng_depth_Callback(hObject, eventdata) %#ok<*INUSD>
        %  depth callback
        rdepth = str2double(get(hRDepthedt, 'String'));
        try
            validateattributes(rdepth, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan','positive', 'nonempty', 'real'});
        catch errObj
            holetxt = get(hRDepthTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:NzPosNumErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hRDepthedt, 'String', getappdata(hRDepthedt, 'rDepth'));
            return;
        end
        % Check if depth is greater than 2m
        depth_valid = rdepth;
        cond1 = ((get(hRDepthunit, 'Value') == 1) && (depth_valid < 2));
        cond2 = ((get(hRDepthunit, 'Value') == 2) && (depth_valid < 0.002));
        cond3 = ((get(hRDepthunit, 'Value') == 3) && (depth_valid < 2*6.2139e-04));
        cond4 = ((get(hRDepthunit, 'Value') == 4) && (depth_valid < 2*5.3996e-04));
        if (cond1||cond2||cond3||cond4)
            holetxt = get(hRDepthTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinimumMeterErr',holetxt(1:end-1),2));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hRDepthedt, 'String', num2str(10000));
            set(hRDepthunit, 'Value', 1);   % Reset it back to m
            Units.RangeDepth = 1;
        end
        setappdata(hRDepthedt, 'rDepth', get(hRDepthedt, 'String'));
        calculation_Callback();
    end

    function transloss_Callback(hObject, eventdata)
        % Transloss callback
        transloss = str2double(get(hTLedt, 'String'));
        try
            validateattributes(transloss, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan', 'positive', 'nonzero'});
        catch errObj
            holetxt = get(hTLTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:NzPosRealErr',holetxt(1:end-1)));
            errordlg(txt,  getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal')
            set(hTLedt, 'String', getappdata(hTLedt, 'transloss'));
            return;
        end
        setappdata(hTLedt, 'transloss', get(hTLedt, 'String'));
        calculation_Callback();
    end

    function snr_Callback(hObject, eventdata)
        % Validate and call calculation for snr
        snr = str2double(get(hSnredt, 'String'));
        try
            validateattributes(snr, {'numeric'}, ...
                {'scalar', 'finite',  'nonnan', 'nonempty', 'real'});
        catch errObj
            holetxt = get(hSnrtxt, 'String');
            txt = getString(message('phased:apps:sonareqapp:RealValueErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hSnredt, 'String', getappdata(hSnredt, 'snr'));
            return;
        end
        setappdata(hSnredt, 'snr', get(hSnredt, 'String'));
        populate_pdpfa();
        calculation_Callback();
    end

    function SourceLevel_Callback(hObject, eventdata)
        % Validate and call calculation for SourceLevel
        peakSourceLevel = str2double(get(hSourceLeveledt, 'String'));
        try
            validateattributes(peakSourceLevel, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan','real'});
        catch errObj
            holetxt = get(hSourceLeveltxt, 'String');
            txt=getString(message('phased:apps:sonareqapp:RealValueErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hSourceLeveledt, 'String', getappdata(hSourceLeveledt, 'peakSourceLevel'));
            return;
        end
        setappdata(hSourceLeveledt, 'peakSourceLevel', get(hSourceLeveledt, 'String'));
        calculation_Callback();
    end

    function snr_det_btn_Callback(hObject, eventdata)
        btn_label = get(hSnrDetbtn, 'String');
        if strcmp(btn_label, DBLRTARROW)
            % Change string
            set(hSnrDetbtn, 'String', DBLLTARROW);
            set(hSnredt, 'Enable', 'off');
            switch_snrdetpanel('on');
            % Populate Pd and Pfa editbox using values from rocsnr()
            populate_pdpfa();
        elseif strcmp(btn_label, DBLLTARROW)
            switch_snrdetpanel('off');
            drawnow();
            % Change String
            set(hSnrDetbtn, 'String', DBLRTARROW);
            set(hSnredt, 'Enable', 'on');
        end
    end

    function tl_btn_Callback(hObject, eventdata)
        btn_label = get(hTLbtn, 'String');
        if strcmp(btn_label, DBLRTARROW)
            % Change string
            set(hTLbtn, 'String', DBLLTARROW);
            set(hTLedt, 'Enable', 'off');
            switch_tlpanel('on');
        elseif strcmp(btn_label, DBLLTARROW)
            switch_tlpanel('off');
            drawnow();
            % Change String
            set(hTLbtn, 'String', DBLRTARROW);
            set(hTLedt, 'Enable', 'on');
        end
    end

    function gen_report_Callback(~, ~)
        % Generate report and open save dialog
        sparms = get_parm_struct();
        gensonareqreport(sparms);
    end

    function gen_mcode_Callback(~, ~)
        sparms = get_parm_struct();
        gensonareqmcode(sparms);
    end

    function close_Callback(~, ~)
        close(hGui);
    end

    function sonareqapp_help_Callback(~,~)
        helpview([docroot, '\phased\helptargets.map'],  'sonar_app');
    end

    function phased_help_Callback(~, ~)
        % Launch phased array system toolbox documentation
        helpview([docroot, '\phased\helptargets.map'],  'phased_doc');
    end

    function about_Callback(~,~)
        aboutphasedtbx;
    end
% -------------------- Transmission Loss Calculation callbacks -----------------------%
    function tgtrangetx_Callback(hObject, eventdata)
        % Validate and call calculation for range parameter
        tgtrangetx = str2double(get(hTgtRngTxedt, 'String'));
        try
            validateattributes(tgtrangetx, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan', 'positive','nonzero'});
        catch errObj
            holetxt = get(hTgtRngTx, 'String');
            txt =getString(message('phased:apps:sonareqapp:NzPosNumErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal')
            set(hTgtRngTxedt, 'String', getappdata(hTgtRngTxedt, 'tgtrangetx'));
            return;
        end
        % Check if range is greater than 1m
        rng_valid = tgtrangetx;
        cond1 = ((get(hTgtRngTxunit, 'Value') == 1) && (rng_valid < 1));
        cond2 = ((get(hTgtRngTxunit, 'Value') == 2) && (rng_valid < 0.001));
        cond3 = ((get(hTgtRngTxunit, 'Value') == 3) && (rng_valid < 6.2139e-04));
        cond4 = ((get(hTgtRngTxunit, 'Value') == 4) && (rng_valid < 5.3996e-04));
        if (cond1||cond2||cond3||cond4)
            holetxt = get(hTgtRngTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinimumMeterErr',holetxt(1:end-1),1));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hTgtRngTxedt, 'String', num2str(10000));
            set(hTgtRngTxunit, 'Value', 1);   % Reset it back to m
            Units.Range = 1;
        end
        update_tl_Callback()
        calculation_Callback();
    end

    function freq_Callback(hObject, eventdata)
        % Validate and call calculation for range parameter
        freq = str2double(get(hFreqedt, 'String'));
        try
            validateattributes(freq, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan', 'positive', 'nonzero'});
        catch errObj
            holetxt = get(hFreqlbl, 'String');
            txt =getString(message('phased:apps:sonareqapp:NzPosNumErr',holetxt(1:end-1)));
            errordlg(txt, getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal')
            set(hFreqedt, 'String', getappdata(hFreqedt, 'freq'));
            return;
        end
        % Check if frequency is less than 2MHz
        freq_valid = freq;
        cond1 = ((get(hFrequnit, 'Value') == 1) && (freq_valid > 2000000));
        cond2 = ((get(hFrequnit, 'Value') == 2) && (freq_valid > 2000));
        cond3 = ((get(hFrequnit, 'Value') == 3) && (freq_valid > 2));
        if (cond1||cond2||cond3)
            holetxt = get(hFreqlbl, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinFreqErr',holetxt(1:end-1)));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hFreqedt, 'String', num2str(2));
            set(hFrequnit, 'Value', 2);   % Reset it back to kHz
            Units.Freq = 1e3;
        end
        update_tl_Callback()
        calculation_Callback();
    end

    function depth_Callback(hObject, eventdata)
        % Validate and call calculation for range parameter
        depth = str2double(get(hDepthedt, 'String'));
        try
            validateattributes(depth, {'numeric'}, ...
                {'scalar', 'finite', 'nonnan', 'positive', 'nonzero'});
        catch errObj
            holetxt = get(hDepthlbl, 'String');
            txt = getString(message('phased:apps:sonareqapp:NzPosNumErr',holetxt(1:end-1)));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal')
            set(hDepthedt, 'String', getappdata(hDepthedt, 'depth'));
            return;
        end
        % Check if depth is greater than 2m
        depth_valid = depth;
        cond1 = ((get(hDepthunit, 'Value') == 1) && (depth_valid < 2));
        cond2 = ((get(hDepthunit, 'Value') == 2) && (depth_valid < 0.002));
        cond3 = ((get(hDepthunit, 'Value') == 3) && (depth_valid < 2*6.2139e-04));
        cond4 = ((get(hDepthunit, 'Value') == 4) && (depth_valid < 2*5.3996e-04));
        if (cond1||cond2||cond3||cond4)
            holetxt = get(hDepthlbl, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinimumMeterErr',holetxt(1:end-1),2));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hDepthedt, 'String', num2str(10000));
            set(hDepthunit, 'Value', 1);   % Reset it back to m
            Units.Depth = 1;
        end
        update_tl_Callback()
        calculation_Callback();
    end

% -------------------- Detection setting callbacks -----------------------%
    function pd_Callback(hObject, eventdata)
        % Prob of detection edit box callback Validate and call calculation
        % for snr parameter
        pd = str2double(get(hProbDetedt, 'String'));
        try
            validateattributes(pd, {'numeric'}, ...
                {'scalar', 'nonnan', 'positive', '>=', 0.1, '<=', 0.99});
        catch errObj
            holetxt = get(hProbDetlbl, 'String');
            errordlg(getString(message('phased:apps:sonareqapp:ProbDetErr', ...
                holetxt(1:end-1), num2str(0.1), num2str(0.99))), getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hProbDetedt, 'String', getappdata(hProbDetedt, 'probdet'));
            return;
        end
        setappdata(hProbDetedt, 'probdet', get(hProbDetedt, 'String'));
        update_snr_Callback();
        calculation_Callback();
    end

    function pfa_Callback(Object, eventdata)
        % Prob of false alarm edit box callback
        pfa = str2double(get(hProbFaedt, 'String'));
        try
            validateattributes(pfa, {'numeric'}, ...
                {'scalar', 'nonnan', 'positive', 'nonzero', '<', 1, '>=', 1e-20});
        catch errObj
            holetxt = get(hProbFalbl, 'String');
            errordlg(getString(message('phased:apps:sonareqapp:ProbFaErr', ...
                holetxt(1:end-1), num2str(1e-20), 1)), getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hProbFaedt, 'String', getappdata(hProbFaedt, 'probfa'));
            return;
        end
        % Make sure Pd is <= 0.99
        if (str2double(get(hProbDetedt, 'String')) > 0.99)
            holetxt = get(hProbDetlbl, 'String');
            errordlg(getString(message('phased:apps:sonareqapp:ProbDetErr', ...
                holetxt(1:end-1), num2str(0.1), num2str(0.99))), getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hProbFaedt, 'String', getappdata(hProbFaedt, 'probfa'));
            return;
        end
        setappdata(hProbFaedt, 'probfa', get(hProbFaedt, 'String'));
        update_snr_Callback();
        calculation_Callback();
    end

    function numpulse_Callback(hObject, eventdata)
        % Number of pulses callback
        numpulse = str2double(get(hNumPulseedt, 'String'));
        try
            validateattributes(numpulse, {'numeric'}, ...
                {'scalar', 'nonnan', 'positive', 'nonzero', 'integer'});
        catch errObj
            holetxt = get(hNumPulselbl, 'String');
            errordlg(getString(message('phased:apps:sonareqapp:NzPosIntErr', ...
                holetxt(1:end-1))), getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hNumPulseedt, 'String', getappdata(hNumPulseedt, 'numpulse'));
            return;
        end
        % Make sure Pd is <= 0.99
        if (str2double(get(hProbDetedt, 'String')) > 0.99)
            holetxt = get(hProbDetlbl, 'String');
            errordlg(getString(message('phased:apps:sonareqapp:ProbDetErr', ...
                holetxt(1:end-1), num2str(0.1), num2str(0.99))), getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hProbFaedt, 'String', getappdata(hProbFaedt, 'probfa'));
            return;
        end
        setappdata(hNumPulseedt, 'numpulse', get(hNumPulseedt, 'String'));
        update_snr_Callback();
        calculation_Callback();
    end

    function swerling_num_Callback(hObject, eventdata)
        % Swerling number popup callback
        update_snr_Callback();
        calculation_Callback();
    end
%
% % -------------------- Unit UI callbacks ---------------------------------%

    function result_range_unit_Callback(hObject, eventdata)
        % Callback for range unit popup in the result panel Convert from m
        % to km, mi, nmi
        switch get(hObject, 'Value')
            case 1  % m
                Units.ResultRange = 1;
            case 2  % m to km
                Units.ResultRange = unitsratio('km', 'm');
            case 3  % m to mi
                Units.ResultRange = unitsratio('mi', 'm');
            case 4  % m to nmi
                Units.ResultRange = unitsratio('nm', 'm');
        end
        calculation_Callback();   % Update calculation
    end

    function range_depth_unit_Callback(hObject, eventdata)
        % Callback for depth unit popup in the Range Calculation. Convert
        % from m to km, mi, nmi
        switch get(hObject, 'Value')
            case 1  % m
                Units.RangeDepth = 1;
            case 2  % m to km
                Units.RangeDepth = unitsratio('m', 'km');
            case 3  % m to mi
                Units.RangeDepth = unitsratio('m', 'mi');
            case 4  % m to nmi
                Units.RangeDepth = unitsratio('m', 'nm');
        end
        % Check if depth is greater than 2m
        depth = str2double(get(hRDepthedt, 'String'));
        depth_unit = get(hRDepthunit, 'Value');
        
        cond1 = ((depth_unit == 1) && (depth < 2));
        cond2 = ((depth_unit == 2) && (depth < 0.002));
        cond3 = ((depth_unit == 3) && (depth < 2*6.2139e-04));
        cond4 = ((depth_unit == 4) && (depth < 2*5.3996e-04));
        if (cond1||cond2||cond3||cond4)
            holetxt = get(hRDepthTx, 'String');
            txt=[holetxt(1:end-1) ' must be >= 2m'];
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hRDepthedt, 'String', num2str(10000));
            set(hRDepthunit, 'Value', 1);   % Reset it back to m
            Units.RangeDepth = 1;
        end
        setappdata(hRDepthedt, 'rDepth', get(hRDepthedt, 'String'));
        calculation_Callback();   % Update calculation
    end

    function range_freq_unit_Callback(hObject, eventdata)
        % Callback for freq unit popup in the Range Calculation. Convert
        % from kHz/MHz to Hz
        switch get(hObject, 'Value')
            case 1  % Hz 
                Units.RangeFreq = 1;
            case 2  % kHz to Hz
                Units.RangeFreq = unitsratio('m', 'km');
            case 3  % MHz to Hz
                Units.RangeFreq = unitsratio('mm', 'km');
        end
        rfreq = str2double(get(hRFreqedt, 'String'));
        rfreq_unit = get(hRFrequnit, 'Value');
        %         rfreq_valid = rfreq*Units.RangeFreq;
        cond1 = ((rfreq_unit == 1) && (rfreq > 2000000));
        cond2 = ((rfreq_unit == 2) && (rfreq > 2000));
        cond3 = ((rfreq_unit == 3) && (rfreq > 2));
        if (cond1||cond2||cond3)
            holetxt = get(hRFreqTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinFreqErr',holetxt(1:end-1)));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hRFreqedt, 'String',num2str(2)); 
            set(hRFrequnit, 'Value', 2);   % Reset it back to kHz
            Units.RangeFreq = 1e3;
        end
        setappdata(hRFreqedt, 'rFreq', get(hRFreqedt, 'String'));
        calculation_Callback();   % Update calculation
    end

    function target_range_unit_Callback(hObject, eventdata)
        % Callback for range unit popup in the TL widget. Convert
        % from m to km, mi, nmi
        switch get(hObject, 'Value')
            case 1  % m
                Units.Range = 1;
            case 2  % m to km
                Units.Range = unitsratio('m', 'km');
            case 3  % m to mi
                Units.Range = unitsratio('m', 'mi');
            case 4  % m to nmi
                Units.Range = unitsratio('m', 'nm');
        end
        % Check if range is greater than 1m
        rng = str2double(get(hTgtRngTxedt, 'String'));
        rng_unit = get(hTgtRngTxunit, 'Value');
        
        cond1 = ((rng_unit == 1) && (rng < 1));
        cond2 = ((rng_unit == 2) && (rng < 0.001));
        cond3 = ((rng_unit == 3) && (rng < 6.2139e-04));
        cond4 = ((rng_unit == 4) && (rng < 5.3996e-04));
        if (cond1||cond2||cond3||cond4)
            holetxt = get(hTgtRngTx, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinimumMeterErr',holetxt(1:end-1),1));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hTgtRngTxedt, 'String', num2str(10000));
            set(hTgtRngTxunit, 'Value', 1);   % Reset it back to m
            Units.Range = 1;
        end
        setappdata(hTgtRngTxedt, 'tgtrangetx', get(hTgtRngTxedt, 'String'));
        update_tl_Callback();
        calculation_Callback();   % Update calculation
    end

    function freq_unit_Callback(hObject, eventdata)
        % Callback for freq unit popup in the Range Calculation. Convert
        % from kHz/MHz to Hz
        switch get(hObject, 'Value')
            case 1  % Hz 
                Units.Freq = 1;
            case 2  % kHz to Hz
                Units.Freq = unitsratio('m', 'km');
            case 3  % MHz to Hz
                Units.Freq = unitsratio('mm', 'km');
        end
        freq = str2double(get(hFreqedt, 'String'));
        freq_unit = get(hFrequnit, 'Value');
        %         rfreq_valid = rfreq*Units.RangeFreq;
        cond1 = ((freq_unit == 1) && (freq > 2000000));
        cond2 = ((freq_unit == 2) && (freq > 2000));
        cond3 = ((freq_unit == 3) && (freq > 2));
        if (cond1||cond2||cond3)
            holetxt = get(hFreqlbl, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinFreqErr',holetxt(1:end-1)));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hFreqedt, 'String', num2str(2));
            set(hFrequnit, 'Value', 2);   % Reset it back to kHz
            Units.Freq = 1e3;
        end
        setappdata(hFreqedt, 'freq', get(hFreqedt, 'String'));
        update_tl_Callback();
        calculation_Callback();   % Update calculation
    end
    function depth_unit_Callback(hObject, eventdata)
        % Callback for depth unit popup in the TL widget. Convert
        % from m to km, mi, nmi
        switch get(hObject, 'Value')
            case 1  % m
                Units.Depth = 1;
            case 2  % m to km
                Units.Depth = unitsratio('m', 'km');
            case 3  % m to mi
                Units.Depth = unitsratio('m', 'mi');
            case 4  % m to nmi
                Units.Depth = unitsratio('m', 'nm');
        end
        % Check if depth is greater than 2m
        depth = str2double(get(hDepthedt, 'String'));
        depth_unit = get(hDepthunit, 'Value');
        
        cond1 = ((depth_unit == 1) && (depth < 2));
        cond2 = ((depth_unit == 2) && (depth < 0.002));
        cond3 = ((depth_unit == 3) && (depth < 2*6.2139e-04));
        cond4 = ((depth_unit == 4) && (depth < 2*5.3996e-04));
        if (cond1||cond2||cond3||cond4)
            holetxt = get(hDepthlbl, 'String');
            txt=getString(message('phased:apps:sonareqapp:MinimumMeterErr',holetxt(1:end-1),2));
            errordlg(txt,getString(message('phased:apps:sonareqapp:InvalidInput')), 'modal');
            set(hDepthedt, 'String', num2str(10000));
            set(hDepthunit, 'Value', 1);   % Reset it back to m
            Units.Depth = 1;
        end
        setappdata(hDepthedt, 'depth', get(hDepthedt, 'String'));
        update_tl_Callback();
        calculation_Callback();   % Update calculation
    end

% --------------------- Helper functions ---------------------------------%
    function render_result(calc_type, result)
        % Renders the result panel and its contents
        range_label = {[getString(message('phased:apps:sonareqapp:Range')) ':']};
        set(hResulttxt, 'String', num2str(result, 4));
        switch calc_type
            case 1  % Range
                set(hResultParmtxt, 'String', range_label);
                set(hResultRangeunit, 'Visible', 'on');
                set(hResultTLunit, 'Visible', 'off')
                set(hResultSourceLevelunit, 'Visible', 'off');
                set(hResultSnrunit, 'Visible', 'off');
                hLmanResult.remove(1,3);
                hLmanResult.add(hResultRangeunit, 1, 3, ...
                    'MinimumWidth', 0.95*UNIT_WIDTH,  'MinimumHeight', 30, ...
                    'TopInset', 0.75,'Anchor', 'Northwest');
                drawnow();
                
            case 2 % Transmission Loss
                set(hResultParmtxt, 'String', [calc_parms{calc_type} ':']);
                set(hResultTLunit, 'Visible', 'on');
                set(hResultRangeunit, 'Visible', 'off');
                set(hResultSourceLevelunit, 'Visible', 'off');
                set(hResultSnrunit, 'Visible', 'off');
                hLmanResult.remove(1,3);
                hLmanResult.add(hResultTLunit, 1, 3, ...
                    'MinimumWidth', 0.95*UNIT_WIDTH,  'MinimumHeight', 30, ...
                    'TopInset',TOP_INSET);
                drawnow();
                
            case 3  % SourceLevel
                set(hResultParmtxt, 'String', [calc_parms{calc_type} ':']);
                set(hResultRangeunit, 'Visible', 'off');
                set(hResultTLunit, 'Visible', 'off')
                set(hResultSourceLevelunit, 'Visible', 'on');
                set(hResultSnrunit, 'Visible', 'off');
                hLmanResult.remove(1,3);
                hLmanResult.add(hResultSourceLevelunit, 1, 3, ...
                    'MinimumWidth',0.95*UNIT_WIDTH, 'MinimumHeight', 30, ...
                    'TopInset',TOP_INSET);
                drawnow();
                
            case 4  % SNR
                set(hResultParmtxt,'String', [calc_parms{calc_type} ':']);
                set(hResultRangeunit, 'Visible', 'off');
                set(hResultSourceLevelunit, 'Visible', 'off');
                set(hResultSnrunit, 'Visible', 'on');
                hLmanResult.remove(1,3);
                hLmanResult.add(hResultSnrunit, 1, 3, ...
                    'MinimumWidth', 0.95*UNIT_WIDTH, 'MinimumHeight', 30, ...
                    'TopInset', TOP_INSET);
                drawnow();
        end
    end

    function render_sonarconfig_ui(calc_type,mode_type)
        if calc_type == 1 && (mode_type == 1||mode_type==2)
            render_sonarconf(1, 2);
            % Range
            toggle_key = {'off', 'off','off','off','off','off' ...
                'off', 'off', 'off','off','off','off' ...
                'off', 'on', 'on','on','on','on','on' ...
                };
            %  Hide Transmission Loss uicontrols when Range is chosen for calculation
            btn_label = get(hTLbtn, 'String');
            % If TL interface setting was open, close it too.
            if strcmp(btn_label, DBLLTARROW)  % If detection setting was open
                switch_tlpanel('off');
                set(hTLedt, 'Enable', 'on');
                set(hTLbtn, 'String', DBLRTARROW);
            end
            toggle_sonarconfig_ui(toggle_key);
            
        elseif calc_type == 2 &&( mode_type == 2||mode_type == 1)
            % Transmission Loss
            render_sonarconf(3, 4);
            toggle_key = {'off', 'off', 'off','off','off','off'...
                'off', 'off', 'off','off','off','off' ...
                'off', 'off', 'off', 'off','off','off','off'...
                };
            
            % Hide Transmission Loss uicontrols when Transmission Loss is chosen for calculation
            btn_label = get(hTLbtn, 'String');
            % If TL interface setting was open, close it too.
            if strcmp(btn_label, DBLLTARROW)  % If detection setting was open
                switch_tlpanel('off');
                set(hTLedt, 'Enable', 'on');
                set(hTLbtn, 'String', DBLRTARROW);
            end
            toggle_sonarconfig_ui(toggle_key);
            
        elseif (calc_type == 3 || calc_type == 4) && (mode_type == 1||mode_type==2)
            render_sonarconf(4, 3);
            % SourceLevel or SNR
            toggle_key = {'on', 'on', 'on','on','on','on'...
                'on', 'on', 'on','on','on','on' ...
                'on', 'off', 'off','off','off','off','off' ...
                };
            toggle_sonarconfig_ui(toggle_key);
            
        end
    end

    function toggle_sonarconfig_ui(toggle_key)
        % Toggles visibility of sonar config ui components
        set(hTLbtn, 'Visible', toggle_key{1});
        set(hTLTx, 'Visible', toggle_key{2});
        set(hTLedt, 'Visible', toggle_key{3});
        set(hTLunit, 'Visible', toggle_key{4});
        set(hTgtRngTx, 'Visible', toggle_key{5});
        set(hTgtRngTxedt, 'Visible', toggle_key{6});
        set(hTgtRngTxunit, 'Visible', toggle_key{7});
        set(hFreqlbl, 'Visible', toggle_key{8});
        set(hFreqedt, 'Visible', toggle_key{9});
        set(hFrequnit, 'Visible', toggle_key{10});
        set(hDepthlbl, 'Visible', toggle_key{11});
        set(hDepthedt, 'Visible', toggle_key{12});
        set(hDepthunit, 'Visible', toggle_key{13});
        set(hRFreqTx, 'Visible', toggle_key{14});
        set(hRFreqedt, 'Visible', toggle_key{15});
        set(hRFrequnit, 'Visible', toggle_key{16});
        set(hRDepthTx, 'Visible', toggle_key{17});
        set(hRDepthedt, 'Visible', toggle_key{18});
        set(hRDepthunit, 'Visible', toggle_key{19});
    end

    function render_SourceLevel_snr(calc_type)
        % Manages dynamic rendering of SourceLevel and SNR related ui controls.
        if (calc_type == 1|| calc_type ==2) % Range/Transmission Loss
            % Show both SourceLevel and snr ui controls
            toggle_key = {'on', 'on', 'on', ...
                'on', 'on', 'on', 'on'};
            toggle_SourceLevel_snr(toggle_key);
            swap_SourceLevel_snr_pos(1, 2);
            
        elseif calc_type == 3  % SourceLevel
            % Hide peak SourceLevel uicontrols when Peak SourceLevel chosen for calc
            % Toggle visibility
            swap_SourceLevel_snr_pos(2, 1);   % Reposition
            % Visibility
            toggle_key = {'off', 'off', 'off', ...
                'on', 'on', 'on', 'on'};
            toggle_SourceLevel_snr(toggle_key);
            
        elseif calc_type == 4   % SNR
            % %             Hide SNR uicontrols when SNR is chosen for calculation
            btn_label = get(hSnrDetbtn, 'String');
            % If detection setting was open, close it too.
            if strcmp(btn_label, DBLLTARROW)  % If detection setting was open
                switch_snrdetpanel('off');
                set(hSnredt, 'Enable', 'on');
                set(hSnrDetbtn, 'String', DBLRTARROW);
            end
            swap_SourceLevel_snr_pos(1, 2);
            toggle_key = {'on', 'on', 'on', ...
                'off', 'off', 'off', 'off'};
            toggle_SourceLevel_snr(toggle_key);
        end
    end

    function swap_SourceLevel_snr_pos(SourceLevel_pos, snr_pos)
        % Swaps position of SourceLevel ans snr uicontrols
        remove_pow_lm();
        remove_snr_lm();
        layout_snr(snr_pos);
        %%drawnow();
        layout_SourceLevel(SourceLevel_pos);
        hLmanPowSnr.HorizontalWeights = [1 1 0 0];
        hLmanPowSnr.update();
    end

    function remove_pow_lm()
        % Remove SourceLevel uicontrols from LM
        hLmanPowSnr.remove(hSourceLeveltxt);
        hLmanPowSnr.remove(hSourceLeveledt);
        hLmanPowSnr.remove(hSourceLevelunit);
    end

    function remove_snr_lm()
        % Remove from LM
        hLmanPowSnr.remove(hSnrtxt);
        hLmanPowSnr.remove(hSnredt);
        hLmanPowSnr.remove(hSnrunit);
        hLmanPowSnr.remove(hSnrDetbtn);
    end

    function layout_SourceLevel(pos)
        % Layout SourceLevel ui controls at row pos
        hLmanPowSnr.add(hSourceLeveltxt, pos, 1:2, 'Fill', 'Horizontal', ...
            'TopInset', 5, 'Anchor', 'NorthWest');
        hLmanPowSnr.add(hSourceLeveledt, pos, 3, 'MinimumWidth', EDIT_WIDTH, ...
            'MinimumHeight', EDIT_HEIGHT, 'Anchor', 'North' );
        hLmanPowSnr.add(hSourceLevelunit, pos, 4, 'MinimumWidth', 0.95*UNIT_WIDTH, ...
            'MinimumHeight', UNIT_HEIGHT, 'TopInset', TOP_INSET);
        hLmanPowSnr.update();
        %drawnow();
    end

    function layout_snr(pos)
        % Layout snr ui controls at row pos
        hLmanPowSnr.add(hSnrtxt, pos, 1, 'MinimumWidth', 32, ...
            'TopInset', 5, 'Anchor', 'NorthWest');
        hLmanPowSnr.add(hSnrDetbtn, pos, 2,'Anchor', 'Northwest',...
            'MinimumWidth', 35, 'MinimumHeight', SNRBTN_HEIGHT, 'LeftInset', -50);
        hLmanPowSnr.add(hSnredt, pos, 3, 'MinimumWidth', EDIT_WIDTH, ...
            'MinimumHeight', EDIT_HEIGHT, 'Anchor', 'Northeast');
        hLmanPowSnr.add(hSnrunit, pos, 4, 'MinimumWidth', 0.95*UNIT_WIDTH , ...
            'TopInset', TOP_INSET);
        hLmanPowSnr.update();
    end

    function toggle_SourceLevel_snr(toggle_key)
        % Makes SourceLevel and snr related ui controls visible or invisible
        % according to the toggle_key
        set(hSourceLeveltxt, 'Visible', toggle_key{1});
        set(hSourceLeveledt, 'Visible', toggle_key{2});
        set(hSourceLevelunit, 'Visible', toggle_key{3});
        
        set(hSnrtxt, 'Visible', toggle_key{4});
        set(hSnredt, 'Visible', toggle_key{5});
        set(hSnrDetbtn, 'Visible', toggle_key{6});
        set(hSnrunit, 'Visible', toggle_key{7});
    end

    function resize_gui(size_state, rszlen)
        % Resize the main GUI figure when SNR detection assistant button is
        % pressed. size_state: Is a string set to 'on' or 'off'. When set
        % to 'on' the size is increased and vice versa when set to 'off'.
        % rszlen: The resize height value
        if strcmp(size_state, 'on')
            set(hGui, 'Position', get(hGui, 'Position') + [0 -rszlen 0 rszlen]);
        elseif strcmp(size_state, 'off')
            set(hGui, 'Position', get(hGui, 'Position') + [0 rszlen 0 -rszlen]);
        end
        drawnow();
    end

    function resize_gui_tl(size_state, rszlen)
        % Resize the main GUI figure when TL calcualtion assistant button is
        % pressed. size_state: Is a string set to 'on' or 'off'. When set
        % to 'on' the size is increased and vice versa when set to 'off'.
        % rszlen: The resize height value
        if strcmp(size_state, 'on')
            set(hGui, 'Position', get(hGui, 'Position') + [0 -rszlen 0 rszlen]);
        elseif strcmp(size_state, 'off')
            set(hGui, 'Position', get(hGui, 'Position') + [0 rszlen 0 -rszlen]);
        end
        drawnow();
    end

    function populate_pdpfa()
        % Populate Pd edit box using SNR value and Pfa ( set to 1e-5 by
        % default). This function is called when user sets the SNR edit
        % box.
        SNRdB = str2double(get(hSnredt, 'String'));
        num_pulses = str2double(get(hNumPulseedt, 'String'));
        tgt_case = get(hSwerlingNumpopup, 'Value');
        pfa = str2double(get(hProbFaedt,'String'));
        pd = get_pdpfa(SNRdB, pfa, num_pulses, tgt_case);
        % Round off to two decimal places, shnidman limits Pd <= 0.99
        %pd = round(pd*100)/100;
        set(hProbDetedt, 'String', num2str(pd));
        set(hProbFaedt, 'String', num2str(pfa, 3));
        setappdata(hSnredt, 'probfa', get(hProbFaedt, 'String'));
        setappdata(hProbDetedt, 'probdet', get(hProbDetedt, 'String'));
        setappdata(hProbFaedt, 'probfa', get(hProbFaedt, 'String'));
    end

    function update_snr_Callback()
        Pfa = str2double(get(hProbFaedt, 'String'));
        Pd = str2double(get(hProbDetedt,'String'));
        N = str2double(get(hNumPulseedt, 'String'));
        
        swerling_num = get(hSwerlingNumpopup, 'Value') - 1;
        try
            SNRdB = shnidman(Pd, Pfa, N, swerling_num);
        catch errObj
            oldPfa = get(hProbFaedt, 'String');
            oldPd = get(hProbDetedt, 'String');
            oldN = get(hNumPulseedt, 'String');
            oldSwerling = get(hSwerlingNumpopup, 'Value') - 1;
            errordlg(errObj.message, 'Invalid Input', 'modal');
            set(hProbFaedt, 'String', oldPfa);
            set(hProbDetedt, 'String', oldPd);
            set(hNumPulseedt, 'String', oldN);
            %  set(hSwerlingNumpopup, 'String', oldSwerling);
            set(hSwerlingNumpopup, 'Value', oldSwerling);
            return;
        end
        % Round off to two decimal places to match with rocsnr results
        %SNRdB = round((SNRdB)*10)/10;
        set(hSnredt, 'String', num2str(SNRdB));
        setappdata(hSnredt, 'snr', get(hSnredt, 'String'));
    end

    function update_tl_Callback()
        freq = str2double(get(hFreqedt, 'String'))*Units.Freq;
        depth = str2double(get(hDepthedt,'String'))* Units.Depth;
        rng = str2double(get(hTgtRngTxedt, 'String'))* Units.Range;
        
        transloss_update = range2tl(rng,freq,depth);
        
        set(hTLedt, 'String', num2str(transloss_update));
        setappdata(hTLedt, 'transloss', get(hTLedt, 'String'));
    end

    function [pd, pfa] = get_pdpfa(SNRdB, pfa, num_pulses, tgt_case)
        % Returns probability of detection for given SNR and Pfa using the
        % builtin function rocsnr() Signal type dictionary
        signaltype = {'NonfluctuatingNoncoherent', ...
            'Swerling1', ...
            'Swerling2', ...
            'Swerling3', ...
            'Swerling4'};
        
        % Get Pd using rocsnr for giving SNR, keep Pfa max
        [pd, pfa] = rocsnr(SNRdB, 'NumPulses', num_pulses, ...
            'MaxPfa', pfa, 'MinPfa', pfa/10, 'SignalType', signaltype{tgt_case}, 'NumPoints', 1);
    end

    function switch_snrdetpanel(state)
        % Handles all the actions required when showing or hiding the
        % detection setting panel
        if strcmp(state, 'on')
            switch_lines('off');
            drawnow();
            % Resize GUI
            resize_gui('on', RSZ_HEIGHT);
            switch_lines('on');
            drawnow();
            % Display SnrDet panel
            calc_type = get(hCalcpopup, 'Value');
            
            hLmanParm.add(hSnrDetpanel, 5, 1, 'MinimumWidth', GUI_WIDTH - 20, ...
                'MinimumHeight', RSZ_HEIGHT, 'Anchor', 'Northwest', 'LeftInset', -5);
            
            set(hSnrDetpanel, 'Visible', 'on');
        elseif strcmp(state, 'off')
            switch_lines('off');
            % Hide SnrDet panel
            set(hSnrDetpanel, 'Visible', 'off');
            hLmanParm.remove(hSnrDetpanel);
            drawnow();
            % Resize GUI
            resize_gui('off', RSZ_HEIGHT);
            switch_lines('on');
            drawnow();
        end
    end

    function switch_tlpanel(state)
        % Handles all the actions required when showing or hiding the
        % TL calculation panel
        
        if strcmp(state, 'on')
            switch_lines('off');
            drawnow();
            % Resize GUI
            resize_gui_tl('on', RSZ_HEIGHT_TL);
            switch_lines('on');
            drawnow();
            % Display SnrDet panel
            hLmanParm.add(hTLpanel, 3, 1, 'MinimumWidth', GUI_WIDTH - 20, ...
                'MinimumHeight', RSZ_HEIGHT_TL, 'Anchor', 'Northwest', 'LeftInset', -5);
            set(hTLpanel, 'Visible', 'on');
        elseif strcmp(state, 'off')
            switch_lines('off');
            % Hide SnrDet panel
            set(hTLpanel, 'Visible', 'off');
            hLmanParm.remove(hTLpanel);
            drawnow();
            % Resize GUI
            resize_gui_tl('off', RSZ_HEIGHT_TL);
            switch_lines('on');
            drawnow();
        end
    end
    function switch_lines(state)
        % Switch visibility of two horizontal lines
        set(hLineWB, 'Visible', state);
        set(hLineBB, 'Visible', state);
        set(hLineWT, 'Visible', state);
        set(hLineBT, 'Visible', state);
    end

    function sonarparm = get_parm_struct()
        
        % Create sonarparm structure to hold all parameters in the GUI
        sonarparm.calc_type = get(hCalcpopup, 'Value');
        sonarparm.mode_type=get(hModepopup,'Value');
        sonarparm.rngfreq = str2double(get(hRFreqedt,'String'))*Units.RangeFreq;
        sonarparm.rngdepth = str2double(get(hRDepthedt,'String'))*Units.RangeDepth;
        sonarparm.mode_string=Mode_parms{sonarparm.mode_type};
        sonarparm.freq = str2double(get(hFreqedt,'String'))* Units.Freq;
        sonarparm.depth = str2double(get(hDepthedt,'String'))*Units.Depth;
        sonarparm.range = str2double(get(hTgtRngTxedt,'String'))*Units.Range;
        sonarparm.NoiseLevel = str2double(get(hNoiseLeveledt, 'String'));
        sonarparm.target_strength = str2double(get(hTSedt, 'String'));
        sonarparm.directivity_index=str2double(get(hDIedt,'String'));
        sonarparm.SL = str2double(get(hSourceLeveledt, 'String'));
        sonarparm.TL = str2double(get(hTLedt, 'String'));
        sonarparm.SNR = str2double(get(hSnredt, 'String'));
        sonarparm.result = guidata(hGui);
        sonarparm.Pd = str2double(get(hProbDetedt, 'String'));
        sonarparm.Pfa = str2double(get(hProbFaedt, 'String'));
        sonarparm.numpulse = str2double(get(hNumPulseedt, 'String'));
        sonarparm.swerlingcase = get(hSwerlingNumpopup, 'Value') - 1;
        
        if strcmp(get(hSnrDetbtn, 'String'), DBLRTARROW)
            sonarparm.detection = 0;    % OFF
        elseif strcmp(get(hSnrDetbtn, 'String'), DBLLTARROW)
            sonarparm.detection = 1;    % ON
        end
        if strcmp(get(hTLbtn, 'String'), DBLRTARROW)
            sonarparm.tldetection = 0;    % OFF
        elseif strcmp(get(hTLbtn, 'String'), DBLLTARROW)
            sonarparm.tldetection = 1;    % ON
        end
    end

    function render_sonarconf(pos1, pos2)
        % Renders sonar configuration group pos1 : Position of Freq in
        % controls pos2 : Position of depth.
        
        hLmanSonarConfig = siglayout.gridbaglayout(hmediumPanel, ...
            'VerticalGap', 5, 'HorizontalGap', 5);
        calc_type = get(hCalcpopup, 'Value');
        
        if calc_type ==1 % In case of Range Calculation
            hLmanSonarConfig.remove(hTLTx);
            hLmanSonarConfig.remove(hTLbtn);
            hLmanSonarConfig.remove(hTLedt);
            hLmanSonarConfig.remove(hTLunit);
            hLmanSonarConfig.add(hRFreqTx, pos1, 1:2, 'Fill', 'Horizontal', ....
                'TopInset',5, 'LeftInset', 2,'Anchor', 'NorthWest');
            hLmanSonarConfig.add(hRDepthTx, pos2, 1:2, 'Fill', 'Horizontal',.....
                'TopInset', 0, 'LeftInset', 2);
            hLmanSonarConfig.add(hRFreqedt, pos1, 3, 'MinimumWidth', EDIT_WIDTH, ...
                'MinimumHeight', 25,'Anchor', 'NorthWest');
            hLmanSonarConfig.add(hRDepthedt,pos2, 3, 'MinimumWidth', EDIT_WIDTH, ...
                'MinimumHeight', 25);
            hLmanSonarConfig.add(hRFrequnit, pos1, 4, 'MinimumWidth', UNIT_WIDTH, ...
                'MinimumHeight', UNIT_HEIGHT,'TopInset', UNIT_INSET);
            hLmanSonarConfig.add(hRDepthunit, pos2, 4, 'MinimumWidth', UNIT_WIDTH, ...
                'MinimumHeight', UNIT_HEIGHT,'TopInset', UNIT_INSET);
        elseif ~( calc_type == 1 || calc_type == 2 ) % Not Range and Transmission Loss
            hLmanSonarConfig.add(hTLTx,1, 1,'MinimumWidth', 120, 'TopInset', 5,....
                'LeftInset', 2,'Anchor', 'NorthWest');
            hLmanSonarConfig.add(hTLbtn, 1, 2,...
                'MinimumWidth', 35, 'MinimumHeight', SNRBTN_HEIGHT, 'TopInset', 0,...
                'LeftInset', -50,'Anchor', 'NorthWest');
            hLmanSonarConfig.add(hTLedt,1, 3, 'MinimumWidth', EDIT_WIDTH, ...
                'MinimumHeight', 25,'Anchor', 'NorthWest');
            hLmanSonarConfig.add(hTLunit, 1, 4, 'MinimumWidth', UNIT_WIDTH, ...
                'MinimumHeight', UNIT_HEIGHT, ...
                'TopInset', UNIT_INSET,'Anchor', 'NorthWest');
        end
        hLmanSonarConfig.HorizontalWeights = [1 0 0];
        drawnow();
    end
end

% [EOF]
