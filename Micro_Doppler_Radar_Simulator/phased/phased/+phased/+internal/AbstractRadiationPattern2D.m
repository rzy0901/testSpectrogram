classdef (Hidden) AbstractRadiationPattern2D < phased.internal.AbstractRespPattern2D
%This class is for internal use only. It may be removed in the future.

%AbstractRadiationPattern2D Class definition for 
%phased.internal.AbstractRadiationPattern2D class
%   This is an abstract class to support basic functionality for response
%   patterns.

%   Copyright 2013-2017 The MathWorks, Inc.
    
    
    properties
        %Angle - Sampling angle vector
        %   Angle is a column vector containing the angles (in degrees) at
        %   which the response pattern are sampled.
        Angle
        %SliceDir - Slicing direction
        %   SliceDir specifies along which direction the 2D response pattern
        %   is formed.  SliceDir can be one of the following: [{'El'} | 'Az'].
        SliceDir
        %CutAngle - A scalar specifying the cut angle that the pattern
        %   was taken along.
        CutAngle
    end
    
    methods
        function obj = AbstractRadiationPattern2D(varargin)
        %RespPattern2D Constructor of phased.internal.RespPattern2D class
            Angle = (-180:180).'; %#ok<*PROP>
            Freq = 1e9;
            Pattern = ones(size(Angle));
            SliceDir = 'El';
            CutAngle = 0;
            sigutils.pvparse(varargin{:});
            obj.Angle = Angle;
            obj.Freq = Freq; %store frequencies corresponding to patterns
            obj.Pattern = Pattern;
            obj.SliceDir = SliceDir;
            obj.Type = '2D Radiation Pattern';
            obj.CutAngle = CutAngle;
        end
        
        function varargout = polar(obj,varargin)
        %POLAR Plot the array response pattern using polar coordinates
        %   POLAR(Hresp) plots the array response pattern Hresp using polar
        %   coordinates.
        %
        %   POLAR(...,'Units',UNITS) specifies the units of the radial axis as
        %   a string ['mag' | 'pow' | {'db'}].
        %
        %   POLAR(...,'NormalizeResp',NFLAG) specifies whether to normalize the
        %   array response.  NFLAG is specified as a logical value [true |
        %   {false}].
        %
        %   H = POLAR(...) returns a handle to the plotted line objects in H.
        %
        %   Example:
        %    % Construct and plot a 2D isotropic array response pattern
        %    hresp = phased.internal.RespPattern2D;
        %    polar(hresp)
                    
                        
            plotoption = getPlotOption(obj,varargin{:});
            plotoption.PlotType = 'polar';
            if isempty(plotoption.Title)
                plotoption.Title = genPlotTitle(obj);
            end
            [angles, response] = chkvectorinput(obj);
            response = privRespPattern(plotoption,response);
            orientation = plotoption.BroadsideOrientation;
            angleMat = repmat(angles,1,size(response,2)); %condition angle data for input to HGPOLAR
            
            % Get existing legend entries to be appended
            [hExistingLines,existingLabels] = obj.getLegendEntriesToAppend;
        
            h = hgpolar(strcmpi(plotoption.Units,'db'), deg2rad(angleMat'),response',orientation);
       
            annotatePlot(obj,plotoption);
            
            % Hide data cursor in polar mode
            hfig = get(get(h(1),'Parent'),'Parent');
            if ishghandle(hfig) && strcmp(get(hfig,'Type'),'figure')
                hdcm = datacursormode(hfig);
                set(hdcm,'Enable','off');
                set(hdcm,'UpdateFcn',@(dummy,event_obj)getString(...
                    message('phased:apps:arrayapp:DataCursorNotSupported')));
                hdcb = findall(hfig,'Tag','Exploration.DataCursor');
                set(hdcb,'Visible','off');
            end
            
            % Create a legend and append existing legend entries to it if
            % needed.
            if isempty(hExistingLines)
                addLegend(obj,h);
            else
                % If appending previous legend entries, pass these
                % entries to ADDLEGEND.
                addLegend(obj,h,hExistingLines,existingLabels);
            end
 
            
            if nargout == 1,
                varargout{1} = h;
            end
        end
        
    end
        
    methods (Access = protected)
        function g = getSliceGroup(obj) 
            % either frequency or cut angle is a vector
            if isscalar(obj.CutAngle)
                g = getSliceGroup@phased.internal.AbstractRespPattern2D(obj);
            else
                if strncmpi(obj.SliceDir,'a',1)
                    g = 'Elevation';
                else
                    g = 'Azimuth';
                end
            end
        end
        
        function angles = getPlotAngles(obj)
            angles = obj.Angle;
        end
        
        function validatePattern(obj,value)
            validateattributes(value,{'numeric'},...
                {'nonnegative',...
                'nrows',numel(obj.Angle)},...
                sprintf('%s.Pattern',class(obj)),'Pattern');   
        end
        
        function xlbl = getXLabel(obj,plotoptionobj)
            %getXLabel Return x-axis label of the plot
            if strcmpi(plotoptionobj.PlotType,'polar')
                xlbl = sprintf('%s, Broadside at %.2f degrees',...
                    getRespLabel(plotoptionobj),plotoptionobj.BroadsideOrientation);
            else
                xlblSuffix = 'Angle (degrees)';
                switch lower(obj.SliceDir)
                    case 'az',
                        xlbl = sprintf('Azimuth %s',xlblSuffix);
                    case 'el',
                        xlbl = sprintf('Elevation %s',xlblSuffix);
                end
            end
            
        end
        
        function plotTitle = genPlotTitle(obj)
            %GENPLOTTITLE Generate a title string for ResponsePattern2D
            %plots.
            
            % either cut angle is scalar or freq is scalar
            if isscalar(obj.CutAngle)
                if strncmpi(obj.SliceDir,'a',1)
                    cut_title = 'Azimuth';
                    cutAngleDir = 'elevation';
                else
                    cut_title = 'Elevation';
                    cutAngleDir = 'azimuth';
                end

                % LaTeX degree symbol for labeling the cut angle
                degSym = '{^\circ}';

                % Built a plot title string containing the cut type and angle
                plotTitle = sprintf('%s Cut (%s angle = %3.1f%s)',...
                    cut_title,cutAngleDir,obj.CutAngle,degSym);
            else  %multiple angle cuts
                if strncmpi(obj.SliceDir,'a',1)
                    cut_title = 'Azimuth';
                else
                    cut_title = 'Elevation';
                end
                % Convert frequency values from Hz units to engineering units
                precision = 5; %number of digits to use after the decimal point.
                [convertedVals,unitPrefix] = convert2engstrs(obj.Freq,precision);
                
                plotTitle = sprintf('%s Cut (frequency = %s %sHz)',...
                    cut_title,convertedVals,unitPrefix);
            end
        end
        
        function addLegend(obj,hlines,varargin)
            %ADDLEGEND Add a legend to label each pattern's frequency
            %   ADDLEGEND(obj,hlines) adds a legend to the current axes
            %   with frequency labels corresponding to each line in handle
            %   array HLINES. If a single pattern is being plotted to an
            %   axes that does not contain a pattern previously labelled in
            %   a legend by AddLegend or saved to the AppData, the legend
            %   label and line handle are saved to the axes' AppData.  The
            %   AppData entries are used to populate a legend if ADDLEGEND
            %   is called again with the same current axes and the axes'
            %   hold state is set to 'on'.
            %
            %   ADDLEGEND(obj,hlines,hExistingLines,existingLabels) adds a
            %   legend with frequency labels corresponding to each line in
            %   HLINES.  Existing line handles in array HEXISTINGLINES and
            %   existing legend labels in EXISTINGLABELS are included in
            %   the legend.
            
            % either cut angle is scalar or freq is scalar
            if isscalar(obj.CutAngle)
                addLegend@phased.internal.AbstractRespPattern2D(obj,hlines,varargin{:});
            else  % multiple angle cuts
            
                if strncmpi(obj.SliceDir,'a',1)
                    cutAngleDir = 'elevation';
                else
                    cutAngleDir = 'azimuth';
                end
                
                % Convert angel values to labels
                labels = arrayfun(@(x)sprintf('%3.1f deg %s',x,cutAngleDir),...
                    obj.CutAngle,'UniformOutput',false);

                % Determine if existing legend entries are to be appended
                appendFlag = nargin>2;

                % Determine if a single line plot is being created in an axes
                % without an existing legend.
                isNewSingleLinePlot =  ~appendFlag && isscalar(hlines);

                % Append each frequency value with its corresponding unit
                % prefix.  For multi-frequency plots, each line is given a
                % descriptive Tag.
                if ~isNewSingleLinePlot
                    for n=1:numel(labels)
                        set(hlines(n),'tag',...
                            sprintf('Angle%s',num2str(obj.CutAngle(n))));
                    end
                end


                % This section is executed when the figure's previous legend
                % entries are going to be included in the new legend.  Previous
                % legend labels and line handles are appended to the newly
                % created ones in the LABELS and LINES cell arrays.
                if appendFlag
                    % Get existing line handles and legend labels
                    hExistingLines = varargin{1};
                    existingLabels = varargin{2};

                    % Append legend labels and line handles to any existing ones
                    labels = [existingLabels(:); labels(:)];
                    hlines = [hExistingLines ; hlines];
                end

                % If the current axes does not contain an existing legend and a
                % single frequency is being plotted, save the line handle and
                % label to AppData for possible later use.  Otherwise, create a
                % legend.
                if isNewSingleLinePlot
                     hLeg = legend(hlines,labels,'AutoUpdate','off');
                     set(hLeg,'Visible','off');
                else
                    legend(hlines,labels,'location','BestOutside','AutoUpdate','off');
                end
            end
        end
        
    end
    
    methods (Access = protected)
        function sortedList = getSortedPropDispList(this)  %#ok<MANU>
            % Get the sorted list of the properties to be displayed.
            sortedList = {'Type','SliceDir','Angle','Pattern','Freq','CutAngle'};
        end
        
        function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
            plotoptionobj = phased.internal.RespPattern2DPlotOption(varargin{:});
        end
        
    end


    methods
        
        function obj = set.Angle(obj,value)
            validateattributes(value,{'double'},{'nonempty','real','column'},...
                sprintf('%s.Angle',class(obj)),'Angle');
            obj.Angle = value;
        end
        
        function obj = set.SliceDir(obj,value)
            value = validatestring(value,{'El','Az'},...
                sprintf('%s.SliceDir',class(obj)),'SliceDir');
            obj.SliceDir = value;
        end
                
    end
    
end

% Local utility functions

function hpol = hgpolar(dbflag, varargin)
%   POLAR(...,'BroadsideOrientation', ANGLE) specifies the angle
%   orientation of the broadside direction in radians measured
%   clockwise with respect to the horizontal.  The default value is 0.
%
%   POLAR(...,'MaxDynamicRange', MAXDYNRGE) specifies the maximum
%   dynamic range in dB.  The default value is 100 dB. This value is
%   ignored when 'Units' is not 'db'.
%
%   BroadsideOrientation and MaxDynamicRange are currently hidden
% Read input arguments

[cax,args,nargs] = axescheck(varargin{:});  %#ok<ASGLU>
[theta,rho] = deal(args{1:2});
line_style = 'auto';
orientation = args{3};

% Get hold state
cax = newplot(cax);
next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% Save the original data
saved_data = getappdata(cax,'saved_data');
if ~isempty(saved_data)
    saved_data = [saved_data; rho];
else
    saved_data = {rho};
end

% Get x-axis text color so grid is in same color
tc = get(cax,'XColor');
ls = get(cax,'GridLineStyle');

% Hold to current Text defaults, reset to Axes' font
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   10, ...  % overwrite FontSize to match annotation
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% Angle offset in degrees
theta_offset = orientation;

% Flag that indicates whether to redraw data on hold or not
redraw_hold = false;
% Transform current data when not in hold state
if ~hold_state
    % clip current data using range limit
    if dbflag
        rho = phased.internal.RespPattern2D.limitDynamicdBRange(rho,50);
    end
    % find the scale of the current data
    minrho = min(rho(:));
    maxrho = max(rho(:));
    % offset current data to make sure all radial data is positive
    if dbflag
        rmin = sqrt(eps);
        rho = rho - minrho + rmin;
    else
        minrho = 0;
        rmin = 0;
    end
end

% If hold is ON, read data from existing plot
if hold_state
    % get data from cax, including original data and clipped data
    retrieve_data_hold = getappdata(cax,'saved_data');
    minrho_hold = getappdata(cax, 'RadialOffset');
    maxrho_hold = getappdata(cax, 'MaxPatternValue');
    rmin_hold = getappdata(cax, 'MinRadialValue');
    theta_offset_hold = getappdata(cax, 'AngleOffset');
    % only manipulate hold polar data if the hold plot has been created by
    % this program
    if ( ~isempty(minrho_hold) && ~isempty(maxrho_hold) && ...
            ~isempty(rmin_hold) && ~isempty(theta_offset_hold) )
        % Delete the existing polar grid and background since it is going
        % to be redrawn anyway Clear all children in axes except for the
        % actual data lines
        hold_lines = get(cax,'Children');
        hh_hold = allchild(cax);
        delete(hh_hold(~ismember(hh_hold, hold_lines)));
        % read the handles to the children of the current axis
        % make sure we only get objects with visible handles on
        fshh = get(0, 'ShowHiddenHandles');
        set(0, 'ShowHiddenHandles', 'off');
        currlines = get(cax, 'children');
        set(0, 'ShowHiddenHandles', fshh);
        
        NumCurrlines = size(currlines,1);
        theta_hold = cell(NumCurrlines,1);
        rho_hold = cell(NumCurrlines,1);
        
        % loop through all children
        for ii = 1:NumCurrlines
            xx_hold = get(currlines(ii), 'XData')';
            yy_hold = get(currlines(ii), 'YData')';
            [theta_hold{ii}, rho_hold{ii}] = cart2pol(xx_hold, yy_hold);
            % transform angle to absolute value (no angle offset)
            theta_hold{ii} = -theta_hold{ii}-theta_offset_hold;
            % bring data to original scale
            rho_hold{ii} = offsetrad2absrad(rho_hold{ii}, minrho_hold, ...
                rmin_hold);
        end
        
        % max rho for hold_state is true
        maxrho = max(rho(:));
        
        % clip the current and retrieved hold data
        if dbflag
            if maxrho-50 >= minrho_hold
                rho = phased.internal.RespPattern2D.limitDynamicdBRange(rho,...
                    maxrho-minrho_hold);
            else
                rho = phased.internal.RespPattern2D.limitDynamicdBRange(rho,50);
                retrieve_data_hold =...
                    phased.internal.RespPattern2D.limitDynamicdBRange(...
                    retrieve_data_hold,...
                    max(cell2mat(cellfun(@max,retrieve_data_hold,'UniformOutput',false)))-(max(rho(:))-50));
                rho_hold = cellfun(@(x)x.',retrieve_data_hold,'UniformOutput',false);
            end
        end
        
        % find the scale of the current data
        minrho = min(rho(:));
        maxrho = max(rho(:));
        
        % adjust min and max radii for the current plot if change of scale
        % is needed
        if ( ( minrho < 0 ) || ( minrho_hold < 0 ) || ...
                ( abs(maxrho_hold - maxrho) > eps ) || ...
                ( abs(minrho_hold - minrho) > eps ) || ...
                ( abs(theta_offset - theta_offset_hold) > eps ) )
            % initialize redrawing flag
            redraw_hold = false;
            
            if ( minrho >= minrho_hold )
                % No need to redraw existing lines
                minrho = minrho_hold;
            else
                % Need to redraw existing lines
                redraw_hold = true;
            end
            
            if ( maxrho_hold >= maxrho )
                % No need to redraw existing lines
                maxrho = maxrho_hold;
            else
                % Need to redraw existing lines
                redraw_hold = true;
            end
            
            if ( abs(theta_offset - theta_offset_hold) > eps )
                redraw_hold = true;
            end
            
            % offset current data by minrho
            rmin = sqrt(eps);
            rho = absrad2offsetrad(rho, minrho, rmin);
            
            % offset hold data by the same amount
            for ii = 1:NumCurrlines
                rho_hold{ii} = absrad2offsetrad(rho_hold{ii}, minrho, rmin);
            end
        else
            rmin = 0;
        end
        
    else % hold plot has not been created using this polar
        error(message('phased:phased:internal:RespPattern2D:InvalidHoldPlot'));
    end
    
end % (if hold_state)

isFlatdBResponse = flatdbresponse(minrho,maxrho,dbflag);

% make a radial grid

% Set 'hold' to 'on' only if the axes' hold state is currently 'off'. This
% ensures that an existing PlotHoldStyle property value is not changed when
% the plot's hold status is set to 'ALL'.
if ~hold_state
    hold(cax,'on');
end

% maximum radius value including the radial offset
maxrho_offset = maxrho - minrho + rmin;
if maxrho_offset == rmin
    maxrho_offset = 1;
end
hhh=line([-maxrho_offset -maxrho_offset maxrho_offset maxrho_offset], ...
    [-maxrho_offset maxrho_offset maxrho_offset -maxrho_offset], ...
    'parent',cax);
%   tighten the axis to maximize plot area
axis tight;

set(cax,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatioMode','auto')
v = [get(cax,'xlim') get(cax,'ylim')];

% Find the inter-tick distance
rinc = get(cax,'YTick');
rinc = rinc(2) - rinc(1);
delete(hhh);

% check radial limits and ticks
rmax = v(4);

% define a circle
th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
inds = 1:(length(th)-1)/4:length(th);
xunit(inds(2:2:4)) = zeros(2,1);
yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
if ~ischar(get(cax,'color')),
    p = patch('XData',xunit*rmax,'YData',yunit*rmax, ...
        'EdgeColor',tc,'FaceColor',get(cax,'color'),...
        'HandleVisibility','off','parent',cax);
    set(p,'ZData',-1*ones(size(get(p, 'XData'))));
end

% draw radial circles

if isFlatdBResponse
    rtickval = minrho;
else
    rtickval = absrad2offsetrad(minrho-mod(minrho,rinc)+rinc:rinc:maxrho, minrho, rmin);
end

% set angle for the circular tick marks
c82 = cos(82*pi/180 - theta_offset);
s82 = sin(82*pi/180 - theta_offset);

% determine the angle of the offset angle in order to change the
% location of the tick marks
if ( sign(c82) <= 0 )
    horizontalalignment = 'right';
else
    horizontalalignment = 'left';
end

if ( sign(s82) <= 0 )
    verticalalignment = 'top';
else
    verticalalignment = 'bottom';
end


% draw radial grid
if isFlatdBResponse
    text((maxrho_offset+rinc/20)*c82,(maxrho_offset+rinc/20)*s82, ...
        ['  ' num2str(rtickval)], ...
        'verticalalignment',verticalalignment,...
        'horizontalalignment',horizontalalignment,...
        'handlevisibility','off','parent',cax,'Tag',...
        sprintf('CircleTick%s',num2str(rtickval)))
else
    for ii = rtickval
        line(xunit*ii,yunit*ii,'linestyle',ls,'color',tc,'linewidth',1,...
            'handlevisibility','off','parent',cax);
        ticklabel = num2str(roundtick(offsetrad2absrad(ii, minrho, rmin)));
        text((ii+rinc/20)*c82,(ii+rinc/20)*s82, ...
            sprintf('  %s',ticklabel), ...
            'verticalalignment',verticalalignment,...
            'horizontalalignment',horizontalalignment,...
            'handlevisibility','off','parent',cax,'Tag',...
            sprintf('CircleTick%s',ticklabel))
    end
end
    

% plot spokes
th = (1:6)*2*pi/12;

% offset selected angle by theta_offset
cst = cos(th-theta_offset); snt = sin(th-theta_offset);
cs = [-cst; cst];
sn = [-snt; snt];
line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
    'handlevisibility','off','parent',cax)

% tick clockwise setting
% -1 for counter clockwise, 1 for clockwise
clockwise_dir = -1;

% annotate spokes in degrees
rt = 1.1*rmax;

for ii = 1:length(th)
    if ii == length(th)
        loc = int2str(ii*30);
    else
        loc = int2str(-clockwise_dir*ii*30);
    end
    text(rt*cst(ii),rt*snt(ii), loc,...
        'horizontalalignment','center',...
        'handlevisibility','off','parent',cax);
    if ii == length(th)
        loc = int2str(0);
    else
        loc = int2str(clockwise_dir*(180-ii*30));
    end
    text(-rt*cst(ii),-rt*snt(ii),loc,'horizontalalignment','center',...
        'handlevisibility','off','parent',cax)
end

% set view to 2-D
view(cax,2);
% set axis limits
axis(cax,rmax*[-1 1 -1.15 1.15]);

% Save user data for use when hold is on
% appdata is [minrho maxrho]
setappdata(cax, 'saved_data', saved_data);
setappdata(cax, 'RadialOffset', minrho);
setappdata(cax, 'MaxPatternValue', maxrho);
setappdata(cax, 'MinRadialValue', rmin);
setappdata(cax, 'AngleOffset', theta_offset);

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% change sign of angle for angle to be measured clockwise
theta = -clockwise_dir*theta-theta_offset;
% transform data to Cartesian coordinates.
if isFlatdBResponse
    [xx, yy] = pol2cart(theta,maxrho_offset);
else
    [xx, yy] = pol2cart(theta, rho);
end

% plot data on top of grid
if strcmp(line_style,'auto')

    q = plot(xx',yy','parent',cax);  
else
    q = plot(xx',yy',line_style,'parent',cax);  
end

% Tagging plot handles
set(q,'Tag','2D polar plot');

if nargout == 1
    hpol = q;
end
% Redraw the plots on hold
if redraw_hold
    for jj = 1:NumCurrlines
        % change sign of angle for angle to be measured clockwise
        theta_hold{jj} = -theta_hold{jj}-theta_offset;
        % transform data to Cartesian coordinates.
        theta_hold_current = theta_hold{jj};
        rho_hold_current = rho_hold{jj};
        [xx_hold, yy_hold] = pol2cart(theta_hold_current(:), rho_hold_current(:));
        % change current lines data
        set(currlines(jj), 'xdata', xx_hold, 'ydata', yy_hold);
    end
end

if ~hold_state
    set(cax,'dataaspectratio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
end
set(get(cax,'xlabel'),'visible','on')
set(get(cax,'ylabel'),'visible','on')
set(get(cax,'title'),'visible','on')

if ~isempty(q) && ~isdeployed
    makemcode('RegisterHandle',cax,'IgnoreHandle',q,'FunctionName','polar');
end

if hold_state
    hold(cax,'on');
else
    hold(cax,'off');
end

end % hgpolar


function rho_offset = absrad2offsetrad(rho_absolute, minrho, rmin)
% Offset radial data
rho_offset = rho_absolute - minrho + rmin;
end

function rho_absolute = offsetrad2absrad(rho_offset, minrho, rmin)
% Bring offset data back to absolute scale
rho_absolute = rho_offset + minrho - rmin;
end

function rt = roundtick(rtickval)
% Round the tick labels to 2 decimal points
rt = round(rtickval*100)/100;
end

function flag = flatdbresponse(minresp,maxresp,dbflag)
flag = (abs(minresp-maxresp) < .01) && dbflag;
end


% [EOF]
