function hpol = setpolaraxes(this)
%   POLAR(...,'BroadsideOrientation', ANGLE) specifies the angle
%   orientation of the broadside direction in radians measured
%   clockwise with respect to the horizontal.  The default value is 0.
%
%   POLAR(...,'MaxDynamicRange', MAXDYNRGE) specifies the maximum
%   dynamic range in dB.  The default value is 100 dB. This value is
%   ignored when 'Units' is not 'db'.
%
%   BroadsideOrientation and MaxDynamicRange are currently hidden

%   Copyright 2011 The MathWorks, Inc.


% Read input arguments
% [cax,args,nargs] = axescheck(varargin{:}); %#ok<NASGU>
% [theta,rho] = deal(args{1:2});
line_style = 'auto';
% orientation = args{3};
orientation = 0;

% Get hold state
%cax = newplot(cax);
cax = this.Axes;
next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% Get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold to current Text defaults, reset to Axes' font
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% Angle offset in degrees
theta_offset = orientation;

% Flag that indicates whether to redraw data on hold or not
redraw_hold = false;
%       Transform current data when not in hold state
  minrho = 0;
      maxrho = 1;
      rmin=0;
if ~hold_state
    % find the scale of the current data
%     minrho = min(rho(:));
%     maxrho = max(rho(:));
      
      minrho = 0;
      maxrho = 1;
% offset current data to make sure all radial data is positive
    if ( minrho < 0 )
        rmin = sqrt(eps);
      %  rho = rho+ abs(minrho)+ rmin;
    else
        minrho = 0;
        rmin = 0;
    end
end

% If hold is ON, read data from existing plot
if  hold_state
    minrho_hold = getappdata(cax, 'RadialOffset');
    maxrho_hold = getappdata(cax, 'MaxPatternValue');
    rmin_hold = getappdata(cax, 'MinRadialValue');
    theta_offset_hold = getappdata(cax, 'AngleOffset');
    % only manipulate hold polar data if the hold plot has been created by this program
    if ( ~isempty(minrho_hold) && ~isempty(maxrho_hold) && ...
            ~isempty(rmin_hold) && ~isempty(theta_offset_hold) )
        % Delete the existing polar grid and background since it is going to be redrawn anyway
        % Clear all children in axes except for the actual data lines
        hold_lines = get(cax,'Children');
        hh_hold = allchild(cax);
        delete(hh_hold(~ismember(hh_hold, hold_lines)));
        % read the handles to the children of the current axis
        % make sure we only get objects with visible handles on
        fshh = get(0, 'ShowHiddenHandles');
        set(0, 'ShowHiddenHandles', 'off');
        currlines = get(cax, 'children');
        set(0, 'ShowHiddenHandles', fshh);
        
        theta_hold = cell(size(currlines,1));
        rho_hold = cell(size(currlines,1));
        
        % loop through all children
        for ii = 1:size(currlines,1)
            xx_hold = get(currlines(ii), 'XData')';
            yy_hold = get(currlines(ii), 'YData')';
            [theta_hold{ii}, rho_hold{ii}] = cart2pol(xx_hold, yy_hold);
            % transform angle to absolute value (no angle offset)
            theta_hold{ii} = -theta_hold{ii}-theta_offset_hold;
            % bring data to original scale
            rho_hold{ii} = offsetrad2absrad(rho_hold{ii}, minrho_hold, rmin_hold);
        end
        
        % find the scale of the current data
%         minrho = min(rho(:));
%         maxrho = max(rho(:));
    
        % adjust min and max radii for the current plot if change of scale is needed
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
%              rho = absrad2offsetrad(rho, minrho, rmin);
            
            % offset hold data by the same amount
            for ii = 1:size(currlines,1)
                rho_hold{ii} = absrad2offsetrad(rho_hold{ii}, minrho, rmin);
            end
        else
            rmin = 0;
        end
        
    end
    
end % (if hold_state)

% make a radial grid
hold(cax,'on');
% maximum radius value including the radial offset
maxrho_offset = maxrho + abs(minrho) + rmin;
if maxrho_offset == 0
    maxrho_offset = 1;
end
%   LINE(X,Y) adds the line in vectors X and Y to the current axes. 
%   If X and Y are matrices the same size, one line per column is added.
%   LINE(X,Y,Z) creates lines in 3-D coordinates.
hhh=line([-maxrho_offset -maxrho_offset maxrho_offset maxrho_offset], ...
    [-maxrho_offset maxrho_offset maxrho_offset -maxrho_offset], ...
    'parent',cax);
%   tighten the axis to maximize plot area
% axis tight;

set(cax,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatioMode','auto')
%v = [get(cax,'xlim') get(cax,'ylim')];
v = [-4 4 -4 4];
% Find the inter-tick distance
rinc = get(cax,'YTick');
% rinc = rinc(2) - rinc(1);
rinc=0.2;
delete(hhh);

% check radial limits and ticks
% rmax = v(4); % original setup
rmax = 1;
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
    %   PATCH(X,Y,C) adds the "patch" or filled 2-D polygon defined by
%   vectors X and Y to the current axes. If X and Y are matrices of
%   the same size, one polygon ("face") per column is added. C
%   specifies the color of the face(s) ("flat" coloring), or the 
%   vertices ("interpolated" coloring), for which bilinear interpolation
%   is used to determine the interior color of the polygon.
    p = patch('XData',xunit*rmax,'YData',yunit*rmax, ...
        'EdgeColor',tc,'FaceColor',get(cax,'color'),...
        'HandleVisibility','off','parent',cax);
    set(p,'ZData',-1*ones(size(get(p, 'XData'))));
end

% draw radial circles

rtickval = absrad2offsetrad(minrho-mod(minrho,rinc)+rinc:rinc:maxrho, minrho, rmin);

% set angle for the circular tick marks
c82 = cos(82*pi/180 - theta_offset);
s82 = sin(82*pi/180 - theta_offset);

% draw radial grid
for ii = rtickval
    line(xunit*ii,yunit*ii,'LineStyle',ls,'color',tc,'LineWidth',1,...
        'HandleVisibility','off','parent',cax);
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
    
    text((ii+rinc/20)*c82,(ii+rinc/20)*s82, ...
        ['  ' num2str(roundtick(offsetrad2absrad(ii, minrho, rmin)))], ...
        'VerticalAlignment',verticalalignment,...
        'HorizontalAlignment',horizontalalignment,...
        'HandleVisibility','off','parent',cax)
end

% plot spokes
th = (1:6)*2*pi/12;

% offset selected angle by theta_offset
cst = cos(th-theta_offset); snt = sin(th-theta_offset);
cs = [-cst; cst];
sn = [-snt; snt];
line(rmax*cs,rmax*sn,'LineStyle',ls,'color',tc,'LineWidth',1,...
    'HandleVisibility','off','parent',cax)

% annotate spokes in degrees
rt = 1.1*rmax;

for ii = 1:length(th)
    if ii == length(th)
        loc = int2str(ii*30);
    else
        loc = int2str(-ii*30);
    end
    text(rt*cst(ii),-rt*snt(ii), loc,...
        'HorizontalAlignment','center',...
        'HandleVisibility','off','parent',cax);
    if ii == length(th)
        loc = int2str(0);
    else
        loc = int2str(180-ii*30);
    end
    text(-rt*cst(ii),rt*snt(ii),loc,'HorizontalAlignment','center',...
        'HandleVisibility','off','parent',cax)
end

% set view to 2-D
% view(cax,2);
% set axis limits
 axis(cax,rmax*[-1 1 -1.15 1.15]);
% axis off;
% Save user data for use when hold is on
% appdata is [minrho maxrho]
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
% theta = -theta-theta_offset;
% transform data to Cartesian coordinates.
% if all(rho == 0)
%     [xx, yy] = pol2cart(theta,1);
% else
%     [xx, yy] = pol2cart(theta, rho');
% end

% plot data on top of grid
% if strcmp(line_style,'auto')
%    set(this.Lines, 'XData', xx', 'YData',yy','parent',cax);
% else
%    set(this.Lines, 'XData', xx', 'YData',yy','LineStyle', line_style,'parent',cax);
% end
% if strcmp(line_style,'auto')
%     q = plot(xx',yy','parent',cax);
% else
%     q = plot(xx',yy',line_style,'parent',cax);
% end
% 
% %Tagging plot handles
% set(q,'Tag','2D polar plot');
% 
% if nargout == 1
%     hpol = q;
% end
%       Redraw the plots on hold
if redraw_hold
    for jj = 1:size(currlines,1)
        % change sign of angle for angle to be measured clockwise
        theta_hold{jj} = -theta_hold{jj}-theta_offset;
        % transform data to Cartesian coordinates.
        [xx_hold, yy_hold] = pol2cart(theta_hold{jj}, rho_hold{jj});
        % change current lines data
        set(currlines(jj), 'XData', xx_hold, 'YData', yy_hold);
    end
end

if ~hold_state
    set(cax,'DataAspectRatio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
end
set(get(cax,'xlabel'),'visible','on')
set(get(cax,'ylabel'),'visible','on')
set(get(cax,'title'),'visible','on')

% if ~isempty(q) && ~isdeployed
%     makemcode('RegisterHandle',cax,'IgnoreHandle',q,'FunctionName','polar');
% end
end % hgpolar

% function varargout = plot(obj,varargin)
%     %PLOT     Plot the response pattern
%     %   plot(Hresp) plots the response pattern Hresp. 
%     %
%     %   plot(Hresp, 'Units', UNIT) plots the response pattern using the
%     %   unit specified in UNIT. UNIT can be any of the following:
%     %   ['mag' | 'power' | {'db'}].
%     %
%     %   plot(..., 'NormalizeResp', NFLAG) plots the normalized response
%     %   pattern if NFLAG is true. NFLAG can be any of the following:
%     %   [true | {false}].
%     %
%     %   plot(..., 'Title', TITLE) uses TITLE as the title of the resulting
%     %   plot.
%     %
%     %   plot returns an handle to the current plot object.
%     %
%     %   Example:
%     %       % Construct a 2D isotropic array response pattern and then plot
%     %       % it.
%     %       hresp = phased.internal.RespPattern2D;
%     %       plot(hresp)
%         plotoption = getPlotOption(obj,varargin{:});
%         plotoption.PlotType = 'plot';
%         if isempty(plotoption.Title)
%             if strncmpi(obj.SliceDir,'a',1)
%                 cut_title = 'Azimuth';
%             else
%                 cut_title = 'Elevation';
%             end
%             plotoption.Title = sprintf('%s Cut',cut_title);
%         end
%         [angles, response] = chkvectorinput(obj);
%         response = privRespPattern(plotoption,response);
%         h = plot(angles,response);
%         annotatePlot(obj,plotoption);
%         
%         limitDynamicPlotRange(plotoption,h);
%         
%         if nargout == 1,
%             varargout{1} = h;
%         end
%     end
function rho_offset = absrad2offsetrad(rho_absolute, minrho, rmin)
% Offset radial data
rho_offset = rho_absolute + abs(minrho) + rmin;
end


function rho_absolute = offsetrad2absrad(rho_offset, minrho, rmin)
% Bring offset data back to absolute scale
rho_absolute = rho_offset - abs(minrho) - rmin;
end

function rt = roundtick(rtickval)
% Round the tick labels to 2 decimal points
rt = round(rtickval*100)/100;
end

