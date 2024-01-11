classdef (Hidden, Sealed) RangeAnglePattern < phased.internal.AbstractRespPattern3D
% This class if for internal use only. It may be removed in the future.

%RangeAnglePattern   Range angle response pattern
%   Hresp =
%   phased.internal.RangeAnglePattern('PropertyName',PropertyValue,...)
%   returns a range angle response pattern object Hresp. See properties
%   list below for valid PropertyNames.
%
%   RangeAnglePattern methods:
%       plot      - Plot the range angle response pattern.
%
%   RangeAnglePattern properties:
%       Type      - 'Range Angle Response Pattern'.  This is read-only
%       Range     - Range angle pattern sampling ranges (meter)
%       Angle     - Range angle pattern sampling angles (deg)
%       PRF       - Pulse Repetition Frequency
%       Pattern   - Magnitude range angle pattern
%
%   Example:
%       % Construct a range angle response pattern and then plot it.
%       hresp = phased.internal.RangeAnglePattern;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern2D,
%   phased.internal.RespPattern3D.

%   Copyright 2018 The MathWorks, Inc.

properties
    %Range   - Range
    %   Range is a row vector containing the ranges (in meters) at which
    %   the response pattern are sampled. Default value of Range is 0:99.
    Range
    %Angle   - Angle
    %   Angle is a row vector containing the azimuth angles (in degrees) at
    %   which the response pattern are sampled. Default value of Angle is
    %   -90:90.
    Angle
    %Pattern - Response pattern
    %   Pattern is a matrix containing the samples of the magnitude
    %   response pattern at each corresponding (Range,Angle) pair. The
    %   number of rows of Pattern must match the number of elements in
    %   Angle and the number of columns of Pattern must match the number
    %   of elements in Range. Default value of Pattern is all ones.
    Pattern
end

methods
    function obj = RangeAnglePattern(varargin)
    %RangeAnglePattern Constructor of phased.internal.RangeAnglePattern class
        Range = 0:99; %#ok<*PROP>
        Angle = -90:90;
        Pattern = ones(numel(Range),numel(Angle));
        sigutils.pvparse(varargin{:});
        obj.Range = Range;
        obj.Angle = Angle;
        obj.Pattern = Pattern;
        obj.Type = 'Range Angle Response Pattern';
    end

end

methods (Access = protected)
    function sortedList = getSortedPropDispList(this)  %#ok<MANU>
        % Get the sorted list of the properties to be displayed. 
        sortedList = {'Type','Range','Angle','Pattern'};
    end
    
    function ylbl = getYLabel(obj,~)  %#ok<INUSD>
    %getYLabel Return y-axis label of the plot
        ylbl = 'Range (meters)';
    end
    
    function xlbl = getXLabel(obj,plotoptionobj)   %#ok<INUSD>
    %getXLabel Return x-axis label of the plot
        xlbl = 'Angle (degrees)';
    end

    function y = getYData(obj,plotoptionobj)
    %getYData Get y-axis data
        if strcmp(plotoptionobj.Style,'polar')
            y = obj.Range.'*cosd(obj.Angle);
        else
            y = obj.Range;
        end
    end
    
    function x = getXData(obj,plotoptionobj) 
    %getYData Get x-axis data
        if strcmp(plotoptionobj.Style,'polar')
            x = obj.Range.'*sind(obj.Angle);
        else
            x = obj.Angle;
        end
    end
    
    function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
        plotoptionobj = phased.internal.RangeAnglePatternPlotOption(varargin{:});
    end
    
    function h = privplot(obj,plotoption)  
        %PRIVPLOT  Choose customized default plot scheme
        x = getXData(obj,plotoption);
        y = getYData(obj,plotoption);
        z = privRespPattern(plotoption,obj.Pattern);
        if strcmp(plotoption.Style,'polar')
            h = surf(x,y,z,'EdgeColor','none','FaceColor','interp','tag','polarsurf');
            view(0,90)
            %sax = get(h,'Parent');
            %axis equal;
            axis off;
            set(get(h,'Parent'),'ActivePositionProperty','position','DataAspect',[1 1 1]);
            %axis tight;
            %title(plotoption.Title);
            
            %[angmin,angmax] = bounds(obj.Angle);
            %[rmin,rmax] = bounds(obj.Range);
            %pax = polaraxes(get(sax,'Parent'),'ThetaZeroLocation','top','ThetaDir','clockwise',...
            %    'ThetaLim',[angmin angmax],'ThetaTick',linspace(angmin,angmax,5),...
            %    'RLim',[rmin rmax],'Color','none','Tag','overlappolarax');
            %linkprop([sax pax],'InnerPosition');
        else
            h = imagesc(x,y,z);
            % image y axis starts from top, so we need to correct this.
            set(get(h,'Parent'),'YDir','normal');
        end
    end
    
    function annotatePlot(obj,plotoption)
    %ANNOTATEPLOT Annotate response pattern plots
        if strcmp(plotoption.Style,'polar')
            hcbar = colorbar('Location','SouthOutside');
            hcbar.Position(2) = 0.1;
            % hcbar.Position = [0.05 0.1 0.3 0.05];
            [fsize,lblColor] = phased.internal.AbstractRespPattern.setAnnotationSizeColor;
            xlabel(hcbar,getRespLabel(plotoption),'fontsize',fsize,'Color',lblColor);
            % pax = findobj(get(hcbar,'Parent'),'Tag','overlappolarax');
            % pax.RAxis.Label.String = 'Range (m)';
            % pax.RAxis.Label.FontSize = fsize;
            % pax.RAxis.Label.Color = lblColor;
            % pax.ThetaAxis.Label.String = 'Angle (deg)';
            % pax.ThetaAxis.Label.FontSize = fsize;
            % pax.ThetaAxis.Label.Color = lblColor;
            
            % plot polar axes
            hsurf = findobj(get(hcbar,'Parent'),'Tag','polarsurf');
            hax = get(hsurf,'Parent');
            [angmin,angmax] = bounds(obj.Angle);
            thetagrid = linspace(angmin,angmax,5);
            [rmin,rmax] = bounds(obj.Range);
            gridcolor = [0.7 0.7 0.7];
            line(hax,[rmin rmax]'*sind(thetagrid),[rmin rmax]'*cosd(thetagrid),'Color',gridcolor);
            rmaxlog = floor(log10(rmax));
            rmaxscale = 10^(rmaxlog);
            rmaxgrid = ceil(rmax/(rmaxscale))*(rmaxscale);  % in the mulitple of nearest 10's scale
            rgrid = linspace(0,rmaxgrid,6);
            rgrid(1) = rmin;
            rgrid(end) = rmax;
            rgridtheta = linspace(angmin,angmax,361);
            line(hax,sind(rgridtheta')*rgrid,cosd(rgridtheta')*rgrid,'Color',gridcolor);
            thetaticklblspacing = 0.05*(rmax-rmin);
            thetagridlabelrng = rmax+thetaticklblspacing;
            for m = 1:numel(thetagrid)
                text(thetagridlabelrng*sind(thetagrid(m)),thetagridlabelrng*cosd(thetagrid(m)),...
                sprintf('%5.1f',thetagrid(m)),'FontSize',fsize,'Color',gridcolor,'Rotation',-thetagrid(m),...
                'HorizontalAlignment','center');
            end
            [~,rexp,runit] = engunits(rgrid(end));
            rticklblspacing = 0.02*(rmax-rmin);
            for m = 1:numel(rgrid)
                if (m>1 && m<numel(rgrid)) || ...
                        (m == numel(rgrid) && (rgrid(end)-rgrid(end-1) >= 0.5*(rgrid(2)-rgrid(1)))) || ...
                        (m == 1 && (rgrid(1) <= 0.5*(rgrid(2)-rgrid(1)))) 
                    text(rgrid(m)*sind(angmax)+rticklblspacing*cosd(angmax),...
                        rgrid(m)*cosd(angmax)-rticklblspacing*sind(angmax),...
                        sprintf('%5.2f',rgrid(m)*rexp),'FontSize',fsize,'Color',gridcolor,'Rotation',-angmax,...
                        'HorizontalAlignment','left');
                end
            end
            thetaaxlblang = (thetagrid(end-1)+thetagrid(end))/2;
            thetaaxlblrng = rmax+0.12*(rmax-rmin);
            text(thetaaxlblrng*sind(thetaaxlblang),thetaaxlblrng*cosd(thetaaxlblang),...
                'Angle (degrees)','Color',gridcolor,'FontSize',fsize,'Rotation',-thetaaxlblang,...
                'HorizontalAlignment','center');
            raxlblrng = rgrid(3);
            text(raxlblrng*sind(angmin)-thetaticklblspacing*cosd(angmin),...
                raxlblrng*cosd(angmin)+thetaticklblspacing*sind(angmin),...
                sprintf('Range (%cm)',runit),'FontSize',fsize,'Color',gridcolor,'Rotation',-(90-angmax),...
                'HorizontalAlignment','center');
           
            title(plotoption.Title);
            
        else
            annotatePlot@phased.internal.AbstractRespPattern3D(obj,plotoption);
            hcbar = colorbar;
            [fsize,lblColor] = phased.internal.AbstractRespPattern.setAnnotationSizeColor;
            ylabel(hcbar,getRespLabel(plotoption),'fontsize',fsize,'Color',lblColor);
        end
    end
end

methods 
    function set.Pattern(obj,value)
        sigdatatypes.checkFiniteNonNegDblMat(obj,'Pattern',value,...
            [numel(obj.Range) numel(obj.Angle)]); %#ok<*MCSUP>
        obj.Pattern = value;
    end
    
    function set.Range(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Range',class(obj)),'Range');
        obj.Range = value;
    end

    function set.Angle(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Angle',class(obj)),'Angle');
        obj.Angle = value;
    end
end

end
% [EOF]
