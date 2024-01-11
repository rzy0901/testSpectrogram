classdef (Hidden, Sealed) RangeSpeedPattern < phased.internal.AbstractRespPattern3D
% This class if for internal use only. It may be removed in the future.

%RangeSpeedPattern   Range speed response pattern
%   Hresp =
%   phased.internal.RangeSpeedPattern('PropertyName',PropertyValue,...)
%   returns a range speed response pattern object Hresp. See properties
%   list below for valid PropertyNames.
%
%   RangeSpeedPattern methods:
%       plot      - Plot the range speed response pattern.
%
%   RangeSpeedPattern properties:
%       Type      - 'Range Speed Response Pattern'.  This is read-only
%       Range     - Range speed pattern sampling ranges (meter)
%       Speed     - Range speed pattern sampling speed (m/s)
%       Pattern   - Magnitude range speed pattern
%
%   Example:
%       % Construct a range speed response pattern and then plot it.
%       hresp = phased.internal.RangeSpeedPattern;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern2D,
%   phased.internal.RespPattern3D.

%   Copyright 2012 The MathWorks, Inc.

properties
    %Range   - Range
    %   Range is a row vector containing the ranges (in meters) at which
    %   the response pattern are sampled. Default value of Range is 0:99.
    Range
    %Speed - Speed
    %   Speed is a row vector containing the speed (in m/s) at which the
    %   response pattern are sampled. Default value of Speed is
    %   -0.5:0.1:0.5.
    Speed
    %Pattern - Response pattern
    %   Pattern is a matrix containing the samples of the magnitude
    %   response pattern at each corresponding (Range,Speed) pair. The
    %   number of rows of Pattern must match the number of elements in
    %   speed and the number of columns of Pattern must match the number
    %   of elements in Range. Default value of Pattern is all ones.
    Pattern
end

methods
    function obj = RangeSpeedPattern(varargin)
    %RangeSpeedPattern Constructor of phased.internal.RangeSpeedPattern class
        Range = 0:99; %#ok<*PROP>
        Speed = -0.5:0.1:0.5;
        Pattern = ones(numel(Range),numel(Speed));
        sigutils.pvparse(varargin{:});
        obj.Range = Range;
        obj.Speed = Speed;
        obj.Pattern = Pattern;
        obj.Type = 'Range Speed Response Pattern';
    end

    function varargout = plot(obj,varargin)
    %PLOT     Plot the response pattern
    %   plot(Hresp) plots the range speed pattern Hresp.
    %
    %   plot(Hresp, 'Units', UNIT) plots the range speed pattern using
    %   the unit specified in UNIT. UNIT can be any of the following:
    %   ['mag' | 'power' | {'db'}].
    %
    %   plot(..., 'NormalizeResp', NFLAG) plots the normalized response
    %   pattern if NFLAG is true. NFLAG can be any of the following: [true
    %   | {false}].
    %
    %   plot(..., 'Title', TITLE) uses TITLE as the title of the resulting
    %   plot.
    %
    %   plot returns an handle to the current plot object.
    %
    %   Example:
    %       % Construct and plot a range speed response pattern.
    %       hresp = phased.internal.RangeSpeedPattern; plot(hresp)
    %
    %   See also phased.internal.RespPattern3D,
    %   phased.internal.RangeSpeedPattern.
        
        h = plot@phased.internal.AbstractRespPattern3D(obj,varargin{:});
        if nargout
            varargout{1} = h;
        end
    end
    
end

methods (Access = protected)
    function sortedList = getSortedPropDispList(this)  %#ok<MANU>
        % Get the sorted list of the properties to be displayed. 
        sortedList = {'Type','Range','Speed','Pattern'};
    end
    
    function ylbl = getYLabel(obj,~)  %#ok<INUSD>
    %getYLabel Return y-axis label of the plot
        ylbl = 'Range (meters)';
    end
    
    function xlbl = getXLabel(obj,~)  %#ok<INUSD>
    %getXLabel Return x-axis label of the plot
        xlbl = 'Speed (m/s)';
    end

    function x = getXData(obj,~)
    %getXData Get x-axis data
        x = obj.Speed;
    end
    
    function y = getYData(obj,~)
    %getYData Get y-axis data
        y = obj.Range;
    end
    
    function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
        plotoptionobj = phased.internal.RangeSpeedPatternPlotOption(varargin{:});
    end
    
    function h = privplot(obj,plotoption)  
        %PRIVPLOT  Choose customized default plot scheme
        x = getXData(obj,plotoption);
        y = getYData(obj,plotoption);
        z = privRespPattern(plotoption,obj.Pattern);
        h = imagesc(x,y,z);
        % image y axis starts from top, so we need to correct this.
        set(get(h,'Parent'),'YDir','normal');
    end
    
    function annotatePlot(obj,plotoption)
    %ANNOTATEPLOT Annotate response pattern plots
        annotatePlot@phased.internal.AbstractRespPattern3D(obj,plotoption);
        hcbar = colorbar;
        [fsize,lblColor] = phased.internal.AbstractRespPattern.setAnnotationSizeColor;
        ylabel(hcbar,getRespLabel(plotoption),'fontsize',fsize,'Color',lblColor);
    end
end

methods 
    function set.Pattern(obj,value)
        sigdatatypes.checkFiniteNonNegDblMat(obj,'Pattern',value,...
            [numel(obj.Range) numel(obj.Speed)],'float'); %#ok<*MCSUP>
        obj.Pattern = value;
    end
    
    function set.Range(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Range',class(obj)),'Range');
        obj.Range = value;
    end

    function set.Speed(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Speed',class(obj)),'Speed');
        obj.Speed = value;
    end
end

end
% [EOF]
