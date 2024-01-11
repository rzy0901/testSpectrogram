classdef (Hidden, Sealed) RangePattern < phased.internal.AbstractRespPattern
% This class if for internal use only. It may be removed in the future.

%RangePattern   Range Doppler response pattern
%   Hresp =
%   phased.internal.RangePattern('PropertyName',PropertyValue,...)
%   returns a range response pattern object Hresp. See properties
%   list below for valid PropertyNames.
%
%   RangePattern methods:
%       plot      - Plot the range response pattern.
%
%   RangePattern properties:
%       Type      - 'Range Response Pattern'.  This is read-only
%       Range     - Range Doppler pattern sampling ranges (meter)
%       Pattern   - Magnitude range pattern
%
%   Example:
%       % Construct a range response pattern and then plot it.
%       hresp = phased.internal.RangePattern;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern2D,
%   phased.internal.RespPattern3D.

%   Copyright 2016 The MathWorks, Inc.

properties
    Pattern
    %Range   - Range
    %   Range is a column vector containing the ranges (in meters) at which
    %   the response pattern are sampled. Default value of Range is
    %   (0:99)'.
    Range
end

methods
    function obj = RangePattern(varargin)
    %RangePattern Constructor of phased.internal.RangePattern class
        Range = (0:99)'; %#ok<*PROP>
        Pattern = ones(numel(Range),1);
        sigutils.pvparse(varargin{:});
        obj.Range = Range;
        obj.Pattern = Pattern;
        obj.Type = 'Range Response Pattern';
    end

    function varargout = plot(obj,varargin)
    %PLOT     Plot the range response pattern
    %   plot(Hresp) plots the range pattern Hresp.
    %
    %   plot(Hresp, 'Units', UNIT) plots the range pattern using the unit
    %   specified in UNIT. UNIT can be any of the following: ['mag' |
    %   'power' | {'db'}].
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
    %       % Construct and plot a range response pattern.
    %       hresp = phased.internal.RangePattern; plot(hresp)
    %
    %   See also phased.internal.RespPattern3D,
    %   phased.internal.RangePattern.
        
    plotoption = getPlotOption(obj,varargin{:});
    h = privplot(obj,plotoption);
    annotatePlot(obj,plotoption);
    
    limitDynamicPlotRange(plotoption,h);
    
    if nargout == 1
        varargout{1} = h;
    end
    end
end

methods (Access = protected)
    function sortedList = getSortedPropDispList(this)  %#ok<MANU>
        % Get the sorted list of the properties to be displayed. 
        sortedList = {'Type','Range','Pattern'};
    end
    
    function xlbl = getXLabel(obj,~)  %#ok<INUSD>
    %getXLabel Return x-axis label of the plot
        xlbl = 'Range (meters)';
    end

    function x = getXData(obj,~)
    %getYData Get x-axis data
        x = obj.Range;
    end
    
    function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
        plotoptionobj = phased.internal.RangePatternPlotOption(varargin{:});
    end
    
    function h = privplot(obj,plotoption)  
        %PRIVPLOT  Choose customized default plot scheme
        x = getXData(obj,plotoption);
        y = privRespPattern(plotoption,obj.Pattern);
        h = plot(x,y);
    end
    
    function annotatePlot(obj,plotoption)
        %ANNOTATEPLOT Annotate response pattern plots
        phased.internal.AbstractRespPattern.annotatePlotSetup;
        ylabel(getRespLabel(plotoption));
        xlabel(getXLabel(obj,plotoption));
        title(plotoption.Title);
    end
end

methods
    function set.Range(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','column'},...
            sprintf('%s.Range',class(obj)),'Range');
        obj.Range = value;
    end

end

end
% [EOF]
