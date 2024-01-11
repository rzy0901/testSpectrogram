classdef (Hidden, Sealed) RangeDopplerPattern < phased.internal.AbstractRespPattern3D
% This class if for internal use only. It may be removed in the future.

%RangeDopplerPattern   Range Doppler response pattern
%   Hresp =
%   phased.internal.RangeDopplerPattern('PropertyName',PropertyValue,...)
%   returns a range Doppler response pattern object Hresp. See properties
%   list below for valid PropertyNames.
%
%   RangeDopplerPattern methods:
%       plot      - Plot the range Doppler response pattern.
%
%   RangeDopplerPattern properties:
%       Type      - 'Range Doppler Response Pattern'.  This is read-only
%       Range     - Range Doppler pattern sampling ranges (meter)
%       Doppler   - Range Doppler pattern sampling Doppler frequencies (Hz)
%       PRF       - Pulse Repetition Frequency
%       Pattern   - Magnitude range Doppler pattern
%
%   Example:
%       % Construct a range Doppler response pattern and then plot it.
%       hresp = phased.internal.RangeDopplerPattern;
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
    %Doppler - Doppler frequency 
    %   Doppler is a row vector containing the Doppler frequencies (in Hz)
    %   at which the response pattern are sampled. Default value of Doppler
    %   is -0.5:0.1:0.5.
    Doppler
    %PRF     - Pulse Repetition Frequency
    %   PRF is the pulse repetition frequency (in Hz) used to generate the
    %   range Doppler response pattern. Default value of PRF is 1.
    PRF
    %Pattern - Response pattern
    %   Pattern is a matrix containing the samples of the magnitude
    %   response pattern at each corresponding (Range,Doppler) pair. The
    %   number of rows of Pattern must match the number of elements in
    %   Doppler and the number of columns of Pattern must match the number
    %   of elements in Range. Default value of Pattern is all ones.
    Pattern
end

methods
    function obj = RangeDopplerPattern(varargin)
    %RangeDopplerPattern Constructor of phased.internal.RangeDopplerPattern class
        Range = 0:99; %#ok<*PROP>
        Doppler = -0.5:0.1:0.5;
        PRF = 1;
        Pattern = ones(numel(Range),numel(Doppler));
        sigutils.pvparse(varargin{:});
        obj.Range = Range;
        obj.Doppler = Doppler;
        obj.PRF = PRF;
        obj.Pattern = Pattern;
        obj.Type = 'Range Doppler Response Pattern';
    end

    function varargout = plot(obj,varargin)
    %PLOT     Plot the response pattern
    %   plot(Hresp) plots the range Doppler pattern Hresp.
    %
    %   plot(Hresp, 'Units', UNIT) plots the range Doppler pattern using
    %   the unit specified in UNIT. UNIT can be any of the following:
    %   ['mag' | 'power' | {'db'}].
    %
    %   plot(..., 'NormalizeResp', NFLAG) plots the normalized response
    %   pattern if NFLAG is true. NFLAG can be any of the following: [true
    %   | {false}].
    %
    %   plot(..., 'NormalizeDoppler', NDFLAG) plots the response pattern in
    %   normalized Doppler frequency scale if NDFLAG is true. NDFLAG can be
    %   any of the following: [{true} | false].
    %
    %   plot(..., 'Title', TITLE) uses TITLE as the title of the resulting
    %   plot.
    %
    %   plot returns an handle to the current plot object.
    %
    %   Example:
    %       % Construct and plot a range Doppler response pattern.
    %       hresp = phased.internal.RangeDopplerPattern; plot(hresp)
    %
    %   See also phased.internal.RespPattern3D,
    %   phased.internal.RangeDopplerPattern.
        
        h = plot@phased.internal.AbstractRespPattern3D(obj,varargin{:});
        if nargout
            varargout{1} = h;
        end
    end
    
end

methods (Access = protected)
    function sortedList = getSortedPropDispList(this)  %#ok<MANU>
        % Get the sorted list of the properties to be displayed. 
        sortedList = {'Type','Range','Doppler','PRF','Pattern'};
    end
    
    function ylbl = getYLabel(obj,~)  %#ok<INUSD>
    %getYLabel Return y-axis label of the plot
        ylbl = 'Range (meters)';
    end
    
    function xlbl = getXLabel(obj,plotoptionobj)  %#ok<INUSL>
    %getXLabel Return x-axis label of the plot
        if plotoptionobj.NormalizeDoppler
            xlblprefix = 'Normalized';
            xlblsuffix = '';
        else
            xlblprefix = '';
            xlblsuffix = '(Hz)';
        end
        xlbl = sprintf('%s Doppler Frequency %s',xlblprefix,xlblsuffix);
        xlbl = strtrim(xlbl);
    end

    function y = getYData(obj,~)
    %getYData Get y-axis data
        y = obj.Range;
    end
    
    function x = getXData(obj,plotoptionobj)
    %getYData Get x-axis data
        x = obj.Doppler;
        if plotoptionobj.NormalizeDoppler
            x = x/obj.PRF;
        end
    end
    
    function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
        plotoptionobj = phased.internal.RangeDopplerPatternPlotOption(varargin{:});
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
            [numel(obj.Range) numel(obj.Doppler)],'float'); %#ok<*MCSUP>
        obj.Pattern = value;
    end
    
    function set.Range(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Range',class(obj)),'Range');
        obj.Range = value;
    end

    function set.Doppler(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Doppler',class(obj)),'Doppler');
        obj.Doppler = value;
    end
end

end
% [EOF]
