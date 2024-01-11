classdef (Hidden, Sealed) AngleDopplerPattern < phased.internal.AbstractRespPattern3D
% This class if for internal use only. It may be removed in the future.

%AngleDopplerPattern   Angle Doppler response pattern
%   Hresp =
%   phased.internal.AngleDopplerPattern('PropertyName',PropertyValue,...)
%   returns an angle Doppler response pattern object Hresp. See properties
%   list below for valid PropertyNames.
%
%   AngleDopplerPattern methods:
%       plot      - Plot the angle Doppler response pattern.
%
%   AngleDopplerPattern properties:
%       Type      - 'Angle Doppler Response Pattern'.  This is read-only
%       Angle     - Angle Doppler pattern sampling angles (degrees)
%       Doppler   - Angle Doppler pattern sampling Doppler frequencies (Hz)
%       PRF       - Pulse repetition frequency (Hz)
%       Pattern   - Magnitude angle Doppler pattern
%
%   Example:
%       % Construct an angle Doppler response pattern and then plot it.
%       hresp = phased.internal.AngleDopplerPattern;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern2D, phased.internal.RespPattern3D.

%   Copyright 2008-2010 The MathWorks, Inc.

properties
    %Angle   - Angle
    %   Angle is a row vector containing the azimuth angles (in degrees) at
    %   which the response pattern are sampled. Default value of Angle is
    %   -90:90.
    Angle
    %Doppler - Doppler frequency
    %   Doppler is a row vector containing the Doppler frequencies (in Hz)
    %   at which the response pattern are sampled. Default value of Doppler
    %   is -0.5:0.1:0.5.
    Doppler
    %PRF     - Pulse Repetition Frequency
    %   PRF is the pulse repetition frequency (in Hz) used to generate the
    %   angle doppler response pattern. Default value of PRF is 1.
    PRF
    %Pattern - Response pattern
    %   Pattern is a matrix containing the samples of the magnitude
    %   response pattern at each corresponding (Angle,Doppler) pair. The
    %   number of rows of Pattern must match the number of elements in
    %   Doppler and the number of columns of Pattern must match the number
    %   of elements in Angle. The default pattern is an isotropic
    %   pattern. Default value of Pattern is all ones.
    Pattern
end

methods
    function obj = AngleDopplerPattern(varargin)
    %AngleDopplerPattern Constructor of phased.internal.AngleDopplerPattern class
        Angle = -90:90; %#ok<*PROP>
        Doppler = -0.5:0.1:0.5;
        Pattern = ones(numel(Doppler),numel(Angle));
        PRF = 1;
        sigutils.pvparse(varargin{:});
        obj.Angle = Angle;
        obj.Doppler = Doppler;
        obj.Pattern = Pattern;
        obj.PRF = PRF;
        obj.Type = 'Angle Doppler Response Pattern';
    end

    function varargout = plot(obj,varargin)
    %PLOT     Plot the response pattern
    %   plot(Hresp) plots the angle Doppler pattern Hresp.
    %
    %   plot(Hresp, 'Units', UNIT) plots the angle Doppler pattern using
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
    %       % Construct and plot an angle Doppler response pattern.
    %       hresp = phased.internal.AngleDopplerPattern; plot(hresp)
    %
    %   See also phased.internal.RespPattern3D, phased.internal.AngleDopplerPattern
        
        h = plot@phased.internal.AbstractRespPattern3D(obj,varargin{:});
        if nargout
            varargout{1} = h;
        end
    end
    
end

methods (Access = protected)
    function sortedList = getSortedPropDispList(this)  %#ok<MANU>
        % Get the sorted list of the properties to be displayed. 
        sortedList = {'Type','Angle','Doppler','PRF','Pattern'};
    end
    
    function xlbl = getXLabel(obj,~) %#ok<MANU>
    %getXLabel Return x-axis label of the plot
        xlbl = 'Angle (degrees)';
    end
    
    function ylbl = getYLabel(obj,plotoptionobj) %#ok<MANU>
    %getYLabel Return y-axis label of the plot
       
        if plotoptionobj.NormalizeDoppler
            ylblprefix = 'Normalized';
            ylblsuffix = '';
        else
            ylblprefix = '';
            ylblsuffix = '(Hz)';
        end
        ylbl = sprintf('%s Doppler Frequency %s',ylblprefix,ylblsuffix);
        ylbl = strtrim(ylbl);
    end

    function x = getXData(obj,~)
    %getXData Get x-axis data
        x = obj.Angle;
    end
    
    function y = getYData(obj,plotoptionobj)
    %getYData Get y-axis data
        y = obj.Doppler;
        if plotoptionobj.NormalizeDoppler
            y = y/obj.PRF;
        end
    end
    
    function plotoptionobj = getPlotOption(obj,varargin) %#ok<MANU>
        plotoptionobj = phased.internal.AngleDopplerPatternPlotOption(varargin{:});
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
            [numel(obj.Doppler) numel(obj.Angle)]); %#ok<*MCSUP>
        obj.Pattern = value;
    end
    
    function set.Angle(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Angle',class(obj)),'Angle');
        obj.Angle = value;
    end

    function set.Doppler(obj,value)
        validateattributes(value,{'numeric'},{'finite','real','nonnan','row'},...
            sprintf('%s.Doppler',class(obj)),'Doppler');
        obj.Doppler = value;
    end
end

end
% [EOF]
