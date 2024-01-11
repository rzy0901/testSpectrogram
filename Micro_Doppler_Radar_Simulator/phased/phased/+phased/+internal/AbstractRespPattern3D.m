classdef (Hidden) AbstractRespPattern3D < phased.internal.AbstractRespPattern
%This class is for internal use only. It may be removed in the future.

%AbstractRespPattern3D   Define the AbstractRespPattern3D class.

%   Copyright 2008-2010 The MathWorks, Inc.

    methods
        function varargout = plot(obj,varargin)
            %PLOT     Plot the response pattern
            %   plot(Hresp) plots the response pattern Hresp.
            %
            %   plot(Hresp, 'Units', UNIT) plots the response pattern using the
            %   unit specified in UNIT. UNIT can be any of the following:
            %   [{'mag'} | 'power' | 'db'].
            %
            %   plot(..., 'NormalizeResp', NFLAG) plots the normalized response
            %   pattern if NFLAG is true. NFLAG can be any of the following:
            %   [true | {false}].
            %
            %   plot(..., 'Title', TITLE) uses TITLE as the title of the resulting
            %   plot.
            %
            %   plot returns an handle to the current plot object.
            %
            %   Example:
            %       % Construct and plot a 3D isotropic array response pattern.
            %       hresp = phased.internal.RespPattern3D;
            %       plot(hresp)
            %
            %   See also phased.internal.RespPattern3D, phased.internal.AngleDopplerPattern
            
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
        function h = privplot(obj,plotoption)  
        %PRIVPLOT  Choose customized default plot scheme
            x = getXData(obj,plotoption);
            y = getYData(obj,plotoption);
            z = privRespPattern(plotoption,obj.Pattern);
            h = surf(x,y,z,'EdgeAlpha',0.1,'LineStyle','none');
        end
                      
        function annotatePlot(obj,plotoptionobj)
            %ANNOTATEPLOT Annotate response pattern plots
            phased.internal.AbstractRespPattern.annotatePlotSetup;
            xlabel(getXLabel(obj,plotoptionobj));
            ylabel(getYLabel(obj,plotoptionobj));
            title(plotoptionobj.Title);
        end
    
    end
    
    methods (Abstract, Access = protected)
        x = getXData(obj,plotoptionobj);
        y = getYData(obj,plotoptionobj);
        xlbl = getXLabel(obj,plotoptionobj);
        ylbl = getYLabel(obj,plotoptionobj);
        plotoptionobj = getPlotOption(obj,varargin);
    end


end

% [EOF]


