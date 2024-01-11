classdef (Hidden) PolarizationResponse < phased.internal.AbstractRespPattern3D
%This class is for internal use only. It may be removed in the future.

%PolarizationResponse   Define the PolarizationResponse class.

%   Copyright 2012 The MathWorks, Inc.

    properties

        OrientationAngle;
        EllipticityAngle;
        Pattern;
        ResponseType;
    end

    methods

        function this = PolarizationResponse(varargin)
            %PolarizationResponse   Construct the PolarizationResponse
            %class.
            OrientationAngle = -90:90;
            EllipticityAngle = -45:45;
            Pattern = zeros(181,91);
            ResponseType = 'Co-Pol';
            sigutils.pvparse(varargin{:});
            this.OrientationAngle = OrientationAngle;
            this.EllipticityAngle = EllipticityAngle;
            this.Pattern = Pattern;
            this.ResponseType = ResponseType;
            this.Type = 'Polarization Response';
        end

        function varargout = plot(this, varargin)
            %method1   Example method
            %   method1(H) Add a complete method description here
            h = plot@phased.internal.AbstractRespPattern3D(this,varargin{:});
            if nargout
                varargout{1} = h;
            end

        end
    end
    
    methods (Access = protected)
        function sortedList = getSortedPropDispList(this) %#ok<MANU>
            sortedList = {'Type','OrientationAngle','EllipticityAngle','Pattern','ResponseType'};
        end
        
        function xlbl = getXLabel(this,~) %#ok<INUSD>
            xlbl = 'Ellipticity Angle (Degrees)';
        end
        
        function ylbl = getYLabel(this,~) %#ok<INUSD>
            ylbl = 'Orientation Angle (Degrees)';
        end
        
        function x = getXData(this,~)
            x = this.EllipticityAngle;
        end
        
        function y = getYData(this,~)
            y = this.OrientationAngle;
        end
        
        function plotoptionobj = getPlotOption(this,varargin)
            plotoptionobj = phased.internal.PolarizationResponsePlotOption(varargin{:});
            if isempty(plotoptionobj.Title)
                if strncmp(this.ResponseType,'Co',2)
                    plotoptionobj.Title = 'Co-Pol Response';
                else
                    plotoptionobj.Title = 'Cross-Pol Response';
                end
            end
        end
        
        function h = privplot(obj,plotoption) 
            x = getXData(obj,plotoption);
            y = getYData(obj,plotoption);
            z = privRespPattern(plotoption,obj.Pattern);
            h = surf(x,y,z,'EdgeColor','none');
        end
        
        function annotatePlot(obj,plotoptionobj)
            annotatePlot@phased.internal.AbstractRespPattern3D(obj,plotoptionobj);
            resplbl = getRespLabel(plotoptionobj);
            zlabel(resplbl);
        end
        
    end
    
    methods
        
        function set.OrientationAngle(this,val)
            sigdatatypes.validateAngle(val,'',...
                'OrientationAngle',{'>=',-90,'<=',90});
            this.OrientationAngle = val;
        end
        
        function set.EllipticityAngle(this,val)
            sigdatatypes.validateAngle(val,'',...
                'EllipticityAngle',{'>=',-45,'<=',45});
            this.EllipticityAngle = val;
        end
               
        function set.ResponseType(this,val)
            val = validatestring(val,{'Co-Pol','Cross-Pol'},...
                '','ResponseType');
            this.ResponseType = val;
        end
    end
end

% [EOF]
