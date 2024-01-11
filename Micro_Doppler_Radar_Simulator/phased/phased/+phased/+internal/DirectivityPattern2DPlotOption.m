classdef (Hidden,Sealed) DirectivityPattern2DPlotOption < phased.internal.AbstractRespPattern2DPlotOption
%This class is for internal use only. It may be removed in the future.

    %DirectivityPattern2DPlotOption   Define the DirectivityPattern2DPlotOption class.

    %   Copyright 2013 The MathWorks, Inc.
    
    properties (Constant)
        BroadsideOrientation = 0;
    end
    
    properties
        PlotType = 'plot';
    end

    methods

        function this = DirectivityPattern2DPlotOption(varargin)
            initPropValuePairs(this,varargin{:});
            this.Units = 'db';
            this.NormalizeResp = false;
        end
        
        function response = privRespPattern(obj,response) %#ok<INUSL>
            %PRIVRESPPATTERN Prepare response pattern
            %   PRIVRESPPATTERN is a protected method which prepares the
            %   response pattern object for plotting.
            
            % no op, Directivity is always unnormalized, in dB
        end
        
        function set.PlotType(obj,val)
            val = validatestring(val,{'plot','polar','waterfall'},...
                sprintf('%s.PlotType',class(obj)),'PlotType');
            obj.PlotType = val;
        end
    end
        
    methods (Access = protected)
        function unitlbl = getUnitLabel(obj) %#ok<MANU>
            unitlbl = 'Directivity (dBi)';
        end
    end

end

% [EOF]
