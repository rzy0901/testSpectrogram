classdef (Hidden,Sealed) RespPattern2DPlotOption < phased.internal.AbstractRespPattern2DPlotOption
%This class is for internal use only. It may be removed in the future.

    %RespPattern2DPlotOption   Define the RespPattern2DPlotOption class.

    %   Copyright 2009-2011 The MathWorks, Inc.
    
    properties (Constant)
        BroadsideOrientation = 0;
    end
    
    properties
        PlotType = 'plot';
    end

    methods

        function this = RespPattern2DPlotOption(varargin)
            initPropValuePairs(this,varargin{:});
        end
        
        function set.PlotType(obj,val)
            val = validatestring(val,{'plot','polar','waterfall'},...
                sprintf('%s.PlotType',class(obj)),'PlotType');
            obj.PlotType = val;
        end
    end
        
end

% [EOF]
