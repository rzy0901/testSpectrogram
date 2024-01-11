classdef (Hidden,Sealed) RespPatternUV2DPlotOption < phased.internal.AbstractRespPattern2DPlotOption
%This class is for internal use only. It may be removed in the future.

    %RespPatternUV2DPlotOption   Define the RespPatternUV2DPlotOption class.

    %   Copyright 2009-2011 The MathWorks, Inc.
    
    properties
        PlotType = 'plot';
    end

    methods

        function this = RespPatternUV2DPlotOption(varargin)
            initPropValuePairs(this,varargin{:});
        end
        
        function set.PlotType(obj,val)
            val = validatestring(val,{'plot','waterfall'},...
                sprintf('%s.PlotType',class(obj)),'PlotType');
            obj.PlotType = val;
        end
    end
        
end

% [EOF]
