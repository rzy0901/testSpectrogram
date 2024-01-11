classdef (Hidden,Sealed) RespPatternUV3DPlotOption < phased.internal.AbstractRespPatternUV3DPlotOption
%This class is for internal use only. It may be removed in the future.

%RespPatternUV3DPlotOption   Define the RespPatternUV3DPlotOption class.

%   Copyright 2009-2013 The MathWorks, Inc.
    
    methods
        
        function this = RespPatternUV3DPlotOption(varargin)
 
            this@phased.internal.AbstractRespPatternUV3DPlotOption(varargin{:})
            this.Title = '3D Response Pattern in u-v space';
            
        end
        
    end
    
end

% [EOF]
