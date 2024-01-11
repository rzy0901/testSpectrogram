classdef (Hidden,Sealed) RespPattern3DPlotOption < phased.internal.AbstractRespPattern3DPlotOption
%This class is for internal use only. It may be removed in the future.

    %RespPattern3DPlotOption   Define the RespPattern3DPlotOption class.
    
    %   Copyright 2009-2013 The MathWorks, Inc.
    %     
    
    methods
        
        function obj = RespPattern3DPlotOption(varargin)
            %RespPattern3DPlotOption   Construct the
            %RespPattern3DPlotOption class.
            obj.Title = '3D Response Pattern';
            initPropValuePairs(obj,varargin{:});
            
        end
        
    end
    
end

% [EOF]
