classdef (Hidden,Sealed) DirectivityPatternUV3DPlotOption < phased.internal.AbstractRespPatternUV3DPlotOption
%This class is for internal use only. It may be removed in the future.

%DirectivityPatternUV3DPlotOption   Define the DirectivityPatternUV3DPlotOption class.

%   Copyright 2009-2013 The MathWorks, Inc.
    
    methods
        
        function this = DirectivityPatternUV3DPlotOption(varargin)
 
            this@phased.internal.AbstractRespPatternUV3DPlotOption(varargin{:})
            this.Units = 'db';          % directivity always in db
            this.NormalizeResp = false; % directivity does not normalize
            this.Title = '3D Directivity Pattern in u-v space';
            
        end
        
        function response = privRespPattern(obj,response) %#ok<INUSL>
            %PRIVRESPPATTERN Prepare response pattern
            %   PRIVRESPPATTERN is a protected method which prepares the
            %   response pattern object for plotting.
            
            % no op, Directivity is always unnormalized, in dB
        end
        
    end
    
    methods (Access = protected)
        function unitlbl = getUnitLabel(obj) %#ok<MANU>
            unitlbl = 'Directivity (dBi)';
        end
    end
end

% [EOF]
