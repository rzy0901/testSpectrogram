classdef (Hidden,Sealed) DirectivityPattern3DPlotOption < phased.internal.AbstractRespPattern3DPlotOption
%This class is for internal use only. It may be removed in the future.

    %DirectivityPattern3DPlotOption   Define the DirectivityPattern3DPlotOption class.
    
    %   Copyright 2013 The MathWorks, Inc.
    %     
    
    methods
        
        function this = DirectivityPattern3DPlotOption(varargin)
            %DirectivityPattern3DPlotOption   Construct the
            %DirectivityPattern3DPlotOption class.
            this@phased.internal.AbstractRespPattern3DPlotOption(varargin{:});
            this.Units = 'db';          % directivity always in db
            this.NormalizeResp = false; % directivity does not normalize
            this.Title = '3D Directivity Pattern';
            
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
