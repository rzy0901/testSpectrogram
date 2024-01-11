classdef (Hidden) AbstractRespPatternUV3DPlotOption < phased.internal.AbstractRespPatternPlotOption
%This class is for internal use only. It may be removed in the future.

%AbstractRespPatternUV3DPlotOption   Define the AbstractRespPatternUV3DPlotOption class.

%   Copyright 2013 The MathWorks, Inc.
    
    properties
        PlotType = 'plot';
    end
    
    methods
        
        function this = AbstractRespPatternUV3DPlotOption(varargin)
 
            initPropValuePairs(this,varargin{:});
            
        end
        
        function hax = limitDynamicPlotRange(this,hplotobj)
            hax = limitDynamicPlotRange@phased.internal.AbstractRespPatternPlotOption(this,hplotobj);
            limlabel = getLimLabel(this);
            if strncmpi(this.Units,'db',2)
                caxis(hax,get(hax,limlabel));
            end
        end
    end
    
    methods (Access = protected)
        function limlabel = getLimLabel(this) %#ok<MANU>
            limlabel = 'ZLim';
        end
        
        function datalim = getDataLimit(this,hplotobj) %#ok<INUSL>
            for m = numel(hplotobj):-1:1
                data = get(hplotobj(m),'ZData');
                templim(m,:) = [min(data(:)) max(data(:))];
            end
            datalim = [min(templim(:,1)),max(templim(:,2))];
        end  
    
    end
    
    methods
        
        function set.PlotType(this,val)
            val = validatestring(val,{'plot'},...
                sprintf('%s.PlotType',class(this)),'PlotType');
            this.PlotType = val;
        end
    end
    
end

% [EOF]
