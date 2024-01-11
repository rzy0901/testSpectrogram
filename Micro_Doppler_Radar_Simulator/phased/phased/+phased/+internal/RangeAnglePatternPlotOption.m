classdef (Hidden,Sealed) RangeAnglePatternPlotOption < phased.internal.AbstractRespPatternPlotOption
%This class is for internal use only. It may be removed in the future.

%RangeDopplerPatternPlotOption   Define the RangeDopplerPatternPlotOption class.

%   Copyright 2018 The MathWorks, Inc.
%     

    properties
        Style = 'rectangular'
    end
    
    methods

        function this = RangeAnglePatternPlotOption(varargin)
            %RangeDopplerPatternPlotOption   Construct the
            %RangeDopplerPatternPlotOption class.
            this.Title = 'Range-Angle Response Pattern';
            initPropValuePairs(this,varargin{:});

        end
        
    end

    methods (Access = protected)
        function limlabel = getLimLabel(this) %#ok<MANU>
            limlabel = 'CLim';
        end
        
        function datalim = getDataLimit(this,hplotobj) %#ok<INUSL>
            for m = numel(hplotobj):-1:1
                data = get(hplotobj(m),'CData');
                templim(m,:) = [min(data(:)) max(data(:))];
            end
            datalim = [min(templim(:,1)),max(templim(:,2))];
        end  
    
    end
    

end

% [EOF]



