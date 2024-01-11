classdef (Hidden,Sealed) RangePatternPlotOption < phased.internal.AbstractRespPatternPlotOption
%This class is for internal use only. It may be removed in the future.

%RangePatternPlotOption   Define the RangePatternPlotOption class.

%   Copyright 2016 The MathWorks, Inc.
%     

    methods

        function this = RangePatternPlotOption(varargin)
            %RangePatternPlotOption   Construct the
            %RangePatternPlotOption class.
            this.Title = 'Range Response Pattern';
            initPropValuePairs(this,varargin{:});

        end
        
    end

    methods (Access = protected)
        function limlabel = getLimLabel(this) %#ok<MANU>
            limlabel = 'YLim';
        end
        
        function datalim = getDataLimit(this,hplotobj) %#ok<INUSL>
            for m = numel(hplotobj):-1:1
                data = get(hplotobj(m),'YData');
                templim(m,:) = [min(data(:)) max(data(:))];
            end
            datalim = [min(templim(:,1)),max(templim(:,2))];
        end  
    
    end
    
end

% [EOF]



