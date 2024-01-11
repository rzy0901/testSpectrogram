classdef (Hidden,Sealed) RangeSpeedPatternPlotOption < phased.internal.AbstractRespPatternPlotOption
%This class is for internal use only. It may be removed in the future.

%RangeSpeedPatternPlotOption   Define the RangeSpeedPatternPlotOption class.

%   Copyright 2012 The MathWorks, Inc.
%     

    methods

        function this = RangeSpeedPatternPlotOption(varargin)
            %RangeSpeedPatternPlotOption   Construct the
            %RangeSpeedPatternPlotOption class.
            this.Title = 'Range-Speed Response Pattern';
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



