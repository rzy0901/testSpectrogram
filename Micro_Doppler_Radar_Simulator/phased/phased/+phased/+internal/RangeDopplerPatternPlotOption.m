classdef (Hidden,Sealed) RangeDopplerPatternPlotOption < phased.internal.AbstractRespPatternPlotOption
%This class is for internal use only. It may be removed in the future.

%RangeDopplerPatternPlotOption   Define the RangeDopplerPatternPlotOption class.

%   Copyright 2012 The MathWorks, Inc.
%     

    properties

       NormalizeDoppler = true;
    end

    methods

        function this = RangeDopplerPatternPlotOption(varargin)
            %RangeDopplerPatternPlotOption   Construct the
            %RangeDopplerPatternPlotOption class.
            this.Title = 'Range-Doppler Response Pattern';
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
    
    methods

        function set.NormalizeDoppler(this,val)
            validateattributes(val,{'logical'},{'scalar'},...
                sprintf('%s.NormalizeDoppler',class(this)),'NormalizeDoppler');
            this.NormalizeDoppler = val;
        end
    end

end

% [EOF]



