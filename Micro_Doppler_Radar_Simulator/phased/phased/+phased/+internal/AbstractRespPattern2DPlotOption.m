classdef (Hidden) AbstractRespPattern2DPlotOption < phased.internal.AbstractRespPatternPlotOption
%This class is for internal use only. It may be removed in the future.

%AbstractRespPattern2DPlotOption   Define the AbstractRespPattern2DPlotOption class.

%   Copyright 2009-2011 The MathWorks, Inc.
    
    properties (Abstract)
        PlotType
    end
    
    properties
        FreqUnit = '';
    end

    methods (Access = protected)

        function this = AbstractRespPattern2DPlotOption(varargin)
        end
        
    end
    
    methods (Access = protected)
        function limlabel = getLimLabel(this) 
            % For line and polar plots, set the limlabel to the Y-axis.
            % For waterfall plots, set the limlabel to the Z-axis.
            if strcmp(this.PlotType,'waterfall')
                limlabel = 'ZLim';
            else
                limlabel = 'YLim';
            end
        end
        
        function datalim = getDataLimit(this,hplotobj)
            % For line and polar plots, set the data is along the Y-axis.
            % For waterfall plots, data is along the Z-axis.
            hax = get(hplotobj(1),'Parent');
            plotobjlist = findobj(hax,'Type',get(hplotobj(1),'Type'));
            for m = numel(plotobjlist):-1:1
                if strcmp(this.PlotType,'waterfall')
                    data = get(plotobjlist(m),'ZData');
                else
                    data = get(plotobjlist(m),'YData');
                end
                templim(m,:) = [min(data(:)) max(data(:))];
            end
            datalim = [min(templim(:,1)),max(templim(:,2))];
        end          
        
        
    end

    methods

        function freq = scalePlotFreqs(this,freq)
            [~,scaleFactor,unit] = engunits(median(freq));
            this.FreqUnit = unit;
            freq = freq * scaleFactor;
        end
        
        function grouplbl = getGroupLabel(this,type)
            % no op for azimuth or elevation since waterfall not supported
            if strcmp(type,'Frequency')
                grouplbl  = sprintf('%s (%sHz)',type,this.FreqUnit); 
            end
        end
        
    end

end

% [EOF]
