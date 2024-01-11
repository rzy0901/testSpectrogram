classdef (Hidden) AbstractRadiationPatternUV2D < phased.internal.AbstractRespPattern2D
%This class is for internal use only. It may be removed in the future.

%   Copyright 2013 The MathWorks, Inc.
    
    
    properties
        %U - Sampling vector for U values
        %   U is a column vector containing the U values at which the
        %   response pattern are sampled.
        U
    end
    
    methods
        function obj = AbstractRadiationPatternUV2D(varargin)
        %AbstractRadiationPatternUV2D Constructor of phased.internal.AbstractRadiationPatternUV2D class
            U = (-1:0.01:1).'; %#ok<*PROP>
            Freq = 1e9;
            Pattern = ones(size(U));
            sigutils.pvparse(varargin{:});
            obj.U = U;
            obj.Freq = Freq; %store frequencies corresponding to patterns
            obj.Pattern = Pattern;
            obj.Type = '2D Response Pattern in U/V space';
        end
        
    end
    
    methods (Access = protected)
        function angles = getPlotAngles(obj)
            angles = obj.U;
        end
        
        function xlbl = getXLabel(obj,~) %#ok<INUSD>
            %getXLabel Return x-axis label of the plot
            xlbl = 'U';            
        end
        
        function plotTitle = genPlotTitle(obj) %#ok<MANU>
            %GENPLOTTITLE Generate a title string for ResponsePattern2D
            %plots.

            % Built a plot title string containing the cut type and angle
            plotTitle = 'Response in U Space';
        end
        
        
    end
    
    methods (Access = protected)
        function sortedList = getSortedPropDispList(this)  %#ok<MANU>
            % Get the sorted list of the properties to be displayed.
            sortedList = {'Type','U','Pattern','Freq'};
        end
        
    end


methods
        
        function obj = set.U(obj,value)
            validateattributes(value,{'double'},...
                {'real','column','<=',1,'>=',-1},...
                sprintf('%s.U',class(obj)),'U');
            obj.U = value;
        end
        
    end
end



% [EOF]
