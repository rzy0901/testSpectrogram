classdef (Hidden) PolarizationResponsePlotOption < ...
phased.internal.AbstractRespPatternPlotOption
%This class is for internal use only. It may be removed in the future.

%PolarizationResponsePlotOption   Define the
%PolarizationResponsePlotOption class.

%   Copyright 2012 The MathWorks, Inc.

    methods

        function this = PolarizationResponsePlotOption(varargin)
            %PolarizationResponsePlotOption   Construct the
            %PolarizationResponsePlotOption class.
            initPropValuePairs(this,varargin{:});

        end
        
    end
    
    methods (Access = protected)

        function limlabel = getLimLabel(this) %#ok<MANU>
            %method1   Example method
            %   method1(H) Add a complete method description here
            limlabel = 'ZLim';

        end
        
        function datalim = getDataLimit(this,hplotobj) %#ok<INUSL>
            for m = numel(hplotobj):-1:1
                data = get(hplotobj(m),'ZData');
                templim(m,:) = [min(data(:)) max(data(:))];
            end
            datalim = [min(templim(:,1)),max(templim(:,2))];
        end  
    
        function unitlbl = getUnitLabel(obj)
            switch obj.Units
                case 'mag',
                    unitlbl = 'Scattering Magnitude';
                case 'power',
                    unitlbl = 'Normalized Cross Section';
                case 'db'
                    unitlbl = 'Normalized Cross Section (dB)';
            end
        end
    end
end

% [EOF]
