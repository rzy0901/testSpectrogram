classdef (Hidden) AbstractAntennaElement < phased.internal.AbstractElement
%This class is for internal use only. It may be removed in the future.

%AbstractAntennaElement   Define the AbstractAntennaElement class.

%   Copyright 2009-2011 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#codegen
    
    methods (Access = protected)

        function obj = AbstractAntennaElement(varargin)
            %AbstractAntenna   Construct the AbstractAntenna class.
            obj@phased.internal.AbstractElement(varargin{:});
        end
    end
    
    methods (Abstract, Hidden)
        %This method is used by the array to compute the response of the
        %antenna element at all combinations of Azimuth, az, and Elevation,
        %el, angles. This is only used in a GPU simulation.
        epat = getgpuElemResponse(obj, az, el, freq);
    end
end




