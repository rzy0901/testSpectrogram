classdef (Hidden) AbstractNarrowbandArrayProcessing < phased.internal.AbstractArrayOperation
%This class is for internal use only. It may be removed in the future.

%AbstractNarrowbandArrayProcessing   Define the AbstractNarrowbandArrayProcessing class.

%   Copyright 2010-2013 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)

        %OperatingFrequency     Operating frequency (Hz)
        %   Specify the operating frequency (in Hz) of the system as a
        %   scalar. The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 3e8;
    end

    methods (Access = protected)

        function obj = AbstractNarrowbandArrayProcessing(varargin)
 
            obj@phased.internal.AbstractArrayOperation(varargin{:});

        end
        
    end
    
    methods

        function set.OperatingFrequency(obj,val)
            
            sigdatatypes.validateFrequency(val,'phased.internal',...
                'OperatingFrequency',{'double','single'},{'scalar'});
            obj.OperatingFrequency = val;
        end
    end
    
    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl(sensorType)
        groups = getPropertyGroupsImpl@phased.internal.AbstractArrayOperation(sensorType);
        props = 'OperatingFrequency';
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
    end
end


