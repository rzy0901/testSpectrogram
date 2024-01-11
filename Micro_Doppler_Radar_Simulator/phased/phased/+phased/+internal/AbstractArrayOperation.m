classdef (Hidden) AbstractArrayOperation < phased.internal.AbstractSampleRateEngine
%This class is for internal use only. It may be removed in the future.

%   Copyright 2010-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %SensorArray    Sensor array
        %   Specify the sensor array as a handle. The sensor array must be
        %   an array object in the phased package. 
        SensorArray;
        %PropagationSpeed   Signal propagation speed (m/s)
        %   Specify the propagation speed (in m/s) of the signal as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('LightSpeed');
    end

    methods
        function set.SensorArray(obj,val)
            privValidateSensorArray(obj,val);
            obj.SensorArray = val;
        end
        
        function set.PropagationSpeed(obj,val)
            sigdatatypes.validateSpeed(val,...
                'phased.internal','PropagationSpeed',...
                {'double','single'},{'scalar','positive'});
            obj.PropagationSpeed = val;
        end
    end
    
    methods (Access = protected)
        
        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end

        function privValidateSensorArray(~,val) 
        %privValidateSensorArray
        %   Certain array operation is limited to certain array geometries.
        %   Each operation can then overload this method to do its own
        %   validation. By default, any array geometry is ok.
        
            validateattributes( val, { 'phased.internal.AbstractArray' }, { 'scalar' }, '', 'SensorArray');
        end
        
        function initializeSensorArray(obj)
        %initializeSensorArray
        %   Certain array operation is not initialized with ULA, so they
        %   could overload this method to initialize SensorArray to their
        %   needs. By default, a ULA is used.

            obj.SensorArray = phased.ULA;
        end
        
    end
    
    methods (Access = protected)

        function obj = AbstractArrayOperation(varargin)
            %AbstractArrayMeasure   Construct the AbstractArrayMeasure
            %class.
            setProperties(obj, nargin, varargin{:});
            if isempty(coder.target)
                if isempty(obj.SensorArray)
                    initializeSensorArray(obj);
                end
            else
                if ~coder.internal.is_defined(obj.SensorArray)
                    initializeSensorArray(obj);
                end
            end            
        end

    end
    
    methods (Access = protected)
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            s.SensorArray = saveobj(obj.SensorArray);
            s.isLocked = isLocked(obj);
        end
        
        function s = loadSubObjects(obj,s)
            if isfield(s.SensorArray,'ClassNameForLoadTimeEval')
                obj.SensorArray = eval(...
                    sprintf('%s.loadobj(s.SensorArray)',s.SensorArray.ClassNameForLoadTimeEval));
            else
                obj.SensorArray = s.SensorArray;
            end
            s = rmfield(s,'SensorArray');
        end
        
        function group = getPropertyGroups(obj)
            % Same display for short and long (no link), but showing long set
            % of properties from getPropertyGroupsImpl.
            group = getPropertyGroupsLongImpl(obj);
        end
        
        function group = getPropertyGroupsLongImpl(obj)
            % Move SensorArray property to top of first tab to get it to be first in
            % the display. Call inherited version to convert
            % getPropertyGroupsImpl to MATLAB property groups. Sensor starts
            % out by itself on the third tab.
            allGroups = getPropertyGroupsLongImpl@phased.internal.AbstractSampleRateEngine(obj);
            group = flattenGroupsAndMoveSensorToTop(obj,'SensorArray',allGroups);
        end
   end
    
    
    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl(sensorType)          
          sensorProp = matlab.system.display.internal.Property('SensorArray', ...
                                'CustomPresenter', 'phased.internal.SensorArrayDialog', ...
                                'CustomPresenterPropertyGroupsArgument', sensorType, ...
                                'Description', 'Sensor array');
          sensorGroup = matlab.system.display.SectionGroup('Title', 'Sensor Array', 'PropertyList', {sensorProp});

                    
        mainGroup = matlab.system.display.SectionGroup('TitleSource', 'Auto', 'PropertyList', {'PropagationSpeed'});
        groups = [mainGroup, sensorGroup];
      end
    end
end
