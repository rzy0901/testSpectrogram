classdef (Hidden) AbstractSensorOperation < phased.internal.AbstractSampleRateEngine
%This class is for internal use only. It may be removed in the future.

%   Copyright 2010-2014 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %Sensor Handle of the sensor
    %   Specify the sensor as a sensor array object or an element object in
    %   the phased package. The default value is the object returned by
    %   phased.ULA.
    Sensor;
    %PropagationSpeed   Propagation speed (m/s)
    %   Specify the propagation speed (in m/s) as a positive scalar. The
    %   default value is the light speed.
    PropagationSpeed = physconst('LightSpeed');
end

properties (Nontunable, Logical) 
    %WeightsInputPort Enable weights input
    %   Set this property to true to input weights. The default value is
    %   false.
    WeightsInputPort = false;
end

properties (Access = protected, Nontunable)
    cSensor;        % handle to the sensor
end

properties (Access = protected, Nontunable, Logical)
    pUseArray = false      % whether to use array
    pNeedSteeringAngle;
end

methods
    function set.Sensor(obj,value)
        validateattributes( value, ...
            {'phased.internal.AbstractArray',...
            'phased.internal.AbstractElement',...
            'phased.internal.AbstractSubarray','em.Antenna'},...
            {'scalar'}, '', 'Sensor');
        obj.Sensor = value;
    end
    function set.PropagationSpeed(obj,value)
        validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'PropagationSpeed');
        obj.PropagationSpeed = value;
    end
end
    
methods (Access = protected)
    
    function obj = AbstractSensorOperation(varargin)
        setProperties(obj, nargin, varargin{:});
        if isempty(coder.target)
            if isempty(obj.Sensor)
                obj.Sensor = phased.ULA;
            end
        else
            if ~coder.internal.is_defined(obj.Sensor)
                obj.Sensor = phased.ULA;
            end
        end
    
    end
    
    function setupImpl(obj,~,~,~,~)

        if isempty(coder.target)
            obj.cSensor = cloneSensor(obj.Sensor);
        else
            if isElementFromAntenna(obj.Sensor)
                coder.internal.errorIf(true, ...
                    'phased:system:element:AntennaToolboxCodegenNotSupported','em.Antenna','phased.CustomAntennaElement');
            else
                obj.cSensor = clonecg(obj.Sensor);
            end
        end

        obj.pUseArray = isa(obj.Sensor,'phased.internal.AbstractArray') || ...
            isa(obj.Sensor,'phased.internal.AbstractSubarray');
        obj.pNeedSteeringAngle = ...
            isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
            ~strncmpi(obj.Sensor.SubarraySteering,'None',1);
    end

    
    function num = getNumInputsImpl(obj)
        num = 2;
        if obj.WeightsInputPort
            num = num+1;
        end
        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmpi(obj.Sensor.SubarraySteering,'None',1);
            num = num+1;
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
        s.Sensor = saveobj(obj.Sensor);
        s.isLocked = isLocked(obj);
        if isLocked(obj)
            s.pUseArray = obj.pUseArray;
            s.cSensor = saveobj(obj.cSensor);
            s.pNeedSteeringAngle = obj.pNeedSteeringAngle;
        end
    end
    
    function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
        if index == 1
            flag = false;
        else
            flag = true;
        end
    end
    
    function s = loadSubObjects(obj,s)
        if isfield(s.Sensor,'ClassNameForLoadTimeEval')
            obj.Sensor = eval(...
                sprintf('%s.loadobj(s.Sensor)',s.Sensor.ClassNameForLoadTimeEval));
        else
            obj.Sensor = s.Sensor;
        end
        s = rmfield(s,'Sensor');
        if isfield(s,'isLocked')
            if s.isLocked
                obj.cSensor = eval(...
                    sprintf('%s.loadobj(s.cSensor)',s.cSensor.ClassNameForLoadTimeEval));
                s = rmfield(s,'cSensor');
            end
        end
    end
    
    function group = getPropertyGroups(obj)
            % Same display for short and long (no link), but showing long set
            % of properties from getPropertyGroupsImpl.
            group = getPropertyGroupsLongImpl(obj);
        end
    
    function group = getPropertyGroupsLongImpl(obj)
        % Move Sensor property to top of first tab to get it to be first in
        % the display. Call inherited version to convert
        % getPropertyGroupsImpl to MATLAB property groups. Sensor starts
        % out by itself on the third tab.
        allGroups = getPropertyGroupsLongImpl@phased.internal.AbstractSampleRateEngine(obj);
        group = flattenGroupsAndMoveSensorToTop(obj,'Sensor',allGroups);
    end
end

methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        sensorProp = matlab.system.display.internal.Property('Sensor', ...
                              'CustomPresenter', 'phased.internal.SensorArrayDialog', ...
                              'Description', 'Sensor array');
        sensorGroup = matlab.system.display.SectionGroup('Title', 'Sensor Array', 'PropertyList', {sensorProp});
        mainGroup = matlab.system.display.SectionGroup('TitleSource', 'Auto', 'PropertyList', {'PropagationSpeed'});
        groups = [mainGroup, sensorGroup];
    end
    
    function sensorObj = getSensorObject(systemHandle)
        sensorObj = [];
        if isa(systemHandle, 'matlab.System')
            sensorObj = systemHandle.Sensor;
        else
            sensor = get_param(systemHandle, 'Sensor');            
            try %#ok<EMTC>
                sensorObj = slResolve(sensor, systemHandle);  
            catch e %#ok<NASGU>
            end
        end
    end
end


end

