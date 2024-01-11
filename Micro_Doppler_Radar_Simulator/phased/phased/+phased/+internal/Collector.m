classdef (Sealed, StrictDefaults) Collector < phased.internal.AbstractSensorOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%Collector Narrowband signal collector
%   H = phased.Collector creates a narrowband signal collector System
%   object, H. The object collects incident narrowband signals from given
%   directions using a sensor array or a single element.
%
%   H = phased.Collector(Name,Value) creates a collector object, H, with
%   the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax when the Wavefront property is 'Plane':
%
%   Y = step(H,X,ANG) collects multiple plane wave signals specified in
%   X with all collecting elements, adds them together, and returns the
%   combined output Y. ANG contains the signal incident directions. For
%   each plane wave signal, the received signal is collected using the
%   phase approximation of the time delays across collecting elements in
%   the far field.
%
%   Each column of X is considered a separate signal. ANG is a 2-row
%   matrix. Each column of ANG specifies the incident direction of the
%   corresponding signal in the form of an [AzimuthAngle; ElevationAngle]
%   pair (in degrees). The number of columns in ANG must equal the number
%   of columns in X.
%
%   Y is a matrix whose number of columns equals the number of subarrays
%   if Sensor contains subarrays, or the number of elements otherwise. Each
%   column of Y contains the output of the corresponding element/subarray
%   in response to all the signals in X.
%
%   Y = step(H,X,ANG,LAXES) specifies the collector's local coordinate
%   system in LAXES when you set the EnablePolarization property to true.
%   LAXES is a 3x3 matrix whose columns specify the local coordinate
%   system's orthonormal x, y, and z axes, respectively. Each axis is
%   specified in [x;y;z] form measured in the global coordinate system.  X
%   is a 1xM struct array where M is the number of entries in ANG. Each
%   struct contains three fields: X, Y, and Z. Each field contains the X,
%   Y, and Z component of the polarized input signal, also measured in the
%   global coordinate system, respectively. Within each field is a column
%   vector representing a separate incoming signal. The signals in X, Y,
%   and Z fields must have the same dimension.
%
%   Y = step(H,X,ANG,W) uses W as the weight vector when
%   the WeightsInputPort is set to true. W must be a column vector of
%   length M, where M is the number of collecting elements or subarrays.
%  
%   Y = step(H,X,ANG,STEER) uses STEER as the subarray
%   steering angle (in degrees). STEER must be a length-2 column
%   vector in the form of [AzimuthAngle; ElevationAngle]. This syntax is
%   only applicable when you use subarrays in the Sensor property and set
%   the SubarraySteering property in the Sensor to either 'Phase' or
%   'Time'.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANG,LAXES,W,STEER)
%
%   Step method syntax when the Wavefront property is 'Unspecified':
%
%   Y = step(H,X,ANG) collects different signals at different
%   collecting elements. The number of columns in X must be the same as
%   the number of collecting elements. Each column of X is assumed to be
%   arriving from the direction specified in the corresponding column of
%   ANG and is only collected by the corresponding collecting element.
%   Note that this syntax does not allow the Sensor to contain subarrays.
%
%   Y = step(H,X,ANG,LAXES) specifies the collector's local coordinate
%   system in LAXES when you set the EnablePolarization property to true.
%   LAXES is a 3x3 matrix whose columns specify the local coordinate
%   system's orthonormal x, y, and z axes, respectively. Each axis is
%   specified in [x;y;z] form measured in the global coordinate system. X
%   is a 1xN struct array where N is the number of elements/subarrays in
%   the sensor array. Each struct contains three fields: X, Y, and Z. Each
%   field contains the X, Y, and Z component of the polarized input signal,
%   also measured in the global coordinate system, respectively. Within
%   each field is a column vector representing a separate signal arriving
%   at the corresponding collecting element.
%
%   Y = step(H,X,ANG,W) uses W as the weight vector when
%   the WeightsInputPort is set to true. W must be a column vector of
%   length M, where M is the number of collecting elements.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANG,LAXES,W)
%
%   Collector methods:
%
%   step     - Collect signals (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a collector object with same property values
%   isLocked - Locked status (logical)
%
%   Collector properties:
%
%   Sensor             - Handle of the sensor
%   PropagationSpeed   - Propagation speed
%   OperatingFrequency - Operating frequency
%   SignalType         - Signal type
%   Wavefront          - Type of incoming wavefront
%   EnablePolarization - Enable polarization
%   WeightsInputPort   - Enable weights input
%
%   % Examples:
%
%   % Example 1: 
%   %   Collect signal with a single antenna.
%
%   ha = phased.IsotropicAntennaElement;
%   hc = phased.Collector('Sensor',ha,'OperatingFrequency',1e9);
%   x = [1;1];
%   incidentAngle = [10 30]';
%   y = step(hc,x,incidentAngle);
%
%   % Example 2: 
%   %   Collect a far field signal with a 5-element array.
%
%   ha = phased.ULA('NumElements',5);
%   hc = phased.Collector('Sensor',ha,'OperatingFrequency',1e9);
%   x = [1;1];
%   incidentAngle = [10 30]';
%   y = step(hc,x,incidentAngle);
%
%   % Example 3: 
%   %   Collect signal with a 3-element antenna array. Each antenna 
%   %   collects a separate input signal from a separate direction.
%
%   ha = phased.ULA('NumElements',3);
%   hc = phased.Collector('Sensor',ha,'OperatingFrequency',1e9,...
%                         'Wavefront','Unspecified');
%   x = rand(10,3);   % Each column is a separate signal for one element
%   incidentAngle = [10 0; 20 5; 45 2]'; % 3 angles for 3 signals
%   y = step(hc,x,incidentAngle);
%
%   See also phased, phased.WidebandCollector.

%   Copyright 2009-2013 The MathWorks, Inc.
%     

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %OperatingFrequency Operating frequency (Hz)
    %   Specify the operating frequency (in Hz) as a positive scalar. The
    %   default value is 3e8.
    OperatingFrequency = 3e8;    
    %SignalType     Signal type
    %   Specify whether the signal should be treated as one of 'Field' |
    %   'Power', where 'Field' is the default. When you set SignalType to
    %   'Field', the collected signal is scaled by the field at the
    %   corresponding direction. When you set SignalType to 'Power', the
    %   collected signal is scaled by the directivity at the corresponding
    %   direction.
    SignalType = 'Field';
    %Wavefront Type of incoming wavefront
    %   Specify the type of incoming wavefront as one of 'Plane' |
    %   'Unspecified', where the default is 'Plane'. When you set the
    %   Wavefront property to 'Plane', the input signals are assumed to be
    %   multiple plane waves impinging on the entire array. Each plane wave
    %   is received by all collecting elements. If you set the Wavefront
    %   property to 'Unspecified', the input signals are assumed to be
    %   individual waves impinging on individual sensors. Wavefront must be
    %   'Plane' if you set SignalType to 'Power'.
    Wavefront = 'Plane';
end

properties (Nontunable, Logical)
    %EnablePolarization  Enable polarization
    %   Set this property to true to enable polarization. Set this property
    %   to false to ignore polarization. The default value of this property
    %   is false. This property applies when the sensor specified in the
    %   Sensor property is capable of simulating polarization.
    EnablePolarization = false;
end

properties(Constant, Hidden)
    WavefrontSet = matlab.system.StringSet({'Plane','Unspecified'});
    SignalTypeSet = matlab.system.StringSet({'Field','Power'});
end

properties (Access = private, Nontunable)
    cSteeringVector;
    cIntegratedPattern;
    pDOF;
end

properties (Access = private, Logical, Nontunable)
    pIsPowerSignal
end

methods
    function set.OperatingFrequency(obj,value)
        validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'OperatingFrequency');
        obj.OperatingFrequency = value;
    end    
    function obj = Collector(varargin)
        obj@phased.internal.AbstractSensorOperation(varargin{:});
    end
end

methods (Access = protected)
    
    function num = getNumInputsImpl(obj)
        num = getNumInputsImpl@phased.internal.AbstractSensorOperation(obj);
        if obj.EnablePolarization
            num = num+1;
        end
    end
    
    function resetImpl(obj)
        if obj.pUseArray
            if (obj.Wavefront(1) == 'P') %'Plane'
                reset(obj.cSteeringVector);
                if obj.pIsPowerSignal
                    reset(obj.cIntegratedPattern);
                end
            else
                reset(obj.cSensor);   
            end
        else
            reset(obj.cSensor);
            if obj.pIsPowerSignal
                reset(obj.cIntegratedPattern);
            end
        end   
    end
    
    function releaseImpl(obj)
        if obj.pUseArray
            if (obj.Wavefront(1) == 'P') %'Plane'
                release(obj.cSteeringVector);
                if obj.pIsPowerSignal
                    release(obj.cIntegratedPattern);
                end
            else
                release(obj.cSensor);
            end
        else
            release(obj.cSensor);
            if obj.pIsPowerSignal
                release(obj.cIntegratedPattern);
            end
        end
    end
    
    function validatePropertiesImpl(obj)
        cond = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                strncmpi(obj.Wavefront,'Unspecified',1);
        if cond
            coder.internal.errorIf(cond,'phased:phased:collector:invalidSettingsForSubarray','Wavefront','Plane','Sensor');
        end
        
        cond = strcmp(obj.SignalType,'Power') && strcmp(obj.Wavefront,'Unspecified');
        if cond
            coder.internal.errorIf(cond,'phased:system:array:IrrelevantSetting','Power','Wavefront','Plane');
        end

        cond = obj.EnablePolarization && ~isPolarizationCapable(obj.Sensor);
        if cond
            coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
        end
    end
    
    function validateInputsImpl(obj,x,ang,lclaxes,wArg,stangArg)
        coder.extrinsic('mat2str');        
        if obj.EnablePolarization
            cond = ~isa(x,'struct');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','struct');
            end
            
            if isempty(coder.target)
                fn = fieldnames(x);
                flag_hasXYZ = isempty(setdiff({'X','Y','Z'},fn));
                cond = ~flag_hasXYZ;
            else
                %fieldnames not supported for codegen
                cond = ~(isfield(x(1),'X') && isfield(x(1),'Y') && isfield(x(1),'Z'));
            end
            if cond
                coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationXYZStruct');
            end
            
            xsize = size(x);
            for m = 1:xsize(2)
                x_x = x(m).X;
                x_y = x(m).Y;
                x_z = x(m).Z;
                cond =  ~isa(x_x,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType',sprintf('X(%d).X',m),'double');
                end
                cond =  ~iscolumn(x_x) || isempty(x_x);
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:inputMustBeColVector',sprintf('X(%d).X',m));
                end
                cond =  ~isa(x_y,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType',sprintf('X(%d).Y',m),'double');
                end
                cond =  ~iscolumn(x_y) || isempty(x_y);
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Y',m));
                end
                cond =  ~isa(x_z,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType',sprintf('X(%d).Z',m),'double');
                end
                cond =  ~iscolumn(x_z) || isempty(x_z);
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Z',m));
                end
                cond = numel(x_x)~=numel(x_y) || numel(x_x)~=numel(x_z);
                if cond
                    coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                        'X,Y,Z',sprintf('X(%d)',m));
                end
            end
        else
            cond =  ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDataType','X','double');
            end
            cond =  ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                      'MATLAB:system:inputMustBeMatrix','X');
            end
            xsize = size(x);
        end

        UseArrayFlag = isa(obj.Sensor,'phased.internal.AbstractArray') || ...
            isa(obj.Sensor,'phased.internal.AbstractSubarray');
        if UseArrayFlag
            M = getNumElements(obj.Sensor);
        else
            M = 1;
        end
        cond = UseArrayFlag && ~strncmpi(obj.Wavefront,'Plane',1) && xsize(2)~=M;
        if cond
            coder.internal.errorIf(cond,'phased:phased:collector:DimensionMismatch');
        end

        angsize = size(ang);
        cond = xsize(2) ~= angsize(2);
        if cond
            coder.internal.errorIf(cond,'phased:phased:collector:AngleDimensionMismatch');
        end
        cond = ~isa(ang,'double');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','Ang','double');
        end
        cond = ~isreal(ang);
        if cond
            coder.internal.errorIf(cond,'phased:step:NeedReal', 'Ang');
        end

        if obj.EnablePolarization
            % check lclAxes is 3x3
            sigdatatypes.validate3DCartCoord(lclaxes,'','LAxes',...
                {'size',[3 3]});
        end

        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmp(obj.Sensor.SubarraySteering,'None',1)
            if ~obj.WeightsInputPort
                if obj.EnablePolarization
                    stang = wArg;
                else
                    stang = lclaxes;
                end
            else
                if ~obj.EnablePolarization
                    stang = wArg;
                else
                    stang = stangArg;
                end
            end
            sz_stang = size(stang);
            cond = ~ismatrix(stang) || isempty(stang);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','Steer');
            end
            cond = sz_stang(1) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:NeedTwoRows','Steer');
            end
            cond = sz_stang(2) > 1;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:NeedOneColumn','Steer');
            end
            cond = ~isreal(stang);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidAngle','Steer');
            end
            cond = ~isa(stang,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Steer','double');
            end

        end
        
        if obj.WeightsInputPort
            if ~obj.EnablePolarization
                w = lclaxes;
            else
                w = wArg;
            end
           
            if UseArrayFlag
                Mw = getDOF(obj.Sensor);
            else
                Mw = 1;
            end
            cond = ~isa(w,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','W','double');
            end
            wsize = size(w);
            cond = ~isequal(wsize, [Mw 1]);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDimensions','W',...
                    coder.const(mat2str([Mw 1])), ...
                    coder.const(mat2str(wsize)));
            end
        end
    end
    
    function setupImpl(obj,~,~,~,~,~,~,~) 
        setupImpl@phased.internal.AbstractSensorOperation(obj);
        
        if obj.pUseArray
            M = getDOF(obj.cSensor);
        else
            M = 1;
        end
        obj.pDOF = M;
        obj.pIsPowerSignal = strcmp(obj.SignalType,'Power');
        if obj.pUseArray 
            if strncmp(obj.Wavefront,'Plane',1)
                obj.cSteeringVector = phased.SteeringVector(...
                    'SensorArray',obj.cSensor,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'IncludeElementResponse',true,...
                    'EnablePolarization',obj.EnablePolarization);
                if obj.pIsPowerSignal
                    sensorElem = getElementHandle(obj.cSensor);
                    if isElementFromAntenna(obj.cSensor) || ...
                            isa(sensorElem,'phased.internal.AntennaAdapter')
                        obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                            'Sensor',obj.cSensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',obj.WeightsInputPort,...
                            'EnablePolarization',obj.EnablePolarization);
                    else
                        obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                            'Sensor',obj.cSensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',obj.WeightsInputPort,...
                            'EnablePolarization',obj.EnablePolarization);
                    end
                end
            end
        elseif obj.pIsPowerSignal
            if isPolarizationEnabled(obj.cSensor)
                if isElementFromAntenna(obj.cSensor) || ...
                        isa(obj.cSensor,'phased.internal.AntennaAdapter')
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                        'Sensor',obj.cSensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',false,... % single element weights is just a scalar, separable
                        'EnablePolarization',true);
                else
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                        'Sensor',obj.cSensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',false,... % single element weights is just a scalar, separable
                        'EnablePolarization',true);
                end
            else
                if isElementFromAntenna(obj.cSensor) || ...
                        isa(obj.cSensor,'phased.internal.AntennaAdapter')
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                        'Sensor',obj.cSensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',false,... % single element weights is just a scalar, separable
                        'EnablePolarization',false);
                else
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                        'Sensor',obj.cSensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',false,... % single element weights is just a scalar, separable
                        'EnablePolarization',false);
                end
            end
        end
    end
    
    function flag = isInputComplexityLockedImpl(obj,index)
        flag = true;
        if index == 1
            flag = false;
        end
        if obj.WeightsInputPort && ...
                ( (obj.EnablePolarization && (index == 4)) || ...
                ( ~obj.EnablePolarization && (index == 3)))
                flag = false;
        end
    end

    function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
        flag = false;
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractSensorOperation(obj);
        if isLocked(obj)
            s.pDOF = obj.pDOF;
            s.cSteeringVector = saveobj(obj.cSteeringVector);
            s.cIntegratedPattern = saveobj(obj.cIntegratedPattern);
            s.pIsPowerSignal = obj.pIsPowerSignal;
        end
    end

    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractSensorOperation(obj,s);
        if isfield(s,'isLocked')
            if s.isLocked
                obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                s = rmfield(s,'cSteeringVector');
                if isfield(s,'cIntegratedPattern') && ...
                            isfield(s.cIntegratedPattern,'ClassNameForLoadTimeEval')
                    obj.cIntegratedPattern = eval(...
                        sprintf('%s.loadobj(s.cIntegratedPattern)',s.cIntegratedPattern.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cIntegratedPattern');
                end
            end
            s = rmfield(s,'isLocked');
        end
    end
    
    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        s = loadSubObjects(obj,s);
        if isfield(s,'pNumSensorArrayElements')
            s = rmfield(s,'pNumSensorArrayElements');
        end
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function yout = stepImpl(obj,x,ang,lclaxes,wArg,stangArg)
        
        if obj.EnablePolarization
            fc = obj.OperatingFrequency;
            N = obj.pDOF;
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end
            if obj.pUseArray
                if strncmp(obj.Wavefront,'Plane',1)
                    if ~obj.WeightsInputPort
                        if obj.pNeedSteeringAngle
                            stang = wArg;
                        end
                        w = ones(N,1);
                    else
                        if obj.pNeedSteeringAngle
                            stang = stangArg;
                        end
                        w = wArg;
                    end
                    if obj.pNeedSteeringAngle
                        arrayEffect = step(...
                            obj.cSteeringVector,fc,ang,stang);
                        if obj.pIsPowerSignal
                            if obj.WeightsInputPort
                                intresp = step(obj.cIntegratedPattern,fc,w,stang);
                            else
                                intresp = step(obj.cIntegratedPattern,fc,stang);
                            end
                            arrayEffect.H = phased.internal.normalizeIntegratedPower(arrayEffect.H,intresp.H,false);
                            arrayEffect.V = phased.internal.normalizeIntegratedPower(arrayEffect.V,intresp.V,false);
                        end
                    else
                        arrayEffect = step(...
                            obj.cSteeringVector,fc,ang);
                        if obj.pIsPowerSignal
                            if obj.WeightsInputPort
                                intresp = step(obj.cIntegratedPattern,fc,w);
                            else
                                intresp = step(obj.cIntegratedPattern,fc);
                            end
                            arrayEffect.H = phased.internal.normalizeIntegratedPower(arrayEffect.H,intresp.H,false);
                            arrayEffect.V = phased.internal.normalizeIntegratedPower(arrayEffect.V,intresp.V,false);
                        end
                    end
                    
                    arrayEffect_x = complex(zeros(N,size(ang,2)));
                    arrayEffect_y = arrayEffect_x;
                    arrayEffect_z = arrayEffect_x;
                    for m = size(ang,2):-1:1
                        arrayEffect_temp = sph2cartvec(...
                            [arrayEffect.H(:,m).';arrayEffect.V(:,m).';zeros(1,N)],...
                            ang(1,m),ang(2,m));
                        arrayEffect_temp = phased.internal.local2globalvec(...
                            arrayEffect_temp,lclaxes);
                        arrayEffect_x(:,m) = arrayEffect_temp(1,:).';
                        arrayEffect_y(:,m) = arrayEffect_temp(2,:).';
                        arrayEffect_z(:,m) = arrayEffect_temp(3,:).';
                    end

                    y = x_X*arrayEffect_x.'+x_Y*arrayEffect_y.'+...
                        x_Z*arrayEffect_z.';
                else  % unspecified wavefront
                    if ~obj.WeightsInputPort
                        w = ones(N,1);
                    else
                        w = wArg;
                    end
                    g = step(obj.cSensor,fc,ang);
                    arrayEffect.H = diag(g.H);
                    arrayEffect.V = diag(g.V);
                    
                    arrayEffect_x = complex(zeros(size(ang,2),1));
                    arrayEffect_y = arrayEffect_x;
                    arrayEffect_z = arrayEffect_x;
                    for m = size(ang,2):-1:1
                        arrayEffect_temp = sph2cartvec(...
                            [arrayEffect.H(m).';arrayEffect.V(m).';0],...
                            ang(1,m),ang(2,m));
                        arrayEffect_temp = phased.internal.local2globalvec(...
                            arrayEffect_temp,lclaxes);
                        arrayEffect_x(m,1) = arrayEffect_temp(1);
                        arrayEffect_y(m,1) = arrayEffect_temp(2);
                        arrayEffect_z(m,1) = arrayEffect_temp(3);
                    end
                    
                    y = bsxfun(@times,x_X,arrayEffect_x.')+...
                        bsxfun(@times,x_Y,arrayEffect_y.')+...
                        bsxfun(@times,x_Z,arrayEffect_z.');
                end

            else
                % Single sensor
                arrayEffect = step(obj.cSensor,fc,ang);
                if ~obj.WeightsInputPort
                    w = 1;
                else
                    w = wArg;
                end
                if obj.pIsPowerSignal
                    intresp = step(obj.cIntegratedPattern,fc);
                    arrayEffect.H = phased.internal.normalizeIntegratedPower(arrayEffect.H,intresp.H,false)*abs(w)^2;
                    arrayEffect.V = phased.internal.normalizeIntegratedPower(arrayEffect.V,intresp.V,false)*abs(w)^2;
                end
                
                arrayEffect_x = complex(zeros(1,size(ang,2)));
                arrayEffect_y = arrayEffect_x;
                arrayEffect_z = arrayEffect_x;
                for m = size(ang,2):-1:1
                    arrayEffect_temp = sph2cartvec(...
                        [arrayEffect.H(m);arrayEffect.V(m);0],...
                        ang(1,m),ang(2,m));
                    arrayEffect_temp = phased.internal.local2globalvec(...
                        arrayEffect_temp,lclaxes);
                    arrayEffect_x(m) = arrayEffect_temp(1).';
                    arrayEffect_y(m) = arrayEffect_temp(2).';
                    arrayEffect_z(m) = arrayEffect_temp(3).';
                end

                y = x_X*arrayEffect_x(:)+x_Y*arrayEffect_y(:)+...
                    x_Z*arrayEffect_z(:);
            end
            yout = bsxfun(@times,y,w.');
            
       
        else  % no polarization
            fc = obj.OperatingFrequency;
            N = obj.pDOF;
            if obj.pUseArray
                if strncmp(obj.Wavefront,'Plane',1)
                    if ~obj.WeightsInputPort
                        if obj.pNeedSteeringAngle
                            stang = lclaxes;
                        end
                        w = ones(N,1);
                    else
                        if obj.pNeedSteeringAngle
                            stang = wArg;
                        end
                        w = lclaxes;
                    end
                    if obj.pNeedSteeringAngle
                        arrayEffect = step(obj.cSteeringVector,fc,ang,stang);
                        if obj.pIsPowerSignal
                            if obj.WeightsInputPort
                                intresp = step(obj.cIntegratedPattern,fc,w,stang);
                            else
                                intresp = step(obj.cIntegratedPattern,fc,stang);
                            end
                            arrayEffect = phased.internal.normalizeIntegratedPower(arrayEffect,intresp,false);
                        end
                    else
                        arrayEffect = step(obj.cSteeringVector,fc,ang);
                        if obj.pIsPowerSignal
                            if obj.WeightsInputPort
                                intresp = step(obj.cIntegratedPattern,fc,w);
                            else
                                intresp = step(obj.cIntegratedPattern,fc);
                            end
                            arrayEffect = phased.internal.normalizeIntegratedPower(arrayEffect,intresp,false);
                        end
                    end
                    y = x*arrayEffect.';
                else
                    g_temp = step(obj.cSensor,fc,ang);
                    if isPolarizationEnabled(obj.cSensor)
                        g = diag(hypot(g_temp.H,g_temp.V));
                    else
                        g = diag(g_temp);
                    end
                    if obj.WeightsInputPort
                        w = lclaxes;
                        w = w.*g;
                    else
                        w = g;
                    end
                    y = x;
                end
            else
                % Single sensor
                g_temp = step(obj.cSensor,fc,ang);
                if ~obj.WeightsInputPort
                    w = 1;
                else
                    w = lclaxes;
                end
                if isPolarizationEnabled(obj.cSensor)
                    if obj.pIsPowerSignal
                        intresp = step(obj.cIntegratedPattern,fc);
                        g_temp.H = phased.internal.normalizeIntegratedPower(g_temp.H,intresp.H,false)*abs(w)^2;
                        g_temp.V = phased.internal.normalizeIntegratedPower(g_temp.V,intresp.V,false)*abs(w)^2;
                    end
                    g = hypot(g_temp.H,g_temp.V);
                else
                    if obj.pIsPowerSignal
                        intresp = step(obj.cIntegratedPattern,fc);
                        g_temp = phased.internal.normalizeIntegratedPower(g_temp,intresp,false);
                    end
                    g = g_temp;
                end
                y = x*g;
            end
            yout = bsxfun(@times,y,w.');            
        end
    
    end
        
end

methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractSensorOperation;
      
            pWavefront = ...
                matlab.system.display.internal.Property('Wavefront', ...
                                                        'IsGraphical', false);

            dEnablePolarization = ...
                matlab.system.display.internal.Property('EnablePolarization', ...
                                                        'IsGraphical', false);
            props = {'OperatingFrequency',...
                        'SignalType',...
                        pWavefront,...
                        dEnablePolarization,...
                        'WeightsInputPort'};
            groups(1).PropertyList = [groups(1).PropertyList props];
            action = matlab.system.display.Action(@phased.Collector.onActionCalled, ...
              'Label', 'Analyze');
          matlab.system.display.internal.setCallbacks(action, ...
          'DialogAppliedFcn', @phased.Collector.onDialogApplied, ...
          'SystemDeletedFcn', @phased.Collector.onSystemDeleted);
            groups(2).Actions = action;
  end
  function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:CollectorTitle')),...
          'Text',getString(message('phased:library:block:CollectorDesc')));
  end
end
methods (Static,Hidden)
    function onActionCalled(actionData, sysObj)
        if ~phased.apps.internal.SensorArrayViewer.SettingsDialog.verifySensorArray(sysObj.Sensor)
            return;
        end
        f = actionData.UserData;
        if isempty(f) || isStale(f);
            actionData.UserData = ...
                phased.apps.internal.SensorArrayViewer.SensorArrayViewer(sysObj);
        end
    end
    function onDialogApplied(actionData, sysObj)
 
        f = actionData.UserData;
        if ~isempty(f)
            if isStale(f)
                actionData.UserData = [];
            else
                if ~phased.apps.internal.SensorArrayViewer.SettingsDialog.verifySensorArray(sysObj.Sensor)
                    return;
                end
                updateSettingsWithSystemObject(f,sysObj);
            end
        end
    end
    function onSystemDeleted(actionData)
      f = actionData.UserData;
      if ~isempty(f) && ~isStale(f)
         close(f);
         actionData.UserData = [];
      end
   end
end
methods (Access = protected) %for Simulink
    function varargout = getOutputNamesImpl(~)
        varargout = {''};
    end
    function varargout = getInputNamesImpl(obj)
        varargout = {'X','Ang','LAxes','W','Steer'};
        lastIdx = 3;
        if obj.EnablePolarization
            varargout(lastIdx) = {'LAxes'};
            lastIdx = lastIdx+1;
        end
        if obj.WeightsInputPort
            varargout(lastIdx) = {'W'};
            lastIdx = lastIdx+1;
        end
        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmpi(obj.Sensor.SubarraySteering,'None',1);
            varargout(lastIdx) = {'Steer'};
        end
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Narrowband\n Rx Array');
    end
    
    function varargout = getOutputSizeImpl(obj)
        szX = propagatedInputSize(obj,1);
        %Output is number of snapshots * number of elements
        if  ~isa(obj.Sensor,'phased.internal.AbstractElement');
            M = getDOF(obj.Sensor);
        else
            M = 1;
        end
        varargout{1} = [szX(1) M];
    end
    function varargout = isOutputFixedSizeImpl(obj)
        varargout{1} = propagatedInputFixedSize(obj, 1);
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj,1);
    end
    function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
        varargout{1} = true;
    end
end

end

