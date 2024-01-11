classdef (Sealed,StrictDefaults) Radiator < phased.internal.AbstractSensorOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%Radiator Narrowband signal radiator
%   H = phased.Radiator creates a narrowband signal radiator System object,
%   H. The object returns radiated narrowband signals for given directions
%   using a sensor array or a single element.
%
%   H = phased.Radiator(Name,Value) creates a radiator object, H, with the
%   specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax when the CombineRadiatedSignals property is true:
%
%   Y = step(H,X,ANG) radiates signal X to multiple directions
%   specified in ANG, and combines the radiated signal for each direction
%   in Y. The combination is performed using the phase approximation of
%   time delays across radiating elements at the far field.
%
%   X can be either a vector or a matrix. If X is a vector, it is
%   radiated through all elements. If X is a matrix, its column number
%   must equal the number of subarrays if Sensor contains subarrays, or the
%   number of elements otherwise. Each column of X is radiated by
%   corresponding element/subarray.
%
%   ANG is a 2-row matrix. Each column of ANG specifies a radiating
%   direction in the form of an [AzimuthAngle; ElevationAngle] pair (in
%   degrees). Y is a matrix whose number of columns equals the number of
%   radiating directions. Each column of Y contains the combined far field
%   output from all elements/subarrays to the corresponding direction.
%
%   Y = step(H,X,ANG,LAXES) specifies the radiator's local coordinate
%   system in LAXES when you set the Polarization property to 'Combined'.
%   LAXES is a 3x3 matrix whose columns specify the local coordinate
%   system's orthonormal x, y, and z axes, respectively. Each axis is
%   specified in [x;y;z] form measured in the global coordinate system. Y
%   is a 1xM struct array where M is the number of entries in ANG. Each
%   struct in the array has three fields: X, Y, and Z. Each field contains
%   the X, Y, and Z component, also measured in the global coordinate
%   system, of the polarized far field output, respectively. Within each
%   field is a column vector containing the combined far field output from
%   all elements/subarrays to the corresponding direction.
%
%   Y = step(H,XH,XV,ANG,LAXES) radiates signal XH from H polarization port
%   and XV from V polarization port to multiple directions specified in
%   ANG, and combines the radiated signal for each direction in Y. The
%   combination is performed using the phase approximation of time delays
%   across radiating elements at the far field. This syntax applies when
%   you set the Polarization property to 'Dual'.
%
%   Both XH and XV can be either a vector or a matrix. If XH or XV is a
%   vector, it is radiated through all elements. If XH or XV is a matrix,
%   its column number must equal the number of subarrays if Sensor contains
%   subarrays, or the number of elements otherwise. Each column of XH or XV
%   is radiated by corresponding element/subarray. The dimensions of XH and
%   XV must be the same.
%
%   Y = step(...,W) uses W as the weight vector when the WeightsInputPort
%   property is true. W must be a column vector whose length is the same as
%   the number of radiating elements or subarrays.
%
%   Y = step(...,STEER) uses STEER as the subarray steering angle (in
%   degrees). STEER must be a length-2 column vector in the form of
%   [AzimuthAngle; ElevationAngle]. This syntax is only applicable when you
%   use subarrays in the Sensor property and set the SubarraySteering
%   property in the Sensor to either 'Phase' or 'Time'.
%
%   Y = step(...,WS) uses WS as the weights applied to each element in
%   the subarray. WS can be either a matrix or a cell array. This syntax is
%   only applicable when you use subarrays in the Sensor property and set
%   the SubarraySteering property in the Sensor to 'Custom'.
%   
%   If the Sensor property is a phased.ReplicatedSubarray, WS must be an
%   NSExN matrix where NSE is the number of elements in each individual
%   subarray and N is the number of subarrays. Each column in WS specifies
%   the weights for the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its individual
%   subarrays have same number of elements, WS must be an NSExN matrix
%   where NSE is the number of elements in each individual subarray and N
%   is the number of subarrays. Each column in WS specifies the weights for
%   the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its subarrays
%   can have different number of elements, WS can be either an NSExN
%   matrix, where NSE indicates the number of elements in the largest
%   subarray and N is the number of subarrays, or a 1xN cell array, where N
%   is the number of subarrays and each cell contains a column vector whose
%   length is the same as the number of elements of the corresponding
%   subarray.  If WS is a matrix, the first K entries in each column, where
%   K is the number of elements in the corresponding subarray, specifies
%   the weights for the elements in the corresponding subarray. If WS is a
%   cell array, each cell in the array is a column vector specifying the
%   weights for the elements in the corresponding subarray. 
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,XH,XV,ANG,LAXES,W,STEER)
%
%   or
%
%   Y = step(H,X,ANG,LAXES,W,WS)
%
%   Step method syntax when the CombineRadiatedSignals property is false:
%
%   Y = step(H,X,ANG) radiates the input signal to different angles for
%   different radiating elements. The number of columns of ANG is the
%   same as the number of radiating elements and each column of ANG
%   specifies the radiating direction for the corresponding input signal
%   through the corresponding radiating element. Note that this syntax does
%   not allow the Sensor to contain subarrays. Each column of Y represents
%   the radiated signal from the corresponding radiating element.
%
%   Y = step(H,X,ANG,LAXES) specifies the radiator's local coordinate
%   system in LAXES when you set the Polarization property to 'Combined'.
%   LAXES is a 3x3 matrix whose columns specify the local coordinate
%   system's orthonormal x, y, and z axes, respectively. Each axis is
%   specified in [x;y;z] form measured in the global coordinate system.  Y
%   is a 1xN struct array where N is the number of elements/subarrays in
%   the sensor array. Each struct in the array has three fields: X, Y, and
%   Z. Each field contains the X, Y, and Z component of the polarized far
%   field output, also measured in the global coordinate system,
%   respectively. Within each field is a column vector representing the
%   radiated signal from the corresponding radiating element.
%  
%   Y = step(H,XH,XV,ANG,LAXES) radiates signal XH from H polarization port
%   and signal XV from V polarization port to different angles for
%   different radiating elements. The dimensions of XH and XV must be the
%   same. This syntax applies when you set the Polarization property
%   to 'Dual'.
%
%   Y = step(...,W) uses W as the weight vector when the WeightsInputPort
%   property is true. W must be a column vector whose length is the same as
%   the number of radiating elements.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,XH,XV,ANG,LAXES,W)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   Radiator methods:
%
%   step     - Radiate signals (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a radiator object with same property values
%   isLocked - Locked status (logical)
%
%   Radiator properties:
%
%   Sensor                      - Handle of the sensor
%   PropagationSpeed            - Propagation speed
%   OperatingFrequency          - Operating frequency
%   CombineRadiatedSignals      - Combine radiated signals
%   SensorGainMeasure           - Sensor gain measure 
%   Polarization                - Polarization configuration
%   WeightsInputPort            - Enable weights input
%
%   % Examples:
%
%   % Example 1: 
%   %   Radiate signal with a single antenna.
%
%   antenna = phased.IsotropicAntennaElement;
%   radiator = phased.Radiator('Sensor',antenna,...
%               'OperatingFrequency',300e6);
%   x = [1;1];
%   radiatingAngle = [30 10]';
%   y = radiator(x,radiatingAngle);
%
%   % Example 2: 
%   %   Radiate a far field signal with a 5-element array.
%
%   array = phased.ULA('NumElements',5);
%   radiator = phased.Radiator('Sensor',array,'OperatingFrequency',300e6);
%   x = [1;1];
%   radiatingAngle = [30 10; 20 0]'; % two directions
%   y = radiator(x,radiatingAngle);
%
%   % Example 3: 
%   %   Radiate signal with a 3-element antenna array. Each antenna 
%   %   radiates a separate signal to a separate direction.
%
%   array = phased.ULA('NumElements',3);
%   radiator = phased.Radiator('Sensor',array,'OperatingFrequency',1e9,...
%                        'CombineRadiatedSignals',false);
%   x = [1 2 3;1 2 3];
%   radiatingAngle = [10 0; 20 5; 45 2]'; % One angle for one antenna
%   y = radiator(x,radiatingAngle);
%
%   See also phased, phased.Collector, phased.WidebandRadiator.

%   Copyright 2009-2017 The MathWorks, Inc.
%     

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] V. N. Bringi and V. Chandrasekar, Polarimetric Doppler Weather
%   Radar: Principles and Applications, Cambridge University Press, 2001


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %OperatingFrequency     Operating frequency (Hz)
    %   Specify the operating frequency (in Hz) as a positive scalar. The
    %   default value is 3e8.
    OperatingFrequency = 3e8;    
end

properties (Nontunable, Logical) 
    %CombineRadiatedSignals Combine radiated signals
    %   Set this property to true to combine radiated signals from all
    %   radiating elements. Set this property to false to obtain the
    %   radiated signal for each radiating element. The default value is
    %   true. 
    CombineRadiatedSignals = true;
end

properties (Nontunable)
    %SensorGainMeasure   Sensor gain measure 
    %   Specify how sensor gain is measured as one of 'dB' | 'dBi', where
    %   'dB' is the default value. This property only applies when you set
    %   the CombineRadiatedSignals property to true. 
    %
    %   When you set this property to 'dB', the input signal power is
    %   scaled by the sensor's power pattern (in dB) at the corresponding
    %   direction and then combined. When you set this property to 'dBi',
    %   the input signal power is scaled by the directivity pattern (in
    %   dBi) at the corresponding direction and then combined. The dBi
    %   option is helpful when you want to compare results with those
    %   predicted by the radar equation that uses dBi to specify the
    %   antenna gain. Note that the computation for the dBi option is
    %   expensive as it requires an integration over all directions to
    %   compute the total radiated power of the sensor.
    SensorGainMeasure = 'dB'
end

properties (Nontunable) 
    %Polarization   Polarization configuration
    %   Specify the radiator's polarization configuration as one of 'None'
    %   | 'Combined' | 'Dual', where the default is 'None'. When you set
    %   this property to 'None', the signal is considered as polarization
    %   independent. When you set this property to 'Combined', the radiated
    %   signal contains the polarization information and the signal
    %   reflects the sensor's native polarization. When you set this
    %   property to 'Dual', the H and V polarization of the sensor can be
    %   fed independently.
    Polarization = 'None'
end

properties (Dependent, Hidden)
    %EnablePolarization  Enable polarization
    %   Set this property to true to enable polarization. Set this property
    %   to false to ignore polarization. The default value of this property
    %   is false. This property applies when the sensor specified in the
    %   Sensor property is capable of simulating polarization.
    EnablePolarization
end

properties (Access = private, Dependent)
    %pEnableDualPolarizationInput  Enable dual polarization input signal
    %   Set this property to true to enable dual polarization input signal
    %   in separate H and V polarization. Set this property to false to use
    %   same signal for both polarizations. The default value of this
    %   property is false. This property applies when you set the
    %   EnablePolarization property to true.
    pEnableDualPolarizationInput 
end

properties(Constant, Hidden)
    PolarizationSet = matlab.system.StringSet({'None','Combined','Dual'});
    SensorGainMeasureSet = matlab.system.StringSet({'dB','dBi'});
end

properties (Access = private, Nontunable)
    cSteeringVector;
    cIntegratedPattern
    pDOF;
    pCommonSignalForChannels;
    pPowerDistributionFactor;
    pApplyDirectivityGain
end

methods
    function set.OperatingFrequency(obj,value)
        validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'OperatingFrequency');
        obj.OperatingFrequency = value;
    end    
    
    function set.EnablePolarization(obj,value)
        if coder.target('MATLAB')
            warning(message('phased:phased:ObsoletePropertyByOneProperties',...
                'EnablePolarization','Polarization'));
        end
        if value
            obj.Polarization = 'Combined';
        else
            obj.Polarization = 'None';
        end
    end
    
    function value = get.EnablePolarization(obj)
        if coder.target('MATLAB')
            warning(message('phased:phased:ObsoletePropertyByOneProperties',...
                'EnablePolarization','Polarization'));
        end
        value = ~strcmp(obj.Polarization,'None');
    end
    
    function value = get.pEnableDualPolarizationInput(obj)
        value = strcmp(obj.Polarization,'Dual');
    end
    
end

methods
    function obj = Radiator(varargin)
        obj@phased.internal.AbstractSensorOperation(varargin{:});
    end
end
    
methods (Access = protected)
      
    function num = getNumInputsImpl(obj)
        num = getNumInputsImpl@phased.internal.AbstractSensorOperation(obj);
        if ~strcmp(obj.Polarization,'None')
            num = num+1;
            if obj.pEnableDualPolarizationInput
                num = num+1;
            end
        end
    end
    
    function validatePropertiesImpl(obj)
        cond = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~obj.CombineRadiatedSignals;
        if cond
            coder.internal.errorIf(cond,'phased:phased:collector:invalidSettingsForSubarray','CombineRadiatedSignals','true','Sensor');
        end
        
        cond = ~strcmp(obj.Polarization,'None') && ~isPolarizationCapable(obj.Sensor);
        if cond
            coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
        end
        
    end

    function flag = isInactivePropertyImpl(obj,prop)
        % Return false if property is visible based on object 
        % configuration, for the command line and System block dialog
        if ~obj.CombineRadiatedSignals && strcmp(prop,'SensorGainMeasure')
            flag = true;
        else
            flag = false;
        end
    end
    
    function validateInputsImpl(obj,x,xvArg,angArg,lclaxesArg,wArg,stangArg)
        coder.extrinsic('mat2str');
        polflag = ~strcmp(obj.Polarization,'None');
        dualpolflag = obj.pEnableDualPolarizationInput;
        UseArrayFlag = isa(obj.Sensor,'phased.internal.AbstractArray') || ....
            isa(obj.Sensor,'phased.internal.AbstractSubarray');
        if UseArrayFlag
            M = getDOF(obj.Sensor);
        else
            M = 1;
        end
        
        % input signal
        cond = ~isa(x,'double');
        if cond
            if dualpolflag
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Xh','double');
            else
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','double');
            end
        end
        cond = ~ismatrix(x) || isempty(x);
        if cond
            if dualpolflag
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','Xh');
            else
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','X');
            end
        end
            
        xsize = size(x);
       
        cond = xsize(2)>1 && xsize(2)~=M;
        if cond
            coder.internal.errorIf(cond,'phased:phased:Radiator:DimMismatch');
        end
        
        validateNumChannels(obj,x);
        
        if dualpolflag
            xv = xvArg;
            cond = ~isa(xv,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Xv','double');
            end
            cond = ~ismatrix(xv) || isempty(xv);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','Xv');
            end

            xvsize = size(xv);
            cond = xvsize(2)>1 && xvsize(2)~=M;
            if cond
                coder.internal.errorIf(cond,'phased:phased:Radiator:DimMismatch');
            end
            
            cond = ~isequal(xvsize,xsize);
            if cond
                coder.internal.errorIf(cond,'phased:phased:sizeMismatch','Xh','Xv');
            end
            
            ang = angArg;
        else
            ang = xvArg;
        end

        % angle
        % For near field, the angle must match the channels
        angsize = size(ang);
        cond = UseArrayFlag && ~obj.CombineRadiatedSignals && angsize(2)~=M;
        if cond
            coder.internal.errorIf(cond,'phased:phased:Radiator:InvalidAngle','CombineRadiatedSignals');
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
        
        if polflag
            if dualpolflag
                lclaxes = lclaxesArg;
            else
                lclaxes = angArg;
            end
            % check lclaxes to be 3x3  
            sigdatatypes.validate3DCartCoord(lclaxes,'','LAxes',...
                {'size',[3 3]});
        end
        
        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmp(obj.Sensor.SubarraySteering,'None',1)
            if ~obj.WeightsInputPort
                if polflag
                    if dualpolflag
                        stang = wArg;
                    else
                        stang = lclaxesArg;
                    end
                else
                    stang = angArg;
                end
            else
                if polflag
                    if dualpolflag
                        stang = stangArg;
                    else
                        stang = wArg;
                    end
                else
                    stang = lclaxesArg;
                end
            end
            if strncmp(obj.Sensor.SubarraySteering,'Phase',1) || ...
                    strncmp(obj.Sensor.SubarraySteering,'Time',1)
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
            else
                % does not support multiple weights for multiple
                % frequency yet because this is still an analog
                % behavior so at any moment, there is only one set of
                % weights can be applied.

                ws = stang;  % weights
                Ns = getNumSubarrays(obj.Sensor);
                cond = (~iscell(ws) && ~ismatrix(ws)) || isempty(ws);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:expectedCellOrMatrix','WS');
                end
                Nse = zeros(1,Ns);
                for m = 1:Ns
                    Nse(m) = getNumElements(obj.Sensor,m);
                end
                if iscell(ws)
                    cond = ~isrow(ws) || (numel(ws)~= Ns);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:phased:expectedMatrixSize','WS',1,Ns);
                    end
                    for m = 1:Ns
                        cond = ~iscolumn(ws{m}) || (numel(ws{m})~=Nse(m));
                        if cond
                            coder.internal.errorIf(cond, ...
                                'phased:system:array:SubarrayElementWeightsSizeMismatch',...
                                m,'WS',Nse(m));
                        end
                        cond = ~isa(ws{m},'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                'phased:system:array:SubarrayElementWeightsInvalidDataType',...
                                m,'WS','double');
                        end
                    end
                else
                    sz_ws = size(ws);
                    Nsemax = max(Nse);
                    cond = ~isequal(sz_ws,[Nsemax Ns]);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:phased:expectedMatrixSize','WS',Nsemax,Ns);
                    end
                    cond = ~isa(ws,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:invalidInputDataType','WS','double');
                    end
                end
            end

        end
        % weights
        if obj.WeightsInputPort
            if polflag
                if dualpolflag
                    w = wArg;
                else
                    w = lclaxesArg;
                end
            else
                w = angArg;
            end
            wsize = size(w);
            cond = ~isequal(wsize, [M 1]);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDimensions','W',...
                    coder.const(mat2str([M 1])), ...
                    coder.const(mat2str(wsize)));
            end
            cond = ~isa(w,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','W','double');
            end
        end
            
            
    end
    
    function setupImpl(obj,x,~,~,~,~,~)
        setupImpl@phased.internal.AbstractSensorOperation(obj);
        if obj.pUseArray
            M = getDOF(obj.cSensor);
        else
            M = 1;            
        end
        obj.pDOF = M;
        
        obj.pNumInputChannels = getNumChannels(obj,x);
        obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        
        if obj.pUseArray && (obj.pValidatedNumInputChannels == 1)
            obj.pCommonSignalForChannels = true;
        else
            obj.pCommonSignalForChannels = false;
        end
        obj.pApplyDirectivityGain = strcmp(obj.SensorGainMeasure,'dBi') && ...
            obj.CombineRadiatedSignals;
        if obj.pUseArray 
            if obj.CombineRadiatedSignals
                obj.cSteeringVector = phased.SteeringVector(...
                    'SensorArray',obj.cSensor,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'IncludeElementResponse',true,...
                    'EnablePolarization',~strcmp(obj.Polarization,'None'));
                if obj.pApplyDirectivityGain
                    sensorElem = getElementHandle(obj.cSensor);
                    if isElementFromAntenna(obj.cSensor) || ...
                            isa(sensorElem,'phased.internal.AntennaAdapter')
                        obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                            'Sensor',obj.cSensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',obj.WeightsInputPort,...
                            'EnablePolarization',~strcmp(obj.Polarization,'None'));
                    else
                        obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                            'Sensor',obj.cSensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',obj.WeightsInputPort,...
                            'EnablePolarization',~strcmp(obj.Polarization,'None'));
                    end
                end
            end
        elseif obj.CombineRadiatedSignals && obj.pApplyDirectivityGain
            if isElementFromAntenna(obj.cSensor) || ...
                    isa(obj.cSensor,'phased.internal.AntennaAdapter')
                obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                    'Sensor',obj.cSensor,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'WeightsInputPort',false,...  % single element weights is just a scalar, separable
                    'EnablePolarization',~strcmp(obj.Polarization,'None'));
            else
                obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                    'Sensor',obj.cSensor,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'WeightsInputPort',false,...  % single element weights is just a scalar, separable
                    'EnablePolarization',~strcmp(obj.Polarization,'None'));
            end
        end
            
        if isa(obj.Sensor,'phased.internal.AbstractSubarray')
            obj.pPowerDistributionFactor = getNumElements(obj.cSensor,...
                1:getNumSubarrays(obj.cSensor));
        else
            obj.pPowerDistributionFactor = 1;
        end
    end
    
    function flag = isInputComplexityLockedImpl(obj,index)
        flag = true; % lock by default
        if obj.pEnableDualPolarizationInput
            if (index == 1) || (index == 2)
                flag = false;
            end
            if obj.WeightsInputPort && (index == 5)
                flag = false;
            end
        elseif ~strcmp(obj.Polarization,'None')
            if index == 1
                flag = false;
            end
            if obj.WeightsInputPort && (index == 4)
                flag = false;
            end
        else
            if index == 1
                flag = false;
            end
            if obj.WeightsInputPort && (index == 3)
                flag = false;
            end
        end
    end

    function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
        flag = false;
    end

    function resetImpl(obj)
        if obj.pUseArray
            if obj.CombineRadiatedSignals
                reset(obj.cSteeringVector);
            else
                reset(obj.cSensor);
            end
        else
            reset(obj.cSensor);
        end
    end
    
    function releaseImpl(obj)
        releaseImpl@phased.internal.AbstractSensorOperation(obj);
        if obj.pUseArray
            if obj.CombineRadiatedSignals
                release(obj.cSteeringVector);
            else
                release(obj.cSensor);
            end
        else
            release(obj.cSensor);
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractSensorOperation(obj);
        if isLocked(obj)
            s.pDOF = obj.pDOF;
            s.cSteeringVector = saveobj(obj.cSteeringVector);
            s.cIntegratedPattern = saveobj(obj.cIntegratedPattern);
            s.pCommonSignalForChannels = obj.pCommonSignalForChannels;
            s.pPowerDistributionFactor = obj.pPowerDistributionFactor;
            s.pApplyDirectivityGain = obj.pApplyDirectivityGain;
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
        if isfield(s,'EnablePolarization')
            if s.EnablePolarization
                obj.Polarization = 'Combined';
            else
                obj.Polarization = 'None';
            end
            s = rmfield(s,'EnablePolarization');
        end
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function y = stepImpl(obj,x,xvArg,angArg,lclaxesArg,wArg,stangArg)
        
        fc = obj.OperatingFrequency;
        if ~strcmp(obj.Polarization,'None')
            if obj.pEnableDualPolarizationInput
                xh = x;
                xv = xvArg;
                ang = angArg;
                lclaxes = lclaxesArg;
            else
                ang = xvArg;
                lclaxes = angArg;
            end
            if obj.pUseArray
                if obj.CombineRadiatedSignals
                    if obj.pNeedSteeringAngle
                        if ~obj.WeightsInputPort
                            if obj.pEnableDualPolarizationInput
                                stang = wArg;
                            else
                                stang = lclaxesArg;
                            end
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc,stang);
                            end
                        else
                            if obj.pEnableDualPolarizationInput
                                stang = stangArg;
                                w = wArg;
                            else
                                stang = wArg;
                                w = lclaxesArg;
                            end
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc,w,stang);
                            end
                        end
                        sv = step(obj.cSteeringVector,fc,ang,stang);
                    else
                        if obj.WeightsInputPort
                            if obj.pEnableDualPolarizationInput
                                w = wArg;
                            else
                                w = lclaxesArg;
                            end
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc,w);
                            end
                        else
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc);
                            end
                        end
                        sv = step(obj.cSteeringVector,fc,ang);
                    end
                    pfactor = 1./sqrt(obj.pPowerDistributionFactor(:));
                    sv_h = bsxfun(@times,pfactor,sv.H);
                    sv_v = bsxfun(@times,pfactor,sv.V);
                    if obj.pApplyDirectivityGain
                        sv_h = phased.internal.normalizeIntegratedPower(sv_h,intresp.H,false);
                        sv_v = phased.internal.normalizeIntegratedPower(sv_v,intresp.V,false);
                    end
                    if obj.WeightsInputPort
                        for m = 1:obj.pDOF
                            sv_h(m,:) = w(m)*sv_h(m,:);
                            sv_v(m,:) = w(m)*sv_v(m,:);
                        end
                    end
                    if obj.pCommonSignalForChannels 
                        svc_h = sum(sv_h,1);
                        svc_v = sum(sv_v,1);
                    else
                        svc_h =  sv_h;
                        svc_v =  sv_v;
                    end
                    if obj.pEnableDualPolarizationInput
                        y_h = xh*svc_h;
                        y_v = xv*svc_v;
                    else
                        y_h = x*svc_h;
                        y_v = x*svc_v;
                    end
                else
                    % Directivity gain doesn't apply here
                    g = step(obj.cSensor,fc,ang);
                    g_h = diag(g.H);
                    g_v = diag(g.V);
                    if obj.WeightsInputPort
                        if obj.pEnableDualPolarizationInput
                            w = wArg;
                        else
                            w = lclaxesArg;
                        end
                        g_h = w.*g_h;
                        g_v = w.*g_v;
                    end
                    if obj.pCommonSignalForChannels
                        if obj.pEnableDualPolarizationInput
                            y_h = xh*g_h.';
                            y_v = xv*g_v.';
                        else
                            y_h = x*g_h.';
                            y_v = x*g_v.';
                        end
                    else
                        if obj.pEnableDualPolarizationInput
                            y_h = complex(xh);
                            y_v = complex(xv);
                        else
                            y_h = complex(x);
                            y_v = complex(x);
                        end
                        for m = 1:obj.pDOF
                            y_h(:,m) = y_h(:,m)*g_h(m);
                            y_v(:,m) = y_v(:,m)*g_v(m);
                        end
                    end
                end
            else
                % Single sensor
                g = step(obj.cSensor,fc,ang);
                g_h = g.H;
                g_v = g.V;
                if obj.WeightsInputPort
                    if obj.pEnableDualPolarizationInput
                        w = wArg;
                    else
                        w = lclaxesArg;
                    end
                    g_h = w*g_h;
                    g_v = w*g_v;
                    if obj.pApplyDirectivityGain
                        intresp = step(obj.cIntegratedPattern,fc,w);
                        g_h = phased.internal.normalizeIntegratedPower(g_h,intresp.H,false);
                        g_v = phased.internal.normalizeIntegratedPower(g_v,intresp.V,false);
                    end
                else
                    if obj.pApplyDirectivityGain
                        intresp = step(obj.cIntegratedPattern,fc);
                        g_h = phased.internal.normalizeIntegratedPower(g_h,intresp.H,false);
                        g_v = phased.internal.normalizeIntegratedPower(g_v,intresp.V,false);
                    end
                end
                if obj.pEnableDualPolarizationInput
                    y_h = xh*g_h.';
                    y_v = xv*g_v.';
                else
                    y_h = x*g_h.';
                    y_v = x*g_v.';
                end
            end
            initYfields = complex(zeros(size(y_h,1),1));
            y = repmat(struct('X',initYfields, 'Y', initYfields,'Z', initYfields), ...
                1, size(ang,2));
                
            for m = size(ang,2):-1:1
                y_temp = sph2cartvec(...
                    [y_h(:,m).';y_v(:,m).';zeros(1,size(x,1))],...
                    ang(1,m),ang(2,m));
                y_temp = phased.internal.local2globalvec(...
                    y_temp,lclaxes);
                y(m).X = y_temp(1,:).';
                y(m).Y = y_temp(2,:).';
                y(m).Z = y_temp(3,:).';
            end
            
        else  % no polarization
            ang = xvArg;
            if obj.pUseArray
                if obj.CombineRadiatedSignals
                    if obj.pNeedSteeringAngle
                        if ~obj.WeightsInputPort
                            stang = angArg;
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc,stang);
                            end
                        else
                            stang = lclaxesArg;
                            w = angArg;
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc,w,stang);
                            end
                        end
                        sv = step(obj.cSteeringVector,fc,ang,stang);
                    else
                        sv = step(obj.cSteeringVector,fc,ang);
                        if obj.WeightsInputPort
                            w = angArg;
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc,w);
                            end
                        else
                            if obj.pApplyDirectivityGain
                                intresp = step(obj.cIntegratedPattern,fc);
                            end
                        end
                    end
                    sv = bsxfun(@times,1./sqrt(obj.pPowerDistributionFactor(:)),sv);
                    if obj.pApplyDirectivityGain
                        sv = phased.internal.normalizeIntegratedPower(sv,intresp,false);
                    end
                    if obj.WeightsInputPort
                        for m = 1:obj.pDOF
                            sv(m,:) = w(m)*sv(m,:);
                        end
                    end
                    if obj.pCommonSignalForChannels 
                       y = x*sum(sv,1);
                    else
                       y = x*sv; 
                    end
                    
                else
                    % Directivity gain does not apply here
                    g_temp = step(obj.cSensor,fc,ang);
                    if isPolarizationEnabled(obj.cSensor)
                        g = diag(hypot(g_temp.H,g_temp.V));
                    else
                        g = diag(g_temp);
                    end
                    if obj.WeightsInputPort
                        w = angArg;
                        g = w.*g;
                    end
                    if obj.pCommonSignalForChannels
                        y = x*g.';
                    else
                        y = complex(zeros(size(x)));
                        for m = 1:obj.pDOF
                            y(:,m) = x(:,m)*g(m);
                        end
                    end
                end
            else
                % Single sensor
                g_temp = step(obj.cSensor,fc,ang);
                if isPolarizationEnabled(obj.cSensor)
                    g = hypot(g_temp.H,g_temp.V);
                else
                    g = g_temp;
                end
                if obj.pApplyDirectivityGain
                    intresp = step(obj.cIntegratedPattern,fc);
                    g = phased.internal.normalizeIntegratedPower(g,intresp,false);
                end
                if obj.WeightsInputPort
                    w = angArg;
                    g = w*g;
                end
                y = x*g.';
            end
        end
        
    end
           
end

methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl
      groups = getPropertyGroupsImpl@phased.internal.AbstractSensorOperation;
      
      pCombineRadiatedSignals = ...
          matlab.system.display.internal.Property('CombineRadiatedSignals', ...
                   'IsGraphical', false);
      dPolarization = ...
                matlab.system.display.internal.Property('Polarization', ...
                                                        'IsGraphical', false);
      props = {'OperatingFrequency',...
                pCombineRadiatedSignals,...
                'SensorGainMeasure',...
                dPolarization,...
                'WeightsInputPort'};
      groups(1).PropertyList = [groups(1).PropertyList props];
      action = matlab.system.display.Action(@phased.Radiator.onActionCalled, ...
          'Label', 'Analyze');
      matlab.system.display.internal.setCallbacks(action, ...
          'DialogAppliedFcn', @phased.Radiator.onDialogApplied, ...
          'SystemDeletedFcn', @phased.Radiator.onSystemDeleted, ...
          'IsEnabledFcn', @phased.Radiator.isActionEnabled);
      groups(2).Actions = action;
  end
  function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:RadiatorTitle')),...
          'Text',getString(message('phased:library:block:RadiatorDesc')));
  end
end
methods (Static,Hidden)
    function onActionCalled(actionData, sysObj)
        if ~phased.apps.internal.SensorArrayViewer.SettingsDialog.verifySensorArray(sysObj.Sensor)
            return;
        end
        f = actionData.UserData;
        if isempty(f) || isStale(f)
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
   
    function isEnabled = isActionEnabled(systemHandle)
        % Return whether action enabled based on Sensor type
        
        sensorObj = phased.internal.AbstractSensorOperation.getSensorObject(systemHandle);
        isEnabled = ~isempty(sensorObj) && ~isa(sensorObj, 'phased.internal.AbstractSubarray');
    end
end 

methods (Access = protected) %for Simulink
    function varargout = getOutputNamesImpl(~)
        varargout = {''};
    end
    function varargout = getInputNamesImpl(obj)
        if obj.pEnableDualPolarizationInput
            varargout = {'Xh','Xv','Ang'};
            lastIdx = 4;
        else
            varargout = {'X','Ang'};
            lastIdx = 3;
        end
        if ~strcmp(obj.Polarization,'None')
            varargout(lastIdx) = {'LAxes'};
            lastIdx = lastIdx+1;
        end
        if obj.WeightsInputPort
            varargout(lastIdx) = {'W'};
            lastIdx = lastIdx+1;
        end
        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmpi(obj.Sensor.SubarraySteering,'None',1)
            if strncmp(obj.Sensor.SubarraySteering,'Phase',1) || ...
                    strncmp(obj.Sensor.SubarraySteering,'Time',1)
                varargout(lastIdx) = {'Steer'};
            else
                varargout(lastIdx) = {'WS'};
            end
        end
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Narrowband\n Tx Array');
    end        
    function varargout = getOutputSizeImpl(obj)
        szX = propagatedInputSize(obj,1);
        szAng = propagatedInputSize(obj,2);
        %Output is number of snapshots * number of angles
        varargout{1} = [szX(1) szAng(2)];
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

