classdef (Sealed,StrictDefaults) WidebandCollector < phased.internal.AbstractSensorOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime
%WidebandCollector Wideband signal collector
%   H = phased.WidebandCollector creates a wideband signal collector System
%   object, H. The object collects incident wideband signals from given
%   directions using a sensor array or a single element.
%
%   H = phased.WidebandCollector(Name,Value) creates a wideband signal
%   collector object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax when the Wavefront property is 'Plane':
%
%   Y = step(H,X,ANG) collects multiple plane wave signals specified in
%   X with all collecting elements, adds them together, and returns the
%   combined output Y. ANG contains the signal incident directions. For
%   each plane wave signal, the received signal is collected by decomposing
%   the signal into multiple subbands and using the phase approximation of
%   the time delays across collecting elements in the far field for each
%   subband.
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
%   system in LAXES when you set the Polarization property to 'Combined'.
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
%   Y = step(H,X,ANG,WS) uses WS as the weights applied to each element in
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
%   [YH,YV] = step(...) outputs YH from H polarization port and YV from V
%   polarization port. YH and YV are matrices whose number of columns
%   equals the number of subarrays if Sensor contains subarrays, or the
%   number of elements otherwise. Each column of YH or YV contains the
%   output of the corresponding element/subarray in response to all the
%   signals in X. This syntax only applies when you set the Polarization
%   property to 'Dual'.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [YH,YV] = step(H,X,ANG,LAXES,W,STEER)
%
%   or
%
%   [YH,YV] = step(H,X,ANG,LAXES,W,WS)
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
%   system in LAXES when you set the Polarization property to 'Combined'.
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
%   [YH,YV] = step(...) outputs YH from H polarization port and YV from V
%   polarization port. YH and YV are matrices whose number of columns
%   equals the number of subarrays if Sensor contains subarrays, or the
%   number of elements otherwise. Each column of YH or YV contains the
%   output of the corresponding element/subarray in response to all the
%   signals in X. This syntax only applies when you set the Polarization
%   property to 'Dual'.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [YH,YV] = step(H,X,ANG,LAXES,W)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   WidebandCollector methods:
%
%   step     - Collect signals (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a wideband collector object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of wideband signal collector object
%
%   WidebandCollector properties:
%
%   Sensor                       - Handle of the sensor
%   PropagationSpeed             - Propagation speed
%   SampleRate                   - Sample rate
%   ModulatedInput               - Assume modulated input
%   CarrierFrequency             - Carrier frequency
%   NumSubbands                  - Number of subbands
%   Wavefront                    - Type of incoming wavefront
%   SensorGainMeasure            - Sensor gain measure 
%   Polarization                 - Polarization configuration
%   WeightsInputPort             - Enable weights input
%
%   % Examples:
%
%   % Example 1: 
%   %   Collect signal with a single antenna.
%
%   array = phased.IsotropicAntennaElement;
%   sigcol = phased.WidebandCollector('Sensor',array);
%   x = [1;1];
%   incidentAngle = [10 30]';
%   y = sigcol(x,incidentAngle);
%
%   % Example 2: 
%   %   Collect a far field signal with a 5-element array.
%
%   array = phased.ULA('NumElements',5);
%   sigcol = phased.WidebandCollector('Sensor',array);
%   x = [1;1];
%   incidentAngle = [10 30]';
%   y = sigcol(x,incidentAngle);
%
%   % Example 3: 
%   %   Collect signal with a 3-element antenna array. Each antenna 
%   %   collects a separate input signal from a separate direction.
%
%   array = phased.ULA('NumElements',3);
%   sigcol = phased.WidebandCollector('Sensor',array,...
%               'Wavefront','Unspecified');
%   x = rand(10,3);   % Each column is a separate signal for one element
%   incidentAngle = [10 0; 20 5; 45 2]'; % 3 angles for 3 signals
%   y = sigcol(x,incidentAngle);
%
%   See also phased, phased.Collector.

%   Copyright 2009-2017 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a positive scalar. The
        %   default value is 1e6 (1 MHz).
        SampleRate = 1e6
        %CarrierFrequency Carrier frequency (Hz)
        %   Specify the carrier frequency (in Hz) as a positive scalar. The
        %   default value is 1e9 (1 GHz). This property applies when the
        %   ModulatedInput property is true.
        CarrierFrequency = 1e9
    end
    
    properties (Nontunable, PositiveInteger)
        %NumSubbands    Number of subbands
        %   Specify the number of subbands used in the subband processing
        %   as a positive integer. The default value of this property is
        %   64.
        NumSubbands = 64
    end
    
    properties (Nontunable)
        %Wavefront Type of incoming wavefront
        %   Specify the type of incoming wavefront as one of 'Plane' |
        %   'Unspecified', where the default is 'Plane'. When you set the
        %   Wavefront property to 'Plane', the input signals are assumed to
        %   be multiple plane waves impinging on the entire array. Each
        %   plane wave is received by all collecting elements. If you set
        %   the Wavefront property to 'Unspecified', the input signals are
        %   assumed to be individual waves impinging on individual sensors.
        Wavefront = 'Plane'
    end

    properties (Nontunable)
        %SensorGainMeasure   Sensor gain measure
        %   Specify how sensor gain is measured as one of 'dB' | 'dBi',
        %   where 'dB' is the default value. This property only applies
        %   when you set the Wavefront property to 'Plane' and the
        %   ModulatedInput property to true.
        %
        %   When you set this property to 'dB', the input signal power is
        %   scaled by the sensor's power pattern (in dB) at the
        %   corresponding direction and then combined. When you set this
        %   property to 'dBi', the input signal power is scaled by the
        %   directivity pattern (in dBi) at the corresponding direction and
        %   then combined. The dBi option is helpful when you want to
        %   compare results with those predicted by the radar equation that
        %   uses dBi to specify the antenna gain. Note that the computation
        %   for the dBi option is expensive as it requires an integration
        %   over all directions to compute the total radiated power of the
        %   sensor.
        SensorGainMeasure = 'dB'
    end
    
    properties (Nontunable, Logical) 
        %ModulatedInput Assume modulated input
        %   Set this property to true to indicate the input signal is
        %   demodulated at a carrier frequency. The default value is true. 
        ModulatedInput = true
    end

    properties (Nontunable)
        %Polarization   Polarization configuration
        %   Specify the collector's polarization configuration as one of
        %   'None' | 'Combined' | 'Dual', where the default is 'None'. When
        %   you set this property to 'None', the signal is considered as
        %   polarization independent. When you set this property to
        %   'Combined', the collected signal contains the polarization
        %   information reflecting the sensor's native polarization. When
        %   you set this property to 'Dual', the H and V polarization of
        %   the sensor can be retrieved independently.
        Polarization = 'None'
    end

    properties (Dependent, Hidden)
        %EnablePolarization  Enable polarization
        %   Set this property to true to enable polarization. Set this
        %   property to false to ignore polarization. The default value of
        %   this property is false. This property applies when the sensor
        %   specified in the Sensor property is capable of simulating
        %   polarization.
        EnablePolarization 
    end
    
    properties (Access = private, Dependent)
        %EnableDualPolarizationOutput  Enable dual polarization output 
        %signal
        %   Set this property to true to enable dual polarization output
        %   signal in separate H and V polarization. Set this property to
        %   false to output one signal with collapsed polarization. The
        %   default value of this property is false. This property applies
        %   when you set the EnablePolarization property to true.
        pEnableDualPolarizationOutput
    end

    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from Simulink time engine. Set SampleRateFromInputCheckbox to
        %   false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true
    end
    
    properties(Constant, Hidden)
        WavefrontSet = matlab.system.StringSet({'Plane','Unspecified'});
        PolarizationSet = matlab.system.StringSet({'None','Combined','Dual'});
        SensorGainMeasureSet = matlab.system.StringSet({'dB','dBi'});
    end

    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end

    properties (Access = private, Nontunable)
        pFFTLength      % number of FFT
        pSubbandFreqs   % center frequency for each FFT bin
        cSteeringVector % handle to steering vector 
        cIntegratedPattern
        pDOF
        pSampleRate
        pApplyDirectivityGain
    end
    
    properties (Access = private)
        pBuffer         % buffer for overlap-add
        pNumOverlap     % number of overlapped samples in overlap-add
    end
    
    properties (Access = private)
        cSubbandDivider
        cSubbandCombiner
    end
    
    methods
        function set.SampleRate(obj,value)
            validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'SampleRate');
            obj.SampleRate = value;
        end
        function set.CarrierFrequency(obj,value)
            validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'CarrierFrequency');
            obj.CarrierFrequency = value;
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

        function value = get.pEnableDualPolarizationOutput(obj)
            value = strcmp(obj.Polarization,'Dual');
        end
    end
    
    methods
        function obj = WidebandCollector(varargin)
            obj@phased.internal.AbstractSensorOperation(varargin{:});
        end
    end

    methods (Access = protected)

        function num = getNumInputsImpl(obj)
            num = getNumInputsImpl@phased.internal.AbstractSensorOperation(obj);
            if ~strcmp(obj.Polarization,'None')
                num = num+1;
            end
        end
    
        function num = getNumOutputsImpl(obj)
            if obj.pEnableDualPolarizationOutput
                num = 2;
            else
                num = 1;
            end
        end

        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractSensorOperation(obj);
            cond = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    (obj.Wavefront(1) == 'U'); %Unspecified
            if cond
                coder.internal.errorIf(cond,'phased:phased:collector:invalidSettingsForSubarray','Wavefront','Plane','Sensor');
            end
            
            cond = ~strcmp(obj.Polarization,'None') && ~isPolarizationCapable(obj.Sensor);
            if cond
                coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
            end
        end

        function validateInputsImpl(obj,x,angle,lclaxes,weightsArg,stangArg)
            coder.extrinsic('mat2str');
            polflag = ~strcmp(obj.Polarization,'None');
            if polflag
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
                
                % realign x size to real signal length
                xsize = [size(x(1).X,1) size(x,2)];
                %xsize = size([x.X]);  
            else
                cond = ~isa(x,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','X','double');
                end
                cond = ~obj.ModulatedInput && ~isreal(x);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:WidebandCollector:InvalidComplexInput');
                end
                xsize = size(x);
            end
            UseArrayFlag = isa(obj.Sensor,'phased.internal.AbstractArray') || ....
                isa(obj.Sensor,'phased.internal.AbstractSubarray');
            if UseArrayFlag
                M = getNumElements(obj.Sensor);
            else
                M = 1;
            end
            cond = UseArrayFlag && ~(obj.Wavefront(1) == 'P') && xsize(2)~=M; %Plane
            if cond
                coder.internal.errorIf(cond,'phased:phased:collector:DimensionMismatch');
            end
            
            validateNumChannels(obj,x);
            
            % angle
            angsize = size(angle);
            cond = xsize(2) ~= angsize(2);
            if cond
                coder.internal.errorIf(cond,'phased:phased:collector:AngleDimensionMismatch');
            end
            cond = ~isa(angle,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Ang','double');
            end
            cond = ~isreal(angle);
            if cond
                coder.internal.errorIf(cond,'phased:step:NeedReal', 'Ang');
            end
            
            % local axes
            if polflag
                % check lclAxes is 3x3
                sigdatatypes.validate3DCartCoord(lclaxes,'','LAxes',...
                    {'size',[3 3]});
            end

            % steering angle for subarray
            if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    ~(obj.Sensor.SubarraySteering(1) == 'N') %None
                if ~obj.WeightsInputPort
                    if polflag
                        stang = weightsArg;
                    else
                        stang = lclaxes;
                    end
                else
                    if ~polflag
                        stang = weightsArg;
                    else
                        stang = stangArg;
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
                if ~polflag
                    weights = lclaxes;
                else
                    weights = weightsArg;
                end
                if UseArrayFlag
                    Mw = getDOF(obj.Sensor);
                else
                    Mw = 1;
                end
                cond = ~isa(weights,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','W','double');
                end
                wsize = size(weights);
                cond = ~isequal(wsize, [Mw 1]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions','W',...
                        coder.const(mat2str([Mw 1])), ...
                        coder.const(mat2str(wsize)));
                end
            end
            
        end

        function setupImpl(obj,x,ang,lclaxes,w,stang) %#ok<INUSD>
            setupImpl@phased.internal.AbstractSensorOperation(obj);
            
            % obj.pSampleRate = getSampleRate(obj,xsize(1),1,obj.SampleRate);
            fs = obj.SampleRate; % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            obj.pSampleRate = fs;
            
            cond = obj.ModulatedInput && ...
                    (obj.pSampleRate > 2*obj.CarrierFrequency);
            if cond
                coder.internal.errorIf(cond,'phased:phased:WidebandCollector:FsTooHigh');
            end
            
            % Determine the FFT length
            nfft = obj.NumSubbands;
            obj.pFFTLength = nfft;
            
            obj.cSubbandDivider = phased.internal.SubbandDivider(...
                'OperatingFrequency',obj.CarrierFrequency,...
                'SampleRate',obj.pSampleRate,...
                'NumSubbands',nfft,'EnableWarning',false);
            obj.cSubbandCombiner = phased.internal.SubbandCombiner(...
                'NumSubbands',nfft,'TimeSignalLengthSource','Inherit',...
                'EnableWarning',false);
            
            % Pre-calculate the center frequency of each FFT bin
            if obj.ModulatedInput
                obj.pSubbandFreqs = phased.internal.subbandCenterFrequency(...
                    obj.CarrierFrequency,obj.pSampleRate,obj.pFFTLength).';
            else
                obj.pSubbandFreqs = phased.internal.subbandCenterFrequency(...
                    0,obj.pSampleRate,obj.pFFTLength).';
            end
            if obj.pUseArray
                obj.pDOF = getDOF(obj.cSensor);
            else
                obj.pDOF = 1;
            end
            
            obj.pApplyDirectivityGain = strcmp(obj.Wavefront,'Plane') && ...
                strcmp(obj.SensorGainMeasure,'dBi') && obj.ModulatedInput;
            if obj.pUseArray 
                if obj.Wavefront(1) == 'P' %Plane
                    if obj.ModulatedInput
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
                    else
                        % Directivity gain doesn't apply here
                        obj.cSteeringVector = phased.SteeringVector(...
                            'SensorArray',obj.cSensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'IncludeElementResponse',false);
                    end
                end
            else
                if obj.ModulatedInput && obj.pApplyDirectivityGain
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
            end
            
            if obj.ModulatedInput
                obj.pBuffer = complex(zeros(nfft,obj.pDOF));
            else
                obj.pBuffer = zeros(nfft,obj.pDOF);
            end
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        end

        function processInputSizeChangeImpl(obj,x,~,~,~,~)
            if ~strcmp(obj.Polarization,'None')
                xsize = [size(x(1).X,1) size(x,2)];
            else
                xsize = size(x);
            end
            obj.pNumOverlap = obj.pFFTLength - xsize(1);
        end

        function flag = isInputComplexityLockedImpl(obj,index)
            flag = true;
            if index == 1
                if obj.ModulatedInput
                    flag = false;
                else
                    flag = true;
                end
            end
            if obj.WeightsInputPort && ...
                    ( (~strcmp(obj.Polarization,'None') && (index == 4)) || ...
                    ( strcmp(obj.Polarization,'None') && (index == 3)))
                    flag = false;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function resetImpl(obj)
            if obj.pUseArray
                if (obj.Wavefront(1) == 'P') %'Plane'
                    reset(obj.cSteeringVector);
                    if ~obj.ModulatedInput
                        reset(obj.cSensor);
                    end
                else
                    reset(obj.cSensor);
                end
            else
                reset(obj.cSensor);
            end

            obj.pBuffer = zeros(size(obj.pBuffer),'like',obj.pBuffer);
            reset(obj.cSubbandDivider);
            reset(obj.cSubbandCombiner);
        end

        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSensorOperation(obj);
            if obj.pUseArray
                if (obj.Wavefront(1) == 'P') %'Plane'
                    release(obj.cSteeringVector);
                    if ~obj.ModulatedInput
                        release(obj.cSensor);
                    end
                else
                    release(obj.cSensor);
                end
            else
                release(obj.cSensor);
            end
            release(obj.cSubbandDivider);
            release(obj.cSubbandCombiner);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSensorOperation(obj);
            if isLocked(obj)
                s.pDOF = obj.pDOF;
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.cSubbandDivider = saveobj(obj.cSubbandDivider);
                s.cSubbandCombiner = saveobj(obj.cSubbandCombiner);
                s.cIntegratedPattern = saveobj(obj.cIntegratedPattern);
                s.pFFTLength = obj.pFFTLength;
                s.pSubbandFreqs = obj.pSubbandFreqs;
                s.pBuffer = obj.pBuffer;
                s.pNumOverlap = obj.pNumOverlap;
                s.pSampleRate = obj.pSampleRate;
                s.pApplyDirectivityGain = obj.pApplyDirectivityGain;
            end
        end

        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSensorOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    if isfield(s,'cSubbandDivider')
                        obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                        s = rmfield(s,'cSteeringVector');
                        obj.cSubbandDivider = phased.internal.SubbandDivider.loadobj(s.cSubbandDivider);
                        s = rmfield(s,'cSubbandDivider');
                        obj.cSubbandCombiner = phased.internal.SubbandDivider.loadobj(s.cSubbandCombiner);
                        s = rmfield(s,'cSubbandCombiner');
                    else
                        obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                        s = rmfield(s,'cSteeringVector');
                        if isfield(s,'cFFT')
                            s = rmfield(s,'cFFT');
                        end
                        if isfield(s,'cIFFT')
                            s = rmfield(s,'cIFFT');
                        end
                        % recover locked sample rate information
                        if isfield(s,'pSampleRate')
                            obj.pSampleRate = s.pSampleRate;
                            s = rmfield(s,'pSampleRate');
                        else
                            obj.pSampleRate = s.SampleRate;
                        end
                        obj.cSubbandDivider = phased.internal.SubbandDivider(...
                            'OperatingFrequency',s.CarrierFrequency,...
                            'SampleRate',obj.pSampleRate,...
                            'NumSubbands',s.pFFTLength,'EnableWarning',false);
                        obj.cSubbandCombiner = phased.internal.SubbandCombiner(...
                            'NumSubbands',s.pFFTLength,...
                            'TimeSignalLength',s.pFFTLength-s.pNumOverlap,'EnableWarning',false);
                    end
                end
                if isfield(s,'cIntegratedPattern') && ...
                            isfield(s.cIntegratedPattern,'ClassNameForLoadTimeEval')
                    obj.cIntegratedPattern = eval(...
                        sprintf('%s.loadobj(s.cIntegratedPattern)',s.cIntegratedPattern.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cIntegratedPattern');
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

        function [y,yv] = stepImpl(obj,x,ang,lclaxes,w,stangArg)

            if ~strcmp(obj.Polarization,'None')
                % Modulate each subband with corresponding phase in frequency
                % domain
                if obj.pUseArray
                    if (obj.Wavefront(1) == 'P') %Plane
                        if obj.pNeedSteeringAngle
                            if ~obj.WeightsInputPort
                                stang = w;
                                if obj.pEnableDualPolarizationOutput
                                    [y_temp_h,y_temp_v] = farFieldFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes,stang);
                                else
                                    y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,stang);
                                end
                            else
                                wt = w;
                                stang = stangArg;
                                if obj.pEnableDualPolarizationOutput
                                    [y_temp_h,y_temp_v] = farFieldFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes,wt,stang);
                                else
                                    y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,wt,stang);
                                end
                            end
                        else
                            if obj.WeightsInputPort
                                wt = w;
                                if obj.pEnableDualPolarizationOutput
                                    [y_temp_h,y_temp_v] = farFieldFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes,wt);
                                else
                                    y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,wt);
                                end
                            else
                                if obj.pEnableDualPolarizationOutput
                                    [y_temp_h,y_temp_v] = farFieldFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes);
                                else
                                    y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes);
                                end
                            end
                        end
                    else
                        % Directivity gain doesn't apply here
                        if obj.WeightsInputPort
                            wt = w;
                        end
                        if obj.pEnableDualPolarizationOutput
                            [y_temp_h,y_temp_v] = independentArrayChannelFreqDomainDualPolarizedModulate(...
                                obj,x,ang,lclaxes);
                        else
                            y_temp = independentArrayChannelFreqDomainPolarizedModulate(...
                                obj,x,ang,lclaxes);
                        end
                    end
                else % single element
                    if obj.WeightsInputPort
                        wt = w;
                        if obj.pEnableDualPolarizationOutput
                            [y_temp_h,y_temp_v] = singleElementFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes,wt);
                        else
                            y_temp = singleElementFreqDomainPolarizedModulate(obj,x,ang,lclaxes,wt);
                        end
                    else
                        if obj.pEnableDualPolarizationOutput
                            [y_temp_h,y_temp_v] = singleElementFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes);
                        else
                            y_temp = singleElementFreqDomainPolarizedModulate(obj,x,ang,lclaxes);
                        end
                    end
                end
                
                % Convert back to time domain
                if obj.pEnableDualPolarizationOutput
                    y_out_h = step(obj.cSubbandCombiner,complex(y_temp_h),complex(x(1).X));
                    y_out_v = step(obj.cSubbandCombiner,complex(y_temp_v),complex(x(1).X));
                else
                    y_out = step(obj.cSubbandCombiner,complex(y_temp),complex(x(1).X));
                end
                xlen = size(x(1).X,1);
                
            else % no polarization
                % Modulate each subband with corresponding phase in frequency
                % domain
                if obj.pUseArray
                    if (obj.Wavefront(1) == 'P') %Plane
                        if obj.pNeedSteeringAngle
                            if ~obj.WeightsInputPort
                                stang = lclaxes;
                                y_temp = farFieldFreqDomainModulate(obj,x,ang,stang);
                            else
                                stang = w;
                                wt = lclaxes;
                                y_temp = farFieldFreqDomainModulate(obj,x,ang,wt,stang);
                            end
                        else
                            if obj.WeightsInputPort
                                wt = lclaxes;
                                y_temp = farFieldFreqDomainModulate(obj,x,ang,wt);
                            else
                                y_temp = farFieldFreqDomainModulate(obj,x,ang);
                            end
                        end
                    else
                        % Directivity gain doesn't apply here
                        if obj.WeightsInputPort
                            wt = lclaxes;
                        end
                        y_temp = independentArrayChannelFreqDomainModulate(...
                            obj,x,ang);
                    end
                else % single element
                    if obj.WeightsInputPort
                        wt = lclaxes;
                        y_temp = singleElementFreqDomainModulate(obj,x,ang,wt);
                    else
                        y_temp = singleElementFreqDomainModulate(obj,x,ang);
                    end
                end

                % Convert back to time domain
                y_out = step(obj.cSubbandCombiner,complex(y_temp),complex(x));
                xlen = size(x,1);
            end
            
            if ~obj.ModulatedInput
                if obj.pEnableDualPolarizationOutput
                    y = real(y_out_h(1:xlen,:));
                    yv = real(y_out_v(1:xlen,:));
                else
                    y = real(y_out(1:xlen,:));
                end
            else
                if obj.pEnableDualPolarizationOutput
                    y = complex(y_out_h(1:xlen,:));
                    yv = complex(y_out_v(1:xlen,:));
                else
                    y = complex(y_out(1:xlen,:));
                end
            end
            %y = overlapadd(obj,y,obj.pNumOverlap);

            % Apply weights
            if obj.WeightsInputPort
                if obj.pUseArray
                    for m = 1:obj.pDOF
                        if obj.pEnableDualPolarizationOutput
                            y(:,m) = y(:,m)*wt(m);
                            yv(:,m) = yv(:,m)*wt(m);
                        else
                            y(:,m) = y(:,m)*wt(m);
                        end
                    end
                elseif ~obj.pApplyDirectivityGain  % single sensor, no directivity gain
                    for m = 1:obj.pDOF
                        if obj.pEnableDualPolarizationOutput
                            y(:,m) = y(:,m)*wt(m);
                            yv(:,m) = yv(:,m)*wt(m);
                        else
                            y(:,m) = y(:,m)*wt(m);
                        end
                    end
                end
            end

        end

        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if ~obj.ModulatedInput && strcmp(prop, 'CarrierFrequency')
                flag = true;
            end
            if (~strcmp(obj.Wavefront,'Plane') || ~obj.ModulatedInput) && ...
                    strcmp(prop,'SensorGainMeasure')
                flag = true;
            end
        end

    end

    methods (Access = private)

        function phaseModulatedXFreq = ...
                farFieldFreqDomainModulate(obj,x,ang,wArg,stangArg)
        % collect far field signals and sum them together

            % Calculate steering vector for different subbands and angles
            % in one step. In case of unmodulated input, because there are
            % negative frequencies, current steering vector does not
            % support it when including the element response, therefore, we
            % get the steering vector without element response first and
            % then add in element response manually.
            subfreq = obj.pSubbandFreqs;
            
            if obj.pNeedSteeringAngle
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;  
                        stang = stangArg;
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,w,stang);
                        end
                    else
                        stang = wArg;
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,stang);
                        end
                    end
                    phaseMtx = step(obj.cSteeringVector,subfreq,ang,stang); % chan x ang x freq
                    if obj.pApplyDirectivityGain
                        for m = 1:numel(subfreq)
                            phaseMtx(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx(:,:,m),intresp(m),false);
                        end
                    end
                else
                    % Directivity gain doesn't apply here
                    negidx = find(subfreq<0);
                    phaseMtx = step(obj.cSteeringVector,abs(subfreq),ang);  % chan x ang x freq
                    phaseMtx(:,:,negidx) = conj(phaseMtx(:,:,negidx));      % chan x ang x freq
                end
            else
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;  
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,w);
                        end
                    else
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq);
                        end
                    end
                    phaseMtx = step(obj.cSteeringVector,subfreq,ang);       % chan x ang x freq
                    if obj.pApplyDirectivityGain
                        for m = 1:numel(subfreq)
                            phaseMtx(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx(:,:,m),intresp(m),false);
                        end
                    end
                else
                    % Directivity gain doesn't apply here
                    negidx = find(subfreq<0);
                    phaseMtx = step(obj.cSteeringVector,abs(subfreq),ang);  % chan x ang x freq
                    phaseMtx(:,:,negidx) = conj(phaseMtx(:,:,negidx));      % chan x ang x freq
                end
            end
            if ~obj.ModulatedInput
                if obj.pNeedSteeringAngle
                    if obj.WeightsInputPort
                        stang = stangArg;
                    else
                        stang = wArg;
                    end
                    if isPolarizationEnabled(obj.cSensor)
                        resp_temp = step(obj.cSensor,abs(subfreq),ang,obj.PropagationSpeed,stang);
                        phaseMtx = phaseMtx.*...
                            hypot(resp_temp.H,resp_temp.V);
                    else
                        phaseMtx = phaseMtx.*...
                            step(obj.cSensor,abs(subfreq),ang,obj.PropagationSpeed,stang);
                    end
                else
                    if isPolarizationEnabled(obj.cSensor)
                        resp_temp = step(obj.cSensor,abs(subfreq),ang);
                        phaseMtx = phaseMtx.*...
                            hypot(resp_temp.H,resp_temp.V);
                    else
                        phaseMtx = phaseMtx.*...
                            step(obj.cSensor,abs(subfreq),ang);
                    end
                end
            end

            % Translate to frequency domain
            xtempfreq = step(obj.cSubbandDivider,complex(x));  % snap x ang x freq
            NumChannel = obj.pDOF;

            % Modulate each channel at different subbands and angles
            phaseModulatedXFreq = complex(zeros(size(xtempfreq,1),NumChannel,size(xtempfreq,3))); % snap x chan x freq
            for m = 1:size(xtempfreq,1)
                phaseModulatedXFreq(m,:,:) = permute(sum(bsxfun(@times,...
                    xtempfreq(m,:,:),phaseMtx),2),[2 1 3]);
            end
            
        end

        function phaseModulatedXFreq = ...
                farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,wArg,stangArg)
        % collect far field signals and sum them together

            % Calculate steering vector for different subbands and angles
            % in one step. In case of unmodulated input, because there are
            % negative frequencies, current steering vector does not
            % support it when including the element response, therefore, we
            % get the steering vector without element response first and
            % then add in element response manually.
            subfreq = obj.pSubbandFreqs;
            
            if obj.pNeedSteeringAngle
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;  
                        stang = stangArg;
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,w,stang);
                        end
                    else
                        stang = wArg;
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,stang);
                        end
                    end
                    phaseMtx_temp = step(obj.cSteeringVector,subfreq,ang,stang); % chan x ang x freq
                    if obj.pApplyDirectivityGain
                        for m = 1:numel(subfreq)
                            phaseMtx_temp.H(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.H(:,:,m),intresp.H(m),false);
                            phaseMtx_temp.V(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.V(:,:,m),intresp.V(m),false);
                        end
                    end
                else
                    % Directivity gain doesn't apply here
                    negidx = find(subfreq<0);
                    phaseMtx_temp = step(obj.cSteeringVector,abs(subfreq),ang);
                    phaseMtx_temp(:,:,negidx) = conj(phaseMtx_temp(:,:,negidx)); % chan x ang x freq
                end
            else
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;  
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,w);
                        end
                    else
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq);
                        end
                    end
                    phaseMtx_temp = step(obj.cSteeringVector,subfreq,ang);       % chan x ang x freq
                    if obj.pApplyDirectivityGain
                        for m = 1:numel(subfreq)
                            phaseMtx_temp.H(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.H(:,:,m),intresp.H(m),false);
                            phaseMtx_temp.V(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.V(:,:,m),intresp.V(m),false);
                        end
                    end
                else
                    % Directivity gain doesn't apply here
                    negidx = find(subfreq<0);
                    phaseMtx_temp = step(obj.cSteeringVector,abs(subfreq),ang);  % chan x ang x freq
                    phaseMtx_temp(:,:,negidx) = conj(phaseMtx_temp(:,:,negidx)); % chan x ang x freq 
                end
            end
            if ~obj.ModulatedInput
                if obj.pNeedSteeringAngle
                    if obj.WeightsInputPort
                        stang = stangArg;
                    else
                        stang = wArg;
                    end
                    sensor_resp = step(obj.cSensor,abs(subfreq),ang,obj.PropagationSpeed,stang);
                else
                    sensor_resp = step(obj.cSensor,abs(subfreq),ang);
                end
                phaseMtx.H = phaseMtx_temp.*sensor_resp.H;                     % chan x ang x freq
                phaseMtx.V = phaseMtx_temp.*sensor_resp.V;
            else
                phaseMtx = phaseMtx_temp;                                      % chan x ang x freq
            end
            
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end

            % Translate to frequency domain
            xtempfreq_x = step(obj.cSubbandDivider,complex(x_X));                       % snap x ang x freq
            xtempfreq_y = step(obj.cSubbandDivider,complex(x_Y));
            xtempfreq_z = step(obj.cSubbandDivider,complex(x_Z));
            NumChannel = obj.pDOF;
            NumAngles = size(ang,2);
            
            nfft = obj.pFFTLength;
            phaseModulatedXFreq = complex(zeros(size(xtempfreq_x,1),NumChannel,nfft));  % snap x ang x freq
            
            for m = 1:NumAngles
                
                temp_phaseMtx_h = squeeze(phaseMtx.H(:,m,:));
                temp_phaseMtx_v = squeeze(phaseMtx.V(:,m,:));
                
                arrayEffect_temp = sph2cartvec(...
                    [temp_phaseMtx_h(:).';temp_phaseMtx_v(:).';zeros(1,NumChannel*nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp = phased.internal.local2globalvec(...
                    arrayEffect_temp,lclaxes);
                
                arrayEffect_x = reshape(arrayEffect_temp(1,:).',NumChannel,[]).';     % freq x chan
                arrayEffect_y = reshape(arrayEffect_temp(2,:).',NumChannel,[]).';
                arrayEffect_z = reshape(arrayEffect_temp(3,:).',NumChannel,[]).';
                
                phaseModulatedXFreq = phaseModulatedXFreq + permute(...
                    bsxfun(@times,permute(arrayEffect_x,[3 1 2]),reshape(xtempfreq_x(:,m,:),[],nfft)) + ... % snap x freq x chan
                    bsxfun(@times,permute(arrayEffect_y,[3 1 2]),reshape(xtempfreq_y(:,m,:),[],nfft)) + ...
                    bsxfun(@times,permute(arrayEffect_z,[3 1 2]),reshape(xtempfreq_z(:,m,:),[],nfft)),...
                    [1 3 2]);       % snap x chan x freq
                
            end
            
        end

        function [phaseModulatedXFreq_h,phaseModulatedXFreq_v] = ...
                farFieldFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes,wArg,stangArg)
        % collect far field signals and sum them together

            % Calculate steering vector for different subbands and angles
            % in one step. In case of unmodulated input, because there are
            % negative frequencies, current steering vector does not
            % support it when including the element response, therefore, we
            % get the steering vector without element response first and
            % then add in element response manually.
            subfreq = obj.pSubbandFreqs;
            
            if obj.pNeedSteeringAngle
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;  
                        stang = stangArg;
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,w,stang);
                        end
                    else
                        stang = wArg;
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,stang);
                        end
                    end
                    phaseMtx_temp = step(obj.cSteeringVector,subfreq,ang,stang); % chan x ang x freq
                    if obj.pApplyDirectivityGain
                        for m = 1:numel(subfreq)
                            phaseMtx_temp.H(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.H(:,:,m),intresp.H(m),false);
                            phaseMtx_temp.V(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.V(:,:,m),intresp.V(m),false);
                        end
                    end
                else
                    % Directivity gain doesn't apply here
                    negidx = find(subfreq<0);
                    phaseMtx_temp = step(obj.cSteeringVector,abs(subfreq),ang);
                    phaseMtx_temp(:,:,negidx) = conj(phaseMtx_temp(:,:,negidx)); % chan x ang x freq
                end
            else
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;  
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq,w);
                        end
                    else
                        if obj.pApplyDirectivityGain
                            intresp = step(obj.cIntegratedPattern,subfreq);
                        end
                    end
                    phaseMtx_temp = step(obj.cSteeringVector,subfreq,ang);       % chan x ang x freq
                    if obj.pApplyDirectivityGain
                        for m = 1:numel(subfreq)
                            phaseMtx_temp.H(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.H(:,:,m),intresp.H(m),false);
                            phaseMtx_temp.V(:,:,m) = phased.internal.normalizeIntegratedPower(phaseMtx_temp.V(:,:,m),intresp.V(m),false);
                        end
                    end
                else
                    % Directivity gain doesn't apply here
                    negidx = find(subfreq<0);
                    phaseMtx_temp = step(obj.cSteeringVector,abs(subfreq),ang);  % chan x ang x freq
                    phaseMtx_temp(:,:,negidx) = conj(phaseMtx_temp(:,:,negidx)); % chan x ang x freq 
                end
            end
            if ~obj.ModulatedInput
                if obj.pNeedSteeringAngle
                    if obj.WeightsInputPort
                        stang = stangArg;
                    else
                        stang = wArg;
                    end
                    sensor_resp = step(obj.cSensor,abs(subfreq),ang,obj.PropagationSpeed,stang);
                else
                    sensor_resp = step(obj.cSensor,abs(subfreq),ang);
                end
                phaseMtx.H = phaseMtx_temp.*sensor_resp.H;                     % chan x ang x freq
                phaseMtx.V = phaseMtx_temp.*sensor_resp.V;
            else
                phaseMtx = phaseMtx_temp;                                      % chan x ang x freq
            end
            
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end

            % Translate to frequency domain
            xtempfreq_x = step(obj.cSubbandDivider,complex(x_X));                       % snap x ang x freq
            xtempfreq_y = step(obj.cSubbandDivider,complex(x_Y));
            xtempfreq_z = step(obj.cSubbandDivider,complex(x_Z));
            NumChannel = obj.pDOF;
            NumAngles = size(ang,2);
            
            nfft = obj.pFFTLength;
            phaseModulatedXFreq_h = complex(zeros(size(xtempfreq_x,1),NumChannel,nfft));  % snap x ang x freq
            phaseModulatedXFreq_v = complex(zeros(size(xtempfreq_x,1),NumChannel,nfft));  % snap x ang x freq
            
            for m = 1:NumAngles
                
                temp_phaseMtx_h = squeeze(phaseMtx.H(:,m,:));
                temp_phaseMtx_v = squeeze(phaseMtx.V(:,m,:));
                
                arrayEffect_temp_h = sph2cartvec(...
                    [temp_phaseMtx_h(:).';zeros(1,NumChannel*nfft);zeros(1,NumChannel*nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp_h = phased.internal.local2globalvec(...
                    arrayEffect_temp_h,lclaxes);
                arrayEffect_temp_v = sph2cartvec(...
                    [zeros(1,NumChannel*nfft);temp_phaseMtx_v(:).';zeros(1,NumChannel*nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp_v = phased.internal.local2globalvec(...
                    arrayEffect_temp_v,lclaxes);
                
                arrayEffect_h_x = reshape(arrayEffect_temp_h(1,:).',NumChannel,[]).';     % freq x chan
                arrayEffect_h_y = reshape(arrayEffect_temp_h(2,:).',NumChannel,[]).';
                arrayEffect_h_z = reshape(arrayEffect_temp_h(3,:).',NumChannel,[]).';
                arrayEffect_v_x = reshape(arrayEffect_temp_v(1,:).',NumChannel,[]).';     % freq x chan
                arrayEffect_v_y = reshape(arrayEffect_temp_v(2,:).',NumChannel,[]).';
                arrayEffect_v_z = reshape(arrayEffect_temp_v(3,:).',NumChannel,[]).';
                
                phaseModulatedXFreq_h = phaseModulatedXFreq_h + permute(...
                    bsxfun(@times,permute(arrayEffect_h_x,[3 1 2]),reshape(xtempfreq_x(:,m,:),[],nfft)) + ... % snap x freq x chan
                    bsxfun(@times,permute(arrayEffect_h_y,[3 1 2]),reshape(xtempfreq_y(:,m,:),[],nfft)) + ...
                    bsxfun(@times,permute(arrayEffect_h_z,[3 1 2]),reshape(xtempfreq_z(:,m,:),[],nfft)),...
                    [1 3 2]);       % snap x chan x freq
                phaseModulatedXFreq_v = phaseModulatedXFreq_v + permute(...
                    bsxfun(@times,permute(arrayEffect_v_x,[3 1 2]),reshape(xtempfreq_x(:,m,:),[],nfft)) + ... % snap x freq x chan
                    bsxfun(@times,permute(arrayEffect_v_y,[3 1 2]),reshape(xtempfreq_y(:,m,:),[],nfft)) + ...
                    bsxfun(@times,permute(arrayEffect_v_z,[3 1 2]),reshape(xtempfreq_z(:,m,:),[],nfft)),...
                    [1 3 2]);       % snap x chan x freq
                
            end
            
        end

        function phaseModulatedXFreq = ...
                independentArrayChannelFreqDomainModulate(obj,x,ang)
        % each element collect its own signal

            % Find response for each element at different subbands for
            % given angles
            phaseMtx_temp = step(obj.cSensor,abs(obj.pSubbandFreqs),ang);    % chan x ang x freq
            if isPolarizationEnabled(obj.cSensor)
                phaseMtx = hypot(phaseMtx_temp.H,phaseMtx_temp.V);
            else
                phaseMtx = phaseMtx_temp; 
            end
            nfft = obj.pFFTLength;
            phaseMtxResh = reshape(phaseMtx,[],nfft);
            % We actually have more responses than we need. We extract the
            % response for nth channel corresponding to nth angle. Since
            % number of channels and number of angles always match in this
            % case, these are the responses we need.
            phaseMtxChan = phaseMtxResh(1:obj.pDOF+1:end,:).';               % freq x chan  (chan = ang)

            % Translate to frequency domain and modulate.
            xtempfreq = step(obj.cSubbandDivider,complex(x));                % snap x chan x freq
            phaseModulatedXFreq = bsxfun(@times,xtempfreq,permute(phaseMtxChan,[3 2 1]));  % snap x chan x freq

        end

        function phaseModulatedXFreq = ...
                independentArrayChannelFreqDomainPolarizedModulate(obj,x,ang,lclaxes)
        % each element collect its own signal

            % Find response for each element at different subbands for
            % given angles
            phaseMtx = step(obj.cSensor,abs(obj.pSubbandFreqs),ang);         % chan x ang x freq
            
            phaseMtx_h = phaseMtx.H;
            phaseMtx_v = phaseMtx.V;
            
            nfft = obj.pFFTLength;
            phaseMtx_h_Resh = reshape(phaseMtx_h,[],nfft);
            phaseMtx_v_Resh = reshape(phaseMtx_v,[],nfft);
            % We actually have more responses than we need. We extract the
            % response for nth channel corresponding to nth angle. Since
            % number of channels and number of angles always match in this
            % case, these are the responses we need.
            temp_phaseMtx_h = phaseMtx_h_Resh(1:obj.pDOF+1:end,:).';         % freq x chan  (chan = ang)
            temp_phaseMtx_v = phaseMtx_v_Resh(1:obj.pDOF+1:end,:).';         % freq x chan  (chan = ang)
           
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end

            % Translate to frequency domain and modulate.
            xtempfreq_x = step(obj.cSubbandDivider,complex(x_X));            % snap x chan x freq
            xtempfreq_y = step(obj.cSubbandDivider,complex(x_Y));            % snap x chan x freq
            xtempfreq_z = step(obj.cSubbandDivider,complex(x_Z));            % snap x chan x freq
            
            NumChannel = obj.pDOF;
            NumAngles = size(ang,2);
            
            phaseModulatedXFreq = complex(zeros(size(xtempfreq_x,1),NumChannel,nfft));  % snap x chan x freq
            
            for m = 1:NumAngles
                
                arrayEffect_temp = sph2cartvec(...
                    [temp_phaseMtx_h(:,m).';temp_phaseMtx_v(:,m).';zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp = phased.internal.local2globalvec(...
                    arrayEffect_temp,lclaxes);
                
                arrayEffect_x = arrayEffect_temp(1,:).';        % freq x chan
                arrayEffect_y = arrayEffect_temp(2,:).';
                arrayEffect_z = arrayEffect_temp(3,:).';
                
                phaseModulatedXFreq(:,m,:) = phaseModulatedXFreq(:,m,:) + ...
                    bsxfun(@times,xtempfreq_x(:,m,:),permute(arrayEffect_x,[3 2 1])) + ...
                    bsxfun(@times,xtempfreq_y(:,m,:),permute(arrayEffect_y,[3 2 1])) + ...
                    bsxfun(@times,xtempfreq_z(:,m,:),permute(arrayEffect_z,[3 2 1]));  % snap x chan x freq
                
            end
            
        end
        
        function [phaseModulatedXFreq_h,phaseModulatedXFreq_v] = ...
                independentArrayChannelFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes)
        % each element collect its own signal

            % Find response for each element at different subbands for
            % given angles
            phaseMtx = step(obj.cSensor,abs(obj.pSubbandFreqs),ang);         % chan x ang x freq
            
            phaseMtx_h = phaseMtx.H;
            phaseMtx_v = phaseMtx.V;
            
            nfft = obj.pFFTLength;
            phaseMtx_h_Resh = reshape(phaseMtx_h,[],nfft);
            phaseMtx_v_Resh = reshape(phaseMtx_v,[],nfft);
            % We actually have more responses than we need. We extract the
            % response for nth channel corresponding to nth angle. Since
            % number of channels and number of angles always match in this
            % case, these are the responses we need.
            temp_phaseMtx_h = phaseMtx_h_Resh(1:obj.pDOF+1:end,:).';         % freq x chan  (chan = ang)
            temp_phaseMtx_v = phaseMtx_v_Resh(1:obj.pDOF+1:end,:).';         % freq x chan  (chan = ang)
           
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end

            % Translate to frequency domain and modulate.
            xtempfreq_x = step(obj.cSubbandDivider,complex(x_X));            % snap x chan x freq
            xtempfreq_y = step(obj.cSubbandDivider,complex(x_Y));            % snap x chan x freq
            xtempfreq_z = step(obj.cSubbandDivider,complex(x_Z));            % snap x chan x freq
            
            NumChannel = obj.pDOF;
            NumAngles = size(ang,2);
            
            phaseModulatedXFreq_h = complex(zeros(size(xtempfreq_x,1),NumChannel,nfft));  % snap x chan x freq
            phaseModulatedXFreq_v = complex(zeros(size(xtempfreq_x,1),NumChannel,nfft));  % snap x chan x freq
            
            for m = 1:NumAngles
                
                arrayEffect_temp_h = sph2cartvec(...
                    [temp_phaseMtx_h(:,m).';zeros(1,nfft);zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp_h = phased.internal.local2globalvec(...
                    arrayEffect_temp_h,lclaxes);
                arrayEffect_temp_v = sph2cartvec(...
                    [zeros(1,nfft);temp_phaseMtx_v(:,m).';zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp_v = phased.internal.local2globalvec(...
                    arrayEffect_temp_v,lclaxes);
                
                arrayEffect_h_x = arrayEffect_temp_h(1,:).';        % freq x chan
                arrayEffect_h_y = arrayEffect_temp_h(2,:).';
                arrayEffect_h_z = arrayEffect_temp_h(3,:).';
                arrayEffect_v_x = arrayEffect_temp_v(1,:).';        % freq x chan
                arrayEffect_v_y = arrayEffect_temp_v(2,:).';
                arrayEffect_v_z = arrayEffect_temp_v(3,:).';
                
                phaseModulatedXFreq_h(:,m,:) = phaseModulatedXFreq_h(:,m,:) + ...
                    bsxfun(@times,xtempfreq_x(:,m,:),permute(arrayEffect_h_x,[3 2 1])) + ...
                    bsxfun(@times,xtempfreq_y(:,m,:),permute(arrayEffect_h_y,[3 2 1])) + ...
                    bsxfun(@times,xtempfreq_z(:,m,:),permute(arrayEffect_h_z,[3 2 1]));  % snap x chan x freq
                phaseModulatedXFreq_v(:,m,:) = phaseModulatedXFreq_v(:,m,:) + ...
                    bsxfun(@times,xtempfreq_x(:,m,:),permute(arrayEffect_v_x,[3 2 1])) + ...
                    bsxfun(@times,xtempfreq_y(:,m,:),permute(arrayEffect_v_y,[3 2 1])) + ...
                    bsxfun(@times,xtempfreq_z(:,m,:),permute(arrayEffect_v_z,[3 2 1]));  % snap x chan x freq               
            end
            
        end
        
        function phaseModulatedXFreq = ...
                singleElementFreqDomainModulate(obj,x,ang,w)   %#ok<INUSD>
        % collect all signals with one element and sum them up

            % Calculate the element response for each subband at given
            % angles
            subfreq = obj.pSubbandFreqs;
            phaseMtx_temp = step(obj.cSensor,abs(subfreq),ang);   % ang x freq (chan = 1)
            if isPolarizationEnabled(obj.cSensor)
                phaseMtx = hypot(phaseMtx_temp.H,phaseMtx_temp.V).';
            else
                phaseMtx = phaseMtx_temp.'; % freq x ang (chan = 1)
            end
            if obj.ModulatedInput && obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);  % 1 x freq
                for m = 1:numel(subfreq)
                    phaseMtx(m,:) = phased.internal.normalizeIntegratedPower(phaseMtx(m,:),intresp(m),false);
                end
            end

            % Modulate to frequency domain and modulate
            % Combine all signals from different directions together
            xtempfreq = step(obj.cSubbandDivider,complex(x));      % snap x ang x freq (chan = 1)
            phaseModulatedXFreq = sum(bsxfun(@times,...
                xtempfreq,permute(phaseMtx,[3 2 1])),2);           % snap x chan x freq

        end

        function phaseModulatedXFreq = ...
                singleElementFreqDomainPolarizedModulate(obj,x,ang,lclaxes,w) %#ok<INUSD>
        % collect all signals with one element and sum them up

            % Calculate the element response for each subband at given
            % angles
            subfreq = obj.pSubbandFreqs;
            phaseMtx = step(obj.cSensor,abs(subfreq),ang);         % ang x freq (chan = 1)
            if obj.ModulatedInput && obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);
                for m = 1:numel(subfreq)
                    phaseMtx.H(:,m) = phased.internal.normalizeIntegratedPower(phaseMtx.H(:,m),intresp.H(m),false);
                    phaseMtx.V(:,m) = phased.internal.normalizeIntegratedPower(phaseMtx.V(:,m),intresp.V(m),false);
                end
            end
            
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end

            % Modulate to frequency domain and modulate
            xtempfreq_x = step(obj.cSubbandDivider,complex(x_X));   % snap x ang x freq
            xtempfreq_y = step(obj.cSubbandDivider,complex(x_Y));
            xtempfreq_z = step(obj.cSubbandDivider,complex(x_Z));
                       
            NumAngles = size(ang,2);
            nfft = obj.pFFTLength;
            phaseModulatedXFreq = complex(zeros(size(xtempfreq_x,1),1,nfft));
            
            for m = 1:NumAngles
                
                temp_phaseMtx_h = phaseMtx.H(m,:);
                temp_phaseMtx_v = phaseMtx.V(m,:);
                
                arrayEffect_temp = sph2cartvec(...
                    [temp_phaseMtx_h(:).';temp_phaseMtx_v(:).';zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp = phased.internal.local2globalvec(...
                    arrayEffect_temp,lclaxes);
                
                arrayEffect_x = arrayEffect_temp(1,:).';        % freq x chan (chan = 1)
                arrayEffect_y = arrayEffect_temp(2,:).';
                arrayEffect_z = arrayEffect_temp(3,:).';
                
                phaseModulatedXFreq = phaseModulatedXFreq + ...
                    bsxfun(@times,permute(arrayEffect_x,[3 2 1]),xtempfreq_x(:,m,:)) + ...
                    bsxfun(@times,permute(arrayEffect_y,[3 2 1]),xtempfreq_y(:,m,:)) + ...
                    bsxfun(@times,permute(arrayEffect_z,[3 2 1]),xtempfreq_z(:,m,:));
                
            end

        end

        function [phaseModulatedXFreq_h,phaseModulatedXFreq_v] = ...
                singleElementFreqDomainDualPolarizedModulate(obj,x,ang,lclaxes,w) %#ok<INUSD>
        % collect all signals with one element and sum them up

            % Calculate the element response for each subband at given
            % angles
            subfreq = obj.pSubbandFreqs;
            phaseMtx = step(obj.cSensor,abs(subfreq),ang);         % ang x freq (chan = 1)
            if obj.ModulatedInput && obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);
                for m = 1:numel(subfreq)
                    phaseMtx.H(:,m) = phased.internal.normalizeIntegratedPower(phaseMtx.H(:,m),intresp.H(m),false);
                    phaseMtx.V(:,m) = phased.internal.normalizeIntegratedPower(phaseMtx.V(:,m),intresp.V(m),false);
                end
            end
            
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end

            % Modulate to frequency domain and modulate
            xtempfreq_x = step(obj.cSubbandDivider,complex(x_X));   % snap x ang x freq
            xtempfreq_y = step(obj.cSubbandDivider,complex(x_Y));
            xtempfreq_z = step(obj.cSubbandDivider,complex(x_Z));
                       
            NumAngles = size(ang,2);
            nfft = obj.pFFTLength;
            phaseModulatedXFreq_h = complex(zeros(size(xtempfreq_x,1),1,nfft));
            phaseModulatedXFreq_v = complex(zeros(size(xtempfreq_x,1),1,nfft));
            
            for m = 1:NumAngles
                
                temp_phaseMtx_h = phaseMtx.H(m,:);
                temp_phaseMtx_v = phaseMtx.V(m,:);
                
                arrayEffect_temp_h = sph2cartvec(...
                    [temp_phaseMtx_h(:).';zeros(1,nfft);zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp_h = phased.internal.local2globalvec(...
                    arrayEffect_temp_h,lclaxes);
                arrayEffect_temp_v = sph2cartvec(...
                    [zeros(1,nfft);temp_phaseMtx_v(:).';zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp_v = phased.internal.local2globalvec(...
                    arrayEffect_temp_v,lclaxes);
                
                arrayEffect_h_x = arrayEffect_temp_h(1,:).';        % freq x chan (chan = 1)
                arrayEffect_h_y = arrayEffect_temp_h(2,:).';
                arrayEffect_h_z = arrayEffect_temp_h(3,:).';
                arrayEffect_v_x = arrayEffect_temp_v(1,:).';        % freq x chan (chan = 1)
                arrayEffect_v_y = arrayEffect_temp_v(2,:).';
                arrayEffect_v_z = arrayEffect_temp_v(3,:).';
                
                phaseModulatedXFreq_h = phaseModulatedXFreq_h + ...
                    bsxfun(@times,permute(arrayEffect_h_x,[3 2 1]),xtempfreq_x(:,m,:)) + ...
                    bsxfun(@times,permute(arrayEffect_h_y,[3 2 1]),xtempfreq_y(:,m,:)) + ...
                    bsxfun(@times,permute(arrayEffect_h_z,[3 2 1]),xtempfreq_z(:,m,:));
                phaseModulatedXFreq_v = phaseModulatedXFreq_v + ...
                    bsxfun(@times,permute(arrayEffect_v_x,[3 2 1]),xtempfreq_x(:,m,:)) + ...
                    bsxfun(@times,permute(arrayEffect_v_y,[3 2 1]),xtempfreq_y(:,m,:)) + ...
                    bsxfun(@times,permute(arrayEffect_v_z,[3 2 1]),xtempfreq_z(:,m,:));
                
            end

        end
        
        function y = overlapadd(obj,x,NumOverlap)
            % add overlap from last block
            x(1:NumOverlap,:) = x(1:NumOverlap,:)+obj.pBuffer(1:NumOverlap,:);
            M = size(x,1) - NumOverlap;
            obj.pBuffer(1:NumOverlap,:) = x(M+(1:NumOverlap),:);
            y = x(1:M,:);
        end

    end

    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractSensorOperation;
            dWavefront = ...
                matlab.system.display.internal.Property('Wavefront', ...
                                                        'IsGraphical', false);
            dPolarization = ...
                matlab.system.display.internal.Property('Polarization', ...
                                                        'IsGraphical', false);
            % dSampleRate = matlab.system.display.internal.Property(...
            %     'SampleRate','IsObjectDisplayOnly',true);
            
            props = {...
                'SampleRateFromInputCheckbox',...
                'SampleRate',...
                'ModulatedInput',...
                'CarrierFrequency',...
                'NumSubbands',...
                dWavefront,...
                'SensorGainMeasure',...
                dPolarization,...
                'WeightsInputPort'};
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:WidebandCollectorTitle')),...
                'Text',getString(message('phased:library:block:WidebandCollectorDesc')));
        end    
    end
    methods (Access = protected) %for Simulink
        function varargout = getOutputNamesImpl(obj)
            if obj.pEnableDualPolarizationOutput
                varargout = {'Yh','Yv'};
            else
                varargout = {'Y'};
            end
        end
        function varargout = getInputNamesImpl(obj)
            varargout = {'X','Ang'};
            lastIdx = 3;
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
            str = sprintf('Wideband\n Rx Array');
        end        
        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            %Output is number of snapshots * number of elements
            if  ~isa(obj.Sensor,'phased.internal.AbstractElement')
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
        function varargout = isOutputComplexImpl(obj)
            if obj.ModulatedInput
                varargout{1} = true;
            else
                varargout{1} = false;
            end
        end
    end    
    
end
