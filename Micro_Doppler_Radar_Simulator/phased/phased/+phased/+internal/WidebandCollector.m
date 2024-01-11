classdef (Sealed, StrictDefaults) WidebandCollector < phased.internal.AbstractSensorOperation & ...
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
%   WidebandCollector methods:
%
%   step     - Collect signals (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a wideband collector object with same property values
%   isLocked - Locked status (logical)
%   reset    - Reset states of wideband signal collector object
%
%   WidebandCollector properties:
%
%   Sensor             - Handle of the sensor
%   PropagationSpeed   - Propagation speed
%   SampleRate         - Sample rate
%   ModulatedInput     - Assume modulated input
%   CarrierFrequency   - Carrier frequency
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
%   hc = phased.WidebandCollector('Sensor',ha);
%   x = [1;1];
%   incidentAngle = [10 30]';
%   y = step(hc,x,incidentAngle);
%
%   % Example 2: 
%   %   Collect a far field signal with a 5-element array.
%
%   ha = phased.ULA('NumElements',5);
%   hc = phased.WidebandCollector('Sensor',ha);
%   x = [1;1];
%   incidentAngle = [10 30]';
%   y = step(hc,x,incidentAngle);
%
%   % Example 3: 
%   %   Collect signal with a 3-element antenna array. Each antenna 
%   %   collects a separate input signal from a separate direction.
%
%   ha = phased.ULA('NumElements',3);
%   hc = phased.WidebandCollector('Sensor',ha,'Wavefront','Plane');
%   x = rand(10,3);   % Each column is a separate signal for one element
%   incidentAngle = [10 0; 20 5; 45 2]'; % 3 angles for 3 signals
%   y = step(hc,x,incidentAngle);
%
%   See also phased, phased.Collector.

%   Copyright 2009-2017 The MathWorks, Inc.
%     

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
        %SignalType     Signal type
        %   Specify whether the signal should be treated as one of 'Field'
        %   | 'Power', where 'Field' is the default. When you set
        %   SignalType to 'Field', the collected signal is scaled by the
        %   field at the corresponding direction. When you set SignalType
        %   to 'Power', the collected signal is scaled by the directivity
        %   at the corresponding direction.
        SignalType = 'Field';
        %Wavefront Type of incoming wavefront
        %   Specify the type of incoming wavefront as one of 'Plane' |
        %   'Unspecified', where the default is 'Plane'. When you set the
        %   Wavefront property to 'Plane', the input signals are assumed to
        %   be multiple plane waves impinging on the entire array. Each
        %   plane wave is received by all collecting elements. If you set
        %   the Wavefront property to 'Unspecified', the input signals are
        %   assumed to be individual waves impinging on individual sensors.
        %   Wavefront must be 'Plane' if you set SignalType to 'Power'.
        Wavefront = 'Plane'
    end

    properties (Nontunable, Logical) 
        %ModulatedInput Assume modulated input
        %   Set this property to true to indicate the input signal is
        %   demodulated at a carrier frequency. The default value is true. 
        %   ModulatedInput must be true if you set SignalType to 'Power'.
        ModulatedInput = true
    end

    properties (Nontunable, Logical)
        %EnablePolarization  Enable polarization
        %   Set this property to true to enable polarization. Set this
        %   property to false to ignore polarization. The default value of
        %   this property is false. This property applies when the sensor
        %   specified in the Sensor property is capable of simulating
        %   polarization.
        EnablePolarization = false
    end

    properties(Constant, Hidden)
        WavefrontSet = matlab.system.StringSet({'Plane','Unspecified'});
        SignalTypeSet = matlab.system.StringSet({'Field','Power'});
    end

    properties (Access = private, Nontunable)
        pFFTLength      % number of FFT
        pSubbandFreqs   % center frequency for each FFT bin
        pNumOverlap     % number of overlapped samples in overlap-add
        cSteeringVector % handle to steering vector 
        cIntegratedPattern
        pDOF
        pSampleRate
    end
    
    properties (Access = private)
        pBuffer         % buffer for overlap-add
    end
    
    properties (Access = private)
        cFFT
        cIFFT
    end

    properties (Access = private, Logical, Nontunable)
        pIsPowerSignal
    end
    
    methods
        function set.SampleRate(obj,value)
            validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'SampleRate');
            obj.SampleRate = value;
        end
        function set.CarrierFrequency(obj,value)
            validateattributes( value, { 'double' }, { 'scalar', 'positive', 'integer' }, '', 'CarrierFrequency');
            obj.CarrierFrequency = value;
        end
        function obj = WidebandCollector(varargin)
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
    
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractSensorOperation(obj);
            cond = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    (obj.Wavefront(1) == 'U'); %Unspecified
            if cond
                coder.internal.errorIf(cond,'phased:phased:collector:invalidSettingsForSubarray','Wavefront','Plane','Sensor');
            end
            
            cond = strcmp(obj.SignalType,'Power') && ~obj.ModulatedInput;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:IrrelevantSetting','Power','ModulatedInput','false');
            end

            cond = strcmp(obj.SignalType,'Power') && strcmp(obj.Wavefront,'Custom');
            if cond
                coder.internal.errorIf(cond,'phased:system:array:IrrelevantSetting','Power','Wavefront','Plane');
            end
            
            cond = obj.EnablePolarization && ~isPolarizationCapable(obj.Sensor);
            if cond
                coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
            end
        end

        function validateInputsImpl(obj,x,angle,lclaxes,weightsArg,stangArg)
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
            if obj.EnablePolarization
                % check lclAxes is 3x3
                sigdatatypes.validate3DCartCoord(lclaxes,'','LAxes',...
                    {'size',[3 3]});
            end

            % steering angle for subarray
            if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    ~(obj.Sensor.SubarraySteering(1) == 'N') %None
                if ~obj.WeightsInputPort
                    if obj.EnablePolarization
                        stang = weightsArg;
                    else
                        stang = lclaxes;
                    end
                else
                    if ~obj.EnablePolarization
                        stang = weightsArg;
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
            
            % weights
            if obj.WeightsInputPort
                if ~obj.EnablePolarization
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
            if obj.EnablePolarization
                xsize = [size(x(1).X,1) size(x,2)];
            else
                xsize = size(x);
            end
            
            obj.pSampleRate = getSampleRate(obj,xsize(1),1,obj.SampleRate);
            
            cond = obj.ModulatedInput && ...
                    (obj.pSampleRate > 2*obj.CarrierFrequency);
            if cond
                coder.internal.errorIf(cond,'phased:phased:WidebandCollector:FsTooHigh');
            end
            
            % Determine the FFT length
            nfft = ceil(xsize(1)*1.5); % 50% overlapping
            nfft = max(nfft,128); % minimum FFT length is 128.
            obj.pFFTLength = nfft;
            
            obj.cFFT = dsp.FFT('FFTLengthSource','Property','FFTLength',nfft);
            obj.cIFFT = dsp.IFFT('FFTLengthSource','Property','FFTLength',nfft);
            
            obj.pNumOverlap = nfft - xsize(1);
            % Pre-calculate the center frequency of each FFT bin
            nfft = obj.pFFTLength;
            binStart = floor(nfft/2);
            tempSubbandFreqs = ifftshift((1:nfft)-binStart-1)*obj.pSampleRate/nfft;
            % Adjust the frequencies for modulated signal.
            if obj.ModulatedInput
                obj.pSubbandFreqs = obj.CarrierFrequency+tempSubbandFreqs;
            else
                obj.pSubbandFreqs = tempSubbandFreqs;
            end
            if obj.pUseArray
                obj.pDOF = getDOF(obj.cSensor);
            else
                obj.pDOF = 1;
            end
            obj.pIsPowerSignal = strcmp(obj.SignalType,'Power');
            
            if obj.pUseArray && (obj.Wavefront(1) == 'P') %Plane
                if obj.ModulatedInput
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
                else
                    obj.cSteeringVector = phased.SteeringVector(...
                        'SensorArray',obj.cSensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'IncludeElementResponse',false);
                end
            elseif ~obj.pUseArray && obj.ModulatedInput && obj.pIsPowerSignal
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
                if obj.ModulatedInput
                    flag = false;
                else
                    flag = true;
                end
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
        
        function resetImpl(obj)
            if obj.pUseArray
                if (obj.Wavefront(1) == 'P') %'Plane'
                    reset(obj.cSteeringVector);
                    if ~obj.ModulatedInput
                        reset(obj.cSensor);
                    elseif obj.pIsPowerSignal
                        reset(obj.cIntegratedPattern);
                    end
                else
                    reset(obj.cSensor);
                end
            else
                reset(obj.cSensor);
                if obj.ModulatedInput && obj.pIsPowerSignal
                    reset(obj.cIntegratedPattern);
                end
            end

            if obj.ModulatedInput
                obj.pBuffer = complex(zeros(obj.pNumOverlap,obj.pDOF));
            else
                obj.pBuffer = zeros(obj.pNumOverlap,obj.pDOF);
            end
            reset(obj.cFFT);
            reset(obj.cIFFT);
        end

        function releaseImpl(obj)
            if obj.pUseArray
                if (obj.Wavefront(1) == 'P') %'Plane'
                    release(obj.cSteeringVector);
                    if ~obj.ModulatedInput
                        release(obj.cSensor);
                    elseif obj.pIsPowerSignal
                        release(obj.cIntegratedPattern);
                    end
                else
                    release(obj.cSensor);
                end
            else
                release(obj.cSensor);
                if obj.ModulatedInput && obj.pIsPowerSignal
                    release(obj.cIntegratedPattern);
                end
            end
            release(obj.cFFT);
            release(obj.cIFFT);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSensorOperation(obj);
            if isLocked(obj)
                s.pDOF = obj.pDOF;
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.cIntegratedPattern = saveobj(obj.cIntegratedPattern);
                s.cFFT = saveobj(obj.cFFT);
                s.cIFFT = saveobj(obj.cIFFT);
                s.pFFTLength = obj.pFFTLength;
                s.pSubbandFreqs = obj.pSubbandFreqs;
                s.pBuffer = obj.pBuffer;
                s.pNumOverlap = obj.pNumOverlap;
                s.pSampleRate = obj.pSampleRate;
                s.pIsPowerSignal = obj.pIsPowerSignal;
            end
        end

        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSensorOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
                    if isfield(s,'cFFT')
                        obj.cFFT = dsp.FFT.loadobj(s.cFFT);
                        s = rmfield(s,'cFFT');
                    end
                    if isfield(s,'cIFFT')
                        obj.cIFFT = dsp.IFFT.loadobj(s.cIFFT);
                        s = rmfield(s,'cIFFT');
                    end
                    if isfield(s,'cIntegratedPattern') && ...
                            isfield(s.cIntegratedPattern,'ClassNameForLoadTimeEval')
                        obj.cIntegratedPattern = eval(...
                            sprintf('%s.loadobj(s.cIntegratedPattern)',s.cIntegratedPattern.ClassNameForLoadTimeEval));
                        s = rmfield(s,'cIntegratedPattern');
                    end
                    % recover locked sample rate information
                    if isfield(s,'pSampleRate')
                        obj.pSampleRate = s.pSampleRate;
                        s = rmfield(s,'pSampleRate');
                    else
                        obj.pSampleRate = s.SampleRate;
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

        function y = stepImpl(obj,x,ang,lclaxes,w,stangArg)

            nfft = obj.pFFTLength;
            if obj.EnablePolarization
                % Modulate each subband with corresponding phase in frequency
                % domain
                if obj.pUseArray
                    if (obj.Wavefront(1) == 'P') %Plane
                        if obj.pNeedSteeringAngle
                            if ~obj.WeightsInputPort
                                stang = w;
                                y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft,stang);
                            else
                                wt = w;
                                stang = stangArg;
                                y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft,wt,stang);
                            end
                        else
                            if obj.WeightsInputPort
                                wt = w;
                                y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft,wt);
                            else
                                y_temp = farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft);
                            end
                        end
                    else
                        if obj.WeightsInputPort
                            wt = w;
                        end
                        y_temp = independentArrayChannelFreqDomainPolarizedModulate(...
                            obj,x,ang,lclaxes,nfft);
                    end
                else
                    if obj.WeightsInputPort
                        wt = w;
                        y_temp = singleElementFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft,wt);
                    else
                        y_temp = singleElementFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft);
                    end
                end
            else % no polarization
                % Modulate each subband with corresponding phase in frequency
                % domain
                if obj.pUseArray
                    if (obj.Wavefront(1) == 'P') %Plane
                        if obj.pNeedSteeringAngle
                            if ~obj.WeightsInputPort
                                stang = lclaxes;
                                y_temp = farFieldFreqDomainModulate(obj,x,ang,nfft,stang);
                            else
                                stang = w;
                                wt = lclaxes;
                                y_temp = farFieldFreqDomainModulate(obj,x,ang,nfft,wt,stang);
                            end
                        else
                            if obj.WeightsInputPort
                                wt = lclaxes;
                                y_temp = farFieldFreqDomainModulate(obj,x,ang,nfft,wt);
                            else
                                y_temp = farFieldFreqDomainModulate(obj,x,ang,nfft);
                            end
                        end
                    else
                        if obj.WeightsInputPort
                            wt = lclaxes;
                        end
                        y_temp = independentArrayChannelFreqDomainModulate(...
                            obj,x,ang,nfft);
                    end
                else
                    if obj.WeightsInputPort
                        wt = lclaxes;
                        y_temp = singleElementFreqDomainModulate(obj,x,ang,nfft,wt);
                    else
                        y_temp = singleElementFreqDomainModulate(obj,x,ang,nfft);
                    end
                end

            end
            % Convert back to time domain
            if isreal(y_temp)
                y_ifft_in = complex(y_temp);
            else
                y_ifft_in = y_temp;
            end
            y_out = step(obj.cIFFT,y_ifft_in);
            if ~obj.ModulatedInput
                y = real(y_out);
            else
                if isreal(y_out)
                    y = complex(y_out);
                else
                    y = y_out;
                end
               % y = y_out;
            end
            y = overlapadd(obj,y,obj.pNumOverlap);

            % Apply weights
            if obj.WeightsInputPort
                for m = 1:obj.pDOF
                    y(:,m) = y(:,m)*wt(m);
                end
            end

        end

        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if ~obj.ModulatedInput && strcmp(prop, 'CarrierFrequency')
                flag = true;
            end
        end

    end

    methods (Access = private)

        function phaseModulatedXFreq = ...
                farFieldFreqDomainModulate(obj,x,ang,nfft,wArg,stangArg)
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
                        if obj.pIsPowerSignal
                            intresp = step(obj.cIntegratedPattern,subfreq,w,stang);
                        end
                    else
                        stang = wArg;
                        if obj.pIsPowerSignal
                            intresp = step(obj.cIntegratedPattern,subfreq,stang);
                        end
                    end
                    phaseMtx = step(obj.cSteeringVector,subfreq,ang,stang);
                    if obj.pIsPowerSignal
                        phaseMtx = phased.internal.normalizeIntegratedPower(phaseMtx,intresp,false);
                    end
                else
                    negidx = find(subfreq<0);
                    phaseMtx = step(obj.cSteeringVector,abs(subfreq),ang);
                    phaseMtx(:,:,negidx) = conj(phaseMtx(:,:,negidx));
                end
            else
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;
                        if obj.pIsPowerSignal
                            intresp = step(obj.cIntegratedPattern,subfreq,w);
                        end
                    elseif obj.pIsPowerSignal
                        intresp = step(obj.cIntegratedPattern,subfreq);
                    end
                    phaseMtx = step(obj.cSteeringVector,subfreq,ang);
                    if obj.pIsPowerSignal
                        phaseMtx = phased.internal.normalizeIntegratedPower(phaseMtx,intresp,false);
                    end
                else
                    negidx = find(subfreq<0);
                    phaseMtx = step(obj.cSteeringVector,abs(subfreq),ang);
                    phaseMtx(:,:,negidx) = conj(phaseMtx(:,:,negidx));
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
            if isreal(x)
                x_fft_in = complex(x);
            else
                x_fft_in = x;
            end
            xtempfreq = step(obj.cFFT,x_fft_in);
            NumChannel = obj.pDOF;

            % Modulate each channel at different subbands and angles
            phaseModulatedXFreq = reshape(phaseMtx,NumChannel,[]).*...
                (ones(NumChannel,1)*reshape(xtempfreq.',1,[]));
            % Sum across all angles
            phaseModulatedXFreq = sum(permute(...
                reshape(phaseModulatedXFreq,NumChannel,[],nfft),...
                [1 3 2]),3).';

        end

        function phaseModulatedXFreq = ...
                farFieldFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft,wArg,stangArg)
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
                        if obj.pIsPowerSignal
                            intresp = step(obj.cIntegratedPattern,subfreq,w,stang);
                        end
                    else
                        stang = wArg;
                        if obj.pIsPowerSignal
                            intresp = step(obj.cIntegratedPattern,subfreq,stang);
                        end
                    end
                    phaseMtx_temp = step(obj.cSteeringVector,subfreq,ang,stang);
                    if obj.pIsPowerSignal
                        phaseMtx_temp.H = phased.internal.normalizeIntegratedPower(phaseMtx_temp.H,intresp.H,false);
                        phaseMtx_temp.H = phased.internal.normalizeIntegratedPower(phaseMtx_temp.V,intresp.V,false);
                    end
                else
                    negidx = find(subfreq<0);
                    phaseMtx_temp = step(obj.cSteeringVector,abs(subfreq),ang);
                    phaseMtx_temp(:,:,negidx) = conj(phaseMtx_temp(:,:,negidx));
                end
            else
                if obj.ModulatedInput
                    if obj.WeightsInputPort
                        w = wArg;
                        if obj.pIsPowerSignal
                            intresp = step(obj.cIntegratedPattern,subfreq,w);
                        end
                    elseif obj.pIsPowerSignal
                        intresp = step(obj.cIntegratedPattern,subfreq);
                    end
                    phaseMtx_temp = step(obj.cSteeringVector,subfreq,ang);
                    if obj.pIsPowerSignal
                        phaseMtx_temp.H = phased.internal.normalizeIntegratedPower(phaseMtx_temp.H,intresp.H,false);
                        phaseMtx_temp.H = phased.internal.normalizeIntegratedPower(phaseMtx_temp.V,intresp.V,false);
                    end
                else
                    negidx = find(subfreq<0);
                    phaseMtx_temp = step(obj.cSteeringVector,abs(subfreq),ang);
                    phaseMtx_temp(:,:,negidx) = conj(phaseMtx_temp(:,:,negidx));
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
                phaseMtx.H = phaseMtx_temp.*sensor_resp.H;
                phaseMtx.V = phaseMtx_temp.*sensor_resp.V;
            else
                phaseMtx = phaseMtx_temp;
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
            x_temp = x_X;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_x = step(obj.cFFT,x_fft_in);

            x_temp = x_Y;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_y = step(obj.cFFT,x_fft_in);

            x_temp = x_Z;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_z = step(obj.cFFT,x_fft_in);
            NumChannel = obj.pDOF;
            NumAngles = size(ang,2);
            
            phaseModulatedXFreq = complex(zeros(nfft,NumChannel));
            
            for m = 1:NumAngles
                
                temp_phaseMtx_h = squeeze(phaseMtx.H(:,m,:));
                temp_phaseMtx_v = squeeze(phaseMtx.V(:,m,:));
                
                arrayEffect_temp = sph2cartvec(...
                    [temp_phaseMtx_h(:).';temp_phaseMtx_v(:).';zeros(1,NumChannel*nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp = phased.internal.local2globalvec(...
                    arrayEffect_temp,lclaxes);
                
                arrayEffect_x = reshape(arrayEffect_temp(1,:).',NumChannel,[]).';
                arrayEffect_y = reshape(arrayEffect_temp(2,:).',NumChannel,[]).';
                arrayEffect_z = reshape(arrayEffect_temp(3,:).',NumChannel,[]).';
                
                phaseModulatedXFreq = phaseModulatedXFreq + ...
                    bsxfun(@times,arrayEffect_x,xtempfreq_x(:,m)) + ...
                    bsxfun(@times,arrayEffect_y,xtempfreq_y(:,m)) + ...
                    bsxfun(@times,arrayEffect_z,xtempfreq_z(:,m));
                
            end

        end

        function phaseModulatedXFreq = ...
                independentArrayChannelFreqDomainModulate(obj,x,ang,nfft)
        % each element collect its own signal

            % Find response for each element at different subbands for
            % given angles
            phaseMtx_temp = step(obj.cSensor,abs(obj.pSubbandFreqs),ang);
            if isPolarizationEnabled(obj.cSensor)
                phaseMtx = hypot(phaseMtx_temp.H,phaseMtx_temp.V);
            else
                phaseMtx = phaseMtx_temp; 
            end
            phaseMtxResh = reshape(phaseMtx,[],nfft);
            % We actually have more responses than we need. We extract the
            % response for nth channel corresponding to nth angle. Since
            % number of channels and number of angles always match in this
            % case, these are the responses we need.
            phaseMtxChan = phaseMtxResh(1:obj.pDOF+1:end,:).';

            % Translate to frequency domain and modulate.
            if isreal(x)
                x_fft_in = complex(x);
            else
                x_fft_in = x;
            end
            xtempfreq = step(obj.cFFT,x_fft_in);
            phaseModulatedXFreq = phaseMtxChan.*xtempfreq;

        end

        function phaseModulatedXFreq = ...
                independentArrayChannelFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft)
        % each element collect its own signal

            % Find response for each element at different subbands for
            % given angles
            phaseMtx = step(obj.cSensor,abs(obj.pSubbandFreqs),ang);
            
            phaseMtx_h = phaseMtx.H;
            phaseMtx_v = phaseMtx.V;
            
            phaseMtx_h_Resh = reshape(phaseMtx_h,[],nfft);
            phaseMtx_v_Resh = reshape(phaseMtx_v,[],nfft);
            % We actually have more responses than we need. We extract the
            % response for nth channel corresponding to nth angle. Since
            % number of channels and number of angles always match in this
            % case, these are the responses we need.
            temp_phaseMtx_h = phaseMtx_h_Resh(1:obj.pDOF+1:end,:).';
            temp_phaseMtx_v = phaseMtx_v_Resh(1:obj.pDOF+1:end,:).';
           
            M = size(x,2);
            x_X = complex(zeros(size(x(1).X,1),M));
            x_Y = x_X; x_Z = x_X;
            for m = 1:M
                x_X(:,m) = x(m).X;
                x_Y(:,m) = x(m).Y;
                x_Z(:,m) = x(m).Z;
            end

            % Translate to frequency domain and modulate.
            x_temp = x_X;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_x = step(obj.cFFT,x_fft_in);

            x_temp = x_Y;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_y = step(obj.cFFT,x_fft_in);
            
            x_temp = x_Z;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_z = step(obj.cFFT,x_fft_in);

            NumChannel = obj.pDOF;
            NumAngles = size(ang,2);
            
            phaseModulatedXFreq = complex(zeros(nfft,NumChannel));
            
            for m = 1:NumAngles
                
                arrayEffect_temp = sph2cartvec(...
                    [temp_phaseMtx_h(:,m).';temp_phaseMtx_v(:,m).';zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp = phased.internal.local2globalvec(...
                    arrayEffect_temp,lclaxes);
                
                arrayEffect_x = arrayEffect_temp(1,:).';
                arrayEffect_y = arrayEffect_temp(2,:).';
                arrayEffect_z = arrayEffect_temp(3,:).';
                
                phaseModulatedXFreq(:,m) = phaseModulatedXFreq(:,m) + ...
                    arrayEffect_x.*xtempfreq_x(:,m) + ...
                    arrayEffect_y.*xtempfreq_y(:,m) + ...
                    arrayEffect_z.*xtempfreq_z(:,m);
                
            end
            
        end
        
        function phaseModulatedXFreq = ...
                singleElementFreqDomainModulate(obj,x,ang,nfft,w)  %#ok<INUSL>
        % collect all signals with one element and sum them up

            % Calculate the element response for each subband at given
            % angles
            phaseMtx_temp = step(obj.cSensor,abs(obj.pSubbandFreqs),ang);
            if isPolarizationEnabled(obj.cSensor)
                if obj.pIsPowerSignal
                    intresp = step(obj.cIntegratedPattern,abs(obj.pSubbandFreqs));
                    if obj.WeightsInputPort
                        intresp.H = intresp.H*abs(w)^2;
                        intresp.V = intresp.V*abs(w)^2;
                    end
                    phaseMtx_temp.H = phased.internal.normalizeIntegratedPower(phaseMtx_temp.H,intresp.H,false);
                    phaseMtx_temp.V = phased.internal.normalizeIntegratedPower(phaseMtx_temp.V,intresp.V,false);
                end
                phaseMtx = hypot(phaseMtx_temp.H,phaseMtx_temp.V).';
            else
                if obj.pIsPowerSignal
                    intresp = step(obj.cIntegratedPattern,abs(obj.pSubbandFreqs));
                    if obj.WeightsInputPort
                        intresp = intresp*abs(w)^2;
                    end
                    phaseMtx_temp = phased.internal.normalizeIntegratedPower(phaseMtx_temp,intresp,false);
                end
                phaseMtx = phaseMtx_temp.';
            end

            % Modulate to frequency domain and modulate
            if isreal(x)
                x_fft_in = complex(x);
            else
                x_fft_in = x;
            end
            xtempfreq = step(obj.cFFT,x_fft_in);
            phaseModulatedXFreq = phaseMtx.*xtempfreq;

            % Combine all signals from different directions together
            phaseModulatedXFreq = sum(phaseModulatedXFreq,2);

        end

        function phaseModulatedXFreq = ...
                singleElementFreqDomainPolarizedModulate(obj,x,ang,lclaxes,nfft,w)
        % collect all signals with one element and sum them up

            % Calculate the element response for each subband at given
            % angles
            phaseMtx = step(obj.cSensor,abs(obj.pSubbandFreqs),ang);
            if obj.pIsPowerSignal
                intresp = step(obj.cIntegratedPattern,abs(obj.pSubbandFreqs));
                if obj.WeightsInputPort
                    intresp.H = intresp.H*abs(w)^2;
                    intresp.V = intresp.V*abs(w)^2;
                end
                phaseMtx.H = phased.internal.normalizeIntegratedPower(phaseMtx.H,intresp.H,false);
                phaseMtx.V = phased.internal.normalizeIntegratedPower(phaseMtx.V,intresp.V,false);
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
            x_temp = x_X;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_x = step(obj.cFFT,x_fft_in);
            
            x_temp = x_Y;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_y = step(obj.cFFT,x_fft_in);
            
            x_temp = x_Z;
            if isreal(x_temp)
                x_fft_in = complex(x_temp);
            else
                x_fft_in = x_temp;
            end
            xtempfreq_z = step(obj.cFFT,x_fft_in);
            
            NumAngles = size(ang,2);
            
            phaseModulatedXFreq = complex(zeros(nfft,1));
            
            for m = 1:NumAngles
                
                temp_phaseMtx_h = phaseMtx.H(m,:);
                temp_phaseMtx_v = phaseMtx.V(m,:);
                
                arrayEffect_temp = sph2cartvec(...
                    [temp_phaseMtx_h(:).';temp_phaseMtx_v(:).';zeros(1,nfft)],...
                    ang(1,m),ang(2,m));
                arrayEffect_temp = phased.internal.local2globalvec(...
                    arrayEffect_temp,lclaxes);
                
                arrayEffect_x = arrayEffect_temp(1,:).';
                arrayEffect_y = arrayEffect_temp(2,:).';
                arrayEffect_z = arrayEffect_temp(3,:).';
                
                phaseModulatedXFreq = phaseModulatedXFreq + ...
                    bsxfun(@times,arrayEffect_x,xtempfreq_x(:,m)) + ...
                    bsxfun(@times,arrayEffect_y,xtempfreq_y(:,m)) + ...
                    bsxfun(@times,arrayEffect_z,xtempfreq_z(:,m));
                
            end

        end

        function y = overlapadd(obj,x,NumOverlap)
            % add overlap from last block
            x(1:NumOverlap,:) = x(1:NumOverlap,:)+obj.pBuffer;
            M = size(x,1) - NumOverlap;
            obj.pBuffer = x(M+(1:NumOverlap),:);
            y = x(1:M,:);
        end

    end

    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractSensorOperation;
            dWavefront = ...
                matlab.system.display.internal.Property('Wavefront', ...
                                                        'IsGraphical', false);
            dEnablePolarization = ...
                matlab.system.display.internal.Property('EnablePolarization', ...
                                                        'IsGraphical', false);
            dSampleRate = matlab.system.display.internal.Property(...
                'SampleRate','IsObjectDisplayOnly',true);
            props = {...
                dSampleRate,...
                'ModulatedInput',...
                'CarrierFrequency',...
                'SignalType',...
                dWavefront,...
                dEnablePolarization,...
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
            str = sprintf('Wideband\n Rx Array');
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
        function varargout = isOutputComplexImpl(obj)
            if obj.ModulatedInput
                varargout{1} = true;
            else
                varargout{1} = false;
            end
        end
    end    
    
end
