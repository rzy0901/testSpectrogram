classdef (Sealed,StrictDefaults) WidebandRadiator < phased.internal.AbstractSensorOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime
%Radiator Wideband signal radiator
%   H = phased.WidebandRadiator creates a wideband signal radiator System
%   object, H. The object returns radiated wideband signals for given
%   directions using a sensor array or a single element.
%
%   H = phased.WidebandRadiator(Name,Value) creates a wideband radiator
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Y = step(H,X,ANG) radiates signal X to multiple directions specified in
%   ANG, and combines the radiated signal for each direction in Y. The
%   signal is divided into subbands and the combination is performed using
%   the phase approximation of time delays across radiating elements at the
%   far field in each subband.
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
%   Y = step(...,W) uses W as the weight vector when the
%   WeightsInputPort property is true. W must be a column vector whose
%   length is the same as the number of radiating elements or subarrays.
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
%   Y = step(H,X,ANG,LAXES,W,STEER)
%
%   or
%
%   Y = step(H,X,ANG,LAXES,W,WS)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   WidebandRadiator methods:
%
%   step     - Radiate signals (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a wideband radiator object with same property values
%   isLocked - Locked status (logical)
%
%   WidebandRadiator properties:
%
%   Sensor                      - Handle of the sensor
%   PropagationSpeed            - Propagation speed
%   SampleRate                  - Sample rate
%   CarrierFrequency            - Carrier frequency
%   NumSubbands                 - Number of subbands
%   SensorGainMeasure           - Sensor gain measure 
%   Polarization                - Polarization configuration
%   WeightsInputPort            - Enable weights input
%
%   % Examples:
%
%   % Example 1: 
%   %   Radiate signal with a single antenna.
%
%   ant = phased.IsotropicAntennaElement;
%   sigrad = phased.WidebandRadiator('Sensor',ant,...
%           'CarrierFrequency',300e6);
%   x = rand(10,1)+1i*rand(10,1);
%   radiatingAngle = [30 10]';
%   y = sigrad(x,radiatingAngle);
%
%   % Example 2: 
%   %   Radiate a far field signal with a 5-element array.
%
%   ant = phased.ULA('NumElements',5);
%   sigrad = phased.WidebandRadiator('Sensor',ant,...
%           'CarrierFrequency',300e6);
%   x = rand(10,1)+1i*rand(10,1);
%   radiatingAngle = [30 10; 20 0]'; % two directions
%   y = sigrad(x,radiatingAngle);
%
%   See also phased, phased.Radiator.

%   Copyright 2015-2017 The MathWorks, Inc.
%     

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] V. N. Bringi and V. Chandrasekar, Polarimetric Doppler Weather
%   Radar: Principles and Applications, Cambridge University Press, 2001


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %SampleRate Sample rate (Hz)
    %   Specify the sample rate (in Hz) as a positive scalar. The default
    %   value is 1e6 (1 MHz).
    SampleRate = 1e6
    %CarrierFrequency Carrier frequency (Hz)
    %   Specify the carrier frequency (in Hz) as a positive scalar. The
    %   default value is 1e9 (1 GHz). 
    CarrierFrequency = 1e9
end

properties (Nontunable, PositiveInteger)
    %NumSubbands    Number of subbands
    %   Specify the number of subbands used in the subband processing
    %   as a positive integer. The default value of this property is
    %   64.
    NumSubbands = 64
end

properties (Access = private, Nontunable, Logical) 
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
    %   'dB' is the default value. 
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
    %EnableDualPolarizationInput  Enable dual polarization input signal
    %   Set this property to true to enable dual polarization input signal
    %   in separate H and V polarization. Set this property to false to use
    %   same signal for both polarizations. The default value of this
    %   property is false. This property applies when you set the
    %   EnablePolarization property to true.
    pEnableDualPolarizationInput 
end

properties (Nontunable, Logical)
    %SampleRateFromInputCheckbox Inherit sample rate 
    %   Set SampleRateFromInputCheckbox to true to derive sample rate from
    %   Simulink time engine. Set SampleRateFromInputCheckbox to false to
    %   specify the sample rate. This property applies when used in
    %   Simulink.
    SampleRateFromInputCheckbox = true
end
    
properties (Constant, Hidden)
    SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
        'SystemBlock', 'SampleRateFromInputCheckbox',...
        'getSampleRateInSimulation',false})
    PolarizationSet = matlab.system.StringSet({'None','Combined','Dual'});
    SensorGainMeasureSet = matlab.system.StringSet({'dB','dBi'});
end


properties (Access = private, Nontunable)
    pFFTLength      % number of FFT
    pSubbandFreqs   % center frequency for each FFT bin
    cSteeringVector;
    cIntegratedPattern
    pDOF;
    pSampleRate;
    pCommonSignalForChannels;
    pPowerDistributionFactor;
    pApplyDirectivityGain
end

properties (Access = private)
    pNumOverlap     % number of overlapped samples in overlap-add
    pBuffer         % buffer for overlap-add
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
    
    function value = get.pEnableDualPolarizationInput(obj)
        value = strcmp(obj.Polarization,'Dual');
    end
    
end

methods
    function obj = WidebandRadiator(varargin)
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
        
        %angle
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
    
    function setupImpl(obj,x,ang,~,~,~,~)
        setupImpl@phased.internal.AbstractSensorOperation(obj);
        if obj.pUseArray
            M = getDOF(obj.cSensor);
        else
            M = 1;            
        end
        obj.pDOF = M;
        
        obj.pNumInputChannels = getNumChannels(obj,x);
        obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        
        % obj.pSampleRate = getSampleRate(obj,xsize(1),1,obj.SampleRate);
        fs = obj.SampleRate; % property/method duality
        cond = ~isscalar(fs) || (fs<=0);
        if cond
            coder.internal.errorIf(cond,...
                 'phased:phased:invalidSampleTime');
        end
        obj.pSampleRate = fs;

        cond = obj.pSampleRate > 2*obj.CarrierFrequency;
        if cond
            coder.internal.errorIf(cond,'phased:phased:WidebandCollector:FsTooHigh');
        end
            
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
        obj.pSubbandFreqs = phased.internal.subbandCenterFrequency(...
            obj.CarrierFrequency,obj.pSampleRate,obj.pFFTLength).';
            
        if isa(obj.Sensor,'phased.internal.AbstractSubarray')
            obj.pPowerDistributionFactor = getNumElements(obj.cSensor,...
                1:getNumSubarrays(obj.cSensor));
        else
            obj.pPowerDistributionFactor = 1;
        end
        
        if ~strcmp(obj.Polarization,'None')
            obj.pBuffer = complex(zeros(nfft,size(ang,2),2));
        else
            obj.pBuffer = complex(zeros(nfft,size(ang,2)));
        end
        
        processInputSizeChangeImpl(obj,x);
    end
    
    function processInputSizeChangeImpl(obj,x,~,~,~,~)
        obj.pNumOverlap = obj.pFFTLength - size(x,1);
        obj.pBuffer = complex(zeros(size(obj.pBuffer)));
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
            end
        else
            reset(obj.cSensor);
        end
        obj.pBuffer = complex(zeros(size(obj.pBuffer)));
        reset(obj.cSubbandDivider);
        reset(obj.cSubbandCombiner);
    end
    
    function releaseImpl(obj)
        releaseImpl@phased.internal.AbstractSensorOperation(obj);
        if obj.pUseArray
            if obj.CombineRadiatedSignals
                release(obj.cSteeringVector);
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
            s.pCommonSignalForChannels = obj.pCommonSignalForChannels;
            s.pPowerDistributionFactor = obj.pPowerDistributionFactor;
            s.CombineRadiatedSignals = obj.CombineRadiatedSignals; % future
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

    function s = loadSubObjects(obj,s,wasLocked)
        s = loadSubObjects@phased.internal.AbstractSensorOperation(obj,s);
        if isfield(s,'isLocked')
            s = rmfield(s,'isLocked');
        end
        if wasLocked
            obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
            s = rmfield(s,'cSteeringVector');
            obj.cSubbandDivider = phased.internal.SubbandDivider.loadobj(s.cSubbandDivider);
            s = rmfield(s,'cSubbandDivider');
            obj.cSubbandCombiner = phased.internal.SubbandDivider.loadobj(s.cSubbandCombiner);
            s = rmfield(s,'cSubbandCombiner');
            if isfield(s,'cIntegratedPattern') && ...
                        isfield(s.cIntegratedPattern,'ClassNameForLoadTimeEval')
                obj.cIntegratedPattern = eval(...
                    sprintf('%s.loadobj(s.cIntegratedPattern)',s.cIntegratedPattern.ClassNameForLoadTimeEval));
                s = rmfield(s,'cIntegratedPattern');
            end
            % recover locked sample rate information
            obj.pSampleRate = s.pSampleRate;
            s = rmfield(s,'pSampleRate');
        end
    end
    
    function loadObjectImpl(obj,s,wasLocked) 
        s = loadSubObjects(obj,s,wasLocked);
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
                                [yf_h, yf_v] = farFieldFreqDomainDualPolarizedModulate(...
                                    obj,xh,xv,ang,stang);
                            else
                                stang = lclaxesArg;
                                [yf_h, yf_v] = farFieldFreqDomainPolarizedModulate(...
                                    obj,x,ang,stang);
                            end
                        else
                            if obj.pEnableDualPolarizationInput
                                stang = stangArg;
                                w = wArg;
                                [yf_h, yf_v] = farFieldFreqDomainDualPolarizedModulate(...
                                    obj,xh,xv,ang,w,stang);
                            else
                                stang = wArg;
                                w = lclaxesArg;
                                [yf_h, yf_v] = farFieldFreqDomainPolarizedModulate(...
                                    obj,x,ang,w,stang);
                            end
                        end
                    else
                        if obj.WeightsInputPort
                            if obj.pEnableDualPolarizationInput
                                w = wArg;
                                [yf_h, yf_v] = farFieldFreqDomainDualPolarizedModulate(...
                                    obj,xh,xv,ang,w);
                            else
                                w = lclaxesArg;
                                [yf_h, yf_v] = farFieldFreqDomainPolarizedModulate(...
                                    obj,x,ang,w);
                            end
                        else
                            if obj.pEnableDualPolarizationInput
                                [yf_h, yf_v] = farFieldFreqDomainDualPolarizedModulate(...
                                    obj,xh,xv,ang);
                            else
                                [yf_h, yf_v] = farFieldFreqDomainPolarizedModulate(...
                                    obj,x,ang);
                            end
                        end
                    end
                end
            else  % single sensor
                if obj.WeightsInputPort
                    if obj.pEnableDualPolarizationInput
                        w = wArg;
                        [yf_h, yf_v] = singleElementFreqDomainDualPolarizedModulate(...
                            obj,xh,xv,ang,w);
                    else
                        w = lclaxesArg;
                        [yf_h, yf_v] = singleElementFreqDomainPolarizedModulate(...
                            obj,x,ang,w);
                    end
                else
                    if obj.pEnableDualPolarizationInput
                        [yf_h, yf_v] = singleElementFreqDomainDualPolarizedModulate(...
                            obj,xh,xv,ang);
                    else
                        [yf_h, yf_v] = singleElementFreqDomainPolarizedModulate(...
                            obj,x,ang);
                    end
                end
            end
            
            if obj.pEnableDualPolarizationInput
                y_h = step(obj.cSubbandCombiner,complex(yf_h),complex(xh));
                y_v = step(obj.cSubbandCombiner,complex(yf_v),complex(xv));
            else
                y_h = step(obj.cSubbandCombiner,complex(yf_h),complex(x));
                y_v = step(obj.cSubbandCombiner,complex(yf_v),complex(x));
                % y_h = overlapadd(obj,y_h,obj.pNumOverlap,'H');
                % y_v = overlapadd(obj,y_v,obj.pNumOverlap,'V');
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
                            yf = farFieldFreqDomainModulate(...
                                obj,x,ang,stang);
                        else
                            stang = lclaxesArg;
                            w = angArg;
                            yf = farFieldFreqDomainModulate(...
                                obj,x,ang,w,stang);
                        end
                    else
                        if obj.WeightsInputPort
                            w = angArg;
                            yf = farFieldFreqDomainModulate(...
                                obj,x,ang,w);
                        else
                            yf = farFieldFreqDomainModulate(...
                                obj,x,ang);
                        end
                    end
                    
                end
            else  % single sensor
                if obj.WeightsInputPort
                    w = angArg;
                    yf = singleElementFreqDomainModulate(...
                        obj,x,ang,w);
                else
                    yf = singleElementFreqDomainModulate(...
                        obj,x,ang);
                end
                
            end
            
            y = step(obj.cSubbandCombiner,complex(yf),complex(x));
            % y = overlapadd(obj,y,obj.pNumOverlap);

        end
        
    end
           
end

methods (Access = private)
    function y = overlapadd(obj,x,NumOverlap,pol)
        % add overlap from last block
        buf = obj.pBuffer;
        if ~strcmp(obj.Polarization,'None')
            if strcmp(pol,'H')
                idx = 1;
            else % strcmp(pol,'V')
                idx = 2;
            end
            x(1:NumOverlap,:) = x(1:NumOverlap,:)+buf(1:NumOverlap,:,idx);
            M = size(x,1) - NumOverlap;
            buf(1:NumOverlap,:,idx) = x(M+(1:NumOverlap),:);
            y = x(1:M,:);
        else
            x(1:NumOverlap,:) = x(1:NumOverlap,:)+buf(1:NumOverlap,:);
            M = size(x,1) - NumOverlap;
            buf(1:NumOverlap,:) = x(M+(1:NumOverlap),:);
            y = x(1:M,:);
        end
        obj.pBuffer = buf;
    end
    
    function [yf_h, yf_v] = farFieldFreqDomainPolarizedModulate(...
            obj,x,ang,wArg,stangArg)
        % Calculate steering vector for different subbands and angles
        % in one step. In case of unmodulated input, because there are
        % negative frequencies, current steering vector does not
        % support it when including the element response, therefore, we
        % get the steering vector without element response first and
        % then add in element response manually.
        subfreq = obj.pSubbandFreqs;

        if obj.pNeedSteeringAngle
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
            sv = step(obj.cSteeringVector,subfreq,ang,stang);
        else
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
            sv = step(obj.cSteeringVector,subfreq,ang);
        end
        pfactor = 1./sqrt(obj.pPowerDistributionFactor(:));
        if obj.WeightsInputPort
            sv_h = bsxfun(@times,pfactor.*w,sv.H);  % chan x ang x freq
            sv_v = bsxfun(@times,pfactor.*w,sv.V);  % chan x ang x freq
        else
            sv_h = bsxfun(@times,pfactor,sv.H);  % chan x ang x freq
            sv_v = bsxfun(@times,pfactor,sv.V);  % chan x ang x freq
        end
        if obj.pApplyDirectivityGain
            for m = 1:numel(subfreq)
                sv_h(:,:,m) = phased.internal.normalizeIntegratedPower(...
                    sv_h(:,:,m),intresp.H(m),false);
                sv_v(:,:,m) = phased.internal.normalizeIntegratedPower(...
                    sv_v(:,:,m),intresp.V(m),false);
            end
        end
        
        xtempfreq = step(obj.cSubbandDivider,complex(x));  % snap x chan x freq
        
        yf_h = complex(zeros(size(xtempfreq,1),size(ang,2),size(xtempfreq,3)));
        yf_v = complex(zeros(size(xtempfreq,1),size(ang,2),size(xtempfreq,3)));
        for m = 1:size(xtempfreq,1)
            yf_h(m,:,:) = sum(bsxfun(...
                @times,sv_h,permute(xtempfreq(m,:,:),[2 1 3])),1);  % snap x ang x freq
            yf_v(m,:,:) = sum(bsxfun(...
                @times,sv_v,permute(xtempfreq(m,:,:),[2 1 3])),1);  % snap x ang x freq
        end
        
    end
    
    function [yf_h, yf_v] = farFieldFreqDomainDualPolarizedModulate(...
            obj,xh,xv,ang,wArg,stangArg)
        % Calculate steering vector for different subbands and angles
        % in one step. In case of unmodulated input, because there are
        % negative frequencies, current steering vector does not
        % support it when including the element response, therefore, we
        % get the steering vector without element response first and
        % then add in element response manually.
        subfreq = obj.pSubbandFreqs;

        if obj.pNeedSteeringAngle
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
            sv = step(obj.cSteeringVector,subfreq,ang,stang);
        else
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
            sv = step(obj.cSteeringVector,subfreq,ang);
        end
        pfactor = 1./sqrt(obj.pPowerDistributionFactor(:));
        if obj.WeightsInputPort
            sv_h = bsxfun(@times,pfactor.*w,sv.H);  % chan x ang x freq
            sv_v = bsxfun(@times,pfactor.*w,sv.V);  % chan x ang x freq
        else
            sv_h = bsxfun(@times,pfactor,sv.H);  % chan x ang x freq
            sv_v = bsxfun(@times,pfactor,sv.V);  % chan x ang x freq
        end
        if obj.pApplyDirectivityGain
            for m = 1:numel(subfreq)
                sv_h(:,:,m) = phased.internal.normalizeIntegratedPower(...
                    sv_h(:,:,m),intresp.H(m),false);
                sv_v(:,:,m) = phased.internal.normalizeIntegratedPower(...
                    sv_v(:,:,m),intresp.V(m),false);
            end
        end
        
        xtempfreq_h = step(obj.cSubbandDivider,complex(xh));  % snap x chan x freq
        xtempfreq_v = step(obj.cSubbandDivider,complex(xv));  % snap x chan x freq
        
        yf_h = complex(zeros(size(xtempfreq_h,1),size(ang,2),size(xtempfreq_h,3)));
        yf_v = complex(zeros(size(xtempfreq_v,1),size(ang,2),size(xtempfreq_v,3)));
        for m = 1:size(xtempfreq_h,1)
            yf_h(m,:,:) = sum(bsxfun(...
                @times,sv_h,permute(xtempfreq_h(m,:,:),[2 1 3])),1);  % snap x ang x freq
            yf_v(m,:,:) = sum(bsxfun(...
                @times,sv_v,permute(xtempfreq_v(m,:,:),[2 1 3])),1);  % snap x ang x freq
        end
        
    end
    
    function [yf_h, yf_v] = singleElementFreqDomainPolarizedModulate(...
            obj,x,ang,w)
        % Single sensor
        subfreq = obj.pSubbandFreqs;
        
        g = step(obj.cSensor,subfreq,ang);     % chan (chan==1) x ang x freq
        g_h = g.H;
        g_v = g.V;
        if obj.WeightsInputPort
            g_h = w*g_h;
            g_v = w*g_v;
            if obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);  % scalar w
                intresp.H = w^2*intresp.H;
                intresp.V = w^2*intresp.V;
                for m = 1:numel(subfreq)
                    g_h(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_h(:,m),intresp.H(m),false);
                    g_v(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_v(:,m),intresp.V(m),false);
                end
            end
        else
            if obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);
                for m = 1:numel(subfreq)
                    g_h(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_h(:,m),intresp.H(m),false);
                    g_v(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_v(:,m),intresp.V(m),false);
                end
            end
        end
        
        xtempfreq = step(obj.cSubbandDivider,complex(x));  % snap x chan (chan==1) x freq
        
        yf_h = bsxfun(@times,permute(g_h,[3 1 2]),xtempfreq); % snap x ang x freq
        yf_v = bsxfun(@times,permute(g_v,[3 1 2]),xtempfreq); % snap x ang x freq
    end
    
    function [yf_h, yf_v] = singleElementFreqDomainDualPolarizedModulate(...
            obj,xh,xv,ang,w)
        % Single sensor
        subfreq = obj.pSubbandFreqs;
        
        g = step(obj.cSensor,subfreq,ang);     % chan (chan==1) x ang x freq
        g_h = g.H;
        g_v = g.V;
        if obj.WeightsInputPort
            g_h = w*g_h;
            g_v = w*g_v;
            if obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);  % scalar w
                intresp.H = w^2*intresp.H;
                intresp.V = w^2*intresp.V;
                for m = 1:numel(subfreq)
                    g_h(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_h(:,m),intresp.H(m),false);
                    g_v(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_v(:,m),intresp.V(m),false);
                end
            end
        else
            if obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);
                for m = 1:numel(subfreq)
                    g_h(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_h(:,m),intresp.H(m),false);
                    g_v(:,m) = phased.internal.normalizeIntegratedPower(...
                        g_v(:,m),intresp.V(m),false);
                end
            end
        end
        
        xtempfreq_h = step(obj.cSubbandDivider,complex(xh));  % snap x chan (chan==1) x freq
        xtempfreq_v = step(obj.cSubbandDivider,complex(xv));  % snap x chan (chan==1) x freq
        
        yf_h = bsxfun(@times,permute(g_h,[3 1 2]),xtempfreq_h); % snap x ang x freq
        yf_v = bsxfun(@times,permute(g_v,[3 1 2]),xtempfreq_v); % snap x ang x freq
    end
    
    function yf = farFieldFreqDomainModulate(...
            obj,x,ang,wArg,stangArg)
        % Calculate steering vector for different subbands and angles
        % in one step. In case of unmodulated input, because there are
        % negative frequencies, current steering vector does not
        % support it when including the element response, therefore, we
        % get the steering vector without element response first and
        % then add in element response manually.
        subfreq = obj.pSubbandFreqs;

        if obj.pNeedSteeringAngle
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
            sv = step(obj.cSteeringVector,subfreq,ang,stang);
        else
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
            sv = step(obj.cSteeringVector,subfreq,ang);
        end
        
        pfactor = 1./sqrt(obj.pPowerDistributionFactor(:));
        if obj.WeightsInputPort
            sv = bsxfun(@times,pfactor.*w,sv);  % chan x ang x freq
        else
            sv = bsxfun(@times,pfactor,sv);  % chan x ang x freq
        end
        if obj.pApplyDirectivityGain
            for m = 1:numel(subfreq)
                sv(:,:,m) = phased.internal.normalizeIntegratedPower(...
                    sv(:,:,m),intresp(m),false);
            end
        end
        
        xtempfreq = step(obj.cSubbandDivider,complex(x));  % snap x chan x freq
        
        yf = complex(zeros(size(xtempfreq,1),size(ang,2),size(xtempfreq,3)));
        for m = 1:size(xtempfreq,1)
            yf(m,:,:) = sum(bsxfun(...
                @times,sv,permute(xtempfreq(m,:,:),[2 1 3])),1);  % snap x ang x freq
        end
        
    end
    
    function yf = singleElementFreqDomainModulate(...
            obj,x,ang,w)
        % Single sensor
        subfreq = obj.pSubbandFreqs;
        
        g_temp = step(obj.cSensor,subfreq,ang);  % chan (chan==1) x ang x freq
        if isPolarizationEnabled(obj.cSensor)
            g = hypot(g_temp.H,g_temp.V);
        else
            g = g_temp;
        end
        
        if obj.WeightsInputPort
            g = w*g;
            if obj.pApplyDirectivityGain
                intresp = w^2*step(obj.cIntegratedPattern,subfreq);  % scalar w
                for m = 1:numel(subfreq)
                    g(:,m) = phased.internal.normalizeIntegratedPower(...
                        g(:,m),intresp(m),false);
                end
            end
        else
            if obj.pApplyDirectivityGain
                intresp = step(obj.cIntegratedPattern,subfreq);
                for m = 1:numel(subfreq)
                    g(:,m) = phased.internal.normalizeIntegratedPower(...
                        g(:,m),intresp(m),false);
                end
            end
        end
        
        xtempfreq = step(obj.cSubbandDivider,complex(x));  % snap x chan (chan==1) x freq
        
        yf = bsxfun(@times,permute(g,[3 1 2]),xtempfreq); % snap x ang x freq

    end
end

methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractSensorOperation;
        
        dPolarization = ...
            matlab.system.display.internal.Property('Polarization', ...
            'IsGraphical', false);
        % dSampleRate = matlab.system.display.internal.Property(...
        %     'SampleRate','IsObjectDisplayOnly',true);
        
        props = {'SampleRateFromInputCheckbox',...
            'SampleRate',...
            'CarrierFrequency',...
            'NumSubbands',...
            'SensorGainMeasure',...
            dPolarization,...
            'WeightsInputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
    end
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:WidebandRadiatorTitle')),...
            'Text',getString(message('phased:library:block:WidebandRadiatorDesc')));
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
        str = sprintf('Wideband\n Tx Array');
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

