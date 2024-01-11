classdef (Hidden) AbstractClutterSimulator < phased.internal.AbstractSampleRateEngine & ...
     matlab.system.mixin.Propagates & ...
     matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%AbstractClutterSimulator   Define the abstract clutter simulator class.
%   This is the abstract class for clutter simulation

%   Copyright 2010-2017 The MathWorks, Inc.

%   Assumptions:
%   1. Monostatic
%   2. Free space, propagation is the same for different directions
%   3. Homogeneous terrain
%   4. Clutter patch stationary from pulse to pulse
%   5. Narrowband
%   6. Do not update altitude, don't know how to convert to altitude
%   7. Constant speed
%   8. Patch never changes, although random reflectivity can be regenerated

%   Reference

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
    %OperatingFrequency     Operating frequency (Hz)
    %   Specify the operating frequency (in Hz) of the system as a positive
    %   scalar. The default value of this property is 3e8, i.e., 300 MHz.
    OperatingFrequency = 3e8;
    %SampleRate Sample rate (Hz)
    %   Specify the sample rate (in Hz) as a positive scalar. The default
    %   value of this property is 1e6 (1 MHz).
    SampleRate = 1e6;
    %PRF    Pulse repetition frequency (Hz)
    %   Specify the pulse repetition frequency (in Hz) as a positive scalar
    %   or a row vector. The default value of this property is 1e4 (10
    %   kHz). When PRF is a vector, it represents the case of staggered
    %   PRF, where the output pulses use elements in the vector as their
    %   PRFs one after another in cycle.
    PRF = 1e4;
    %EarthModel     Earth model
    %   Specify the earth model used in clutter simulation as one of 'Flat'
    %   | 'Curved', where the default is 'Flat'. When you set this property
    %   to 'Flat', the earth is assumed to be a flat plane. When you set
    %   this property to 'Curved', the earth is assumed to be a sphere.
    EarthModel = 'Flat';
    %PlatformHeight     Radar height (m)
    %   Specify the radar platform height (in meters) measured upward from
    %   the surface as a nonnegative scalar. The default value of this
    %   property is 300.
    PlatformHeight = 300;
    %PlatformSpeed  Radar speed (m/s)
    %   Specify the radar platform's speed as a nonnegative scalar (in
    %   m/s). The default value of this property is 300.
    PlatformSpeed = 300;  
    %PlatformDirection   Radar motion direction (deg)
    %   Specify the direction of radar platform motion as a 2x1 vector in
    %   the form of [AzimuthAngle; ElevationAngle] (in degrees). The
    %   default value of this property is [90;0], indicating that the
    %   platform moves perpendicular to the radar antenna array's
    %   broadside.
    %
    %   Both azimuth and elevation angle are measured in the local
    %   coordinate system of the radar antenna or antenna array. Azimuth
    %   angle must be between -180 and 180 degrees and elevation angle must
    %   be between -90 and 90 degrees.
    PlatformDirection = [90;0]; % in local coordinate
    %BroadsideDepressionAngle    Broadside depression angle (deg)
    %   Specify the depression angle (in degrees) of the radar antenna
    %   array's broadside as a scalar. The broadside is defined as zero
    %   degrees azimuth and zero degrees elevation. The depression angle is
    %   measured downward from horizontal. The default value of this
    %   property is 0.
    BroadsideDepressionAngle = 0;
    %MaximumRange   Maximum range (m)
    %   Specify the maximum range (in meters) for the clutter simulation as
    %   a positive scalar. The maximum range must be greater than the value
    %   specified in the PlatformHeight property. The default value of this
    %   property is 5000.
    MaximumRange = 5000;
    %AzimuthCoverage    Azimuth coverage (deg)
    %   Specify the azimuth coverage (in degrees) as a positive scalar. 
    %   The clutter simulation covers a region having the specified azimuth
    %   span, symmetric to 0 degrees azimuth. Typically, all clutter 
    %   patches have their azimuth centers within the region, but the 
    %   PatchAzimuthWidth value can cause some patches to extend beyond the 
    %   region. The default value of this property is 60.
    AzimuthCoverage = 60;
    %PatchAzimuthWidth  Clutter patch azimuth span (deg)
    %   Specify the azimuth span (in degrees) of each clutter patch as a
    %   positive scalar. The default value of this property is 1.
    PatchAzimuthWidth = 1;
    %CoherenceTime  Clutter coherence time (s)
    %   Specify the coherence time (in seconds) for the clutter simulation
    %   as a positive scalar. After coherence time elapse, the random
    %   numbers used for the clutter simulation are updated at the next
    %   pulse. The default value of this property is inf, which means that
    %   the random numbers are never updated.
    CoherenceTime = inf;
    %TransmitERP     Effective transmitted power (W)
    %   Specify the transmitted effective radiated power (ERP) (in Watts)
    %   of the radar system as a positive scalar. This property applies
    %   when you set the TransmitSignalInputPort property to false. The
    %   default value of this property is 5000.
    TransmitERP = 5000;
    %OutputFormat     Output signal format
    %   Specify the format of the output signal as one of 'Pulses' |
    %   'Samples', where the default is 'Pulses'. When you set the
    %   OutputFormat property to 'Pulses', the output of the step method is
    %   in the form of multiple pulses, where the number of pulses is the
    %   value of the NumPulses property. When you set the OutputFormat
    %   property to 'Samples', the output is in the form of multiple
    %   samples, where the number of samples is the value of the NumSamples
    %   property.
    OutputFormat = 'Pulses';
    %SeedSource   Source of seed for random number generator
    %   Specify how the random numbers are generated as one of 'Auto' |
    %   'Property', where the default is 'Auto'. When you set this property
    %   to 'Auto', the random numbers are generated using the default
    %   MATLAB random number generator. When you set this property to
    %   'Property', a private random number generator is used with a seed
    %   specified by the value of the Seed property.
    %
    %   To use this object with Parallel Computing Toolbox software, set
    %   this property to 'Auto'.
    SeedSource = 'Auto';
    %Seed     Seed for random number generator
    %   Specify the seed for the random number generator as a nonnegative
    %   integer. The integer must be less than 2^32. This property applies
    %   when you set the SeedSource property to 'Property'. The default
    %   value of this property is 0.
    Seed = 0;
end

properties (Nontunable,PositiveInteger)
    %NumPulses  Number of pulses in output
    %   Specify the number of pulses in the output of the step method as a
    %   positive integer. This property applies only when you set the
    %   OutputFormat property to 'Pulses'. The default value of this
    %   property is 1.
    NumPulses = 1;
    %NumSamples     Number of samples in output
    %   Specify the number of samples in the output of the step method as a
    %   positive integer. This property applies only when you set the
    %   OutputFormat property to 'Samples'. The default value of this
    %   property is 100.
    NumSamples = 100;
end

properties (Logical, Nontunable)
    %TransmitSignalInputPort    Enable transmit signal input
    %   Set this property to true to add input to specify the transmit
    %   signal in the step method. Set this property to false to not
    %   specify the transmit signal. The default value of this property is
    %   false.
    TransmitSignalInputPort = false;
    %PRFSelectionInputPort  Enable PRF selection input
    %   Set this property to true to select which predefined PRF to use
    %   during the simulation via input. Set this property to false to use
    %   the PRF property to define the PRF sequence used in the simulation.
    %   The default value of this property is false.
    PRFSelectionInputPort = false;
end

properties(Constant, Hidden)
    EarthModelSet = matlab.system.StringSet({'Flat','Curved'});
    OutputFormatSet = matlab.system.StringSet({'Pulses','Samples'});
    SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
end

properties (Access = protected, Nontunable)
    cArrayGain;
    pWavelength;
    pRanges;
    pGrazingAngles;
    pDepressionAngles;
    pElAngles;
    pAzAngles;
    pPropagationLoss;
    pPatchAreas;
    pRangeDelta;
    pEffRangeIdx;
    pNumSamplesPerDistinctPulse;
    pRelativeRadialSpeed;
    pDopplerShift;
    cCollector;
    cNoiseSource;
    pPlatformVelocity;
end

properties (Access = protected, Nontunable, Logical)
    pUseArray;
    pNeverUpdateRandomSource;
    pOutputBySamples;
    pNeedSteeringAngle;
    pNeedCustomWeights;
end

properties (Access = protected, Nontunable, PositiveInteger)
    pNumRanges;
    pDOF;
    pBufferLength;
    pNumDistinctPulses;
    %pNumAzimuthPatches  Number of clutter patches in azimuth
    %   Specify the number of clutter patches in azimuth direction as a
    %   positive, odd integer. The default value of this property is 11.
    pNumAzimuthPatches; 
end

properties (Access = protected)
    pTimeOffset;
    pNRCS;
    pTxGain;
    pClutterDeterministicReturn;
    pClutterReturn;
    pClutterAtArrayOutputBuffer;
    pLastRandUpdateTime;
    pRemainingSamplesFromLastPulse;
    pProcessedPulseCount;
end

methods (Access = protected, Abstract)
    nrcs = getNRCS(obj);    % calculate NRCS
end

methods
    function set.Sensor(obj,value)
        validateattributes(value,{'phased.internal.AbstractElement',...
            'phased.internal.AbstractArray','em.Antenna'...
            'phased.internal.AbstractSubarray'},{'scalar'},'','Sensor');
        
        if isa(value,'phased.internal.AbstractElement') || isa(value,'em.Antenna')
            element = value;
        elseif isa(value,'phased.internal.AbstractArray')
            element = getElementHandle(value);
        elseif isa(value,'phased.ReplicatedSubarray')
            element = getElementHandle(value.Subarray);
        else      % isa(value,'phased.PartitionedArray')
            element = getElementHandle(value.Array);
        end
        if iscell(element)
            element = element{1};
        end
        validateattributes(element,...
            {'phased.internal.AbstractAntennaElement','em.Antenna'},...
            {'scalar'},...
            '','Sensor array element');
        obj.Sensor = value;
    end
    
    function set.PropagationSpeed(obj,value)
        sigdatatypes.validateSpeed(value,'','PropagationSpeed',...
            {'scalar','positive'});
        obj.PropagationSpeed = value;
    end
    
    function set.SampleRate(obj, value)
        sigdatatypes.validateFrequency(value,'','SampleRate',...
            {'scalar'});
        obj.SampleRate = value;
    end
    
    function set.PRF(obj,value)
        sigdatatypes.validateFrequency(value,'','PRF',{'row'});
        obj.PRF = value;
    end
    
    function set.OperatingFrequency(obj,value)
        sigdatatypes.validateFrequency(value,'','OperatingFrequency',...
            {'scalar'});
        obj.OperatingFrequency = value;
    end    
    
    function set.PlatformHeight(obj,value)
        sigdatatypes.validateDistance(value,'','PlatformHeight',...
            {'scalar'});
        obj.PlatformHeight = value;
    end
    
    function set.PlatformSpeed(obj,value)
        sigdatatypes.validateSpeed(value,'','PlatformSpeed',...
            {'scalar'});
        obj.PlatformSpeed = value;
    end
    
    function set.PlatformDirection(obj,value)
        sigdatatypes.validateAzElAngle(value,'','PlatformDirection',...
            {'size',[2 1]});
        obj.PlatformDirection = value;
    end
    
    function set.BroadsideDepressionAngle(obj,value)
        sigdatatypes.validateAngle(value,'','BroadsideDepressionAngle',...
            {'scalar'});
        obj.BroadsideDepressionAngle = value;
    end
    
    function set.MaximumRange(obj,value)
        sigdatatypes.validateDistance(value,'','MaximumRange',...
            {'scalar','positive'});
        obj.MaximumRange = value;
    end
    
    function set.AzimuthCoverage(obj,value)
        sigdatatypes.validateAngle(value,'','AzimuthCoverage',...
            {'scalar','positive','<=',360});
        obj.AzimuthCoverage = value;
    end
    
    function set.PatchAzimuthWidth(obj,value)
        sigdatatypes.validateAngle(value,'','PatchAzimuthWidth',...
            {'scalar','positive','<=',360});
        obj.PatchAzimuthWidth = value;
    end
    
    function set.CoherenceTime(obj,value)
        validateattributes(value,{'double'},{'nonnan','nonempty',...
            'positive','scalar'},'','CoherenceTime');
        obj.CoherenceTime = value;
    end

    function set.TransmitERP(obj, value)
        sigdatatypes.validatePower(value,'',...
            'TransmitERP',{'scalar','positive'});
        obj.TransmitERP = value;             
    end
    
    function set.Seed(obj,value)
        validateattributes(value,{'double'},{'scalar','nonnegative',...
            'finite','nonnan','nonempty'},'',...
            'Seed');
        obj.Seed = value;
    end
end

methods (Access = protected)

    function obj = AbstractClutterSimulator(varargin)        
        setProperties(obj,nargin,varargin{:});
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

end

methods (Access = protected)
    function num = getNumInputsImpl(obj)
        num = 0;
        if obj.TransmitSignalInputPort
            num = num+1;
        end
        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmp(obj.Sensor.SubarraySteering,'None',1)
            num = num+1;
        end
        if obj.PRFSelectionInputPort
            num = num+1;
        end    
    end
    
    function validatePropertiesImpl(obj)

        if isa(obj.Sensor,'phased.internal.AbstractArray')
            element = getElementHandle(obj.Sensor);
        elseif isa(obj.Sensor,'phased.ReplicatedSubarray')
            element = getElementHandle(obj.Sensor.Subarray);
        elseif isa(obj.Sensor,'phased.PartitionedArray')
            element = getElementHandle(obj.Sensor.Array);
        else
            element = [];
        end
        
        % if isempty implies sensor is an element
        % no need to re-validate
        if ~isempty(element)
            if iscell(element)
                element = element{1};
            end
            validateattributes(element,...
                               {'phased.internal.AbstractAntennaElement','em.Antenna'},{'scalar'},...
                               '','Sensor array element');
        end

        cond = any(rem(obj.SampleRate, obj.PRF));
        if cond
            coder.internal.errorIf(cond,'phased:Waveform:NeedRatioInteger',...
                'SampleRate', 'PRF');
        end
        
        cond = obj.MaximumRange <= obj.PlatformHeight;
        if cond
            coder.internal.errorIf(cond,'phased:phased:expectedGreaterThan',...
                'MaximumRange',sprintf('%5.2f',obj.PlatformHeight));
        end
        
        cond = obj.PatchAzimuthWidth > obj.AzimuthCoverage;
        if cond
            coder.internal.errorIf(cond,'phased:phased:expectedLessThanOrEqualTo',...
                'PatchAzimuthWidth','AzimuthCoverage');
        end
        
        if rem(obj.AzimuthCoverage, obj.PatchAzimuthWidth) && ...
                rem(floor(obj.AzimuthCoverage/obj.PatchAzimuthWidth),2) && ...
                isempty(coder.target)
            warning(message('phased:clutter:patchCenterOutOfBound',...
                'AzimuthCoverage','PatchAzimuthWidth'));
        end
    end

    function validateInputsImpl(obj,varargin)
        if obj.TransmitSignalInputPort
            x = varargin{1};
            cond = ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','double');
            end
            cond = ~iscolumn(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeColVector','X');
            end
        end

        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmp(obj.Sensor.SubarraySteering,'None',1)
            if obj.TransmitSignalInputPort
                stang = varargin{2};
            else
                stang = varargin{1};
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
            
        if obj.PRFSelectionInputPort
            if obj.TransmitSignalInputPort
                if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    ~strncmp(obj.Sensor.SubarraySteering,'None',1)
                    prfidx = varargin{3};
                else
                    prfidx = varargin{2};
                end
            else
                if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    ~strncmp(obj.Sensor.SubarraySteering,'None',1)
                    prfidx = varargin{2};
                else
                    prfidx = varargin{1};
                end
            end
            cond = ~isa(prfidx,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','PRFIdx','double');
            end
            cond = ~isscalar(prfidx) || isempty(prfidx);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeScalar','PRFIdx');
            end
        end
    end

    function setupImpl(obj,varargin)
        coder.extrinsic('phased.internal.AbstractClutterSimulator.initializeProperties');
        % initialize random source
        obj.cNoiseSource = phased.internal.NoiseSource(...
            'Distribution','Gaussian','OutputComplex',true,...
            'SeedSource',obj.SeedSource);
        
        if (obj.SeedSource(1) == 'P') %Property
            obj.cNoiseSource.Seed = obj.Seed;
        end


        % Initialize properties
        

        initProps  = coder.internal.const(obj.initializeProperties( ...
            obj.PropagationSpeed, ...
            obj.OperatingFrequency, ...
            obj.SampleRate , ...
            obj.PRF, ...
            obj.MaximumRange, ... 
            obj.EarthModel, ...  
            obj.PlatformHeight, ...
            obj.AzimuthCoverage, ...
            obj.PatchAzimuthWidth, ...
            obj.BroadsideDepressionAngle,...
            obj.PlatformDirection, ...
            obj.PlatformSpeed, ...
            obj.CoherenceTime, ...
            obj.OutputFormat));        

         obj.pWavelength = initProps.wavelength;
         obj.pNumSamplesPerDistinctPulse = initProps.numSamplesPerDistinctPulse;
         obj.pNumDistinctPulses = initProps.numDistinctPulses;
         obj.pRangeDelta = initProps.rangeDelta;
         obj.pRanges = initProps.ranges;
         obj.pEffRangeIdx = initProps.effRangeIdx;
         obj.pNumRanges = initProps.numRanges;
         obj.pNumAzimuthPatches = initProps.numAzimuthPatches;
         obj.pAzAngles = initProps.azAngles;
         obj.pGrazingAngles = initProps.grazingAngles;
         obj.pDepressionAngles = initProps.depressionAngles;
         obj.pElAngles = initProps.elAngles;
         obj.pPropagationLoss = initProps.propagationLoss;
         obj.pPatchAreas = initProps.patchAreas;
         obj.pRelativeRadialSpeed = initProps.relativeRadialSpeed;
         obj.pPlatformVelocity = initProps.platformVelocity;
         obj.pDopplerShift = initProps.dopplerShift;
         obj.pNeverUpdateRandomSource = initProps.neverUpdateRandomSource;
         obj.pOutputBySamples = initProps.outputBySamples;
           
         obj.pClutterReturn = complex(zeros(obj.pNumRanges,...
                                         obj.pNumAzimuthPatches));
        
        % if uses an array, instantiate ArrayGain object and specify
        % the number of radiating/collecting elements. Otherwise, it is
        % a single element.
        obj.pUseArray = isa(obj.Sensor,'phased.internal.AbstractArray') || ...
            isa(obj.Sensor,'phased.internal.AbstractSubarray');
        if obj.pUseArray
            obj.cArrayGain = phased.ArrayGain(...
                'SensorArray',obj.Sensor,...
                'PropagationSpeed',obj.PropagationSpeed);
            obj.pDOF = getDOF(obj.Sensor);
        else
            obj.pDOF = 1;
        end

        obj.pNeedSteeringAngle = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
            (strncmp(obj.Sensor.SubarraySteering,'Phase',1) || ...
            strncmp(obj.Sensor.SubarraySteering,'Time',1));

        obj.pNeedCustomWeights = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
            strncmp(obj.Sensor.SubarraySteering,'Cusotm',1);


        % Transmit gain for each clutter patch
        obj.pTxGain = zeros(obj.pNumRanges,obj.pNumAzimuthPatches);
        if obj.pUseArray
            if ~(obj.pNeedSteeringAngle || obj.pNeedCustomWeights)
                for m = 1:obj.pNumAzimuthPatches
                    ang = phased.internal.arbazel2azel(...
                        [obj.pAzAngles(m)*ones(1,obj.pNumRanges);...
                        obj.pElAngles]);
                    obj.pTxGain(:,m) = step(obj.cArrayGain,...
                        obj.OperatingFrequency,ang);
                end
            end
        else
            if isempty(coder.target)
                hac = cloneSensor(obj.Sensor);
            else
                if isElementFromAntenna(obj.Sensor)
                    coder.internal.errorIf(true, ...
                        'phased:system:element:AntennaToolboxCodegenNotSupported','em.Antenna','phased.CustomAntennaElement');
                end
                hac = clonecg(obj.Sensor);
            end
            for m = 1:obj.pNumAzimuthPatches
                ang = phased.internal.arbazel2azel(...
                    [obj.pAzAngles(m)*ones(1,obj.pNumRanges);...
                    obj.pElAngles]);
                if isPolarizationEnabled(hac)
                    resp_temp = step(hac,obj.OperatingFrequency,ang);
                    obj.pTxGain(:,m) = mag2db(hypot(...
                        resp_temp.H,resp_temp.V));
                else
                    obj.pTxGain(:,m) = mag2db(abs(step(hac,...
                        obj.OperatingFrequency,ang)));
                end
            end
        end

        obj.pNRCS = getNRCS(obj);
                
        rcseffect = zeros(max(obj.pEffRangeIdx),1);
        % Combine to get the total reflect gain from each range gate,
        % i.e., 4*pi*A/lambda^2
        rcseffect(obj.pEffRangeIdx,1) = db2mag(...
            aperture2gain(obj.pNRCS(obj.pEffRangeIdx).*...
            obj.pPatchAreas(obj.pEffRangeIdx),...
            obj.pWavelength));

        % clutter return for each patch
        obj.pClutterDeterministicReturn = ...
            zeros(obj.pNumRanges,obj.pNumAzimuthPatches);
        for m = 1:obj.pNumAzimuthPatches
            obj.pClutterDeterministicReturn(obj.pEffRangeIdx,m) = ...
                db2mag(obj.pTxGain(obj.pEffRangeIdx,m)).*...
                rcseffect(obj.pEffRangeIdx)./...
                db2mag(obj.pPropagationLoss(obj.pEffRangeIdx));
        end
        if ~obj.TransmitSignalInputPort
            obj.pClutterDeterministicReturn = ...
                obj.pClutterDeterministicReturn * sqrt(obj.TransmitERP);
        end
        
        % initialize clutter return of each range gate at the
        % receiving sensor/array
        
        if obj.TransmitSignalInputPort
            x = varargin{1};
            sz = size(x);
            obj.pBufferLength = obj.pNumRanges+sz(1)-1;
        else
            obj.pBufferLength = obj.pNumRanges;    
        end
        
        % collector object for the receiving sensor/array
        obj.cCollector = phased.Collector(...
            'Sensor',obj.Sensor,...
            'OperatingFrequency',obj.OperatingFrequency,...
            'PropagationSpeed',obj.PropagationSpeed);

    end
    
    function flag = isInputComplexityLockedImpl(obj,index)
      if obj.TransmitSignalInputPort && index==1
        flag = false;
      else
        flag = true;
      end
    end

    function flag = isOutputComplexityLockedImpl(~,~)
        flag = false;
    end
    
    function resetImpl(obj)
        reset(obj.cCollector);
        if obj.pUseArray
            reset(obj.cArrayGain);
        end
        reset(obj.cNoiseSource);
        
        obj.pClutterAtArrayOutputBuffer = complex(zeros(obj.pBufferLength,...
            obj.pDOF));

        % time offset for the starting time at each simulation step
        obj.pTimeOffset = 0;

        % initialize the time stamp for last randomness update
        obj.pLastRandUpdateTime = 0;
        
        % initialize pulse count
        obj.pRemainingSamplesFromLastPulse = 0;
        obj.pProcessedPulseCount = 0;

    end

    function releaseImpl(obj)
        release(obj.cCollector);
        if obj.pUseArray
            release(obj.cArrayGain);
        end
        release(obj.cNoiseSource);
    end

    function flag = isInactivePropertyImpl(obj,prop)
        if strcmp(prop,'Seed') && strcmp(obj.SeedSource,'Auto')
            flag = true;
        elseif strcmp(prop,'TransmitERP') && obj.TransmitSignalInputPort
            flag = true;
        elseif strcmp(prop,'NumSamples') && ...
                strcmp(obj.OutputFormat,'Pulses')
            flag = true;
        elseif strcmp(prop,'NumPulses') && ...
                strcmp(obj.OutputFormat,'Samples')
            flag = true;
        else
            flag = false;
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        s.Sensor = saveobj(obj.Sensor);
        s.isLocked = isLocked(obj);
        if isLocked(obj)
            s.pUseArray = obj.pUseArray;
            s.cNoiseSource = saveobj(obj.cNoiseSource);
            s.cCollector = saveobj(obj.cCollector);
            if obj.pUseArray
                s.cArrayGain = saveobj(obj.cArrayGain);
            end
            s.pWavelength = obj.pWavelength;
            s.pRanges = obj.pRanges;
            s.pGrazingAngles = obj.pGrazingAngles;
            s.pDepressionAngles = obj.pDepressionAngles;
            s.pElAngles = obj.pElAngles;
            s.pAzAngles = obj.pAzAngles;
            s.pTxGain = obj.pTxGain;
            s.pPropagationLoss = obj.pPropagationLoss;
            s.pNRCS = obj.pNRCS;
            s.pPatchAreas = obj.pPatchAreas;
            s.pRangeDelta = obj.pRangeDelta;
            s.pEffRangeIdx = obj.pEffRangeIdx;
            s.pNumSamplesPerDistinctPulse = obj.pNumSamplesPerDistinctPulse;
            s.pRelativeRadialSpeed = obj.pRelativeRadialSpeed;
            s.pDopplerShift = obj.pDopplerShift;
            s.pClutterDeterministicReturn = obj.pClutterDeterministicReturn;
            s.pPlatformVelocity = obj.pPlatformVelocity;
            s.pNeverUpdateRandomSource = obj.pNeverUpdateRandomSource;
            s.pOutputBySamples = obj.pOutputBySamples;
            s.pNumRanges = obj.pNumRanges;
            s.pDOF = obj.pDOF;
            s.pBufferLength = obj.pBufferLength;
            s.pNumDistinctPulses = obj.pNumDistinctPulses;
            s.pNumAzimuthPatches = obj.pNumAzimuthPatches; 
            s.pTimeOffset = obj.pTimeOffset;
            s.pClutterReturn = obj.pClutterReturn;
            s.pClutterAtArrayOutputBuffer = obj.pClutterAtArrayOutputBuffer;
            s.pLastRandUpdateTime = obj.pLastRandUpdateTime;
            s.pRemainingSamplesFromLastPulse = obj.pRemainingSamplesFromLastPulse;
            s.pProcessedPulseCount = obj.pProcessedPulseCount;
            s.pNeedSteeringAngle = obj.pNeedSteeringAngle;
            s.pNeedCustomWeights = obj.pNeedCustomWeights;
        end
    end

    function s = loadSubObjects(obj,s)
        if isa(s.Sensor,'em.Antena')
            obj.Sensor = em.EmStructures.loadobj(s.Sensor);
        else
            obj.Sensor = eval(...
                sprintf('%s.loadobj(s.Sensor)',s.Sensor.ClassNameForLoadTimeEval));
        end
        s = rmfield(s,'Sensor');
        if isfield(s,'isLocked')
            if s.isLocked
                obj.cNoiseSource = eval(...
                    sprintf('%s.loadobj(s.cNoiseSource)',s.cNoiseSource.ClassNameForLoadTimeEval));
                s = rmfield(s,'cNoiseSource');
                obj.cCollector = phased.Collector.loadobj(s.cCollector);
                s = rmfield(s,'cCollector');
                if s.pUseArray
                    obj.cArrayGain = phased.ArrayGain.loadobj(s.cArrayGain);
                    s = rmfield(s,'cArrayGain');
                end
                if isfield(s,'pNumSensorElements')
                    s = rmfield(s,'pNumSensorElements');
                end
            end
            s = rmfield(s,'isLocked');
        end
    end
    
    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function updateClutterRandomReturn(obj)
        if (obj.pTimeOffset-obj.pLastRandUpdateTime) >= obj.CoherenceTime
            % update based current time offset, not coherent time after
            % last time offset
            calcRandomClutter(obj)
            obj.pLastRandUpdateTime = obj.pTimeOffset;
        end
    end
    
    function N = calcTotalOutputSamples(obj,prfidx)
        if obj.pOutputBySamples
            N = obj.NumSamples;
        elseif obj.PRFSelectionInputPort
            outputPulseIdx = ones(1,obj.NumPulses)*prfidx;
            N = sum(obj.pNumSamplesPerDistinctPulse(outputPulseIdx));
        else
            outputPulseIdx = (1:obj.NumPulses)+obj.pProcessedPulseCount;
            outputDistinctPulseIdx = mod(outputPulseIdx-1,...
                obj.pNumDistinctPulses)+1;
            N = sum(obj.pNumSamplesPerDistinctPulse(outputDistinctPulseIdx));
        end
    end
    
    function N = getNumSamplesForCurrentPulse(obj,prfidx)
        if obj.PRFSelectionInputPort
            currentPulseIdx = prfidx;
        else
            currentPulseIdx = mod(obj.pProcessedPulseCount,...
                obj.pNumDistinctPulses)+1;
        end
        N = obj.pNumSamplesPerDistinctPulse(currentPulseIdx);
    end
end

methods (Access = protected)
    function calcRandomClutter(obj)
        % generate random reflectivity and update random clutter for each
        % patch
        obj.pClutterReturn = complex(zeros(obj.pNumRanges,...
            obj.pNumAzimuthPatches));
        randreflectivity = step(obj.cNoiseSource,1,...
            [numel(obj.pEffRangeIdx),obj.pNumAzimuthPatches]);
        obj.pClutterReturn(obj.pEffRangeIdx,:) = ...
            obj.pClutterDeterministicReturn(obj.pEffRangeIdx,:).*...
            randreflectivity;
    end
    
end

methods(Static, Hidden)
    
    function initProps = initializeProperties( ...
                PropagationSpeed, ...
                OperatingFrequency, ...
                SampleRate , ...
                PRF, ...
                MaximumRange, ... 
                EarthModel, ...  
                PlatformHeight, ...
                AzimuthCoverage, ...
                PatchAzimuthWidth, ...
                BroadsideDepressionAngle,...
                PlatformDirection, ...
                PlatformSpeed, ...
                CoherenceTime, ...
                OutputFormat)
        %This function is for internal use only. It may be removed in the future.

        % calculate wavelength
        wavelength = PropagationSpeed/OperatingFrequency;

        % calculate number of samples per distinct pulse
        numSamplesPerDistinctPulse = round(SampleRate./...
            PRF);
        numDistinctPulses = numel(PRF);

        % Range resolution and possible range gates
        % The effective range is any range that is beyond the height
        rangeDelta = PropagationSpeed/(2*SampleRate);
        ranges = (0:rangeDelta:MaximumRange).';
        if strcmp(EarthModel,'Curved')
            effRangeIdx = find(...
                (ranges >= PlatformHeight) & ...
                (ranges <= horizonrange(PlatformHeight)));
        else
            effRangeIdx = find(ranges >= PlatformHeight);
        end
        numRanges = numel(ranges);

        % Azimuth angle for each clutter patch in the iso-range ring
        if ~rem(AzimuthCoverage,PatchAzimuthWidth)
            numAzimuthPatches = ...
                AzimuthCoverage/PatchAzimuthWidth;
            if ~rem(numAzimuthPatches,2)
                numAzimuthPatches = numAzimuthPatches+1;
            end
        else
            numAzimuthPatches = floor(...
                AzimuthCoverage/PatchAzimuthWidth);
            if ~rem(numAzimuthPatches,2)
                numAzimuthPatches = numAzimuthPatches+1;
            else
                % this is when the patch center resides outside the
                % coverage.
                numAzimuthPatches = numAzimuthPatches+2;
            end
        end
        halfnum = (numAzimuthPatches-1)/2;
        azAngles = (-halfnum:halfnum)*PatchAzimuthWidth;

        % Grazing angle, depression angle, and elevation angle for each
        % range gate
        grazingAngles = zeros(numRanges,1);
        grazingAngles(effRangeIdx,1) = grazingang(...
            PlatformHeight,ranges(effRangeIdx),...
            EarthModel);
        depressionAngles = zeros(numRanges,1);
        depressionAngles(effRangeIdx,1) = depressionang(...
            PlatformHeight,ranges(effRangeIdx),...
            EarthModel);
        elAngles = -(depressionAngles - ...
            BroadsideDepressionAngle).';

        % Propagation loss for each range gate
        propagationLoss = 2*fspl(ranges,wavelength);
        patchAreas = zeros(max(effRangeIdx),1);
        % Azimuth patch area for each range
        patchAreas(effRangeIdx,1) = ...
            phased.internal.AbstractClutterSimulator.patchArea(...
            ranges(effRangeIdx),...
            grazingAngles(effRangeIdx),rangeDelta,...
            PatchAzimuthWidth,...
            rad2deg(rangeDelta/PlatformHeight));

        % Relative radial speed for each clutter patch in respect to
        % the radar system and then convert to Doppler shift
        relativeRadialSpeed = zeros(numRanges,numAzimuthPatches);
        [velx,vely,velz] = sph2cart(...
            deg2rad(PlatformDirection(1)),...
            deg2rad(PlatformDirection(2)),...
            PlatformSpeed);
        platformVelocity = [velx;vely;velz];
        for m = 1:numAzimuthPatches
            PatchPos = [ranges.*cosd(depressionAngles).*cosd(azAngles(m)) ...
                ranges.*cosd(depressionAngles).*sind(azAngles(m)) ...
                ranges.*sind(depressionAngles)].';
            relativeRadialSpeed(:,m) = radialspeed(PatchPos,...
                zeros(3,numRanges),[0;0;0],platformVelocity);
        end
        dopplerShift = 2*relativeRadialSpeed/wavelength;
     
        % setup randomness update flag
        neverUpdateRandomSource = isinf(CoherenceTime);

        % setup output format flag
        outputBySamples = strcmp(OutputFormat,'Samples');        
    
        initProps.wavelength = wavelength;
        initProps.numSamplesPerDistinctPulse = numSamplesPerDistinctPulse;
        initProps.numDistinctPulses = numDistinctPulses;
        initProps.rangeDelta = rangeDelta;
        initProps.ranges = ranges;
        initProps.effRangeIdx = effRangeIdx;
        initProps.numRanges = numRanges;
        initProps.numAzimuthPatches = numAzimuthPatches;
        initProps.azAngles = azAngles;
        initProps.grazingAngles = grazingAngles;
        initProps.depressionAngles = depressionAngles;
        initProps.elAngles = elAngles;
        initProps.propagationLoss = propagationLoss;
        initProps.patchAreas = patchAreas;
        initProps.relativeRadialSpeed = relativeRadialSpeed;
        initProps.platformVelocity = platformVelocity;
        initProps.dopplerShift = dopplerShift;
        initProps.neverUpdateRandomSource = neverUpdateRandomSource;
        initProps.outputBySamples = outputBySamples;        
    
    end
        
    function cparea = patchArea(rng,grazingang,rngdelta,azdelta,eldelta)
    %patcharea Clutter patch area
    %   CPArea = patcharea(RNG,GrazANG,ThetaAZ,DeltaR,ThetaEl) calculates the
    %   clutter patch area based on the clutter patches' ranges RNG (in meters)
    %   and grazing angles GrazANG (in degrees). DeltaR is the width for each
    %   range gate (in meters), i.e., the range between adjacent range samples,
    %   ThetaAz is the azimuth span (in degrees) for each clutter patch,  and
    %   ThetaEl is the minimum elevation span for the patch.
    %   
    %   RNG and GrazANG can be vectors but their dimensions must match.
    %
    %   In classical literature, the area of clutter patch is computed using
    %   both grazing angle, azimuth span. The patch can either be pulse
    %   limited, where the pulse width is also used to calculate the area; or
    %   beam limited, where the elevation beamwidth is used to calculate the
    %   area. In our program, because we are simulating the return at each
    %   range gate, we use the width of range gate, DeltaR, as the basic
    %   calculation. Hence, let's define the corresponding ThetaEl as DeltaR/H
    %   where H is the platform height.
    %
    %   The area then is calculated as 
    %   
    %   A = RNG*ThetaAz*min(DeltaR,RNG*ThetaEl/sin(GrazANG)).
    %
    %   % Example:
    %   %   Calculate the area of the clutter patch if the platform is resided 
    %   %   300 meters above the ground and the clutter patch is 1000 meters 
    %   %   away from the sensor array. Assume each clutter patch occupies 1 
    %   %   degree in azimuth and the minimum elevation span for the patch is
    %   %   14 degrees.
    %
    %   cparea = patchArea(1000,grazingang(300,1000),50,1,14)
    
    AzLength = rng.*deg2rad(azdelta);
    ElLength = min(rngdelta, ...
        rng.*deg2rad(eldelta)./sind(grazingang));
    cparea = AzLength.*ElLength;
    
    end



end

end



% [EOF]
