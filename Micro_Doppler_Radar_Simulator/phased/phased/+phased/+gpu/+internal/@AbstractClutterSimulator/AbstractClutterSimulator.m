classdef (Hidden) AbstractClutterSimulator < matlab.system.internal.gpu.GPUSystem & ...
     matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%AbstractClutterSimulator   Define the abstract clutter simulator class.
%   This is the abstract class for clutter simulation running on the GPU
%
%Use of this object requires a Parallel Computing Toolbox license.

%   Copyright 2012-2017 The MathWorks, Inc.

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
    %   MATLAB GPU random number generator. When you set this property to
    %   'Property', a private random number generator is used with a seed
    %   specified by the value of the Seed property.
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
    %   PRFSelectionInputPort  Enable PRF selection input
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
    pWavelength;
    pRanges;
    pGrazingAngles;
    pDepressionAngles;
    pElAngles;
    pAzAngles;
    pTxGain;
    pPropagationLoss;
    pNRCS;
    pPatchAreas;
    pRangeDelta;
    pEffRangeIdx;
    pNumSamplesPerDistinctPulse;
    pRelativeRadialSpeed;
    pgDopplerShift;                   %on GPU
    pgClutterDeterministicReturn;     %on GPU
    cNoiseSource;
    pPlatformVelocity;
    pgAzAngles      %on GPU - copy of pAzAngles
    pgElAngles      %on GPU - copy of pElAngles
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
    pgClutterReturn;   %on GPU.
    pgClutterAtArrayOutputBuffer;  %on GPU.
    pLastRandUpdateTime;
    pRemainingSamplesFromLastPulse;
    pProcessedPulseCount;
end

methods (Access = protected, Abstract)
    nrcs = getNRCS(obj);    % calculate NRCS
end

methods
    function set.Sensor(obj,value)
        %Validate the Sensor
        validateattributes(value,{'phased.internal.AbstractElement',...
            'phased.internal.AbstractArray','em.Antenna',...
            'phased.internal.AbstractSubarray'},{'scalar'},'','Sensor');

        %Validate the elements
        elem = obj.getElementHandle(value);
        if ~iscell(elem)
          elem = {elem};
        end
        for ii=1:numel(elem) %check all elements
          validateattributes(elem{ii},...
              {'phased.internal.AbstractAntennaElement','em.Antenna'},{'scalar'},...
              '','Sensor array element');
        end
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
        obj.Sensor = phased.ULA;
        setProperties(obj,nargin,varargin{:});
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
        if any(rem(obj.SampleRate, obj.PRF))
            error(message('phased:Waveform:NeedRatioInteger',...
                'SampleRate', 'PRF'));
        end
        
        if obj.MaximumRange <= obj.PlatformHeight
            error(message('phased:phased:expectedGreaterThan',...
                'MaximumRange',sprintf('%5.2f',obj.PlatformHeight)));
        end
        
        if obj.PatchAzimuthWidth > obj.AzimuthCoverage
            error(message('phased:phased:expectedLessThanOrEqualTo',...
                'PatchAzimuthWidth','AzimuthCoverage'));
        end
        
        if rem(obj.AzimuthCoverage, obj.PatchAzimuthWidth) && ...
                rem(floor(obj.AzimuthCoverage/obj.PatchAzimuthWidth),2)
            warning(message('phased:clutter:patchCenterOutOfBound',...
                'AzimuthCoverage','PatchAzimuthWidth'));
        end
        
        %Validate contained objects here since they never have step called.
        hsens = obj.Sensor;

        hsubarr = obj.getSubarrayHandle(hsens);
        if ~isempty(hsubarr)
          validateProperties(hsubarr);
        end

        harr = obj.getArrayHandle(hsens);
        if ~isempty(harr)
          validateProperties(harr);
        end

        helem = obj.getElementHandle(hsens);
        if ~isempty(helem)
          if ~iscell(helem)
            helem = {helem};
          end
          for ii=1:numel(helem)
              if ~isElementFromAntenna(helem{ii})
                  validateProperties(helem{ii});
              end
          end
        end
    end

    function validateInputsImpl(obj,varargin)
        if obj.TransmitSignalInputPort
            x = varargin{1};
            if ~obj.isdouble(x)
                matlab.system.internal.error(...
                    'MATLAB:system:invalidInputDataType','X','double');
            end
            if ~iscolumn(x) || isempty(x)
                matlab.system.internal.error(...
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
                if ~ismatrix(stang) || isempty(stang)
                    matlab.system.internal.error(...
                        'MATLAB:system:inputMustBeMatrix','Steer');
                end
                if sz_stang(1) > 2
                    error(message('phased:system:array:NeedTwoRows','Steer'));
                end
                if sz_stang(2) > 1
                    error(message('phased:system:array:NeedOneColumn','Steer'));
                end
                if ~isreal(stang)
                    error(message('phased:system:array:InvalidAngle','Steer'));
                end
                if ~obj.isdouble(stang)
                    matlab.system.internal.error(...
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
        if obj.PRFSelectionInputPort
            if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    ~strncmp(obj.Sensor.SubarraySteering,'None',1)
                if obj.TransmitSignalInputPort
                    prfidx = varargin{3};
                else
                    prfidx = varargin{2};
                end
            else
                if obj.TransmitSignalInputPort
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

    function setupGPUImpl(obj,varargin)
        % initialize random source
        obj.cNoiseSource = phased.gpu.internal.ComplexWGNSource(...
            'SeedSource',obj.SeedSource);
        if strncmpi(obj.SeedSource,'p',1)
            obj.cNoiseSource.Seed = obj.Seed;
        end

        % calculate wavelength
        obj.pWavelength = obj.PropagationSpeed/obj.OperatingFrequency;

        % calculate number of samples per distinct pulse
        obj.pNumSamplesPerDistinctPulse = round(obj.SampleRate./...
            obj.PRF);
        obj.pNumDistinctPulses = numel(obj.PRF);

        obj.pUseArray = getUseArray(obj);
        obj.pDOF = getClutterDOF(obj);

        obj.pNeedSteeringAngle = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
            (strncmp(obj.Sensor.SubarraySteering,'Phase',1) || ...
            strncmp(obj.Sensor.SubarraySteering,'Time',1));

        obj.pNeedCustomWeights = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
            strncmp(obj.Sensor.SubarraySteering,'Cusotm',1);
        
        % Range resolution and possible range gates
        % The effective range is any range that is beyond the height
        obj.pRangeDelta = obj.PropagationSpeed/(2*obj.SampleRate);
        obj.pRanges = (0:obj.pRangeDelta:obj.MaximumRange).';
        if strcmp(obj.EarthModel,'Curved')
            obj.pEffRangeIdx = find(...
                (obj.pRanges >= obj.PlatformHeight) & ...
                (obj.pRanges <= horizonrange(obj.PlatformHeight)));
        else
            obj.pEffRangeIdx = find(obj.pRanges >= obj.PlatformHeight);
        end
        obj.pNumRanges = numel(obj.pRanges);

        % Azimuth angle for each clutter patch in the iso-range ring
        if ~rem(obj.AzimuthCoverage,obj.PatchAzimuthWidth)
            obj.pNumAzimuthPatches = ...
                obj.AzimuthCoverage/obj.PatchAzimuthWidth;
            if ~rem(obj.pNumAzimuthPatches,2)
                obj.pNumAzimuthPatches = obj.pNumAzimuthPatches+1;
            end
        else
            obj.pNumAzimuthPatches = floor(...
                obj.AzimuthCoverage/obj.PatchAzimuthWidth);
            if ~rem(obj.pNumAzimuthPatches,2)
                obj.pNumAzimuthPatches = obj.pNumAzimuthPatches+1;
            else
                % this is when the patch center resides outside the
                % coverage.
                obj.pNumAzimuthPatches = obj.pNumAzimuthPatches+2;
            end
        end
        halfnum = (obj.pNumAzimuthPatches-1)/2;
        obj.pAzAngles = (-halfnum:halfnum)*obj.PatchAzimuthWidth;

        % Grazing angle, depression angle, and elevation angle for each
        % range gate
        obj.pGrazingAngles = zeros(obj.pNumRanges,1);
        obj.pGrazingAngles(obj.pEffRangeIdx,1) = grazingang(...
            obj.PlatformHeight,obj.pRanges(obj.pEffRangeIdx),...
            obj.EarthModel);
        obj.pDepressionAngles = zeros(obj.pNumRanges,1);
        obj.pDepressionAngles(obj.pEffRangeIdx,1) = depressionang(...
            obj.PlatformHeight,obj.pRanges(obj.pEffRangeIdx),...
            obj.EarthModel);
        obj.pElAngles = -(obj.pDepressionAngles - ...
            obj.BroadsideDepressionAngle).';

        %GPU copies of the angles - convert to radians first
        obj.pgAzAngles = gpuArray(pi* obj.pAzAngles/180); 
        obj.pgElAngles = gpuArray(pi* obj.pElAngles/180);

        % Transmit gain for each clutter patch
        obj.pTxGain = zeros(obj.pNumRanges,obj.pNumAzimuthPatches);
        if obj.pUseArray
            if ~(obj.pNeedSteeringAngle || obj.pNeedCustomWeights)

                if isa(obj.Sensor , 'phased.internal.AbstractSubarray')
                    tempresp = gpuArray.zeros(obj.pNumRanges, ...
                      obj.pNumAzimuthPatches);
                    for ii=1:getNumSubarrays(obj.Sensor)
                        resp = runSubArray(obj, obj.Sensor, ii, ...
                        obj.pgAzAngles, obj.pgElAngles, []);
                        tempresp = tempresp + squeeze((sum(resp,2)));
                    end
                    noiseMtxGain = getNoiseGainMatrix(obj.Sensor, ...
                      obj.OperatingFrequency, obj.PropagationSpeed);
                    gainFactor = real(sum(sum(noiseMtxGain)));

                else %Regular array
                    pos = getElementPosition(obj.Sensor);
                    tempresp  = runArray(obj, obj.OperatingFrequency, ...
                        obj.PropagationSpeed, pos, obj.pgAzAngles, ...
                        obj.pgElAngles, []);
                    gainFactor = sum(obj.pDOF);
                end
                txgfh = @(x) 10*log10(abs(x)*abs(x));
                obj.pTxGain = gather(arrayfun(txgfh, tempresp)) ...
                  - 10*log10(gainFactor);
            end
        else
            hac = cloneSensor(obj.Sensor);
            for m = 1:obj.pNumAzimuthPatches
                ang = [obj.pAzAngles(m)*ones(1,obj.pNumRanges);...
                    obj.pElAngles];
                antOut = step(hac, obj.OperatingFrequency,ang);
                if isPolarizationEnabled(hac)
                    resp = hypot(antOut.H, antOut.V);
                else
                    resp = antOut;
                end
                obj.pTxGain(:,m) = mag2db(abs(resp));
            end
        end 

        % Propagation loss for each range gate
        obj.pPropagationLoss = 2*fspl(obj.pRanges,obj.pWavelength);
        obj.pNRCS = getNRCS(obj);

        % Azimuth patch area for each range
        obj.pPatchAreas(obj.pEffRangeIdx,1) = ...
            localpatcharea(obj.pRanges(obj.pEffRangeIdx),...
            obj.pGrazingAngles(obj.pEffRangeIdx),obj.pRangeDelta,...
            obj.PatchAzimuthWidth,...
            rad2deg(obj.pRangeDelta/obj.PlatformHeight));

        % Combine to get the total reflect gain from each range gate,
        % i.e., 4*pi*A/lambda^2
        rcseffect(obj.pEffRangeIdx,1) = db2mag(...
            aperture2gain(obj.pNRCS(obj.pEffRangeIdx).*...
            obj.pPatchAreas(obj.pEffRangeIdx),...
            obj.pWavelength));

        % Relative radial speed for each clutter patch in respect to
        % the radar system and then convert to Doppler shift
        obj.pRelativeRadialSpeed = zeros(obj.pNumRanges,obj.pNumAzimuthPatches);
        [velx,vely,velz] = sph2cart(...
            deg2rad(obj.PlatformDirection(1)),...
            deg2rad(obj.PlatformDirection(2)),...
            obj.PlatformSpeed);
        obj.pPlatformVelocity = [velx;vely;velz];
        for m = 1:obj.pNumAzimuthPatches
            PatchPos = [obj.pRanges.*cosd(obj.pDepressionAngles).*cosd(obj.pAzAngles(m)) ...
                obj.pRanges.*cosd(obj.pDepressionAngles).*sind(obj.pAzAngles(m)) ...
                obj.pRanges.*sind(obj.pDepressionAngles)].';
            obj.pRelativeRadialSpeed(:,m) = radialspeed(PatchPos,...
                zeros(3,obj.pNumRanges),[0;0;0],obj.pPlatformVelocity);
        end
        obj.pgDopplerShift = ...
            gpuArray(2*obj.pRelativeRadialSpeed/obj.pWavelength);

        % clutter return for each patch. Compute on CPU then ship to GPU.
        tmpgClutterDeterministicReturn = ...
            zeros(obj.pNumRanges,obj.pNumAzimuthPatches);
        for m = 1:obj.pNumAzimuthPatches
            tmpgClutterDeterministicReturn(obj.pEffRangeIdx,m) = ...
                db2mag(obj.pTxGain(obj.pEffRangeIdx,m)).*...
                rcseffect(obj.pEffRangeIdx)./...
                db2mag(obj.pPropagationLoss(obj.pEffRangeIdx));
        end
        if ~obj.TransmitSignalInputPort
            tmpgClutterDeterministicReturn = ...
                tmpgClutterDeterministicReturn * sqrt(obj.TransmitERP);
        end
        %put on the GPU
        obj.pgClutterDeterministicReturn = ...
                  gpuArray(tmpgClutterDeterministicReturn);


        % initialize clutter return of each range gate at the
        % receiving sensor/array
        obj.pBufferLength = obj.pNumRanges;
        if obj.TransmitSignalInputPort
            x = varargin{1};
            sz = size(x);
            obj.pBufferLength = obj.pBufferLength+sz(1)-1;
        end
        
        % setup randomness update flag
        obj.pNeverUpdateRandomSource = isinf(obj.CoherenceTime);

        % setup output format flag
        obj.pOutputBySamples = strcmp(obj.OutputFormat,'Samples');        
    end

    function flag = isInputSizeLockedImpl(~,~)
        %Override the base class implementation
        flag = false;
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
        reset(obj.cNoiseSource);
        
        obj.pgClutterAtArrayOutputBuffer = complex(...
              gpuArray.zeros(obj.pBufferLength,obj.pDOF));
            
        % time offset for the starting time at each simulation step
        obj.pTimeOffset = 0;

        % initialize the time stamp for last randomness update
        obj.pLastRandUpdateTime = 0;
        
        % initialize pulse count
        obj.pRemainingSamplesFromLastPulse = 0;
        obj.pProcessedPulseCount = 0;

    end

    function releaseImpl(obj)
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
        s.base = saveObjectImpl@matlab.system.internal.gpu.GPUSystem(obj);
        s.Sensor = saveobj(obj.Sensor);
        if isLocked(obj)
            s.pUseArray = obj.pUseArray;
            s.cNoiseSource = saveobj(obj.cNoiseSource);
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
            s.pgDopplerShift = obj.pgDopplerShift;
            s.pgClutterDeterministicReturn = ...
                obj.pgClutterDeterministicReturn;
            s.pPlatformVelocity = obj.pPlatformVelocity;
            s.pNeverUpdateRandomSource = obj.pNeverUpdateRandomSource;
            s.pOutputBySamples = obj.pOutputBySamples;
            s.pNumRanges = obj.pNumRanges;
            s.pDOF = obj.pDOF;
            s.pBufferLength = obj.pBufferLength;
            s.pNumDistinctPulses = obj.pNumDistinctPulses;
            s.pNumAzimuthPatches = obj.pNumAzimuthPatches; 
            s.pTimeOffset = obj.pTimeOffset;
            s.pgClutterReturn = obj.pgClutterReturn;
            s.pgClutterAtArrayOutputBuffer = obj.pgClutterAtArrayOutputBuffer;
            s.pLastRandUpdateTime = obj.pLastRandUpdateTime;
            s.pRemainingSamplesFromLastPulse = obj.pRemainingSamplesFromLastPulse;
            s.pProcessedPulseCount = obj.pProcessedPulseCount;
            s.pNeedSteeringAngle = obj.pNeedSteeringAngle;
            s.pNeedCustomWeights = obj.pNeedCustomWeights;
            s.pgAzAngles = obj.pgAzAngles;
            s.pgElAngles = obj.pgElAngles;
        end
    end

    function s = loadSubObjects(obj,s, wasLocked) 
        if isa(s.Sensor,'em.Antenna')
            obj.Sensor = em.EmStructures.loadobj(s.Sensor);
        else
            obj.Sensor = eval(...
                sprintf('%s.loadobj(s.Sensor)',s.Sensor.ClassNameForLoadTimeEval));
        end
        s = rmfield(s,'Sensor');
        if wasLocked
            obj.cNoiseSource = ...
                phased.gpu.internal.ComplexWGNSource.loadobj(...
                s.cNoiseSource);
            s = rmfield(s,'cNoiseSource');
        end
        
    end
    
    function loadObjectImpl(obj,s,wasLocked)
        %Load base class props and then remove the base field
        loadObjectImpl@matlab.system.internal.gpu.GPUSystem(obj, s.base, wasLocked); 
        s = rmfield(s, 'base');
        
        s = loadSubObjects(obj,s, wasLocked);
        
        %set props for this class
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
        obj.pgClutterReturn = complex(...
            gpuArray.zeros(obj.pNumRanges,obj.pNumAzimuthPatches));

        randreflectivity = step(obj.cNoiseSource,1,...
            [numel(obj.pEffRangeIdx),obj.pNumAzimuthPatches]);
        
        obj.pgClutterReturn(obj.pEffRangeIdx,:) = ...
            obj.pgClutterDeterministicReturn(obj.pEffRangeIdx,:).*...
            randreflectivity;
    end
    
    function arrayEffect = runSubArray(obj,arrh, subArrIdx, gAz, gEl, stang)
        %Return the response of the subarray at index subArrIdx 
        %arrayEffect is el x subarrayElemPos x az
        
        freq = obj.OperatingFrequency;
        propSpeed = obj.PropagationSpeed;
        
        %Transform the subarray element positions, az angles, and elev
        %angles to the coordinate system and orientation of the subarray.
        [pos, az, el] = getSubarrayPosAzEl(arrh, subArrIdx, gAz, gEl);

        %Using the original global angles and subarray feed
        %position, calculate path delays. 
        subArrPos = getSubarrayPosition(arrh);
        gAz = reshape(gAz, 1, 1, []);
        gEl = reshape(gEl, [], 1, 1);

        %Calculate the steering vector 
        arrsv = obj.computeSteeringVector(subArrPos(:,subArrIdx), freq, ...
          propSpeed, gAz, gEl);
       

        %compute weights for this subarray if needed
        if obj.pNeedSteeringAngle || obj.pNeedCustomWeights
            subArrSV = getSubArraySteeringVec(arrh, subArrIdx, freq, ...
              propSpeed,  stang);
        else
            subArrSV = [];
        end
       
        if isa(arrh, 'phased.PartitionedArray')
            %Partitioned Array gets weights scaled by subarray feed
            %steering vector
            fh = @(x,y) x*conj(y);
            subArrSV = arrayfun(fh, reshape(subArrSV, 1, [], 1), arrsv);
        else
            subArrSV = reshape(subArrSV, 1, [], 1);
        end
        
        %compute the array response based on the angles relative to the
        %subarray feed position and orientation.
        tempresp = runArray(obj, freq, propSpeed, pos, az, el, subArrSV);
 
        %Multiply the subarray feed delays by the subarray element
        %response.
        arrayEffect = arrayfun(@times, arrsv, tempresp);
    end
        
    function resp = runArray(obj, freq, propSpeed, pos, az, el, weights)
        %Computes the response for a single array or subarray at specified
        %positions, azimuth, elevation, operating frequency. A weight vector
        %is applied if the array is steered.

        arrh = obj.getArrayHandle(obj.Sensor);
        if isvector(el)
            idx = phased.gpu.internal.decomposeCube(numel(el), ...
                numel(az), size(pos, 2) );
            resp = gpuArray.zeros(numel(el), 1, numel(az));
        else
            idx = phased.gpu.internal.decomposeCube(size(el,1), ...
                size(az,3), size(pos, 2) );
            resp = gpuArray.zeros(size(el,1), 1, size(az,3));
        end
        for ii=1:numel(idx)
            %iterate over the positions
            arrsv = obj.computeSteeringVector(pos(:, idx{ii}), freq, ...
              propSpeed, az, el);
          
            arrsv = scaleByElemResponse(arrh, az, el, freq, arrsv, idx{ii});
            
            if obj.pNeedSteeringAngle || obj.pNeedCustomWeights
                %If steering is used, get the steering vector for the subarray
                %elements and do a weighted sum of the subarray response.
                weightedsv = bsxfun(@times, arrsv, weights(:,idx{ii}, :) );
                arrResp = sum(weightedsv, 2);
            else
                
                % For non-steered arrays, weights are equal. 
                %So just sum across elements
                arrResp = sum(arrsv, 2);
            end
            resp = resp + arrResp;
        end
    end

end

methods (Access = private)
    function arrh = getSubarrayHandle(obj, sensor) %#ok<INUSL>
        %Get the handle to array object.
        if isa(sensor, 'phased.internal.AbstractSubarray')
            arrh = sensor;
        else
            arrh = [];
        end
    end
    function arrh = getArrayHandle(obj, sensor) %#ok<INUSL>
        %Get the handle to array object.
        if isa(sensor, 'phased.internal.AbstractArray')
            arrh = sensor;
        elseif isa(sensor, 'phased.ReplicatedSubarray')
            arrh = sensor.Subarray;
        elseif isa(sensor, 'phased.PartitionedArray')
            arrh = sensor.Array;
        else
            arrh = [];
        end
    end
         
    function elemh = getElementHandle(obj, sensor)
        %Get the handle to the element object
        arrh= getArrayHandle(obj, sensor);
        if isempty(arrh)
          elemh = sensor;   %single element
        else
          elemh = getElementHandle(arrh);  %call array method
        end
    end

end

methods (Access = protected)
% Helper methods to assign properties or variables. 
% Used by setupGPUImpl and propagators

  function ua = getUseArray(obj) 
    ua = isa(obj.Sensor,'phased.internal.AbstractArray') || ...
              isa(obj.Sensor,'phased.internal.AbstractSubarray');
  end

  function cdof = getClutterDOF(obj)
    if getUseArray(obj) 
        cdof = getDOF(obj.Sensor);
    else
        cdof = 1;
    end
  end
end
methods (Static, Hidden)
    function sv = computeSteeringVector(pos,freq, propSpeed, az, el)
        %Returns the steering vector for a given set of antenna positions,
        %a frequency, a propagation speed, and azimuth and elevation
        %angles. The computation takes places on a GPU
        
        %The output sv is a cube of dimensions 
        %   numel(el) x numel(pos) x numel(az)
        
        gpos = gpuArray(pos);
        gposx = gpos(1,:);
        gposy = gpos(2,:);
        gposz = gpos(3,:);
        
        %If the input az, el are vectors then reshape for scalar expansion.
        if isvector(az)
            az = reshape( az, 1, 1, []);
            el = reshape( el, [], 1, 1);
        end
        
        sv = arrayfun(@privSteeringVec, gposx, gposy, gposz, ...
            freq, propSpeed, ...
            az, el);
    end
end

end


function cparea = localpatcharea(rng,grazingang,rngdelta,azdelta,eldelta)
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
%   cparea = localpatcharea(1000,grazingang(300,1000),50,1,14)

AzLength = rng.*deg2rad(azdelta);
ElLength = min(rngdelta, ...
    rng.*deg2rad(eldelta)./sind(grazingang));
cparea = AzLength.*ElLength;

end

% [EOF]
