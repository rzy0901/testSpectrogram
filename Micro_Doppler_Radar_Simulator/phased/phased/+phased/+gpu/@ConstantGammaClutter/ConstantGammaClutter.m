classdef (Sealed,StrictDefaults) ConstantGammaClutter < phased.gpu.internal.AbstractClutterSimulator & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%ConstantGammaClutter   Constant gamma clutter simulation on GPU
%   H = phased.gpu.ConstantGammaClutter creates a constant gamma clutter
%   simulation System object, H. This object simulates the clutter return
%   of a monostatic radar system using constant gamma model.
%
%   Use of this object requires a Parallel Computing Toolbox license.
%
%   H = phased.gpu.ConstantGammaClutter(Name,Value) creates a constant
%   gamma clutter object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H) returns an NxM matrix, Y, that contains the collected
%   clutter return. M is the number of radiating/collecting subarrays in
%   the radar system if Sensor contains subarrays, or the number of
%   radiating/collecting sensors otherwise. When you set the OutputFormat
%   property to 'Samples', N is specified in the NumSamples property. When
%   you set the OutputFormat property to 'Pulses', N is the total number of
%   samples in the next L pulses where L is specified in the NumPulses
%   property. This syntax is available when you set the
%   TransmitSignalInputPort property to false.
%
%   The clutter simulation makes the following assumptions:
%       * The radar system is monostatic
%       * The propagation is in free space
%       * The terrain is homogeneous
%       * The clutter patch is stationary during coherence time
%       * The signal is narrowband
%       * The radar system maintains a constant height during simulation
%       * The radar system maintains a constant speed during simulation
%
%   Y = step(H,X) specifies the transmit signal in X as a column vector. In
%   general, transmit signal refers to the signal at the output of the
%   transmitter when the transmitter is turned on during each pulse. This
%   syntax is available when you set the TransmitSignalInputPort property
%   to true.
%
%   Y = step(H,STEER) uses STEER as the subarray steering angle (in
%   degrees). STEER can be a scalar or a length-2 column vector. If STEER
%   is a vector, it is in the form of [AzimuthAngle; ElevationAngle]. If
%   STEER is a scalar, it represents the azimuth angle and the elevation
%   angle is assumed to be 0. This syntax is only applicable when you use
%   subarrays in the Sensor property and set the SubarraySteering property
%   in the Sensor to either 'Phase' or 'Time'.
%
%   Y = step(H,WS) uses WS as the weights applied to each element in the
%   subarray. WS can be either a matrix or a cell array. This syntax is
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
%   can have different number of elements, WS is an NSExN matrix, where NSE
%   indicates the number of elements in the largest subarray and N is the
%   number of subarrays.  The first K entries in each column, where K is
%   the number of elements in the corresponding subarray, specifies the
%   weights for the elements in the corresponding subarray.
%
%   Y = step(H,PRFIDX) specifies the index of pulse repetition frequency
%   (PRF), PRFIDX, as a positive integer. The index is used to identify the
%   entries specified in the PRF property. This syntax applies when you set
%   the PRFSelectionInputPort property to true. Use this syntax for the
%   cases where the transmit pulse needs to be dynamically selected. Under
%   such situations, the PRF property includes a list of predetermined
%   choices of PRFs. During the simulation, based on PRF index input, one
%   of the PRFs is selected as the PRF for the next transmission.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,STEER,PRFIDX)
%
%   or
%
%   Y = step(H,X,WS,PRFIDX)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H) and y = H() are
%   equivalent.
%
%   ConstantGammaClutter methods:
%       step     - Simulate clutter using constant gamma model (see above)
%       release  - Allow property value and input characteristics changes
%       clone    - Create a constant gamma clutter object with same
%                  property values
%       isLocked - Locked status (logical)
%       <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset random number and time count for clutter
%                  simulation
%
%   ConstantGammaClutter properties:
%       Sensor                   - Handle of the sensor
%       Gamma                    - Terrain gamma value
%       EarthModel               - Earth model
%     	MaximumRange             - Maximum range
%     	AzimuthCoverage          - Azimuth coverage
%      	PatchAzimuthWidth        - Clutter patch azimuth span
%       CoherenceTime            - Clutter coherence time
%       PropagationSpeed         - Propagation speed
%       SampleRate               - Sample rate
%       PRF                      - Pulse repetition frequency
%       PRFSelectionInputPort    - Enable PRF selection input
%    	OutputFormat             - Output signal format
%      	NumPulses                - Number of pulses in output
%     	NumSamples               - Number of samples in output
%       OperatingFrequency       - Operating frequency
%    	TransmitSignalInputPort  - Enable transmit signal input
%     	TransmitERP              - Effective transmitted power
%       PlatformHeight           - Radar height
%     	PlatformSpeed            - Radar speed
%    	PlatformDirection        - Radar motion direction
%    	BroadsideDepressionAngle - Broadside depression angle
%      	SeedSource               - Source of seed for random number
%                                  generator
%      	Seed                     - Seed for random number generator
%
%   % Examples:
%
%   % Example 1:
%   %   Assume a radar system with a 4 element ULA. The sample rate is 1
%   %   MHz. The PRF is 10 kHz. The propagation speed is 300,000 km/s. The
%   %   operating frequency is 300 MHz. The radar platform is flying 1 km
%   %   above the ground with a path parallel to the ground along the array
%   %   axis. The platform's speed is 2000 m/s. The mainlobe has a
%   %   depression angle of 30 degrees. The terrain has a gamma value of 0
%   %   dB.
%   %
%   %   Simulate the clutter return for 10 pulses, assuming the earth is
%   %   flat. The maximum clutter range of interest is 5 km and the maximum
%   %   azimuth coverage is +/- 60 degrees. The effective transmitted power
%   %   is 5 kw.
%
%   c = 3e8; prf = 10e3; fs = 1e6; fc = 3e8; height = 1000; speed = 2000;
%   tergamma = 0; Rmax = 5000; Npulse = 10; Azcov = 120; lambda = c/fc;
%   tpower = 5000; Nsamp = fs/prf; Nele = 4; depang = 30;
%
%   % create clutter object
%   ha = phased.ULA('NumElements',Nele,'ElementSpacing',lambda/2);
%   hclutter = phased.gpu.ConstantGammaClutter('Sensor',ha,...
%           'PropagationSpeed',c,'OperatingFrequency',fc,'PRF',prf,...
%           'SampleRate',fs,'Gamma',tergamma,'EarthModel','Flat',...
%           'TransmitERP',tpower,'PlatformHeight',height,...
%           'PlatformSpeed',speed,'PlatformDirection',[90;0],...
%           'BroadsideDepressionAngle',depang,'MaximumRange',Rmax,...
%           'AzimuthCoverage',Azcov);
%
%   % simulate clutter
%   csig = zeros(Nsamp,Nele,Npulse);
%   for m = 1:Npulse
%       csig(:,:,m) = step(hclutter);
%   end
%
%   % plot angle-Doppler response of the clutter at the 20th range bin
%   hresp = phased.AngleDopplerResponse('SensorArray',ha,...
%               'OperatingFrequency',fc,'PropagationSpeed',c,'PRF',prf);
%   plotResponse(hresp,shiftdim(csig(20,:,:)),'NormalizeDoppler',true);
%
%   % Example 2:
%   %   Use the same radar system as the previous example but use a
%   %   rectangular waveform as an additional input. The pulse width of the
%   %   waveform is 2 microseconds.
%
%   c = 3e8; prf = 10e3; fs = 1e6; fc = 3e8; height = 1000; speed = 2000;
%   tergamma = 0; Rmax = 5000; Npulse = 10; Azcov = 120; lambda = c/fc;
%   tpower = 5000; Nsamp = fs/prf; Nele = 4; depang = 30; pw = 2e-6;
%
%   % create clutter object
%   ha = phased.ULA('NumElements',Nele,'ElementSpacing',lambda/2);
%   hclutter = phased.gpu.ConstantGammaClutter('Sensor',ha,...
%           'PropagationSpeed',c,'OperatingFrequency',fc,'PRF',prf,...
%           'SampleRate',fs,'Gamma',tergamma,'EarthModel','Flat',...
%           'PlatformHeight',height,'TransmitSignalInputPort',true,...
%           'PlatformSpeed',speed,'PlatformDirection',[90;0],...
%           'BroadsideDepressionAngle',depang,'MaximumRange',Rmax,...
%           'AzimuthCoverage',Azcov);
%
%   % simulate clutter
%   tx_sig = tpower*ones(floor(pw*fs),1);
%   csig = zeros(Nsamp,Nele,Npulse);
%   for m = 1:Npulse
%       csig(:,:,m) = step(hclutter,tx_sig);
%   end
%
%   % plot angle-Doppler response of the clutter at the 20th range bin
%   hresp = phased.AngleDopplerResponse('SensorArray',ha,...
%               'OperatingFrequency',fc,'PropagationSpeed',c,'PRF',prf);
%   plotResponse(hresp,shiftdim(csig(20,:,:)),'NormalizeDoppler',true);
%
%   See also phased, phased.ConstantGammaClutter, surfacegamma.


%   Copyright 2012-2017 The MathWorks, Inc.

%   Reference
%   [1] Maurice Long, Radar Reflectivity of Land and Sea, 3rd Ed. Artech
%       House, 2001
%   [2] J. R. Guerci, Space Time Adaptive Processing for Radar, Artech
%       House, 2003
%   [3] James Ward, Space-time Adaptive Processing for Airborne Radar,
%       Lincoln Lab Technical Report, 1994
%   [4] David Greig, Time Domain Clutter Model for Airborne Pulse Doppler
%       Radar, IEE Colloquium on Radar System Modeling, 1998
    
    
    properties (Nontunable)
        %Gamma  Terrain gamma value (dB)
        %   Specify the gamma value (in dB) used in the constant gamma
        %   clutter model as a scalar. The default value of this property
        %   is 0. The gamma value depends on both terrain type and the
        %   operating frequency.
        Gamma = 0;
    end
    
    methods
        function set.Gamma(obj,value)
            validateattributes(value,{'double'},...
                {'real','scalar','finite','nonempty'},'','Gamma');
            obj.Gamma = value;
        end
    end
    
    methods
        function obj = ConstantGammaClutter(varargin)
            %ConstantGammaClutter   Construct the ConstantGammaClutter
            %class.
            obj@phased.gpu.internal.AbstractClutterSimulator(varargin{:});
            
        end
    end
    
    methods (Access = protected)
        function nrcs = getNRCS(obj)
            nrcs = db2pow(obj.Gamma).*sind(obj.pGrazingAngles) + ...
                10.*exp(-0.2.*(90-obj.pGrazingAngles));  % [Long03]
        end
        
        function setupGPUImpl(obj,varargin)
            setupGPUImpl@phased.gpu.internal.AbstractClutterSimulator(obj,varargin{:});
            if getNumInputs(obj)>0
                N = obj.getPulseOutputSize(...
                    obj.SampleRate, obj.PRF,obj.OutputFormat, ...
                    obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
                validateSampleRate(obj,N);
            end
        end
        
        function y = stepGPUImpl(obj,xArg, stangArg, prfidxArg)
            % initialize random clutter
            if ~obj.pTimeOffset
                calcRandomClutter(obj);
            end
            
            if obj.pNeedSteeringAngle || obj.pNeedCustomWeights
                if ~obj.TransmitSignalInputPort 
                    x = [];
                    stang_in = xArg;
                    if obj.PRFSelectionInputPort
                        prfidx = stangArg;
                    else
                        prfidx = [];
                    end
                else
                    x = xArg;
                    stang_in = stangArg;
                    if obj.PRFSelectionInputPort
                        prfidx = prfidxArg;
                    else
                        prfidx = [];
                    end
                end
            else
                if ~obj.TransmitSignalInputPort
                    x = [];
                    if obj.PRFSelectionInputPort
                        prfidx = xArg;
                    else
                        prfidx = [];
                    end
                else
                    x = xArg;
                    if obj.PRFSelectionInputPort
                        prfidx = stangArg;
                    else
                        prfidx = [];
                    end
                end
                stang = [];
            end
            
            if obj.pNeedSteeringAngle
                if isscalar(stang_in)
                    stang = [stang_in; 0];
                else
                    stang = stang_in;
                end

                %Validate steering angle
                if ~isempty(stang)
                    azstang = stang(1,:);
                    elstang = stang(2,:);
                    if any(azstang>180) || any(elstang>90)
                        error(message('phased:step:AzElNotLessEqual'));
                    end
                    if any(azstang<-180) || any(elstang<-90)
                        error(message('phased:step:AzElNotGreaterEqual'));
                    end
                end
            elseif obj.pNeedCustomWeights
                stang = stang_in;
            end
            
            stang = gather(stang); %make sure stang is on the CPU.

            prfidx = gather(prfidx); %make sure stang is on the CPU.

            if obj.PRFSelectionInputPort
                N = calcTotalOutputSamples(obj,prfidx);
            else
                N = calcTotalOutputSamples(obj);
            end
            
            y = complex(gpuArray.zeros(N,obj.pDOF));
            outputIdxOffset = 0;
            
            while N > 0
                if obj.pRemainingSamplesFromLastPulse > 0
                    % there are remaining samples from last pulse
                    if N > obj.pRemainingSamplesFromLastPulse
                        num_output_samples = obj.pRemainingSamplesFromLastPulse;
                    else
                        num_output_samples = N;
                    end
                    [y,outputIdxOffset] = outputSamplesFromBuffer(...
                        obj,y,outputIdxOffset,num_output_samples);
                    obj.pRemainingSamplesFromLastPulse = ...
                        obj.pRemainingSamplesFromLastPulse-num_output_samples;
                else
                    % process next pulse
                    gcollectedclutter = getClutterFromCurrentTransmit(obj, x, stang);
                    % cache the clutter from the current simulation step
                    % into buffer, which holds return from each clutter
                    % patch up to current simulation step and not in the
                    % output yet.
                    obj.pgClutterAtArrayOutputBuffer = ...
                        obj.pgClutterAtArrayOutputBuffer + gcollectedclutter;
                    
                    if obj.PRFSelectionInputPort
                        num_samples_current_pulse = ...
                            getNumSamplesForCurrentPulse(obj,prfidx);
                    else
                        num_samples_current_pulse = ...
                            getNumSamplesForCurrentPulse(obj);
                    end
                    
                    obj.pProcessedPulseCount = obj.pProcessedPulseCount + 1;
                    
                    if N < num_samples_current_pulse
                        num_output_samples = N;
                        [y,outputIdxOffset] = outputSamplesFromBuffer(...
                            obj,y,outputIdxOffset,num_output_samples);
                    else
                        num_output_samples = num_samples_current_pulse;
                        [y,outputIdxOffset] = outputSamplesFromBuffer(...
                            obj,y,outputIdxOffset,num_output_samples);
                    end
                    obj.pRemainingSamplesFromLastPulse = ...
                        num_samples_current_pulse - num_output_samples;
                    
                    obj.pTimeOffset = obj.pTimeOffset + ...
                        num_samples_current_pulse/obj.SampleRate;
                    
                    if ~obj.pNeverUpdateRandomSource
                        updateClutterRandomReturn(obj);
                    end
                end
                N = N-num_output_samples;
            end
        end
    end
    
    methods (Access = private)
        function updateClutterAtArrayBuffer(obj,N)
            obj.pgClutterAtArrayOutputBuffer(1:obj.pBufferLength-N,:) = ...
                obj.pgClutterAtArrayOutputBuffer(N+1:obj.pBufferLength,:);
            obj.pgClutterAtArrayOutputBuffer(obj.pBufferLength-N+1:end,:) = 0;
        end
        
       function collectedclutter = getClutterFromCurrentTransmit(obj,x, stang)
           % time span for current simulation step
           timespan = gpuArray.colon(0,obj.pNumRanges-1);
           timespan = reshape(timespan, obj.pNumRanges, 1);
           timespanMat = repmat(timespan, 1, obj.pNumAzimuthPatches);
           
           % at each step, for clutter echo from each patch, modulate
           % with Doppler shift
           currentclutter = obj.pgClutterReturn; % el x az
           currentclutter = arrayfun(@privDopplerShift, timespanMat, ...
               obj.pTimeOffset, obj.SampleRate, ...
               obj.pgDopplerShift, currentclutter);
           clear timespan timespanMat;
           numelang = numel(obj.pgElAngles);
            
           %Convention for subsequent computations: 
           %  Data cube is :elevation x position x azimuth

           %Transmit
           if obj.pNeedSteeringAngle || obj.pNeedCustomWeights
               tempresp = gpuArray.zeros(obj.pNumRanges, obj.pNumAzimuthPatches);
               for ii=1:getNumSubarrays(obj.Sensor)
                   resp = runSubArray(obj, obj.Sensor, ii, obj.pgAzAngles, obj.pgElAngles, stang);
                   tempresp = tempresp + squeeze((sum(resp,2)));
               end

               noiseMtxGain = getNoiseGainMatrix(obj.Sensor, ...
                   obj.OperatingFrequency, obj.PropagationSpeed, stang);
               gainFactor = real(sum(sum(noiseMtxGain)));
               
               txgfh = @(x) 10*log10(abs(x)*abs(x));
               arrayGainScaled = arrayfun(txgfh, tempresp) - 10*log10(gainFactor);%10*log10(sum(obj.pDOF)* gainFactor));
               currentclutter = currentclutter .* arrayfun(@db2mag, arrayGainScaled);
           end
           collectedclutter = complex(gpuArray.zeros(obj.pBufferLength,...
               obj.pDOF));
           clutterPerm = permute(currentclutter, [1 3 2]);
           
           %Receive
           if obj.pUseArray
               if isa(obj.Sensor, 'phased.internal.AbstractArray')
                   %its an array
                   pos = getElementPosition(obj.Sensor);
                   idx = phased.gpu.internal.decomposeCube(...
                       numel(obj.pElAngles), numel(obj.pAzAngles), ...
                       size(pos, 2) );
                   
                   for ii=1:numel(idx)
                       %returns sv which is el x pos x az
                       sv = obj.computeSteeringVector(pos(:,idx{ii}), ...
                           obj.OperatingFrequency, obj.PropagationSpeed, ...
                           obj.pgAzAngles, obj.pgElAngles);
                        
                       sv = scaleByElemResponse(obj.Sensor, ...
                           obj.pgAzAngles, obj.pgElAngles, ...
                           obj.OperatingFrequency, sv, idx{ii});
                       
                       %scale by the permuted clutter
                       sv = bsxfun(@times,sv, clutterPerm);
                        
                       %sum over the azimuth angles
                       collectedclutter(1:numelang, idx{ii}) = sum(sv, 3);
                   end
               else  %it's a subarray
                   for ii=1:getNumSubarrays(obj.Sensor)
                       tempresp = runSubArray(obj, ...
                           obj.Sensor, ii, obj.pgAzAngles, ...
                           obj.pgElAngles, stang);
                       collectedclutter(1:numelang,ii) = sum(...
                           bsxfun(@times, tempresp, clutterPerm),3);
                   end
               end
           else %single element
               if isElementFromAntenna(obj.Sensor)
                   hac = cloneSensor(obj.Sensor);
               else
                   hac = obj.Sensor;
               end
               elPat = getgpuElemResponse(hac, obj.pgAzAngles, ...
                   obj.pgElAngles, obj.OperatingFrequency);
               
               currentclutter = bsxfun(@times, currentclutter, squeeze(elPat));
               collectedclutter(1:numelang) = sum(currentclutter,2);
           end
            
           if isreal(collectedclutter)
               collectedclutter = complex(collectedclutter);
           end
            
           % filter using the transmit signal
           if obj.TransmitSignalInputPort
               collectedclutter = filter(x,1,collectedclutter);
           end
       end
        
        function [y, outputIdxOffset] = outputSamplesFromBuffer(...
                obj,y,outputIdxOffset,num_output_samples)
            
            output_sample_idx_increase = num_output_samples;
            buffer_len = obj.pBufferLength;
            if num_output_samples > buffer_len
                num_output_samples = buffer_len;
            end
            
            temp_buffer_idx = 1:num_output_samples;
            temp_output_idx = temp_buffer_idx+outputIdxOffset;
            y(temp_output_idx,:) = ...
                obj.pgClutterAtArrayOutputBuffer(temp_buffer_idx,:);
            updateClutterAtArrayBuffer(obj,num_output_samples);
            outputIdxOffset = outputIdxOffset + ...
                output_sample_idx_increase;
        end
    end
    
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            
            sensorProp = matlab.system.display.internal.Property('Sensor', ...
                'CustomPresenter', 'phased.internal.SensorArrayDialog', ...
                'CustomPresenterPropertyGroupsArgument', 'antenna', ...
                'Description', 'Sensor array');
            sensorTab = matlab.system.display.SectionGroup('Title', 'Sensor Array', 'PropertyList', {sensorProp});
            
            
            pSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false);
            pSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false);
            
            pTransmitSignalInputPort = matlab.system.display.internal.Property(...
                'TransmitSignalInputPort','IsGraphical',false);
            
            radarProps = {...
                'OperatingFrequency',...
                pTransmitSignalInputPort,...
                'TransmitERP',...
                'PlatformHeight',...
                'PlatformSpeed',...
                'PlatformDirection',...
                'BroadsideDepressionAngle',...
                pSeedSource,...
                pSeed};
            
            radarTab = matlab.system.display.SectionGroup('Title', 'Radar', ...
                'PropertyList', radarProps);
            
            waveformProps = {...
                'SampleRate',...
                'PRF', ...
                'PRFSelectionInputPort',...
                'OutputFormat',...
                'NumPulses',...
                'NumSamples'};
            
            waveformSection = matlab.system.display.Section('Title', 'Reflected signal', ...
                'PropertyList', waveformProps);
            clutterProps = {...
                'Gamma',...
                'EarthModel',...
                'MaximumRange',...
                'AzimuthCoverage',...
                'PatchAzimuthWidth',...
                'CoherenceTime',...
                'PropagationSpeed'};
            clutterSection = matlab.system.display.Section('Title', 'Clutter', ...
                'PropertyList', clutterProps);
            
            mainTab = matlab.system.display.SectionGroup('TitleSource', 'Auto', ...
                'Sections', [clutterSection waveformSection]);
            
            groups = [mainTab,radarTab, sensorTab];
        end
        
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title','Constant Gamma Clutter for GPU',...
                'Text',['Generate, using the GPU, constant gamma clutter reflected from homogeneous', ...
                ' terrain for a monostatic radar transmitting a narrowband', ...
                ' signal into free space. The radar is assumed to be at', ...
                ' constant altitude travelling at a constant speed. Use of this object', ...
                ' requires a Parallel Computing Toolbox license.']);
        end
        
        function retSz = getPulseOutputSize(sampleRate,PRF,outputFormat,...
                numPulses,numSamples,flag)
            %num of samples for each pulse in PRF vector
            pulseLength = round(sampleRate./PRF);
            if strcmp(outputFormat,'Pulses')
                numPRF = numel(PRF);
                if numPRF == 1
                    retSz = pulseLength*numPulses;
                elseif flag
                    retSz = max(pulseLength)*numPulses;
                else
                    %Get total number of samples of all pulse combinations
                    staggeredIndex = bsxfun(@plus,(1:numPRF)'-1,1:numPulses);
                    staggeredIndex = ...
                        phased.internal.AbstractPulseWaveform.getCircularIndex(staggeredIndex,numPRF);
                    currentSizes = sum(pulseLength(staggeredIndex),2);
                    %Get the upperbound
                    retSz = max(currentSizes);
                end
            else
                retSz = numSamples;
            end
        end
        
    end
    
    methods(Access=protected)
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('GPU Constant\nGamma Clutter');
        end
        
        function varargout = getOutputNamesImpl(~)
            varargout = {''};
        end
        
        function varargout = getInputNamesImpl(obj)
            needSteerInput = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmpi(obj.Sensor.SubarraySteering,'None',1);
            needTransmitInput = obj.TransmitSignalInputPort;
            needPRFSelectionInputPort = obj.PRFSelectionInputPort;
            varargout = {};
            if needTransmitInput
                varargout{end+1} = 'X';
            end
            if needSteerInput
                if strncmpi(obj.Sensor.SubarraySteering,'Custom',1)
                    varargout{end+1} = 'WS';
                else
                    varargout{end+1} = 'Steer';
                end
            end
            if needPRFSelectionInputPort
                varargout{end+1} = 'PRFIdx';
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
            allGroups = getPropertyGroupsLongImpl@phased.gpu.internal.AbstractClutterSimulator(obj);
            group = flattenGroupsAndMoveSensorToTop('Sensor',allGroups);
        end
        
        function sts = getSampleTimeImpl(obj)
            if getNumInputs(obj) == 0
                N = obj.getPulseOutputSize(...
                    obj.SampleRate, obj.PRF,obj.OutputFormat, ...
                    obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
                st = N/obj.SampleRate;
                sts = createSampleTime(obj,'Type','Discrete',...
                    'SampleTime',st);
            else
                sts = createSampleTime(obj,'Type','Inherited');
            end
        end
        
        function varargout = getOutputSizeImpl(obj)
            maxNumSample =  obj.getPulseOutputSize(...
                obj.SampleRate, obj.PRF,obj.OutputFormat, ...
                obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
            varargout{1} = [maxNumSample getDOF(obj.Sensor)];
        end
        
        function outtype = getOutputDataTypeImpl(obj) %#ok<MANU>
            outtype = 'double';
        end
        
        function ocplx = isOutputComplexImpl(obj) %#ok
            ocplx = true;
        end
        
        function varargout = isOutputFixedSizeImpl(obj)
            if (strcmp(obj.OutputFormat,'Pulses') && numel(obj.PRF) ~= 1)...
                    || obj.PRFSelectionInputPort
                varargout{1} = false;
            else
                varargout{1} = true;
            end
        end
    end
    methods (Access=private)
        validateSampleRate(obj,N)
    end
end


% [EOF]
