classdef (Hidden) AbstractConstantGammaClutter < phased.internal.AbstractClutterSimulator 
%This class is for internal use only. It may be removed in the future.

%   Copyright 2017 The MathWorks, Inc.

%   Reference
%   [1] Maurice Long, Radar Reflectivity of Land and Sea, 3rd Ed. Artech
%       House, 2001
%   [2] J. R. Guerci, Space Time Adaptive Processing for Radar, Artech
%       House, 2003
%   [3] James Ward, Space-time Adaptive Processing for Airborne Radar,
%       Lincoln Lab Technical Report, 1994
%   [4] David Greig, Time Domain Clutter Model for Airborne Pulse Doppler
%       Radar, IEE Colloquium on Radar System Modeling, 1998

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
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
    
    methods (Access = protected)
        function obj = AbstractConstantGammaClutter(varargin)
            %ConstantGammaClutter   Construct the ConstantGammaClutter
            %class.
            obj@phased.internal.AbstractClutterSimulator(varargin{:});
        end   
    end
    
    methods (Access = protected)
        function nrcs = getNRCS(obj)
            nrcs = db2pow(obj.Gamma).*sind(obj.pGrazingAngles) + ...
                10.*exp(-0.2.*(90-obj.pGrazingAngles));  % [Long03]
        end
        
        function setupImpl(obj,varargin)
            setupImpl@phased.internal.AbstractClutterSimulator(obj,varargin{:});
            if getNumInputs(obj)>0
                N = obj.getPulseOutputSize(...
                    obj.SampleRate, obj.PRF,obj.OutputFormat, ...
                    obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
                validateSampleRate(obj,N,obj.SampleRate);
            end
        end
        
        function flag = isInDVRMode(obj) %#ok<MANU>
            flag = false;
        end

        function y = stepImpl(obj,x,stang,prfidx)
            
            % initialize random clutter
            if ~obj.pTimeOffset
                calcRandomClutter(obj);
            end
            
            if obj.PRFSelectionInputPort
                if obj.TransmitSignalInputPort 
                    if ~obj.pNeedSteeringAngle
                        prfidx = stang;
                    end
                else
                    if obj.pNeedSteeringAngle
                        prfidx = stang;
                    else
                        prfidx = x;
                    end
                end
                N = calcTotalOutputSamples(obj,prfidx);
                if ~obj.pOutputBySamples && isInDVRMode(obj)
                    setNumTicksUntilNextHit(obj,N);
                end
            else
                N = calcTotalOutputSamples(obj);
                if ~obj.pOutputBySamples && isInDVRMode(obj)
                    setNumTicksUntilNextHit(obj,N);
                end
            end
            y = complex(zeros(N,obj.pDOF));
            outputIdxOffset = 0;
            
            while N > 0
                if obj.pRemainingSamplesFromLastPulse > 0
                    % there are remaining samples from last pulse
                    
                    % this can only happen in samples mode since if output
                    % by pulse, there will be no remaining samples
                    
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
                    if obj.TransmitSignalInputPort
                        if obj.pNeedSteeringAngle || obj.pNeedCustomWeights
                            collectedclutter = getClutterFromCurrentTransmit(obj,x,stang);
                        else
                            collectedclutter = getClutterFromCurrentTransmit(obj,x);
                        end
                    else
                        if obj.pNeedSteeringAngle || obj.pNeedCustomWeights 
                            stang = x;
                            collectedclutter = getClutterFromCurrentTransmit(obj,stang);
                        else
                            collectedclutter = getClutterFromCurrentTransmit(obj);
                        end
                    end
                    
                    % cache the clutter from the current simulation step
                    % into buffer, which holds return from each clutter
                    % patch up to current simulation step and not in the
                    % output yet.
                    obj.pClutterAtArrayOutputBuffer = ...
                        obj.pClutterAtArrayOutputBuffer + collectedclutter;
                    
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
            obj.pClutterAtArrayOutputBuffer(1:obj.pBufferLength-N,:) = ...
                obj.pClutterAtArrayOutputBuffer(N+1:obj.pBufferLength,:);
            obj.pClutterAtArrayOutputBuffer(obj.pBufferLength-N+1:end,:) = 0;
        end
        
        function collectedclutter = getClutterFromCurrentTransmit(obj,x,stangArg)
            % time span for current simulation step
            timespan = (0:obj.pNumRanges-1).'/obj.SampleRate + ...
                obj.pTimeOffset;
            
            % at each step, for clutter echo from each patch, modulate
            % with Doppler shift
            currentclutter = obj.pClutterReturn;
            for m = 1:obj.pNumAzimuthPatches
                currentclutter(:,m) = currentclutter(:,m).*...
                    exp(1i*2*pi*obj.pDopplerShift(:,m).*timespan);
            end
            
            % if the subarray has steering capability, modulate the signal
            % with appropriate gain
            if obj.pNeedSteeringAngle
                if ~obj.TransmitSignalInputPort
                    stangx = x;
                else
                    stangx = stangArg;
                end
                if isscalar(stangx)
                    stang = [stangx; 0];
                else
                    stang = stangx;
                end
                azstang = stang(1,:);
                elstang = stang(2,:);
                cond = any(azstang>180) || any(elstang>90);
                if cond
                    coder.internal.errorIf(cond,'phased:step:AzElNotLessEqual');
                end
                cond = any(azstang<-180) || any(elstang<-90);
                if cond
                    coder.internal.errorIf(cond,'phased:step:AzElNotGreaterEqual');
                end
                for m = 1:obj.pNumAzimuthPatches
                    ang = [obj.pAzAngles(m)*ones(1,obj.pNumRanges);...
                        obj.pElAngles];
                    currentclutter(:,m) = currentclutter(:,m) .* ...
                        db2mag(step(obj.cArrayGain,obj.OperatingFrequency,ang,stang));
                end
            elseif obj.pNeedCustomWeights
                if ~obj.TransmitSignalInputPort
                    ws = x;
                else
                    ws = stangArg;
                end
                for m = 1:obj.pNumAzimuthPatches
                    ang = [obj.pAzAngles(m)*ones(1,obj.pNumRanges);...
                        obj.pElAngles];
                    currentclutter(:,m) = currentclutter(:,m) .* ...
                        db2mag(step(obj.cArrayGain,obj.OperatingFrequency,ang,ws));
                end
            end
            
            % collect clutter from each patch for current simulation
            % step
            collectedclutter = complex(zeros(obj.pBufferLength,...
                obj.pDOF));
            if obj.pNeedSteeringAngle
                for m = 1:obj.pNumRanges
                    ang = phased.internal.arbazel2azel([obj.pAzAngles;...
                        obj.pElAngles(m)*ones(1,obj.pNumAzimuthPatches)]);
                    collectedclutter(m,:) = step(obj.cCollector,...
                        currentclutter(m,:),ang,stang);
                end
            elseif obj.pNeedCustomWeights
                for m = 1:obj.pNumRanges
                    ang = phased.internal.arbazel2azel([obj.pAzAngles;...
                        obj.pElAngles(m)*ones(1,obj.pNumAzimuthPatches)]);
                    collectedclutter(m,:) = step(obj.cCollector,...
                        currentclutter(m,:),ang,ws);
                end
            else
                for m = 1:obj.pNumRanges
                    ang = phased.internal.arbazel2azel([obj.pAzAngles;...
                        obj.pElAngles(m)*ones(1,obj.pNumAzimuthPatches)]);
                    collectedclutter(m,:) = step(obj.cCollector,...
                        currentclutter(m,:),ang);
                end
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
                obj.pClutterAtArrayOutputBuffer(temp_buffer_idx,:);
            updateClutterAtArrayBuffer(obj,num_output_samples);
            outputIdxOffset = outputIdxOffset + ...
                output_sample_idx_increase;
        end
        
    end
    
    methods (Access = protected) %For Simulink propagation and mask
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
        function varargout = getOutputSizeImpl(obj)
            maxNumSample =  obj.getPulseOutputSize(...
                obj.SampleRate, obj.PRF,obj.OutputFormat, ...
                obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
            varargout{1} = [maxNumSample getDOF(obj.Sensor)];
        end
        function varargout = isOutputFixedSizeImpl(obj)
            if (strcmp(obj.OutputFormat,'Pulses') && numel(obj.PRF) ~= 1) ...
                    || obj.PRFSelectionInputPort
                varargout{1} = false;
            else
                varargout{1} = true;
            end
        end
        function varargout = getOutputDataTypeImpl(obj)  %#ok<MANU>
            varargout{1} = 'double';
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout{1} = true;
        end
        
        function sts = getSampleTimeImpl(obj)
            if obj.PRFSelectionInputPort && ...
                    strcmp(obj.SimulationTimeSource,'Inherit from Simulink engine')
                sts = createSampleTime(obj,'Type','Inherited');
            else
                if strcmp(obj.OutputFormat, 'Samples') || numel(obj.PRF) == 1
                    N = obj.getPulseOutputSize(...
                        obj.SampleRate, obj.PRF,obj.OutputFormat, ...
                        obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
                    st = phased.internal.samprate2time(obj.SampleRate,N);
                    sts = createSampleTime(obj,'Type','Discrete',...
                        'SampleTime',st);
                else
                    sts = createSampleTime(obj,'Type','Controllable',...
                        'TickTime',1/obj.SampleRate);
                end
            end
        end
        
    end
    
    methods (Access = protected) 
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
            allGroups = getPropertyGroupsLongImpl@phased.internal.AbstractClutterSimulator(obj);
            group = flattenGroupsAndMoveSensorToTop(obj,'Sensor',allGroups);
        end
    end
    
    methods (Static,Hidden,Access=protected)
        function retSz = getPulseOutputSize(sampleRate,PRF,outputFormat,...
                numPulses,numSamples,flag)
            %num of samples for each pulse in PRF vector
            pulseLength = round(sampleRate./PRF(:));
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
end


% [EOF]
