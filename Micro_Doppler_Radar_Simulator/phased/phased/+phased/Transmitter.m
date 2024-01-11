classdef (Sealed,StrictDefaults) Transmitter < phased.internal.AbstractVarSizeEngine & ...
        matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates 
%Transmitter Pulse transmitter
%   H = phased.Transmitter creates a transmitter System object, H. This
%   object transmits the input waveform samples with specified peak power.
%
%   H = phased.Transmitter(Name,Value) creates a transmitter object, H,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax:
%
%   Y = step(H,X) returns the transmitted signal Y, based on the input
%   waveform X. Y is the amplified X where the amplification is based on
%   the characteristics of the transmitter, such as the peak power and the
%   gain.
%
%   [Y,TR] = step(H,X) returns additional output TR as the on/off
%   status of the transmitter when the InUseOutputPort property is true.
%   TR is a logical vector where true indicates the transmitter is on
%   for the corresponding sample time, and false indicates the transmitter
%   is off. The transmitter is on only when the samples in X has non-zero
%   envelope.
%
%   [Y,PH] = step(H,X) returns additional output PH as the random
%   phase noise added to each transmitted sample when the
%   CoherentOnTransmit property is false and the PhaseNoiseOutputPort
%   property is true. PH is a vector which has the same dimension as
%   Y. Each element in PH contains the random phase, between 0 and
%   2*pi, added to the corresponding sample in Y by the transmitter.
%
%   You can combine optional output arguments when their enabling
%   properties are set. Optional outputs must be listed in the same order
%   as the order of the enabling properties. For example,
%
%   [Y,TR,PH] = step(H,X)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   Transmitter methods:
%
%   step     - Transmit pulses (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a transmitter object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of transmitter object
%
%   Transmitter properties:
%
%   PeakPower            - Peak power
%   Gain                 - Gain
%   LossFactor           - Loss factor
%   InUseOutputPort      - Enable transmitter status output
%   CoherentOnTransmit   - Preserve coherence among pulses
%   PhaseNoiseOutputPort - Enable pulse phase noise output
%   SeedSource           - Source of seed for random number generator
%   Seed                 - Seed for random number generator
%
%   % Example:
%   %   Transmit a pulse containing a linear FM waveform. The sample rate 
%   %   is 10 MHz and the pulse repetition frequency is 50 kHz.
%
%   fs = 1e7;
%   waveform = phased.LinearFMWaveform('SampleRate',fs,...
%                   'PulseWidth',1e-5,'SweepBandwidth',5e6);
%   x = waveform(); 
%   tx = phased.Transmitter;
%   y = tx(x); 
%
%   See also phased, phased.ReceiverPreamp, phased.Radiator.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [2] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001
%   [3] Byron Edde, Radar: Principles, Technology, Applications, Prentice
%       Hall, 1993


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %PeakPower Peak power (W)
    %   Specify the transmit peak power (in watts) as a positive scalar.
    %   The default value of this property is 5000.    
    PeakPower = 5000;
    %Gain Gain (dB)
    %   Specify the transmit gain (in dB) as a real scalar. The default
    %   value of this property is 20.
    Gain = 20;
    %LossFactor Loss factor (dB)
    %   Specify the transmit loss factor (in dB) as a non-negative scalar.
    %   The default value of this property is 0.
    LossFactor = 0;    
end
properties (Nontunable, Logical) 
    %InUseOutputPort Enable transmitter status output
    %   Set this property to true to output the transmitter in-use status
    %   for each output sample. 1's indicate the transmitter is on, and 0's
    %   indicate the transmitter is off. Set this property to false to not
    %   output the transmitter status. The default value of this property
    %   is false.
    InUseOutputPort = false;
    %CoherentOnTransmit     Preserve coherence among pulses
    %   Specify whether to preserve coherence among transmitted pulses.
    %   When you set this property to true, the transmitter does not
    %   introduce any random phase to the output pulses. When you set this
    %   property to false, transmitter adds a random phase noise to each
    %   transmitted pulse. The default value of this property is true.
    CoherentOnTransmit = true;
    %PhaseNoiseOutputPort   Enable pulse phase noise output
    %   Set this property to true to output the introduced transmitter
    %   random phase noise for each output sample. This information could
    %   be used in the receiver to simulate coherent on receive systems.
    %   Set this property to false to not output the random phase noise.
    %   The default value of this property is false. This property applies
    %   when you set the CoherentOnTransmit property to false.
    PhaseNoiseOutputPort = false;    
end
properties (Nontunable)
    %SeedSource   Source of seed for random number generator
    %   Specify how the random numbers are generated as one of 'Auto' |
    %   'Property', where the default is 'Auto'. When you set this property
    %   to 'Auto', the random numbers are generated using the default
    %   MATLAB random number generator. When you set this property to
    %   'Property', a private random number generator is used with a seed
    %   specified by the value of the Seed property. This property applies
    %   when you set the CoherentOnTransmit property to false.
    %
    %   To use this object with Parallel Computing Toolbox software, set
    %   this property to 'Auto'.
    SeedSource = 'Auto';
    %Seed     Seed for random number generator
    %   Specify the seed for the random number generator as a non-negative
    %   integer. The integer must be less than 2^32. This property applies
    %   when you set the SeedSource property to 'Property'. The default
    %   value of this property is 0. This property applies when you set the
    %   CoherentOnTransmit property to false and the SeedSource property to
    %   'Property'.
    Seed = 0;    
end

properties(Constant, Hidden)
    SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
end

properties (Access=private, Nontunable)
    pAmpCoeff  % private property to hold amplification coefficient
    cNoiseSource % private random stream
end

properties (Access=private)
    pOutputSize % hold output size,
    pPreviousStatus % private status for last sample in previous sim step
    pPreviousPhaseNoise % private phase noise sample for previous sim step
    
end

methods
    function set.Gain(obj,value)
        validateattributes( value, { 'double' }, { 'scalar', 'real', 'finite' }, '', 'Gain');
        obj.Gain = value;
    end
    function set.LossFactor(obj,value)
        validateattributes( value, { 'double' }, { 'scalar', 'nonnegative', 'finite' }, '', 'LossFactor');
        obj.LossFactor = value;
    end
    function set.PeakPower(obj,value)
        validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'PeakPower');
        obj.PeakPower = value;
    end
    function set.Seed(obj,val)
        validateattributes(val,{'double'},{'scalar','nonnegative',...
            'finite','nonnan','nonempty'},'phased.Transmitter',...
            'Seed');
        obj.Seed = val;
    end
end

methods
    function obj = Transmitter(varargin)
        setProperties(obj, nargin, varargin{:});
    end    
end

methods (Access = protected)
    
    function validateInputsImpl(obj,x)
        cond =  ~isa(x,'double');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','double');
        end
        cond =  ~iscolumn(x) || isempty(x);
        if cond
            coder.internal.errorIf(cond, ...
                                   'MATLAB:system:inputMustBeColVector','X');
        end
        validateNumChannels(obj,x);
    end
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        if obj.CoherentOnTransmit
            if strcmp(prop,'PhaseNoiseOutputPort') || ...
                    strcmp(prop,'SeedSource') || ...
                    strcmp(prop,'Seed')
                flag = true;
            end
        elseif (obj.SeedSource(1)  == 'A') && ... %Auto
                strcmp(prop,'Seed')
            flag = true;
        end
    end
    
    function flag = isInputSizeLockedImpl(~,~)
        flag = false;
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
        if isLocked(obj)
            if ~obj.CoherentOnTransmit
                s.cNoiseSource = saveobj(obj.cNoiseSource);
            end
            s.pAmpCoeff = obj.pAmpCoeff;
            s.pPreviousStatus = obj.pPreviousStatus;
            s.pPreviousPhaseNoise = obj.pPreviousPhaseNoise;
            s.pOutputSize = obj.pOutputSize;
        end
    end

    function s = loadSubObjects(obj,s,wasLocked)
        if isfield(s,'isLocked')    % for backwards compatibility
            if s.isLocked   
                obj.cNoiseSource = phased.internal.NoiseSource(...
                    'SeedSource','Property','Seed',s.pSeed);
                s = rmfield(s,'pRandState'); 
                s = rmfield(s,'pSeed');
            end
            s = rmfield(s,'isLocked');
        elseif wasLocked
            if ~s.CoherentOnTransmit
                obj.cNoiseSource = phased.internal.NoiseSource.loadobj(s.cNoiseSource);
                s = rmfield(s,'cNoiseSource');
            end
        end
    end

    function loadObjectImpl(obj,s,wasLocked)
        s = loadSubObjects(obj,s,wasLocked);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function [yout,TxOn_or_noise, pnoise] = stepImpl(obj,x)
        y = obj.pAmpCoeff*x;
        Tx_on = x~=0;
        nsamp = size(x,1);
        % if obj.pOutputSize == -1
        %     nsamp = size(x,1);
        %     obj.pOutputSize = nsamp;
        % else
        %     nsamp = obj.pOutputSize;
        % end 
        if ~obj.CoherentOnTransmit
            % Get indexes of all rising edges
            pulseStartIdx = find(diff([obj.pPreviousStatus; Tx_on]) > 0);
            if ~isempty(pulseStartIdx)
                NumSegments = numel(pulseStartIdx)+1;
                phaseNoisePulseStart = [1;pulseStartIdx(1:end)];
                phaseNoisePulseEnd = [pulseStartIdx(1:end)-1;nsamp];
                phaseNoiseRandomSample = [obj.pPreviousPhaseNoise;...
                    step(obj.cNoiseSource,1,[NumSegments-1 1])];
            else
                NumSegments = 1;
                phaseNoisePulseStart = 1;
                phaseNoisePulseEnd = nsamp;
                phaseNoiseRandomSample = obj.pPreviousPhaseNoise;
            end
            
            pnoise = complex(zeros(nsamp,1));

            for m = 1:NumSegments
                pnoise(...
                    phaseNoisePulseStart(m):phaseNoisePulseEnd(m)) = ...
                    phaseNoiseRandomSample(m);
            end
            obj.pPreviousStatus = Tx_on(end);
            obj.pPreviousPhaseNoise = pnoise(end);
            pnoise = 2*pi*pnoise;
            yout = addPhaseNoise(y,pnoise);
            if ~obj.InUseOutputPort && obj.PhaseNoiseOutputPort
                TxOn_or_noise = pnoise; % phase noise is second output
            else
                TxOn_or_noise = Tx_on;
            end
        else
            yout = y;
            TxOn_or_noise = Tx_on;
        end
    end
    
    function num = getNumOutputsImpl(obj)
        num = 1;
        if obj.InUseOutputPort
            num = num+1;
        end
        if ~obj.CoherentOnTransmit && obj.PhaseNoiseOutputPort
            num = num+1;
        end
    end
    
    function resetImpl(obj)
        obj.pPreviousStatus = false;
        obj.pPreviousPhaseNoise = complex(0);
        if ~obj.CoherentOnTransmit
            reset(obj.cNoiseSource);
        end
    end
    
    function releaseImpl(obj)
        releaseImpl@phased.internal.AbstractVarSizeEngine(obj);
        if ~obj.CoherentOnTransmit
            release(obj.cNoiseSource);
        end
    end
    
    function setupImpl(obj,x)
        obj.pAmpCoeff = sqrt(obj.PeakPower*db2pow(obj.Gain-obj.LossFactor));
        obj.pNumInputChannels = getNumChannels(obj,x);
        obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        if ~obj.CoherentOnTransmit
            obj.cNoiseSource = phased.internal.NoiseSource(...
                'SeedSource',obj.SeedSource);
            if (obj.SeedSource(1) == 'P') %Property
                obj.cNoiseSource.Seed = obj.Seed;
            end
        end
        obj.pOutputSize = -1; % Initialize to -1
    end   
    
    function processInputSizeChangeImpl(obj,x)
        obj.pOutputSize = size(x,1);
    end
    
    function flag = isOutputComplexityLockedImpl(obj,index)
        flag = false;  % index == 1
        if obj.InUseOutputPort && (index == 2)
            flag = true;
        end
        if obj.PhaseNoiseOutputPort
            if obj.InUseOutputPort && (index == 3)
                flag = true;
            else  
                if index == 2
                    flag = true;
                end
            end
        end
    end

    function flag = isInputComplexityLockedImpl(obj,~)  %#ok<INUSD>
        flag = false;  % index == 1
    end
end

methods (Access = protected, Static, Hidden)
    function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:TransmitterTitle')),...
          'Text',getString(message('phased:library:block:TransmitterDesc')));
    end

    function groups = getPropertyGroupsImpl
        groups = matlab.system.display.Section(...
            'phased.Transmitter');
        dSeedSource = matlab.system.display.internal.Property(...
            'SeedSource','IsGraphical',false,'UseClassDefault',false,...
            'Default','Property');
        dSeed = matlab.system.display.internal.Property(...
            'Seed','IsGraphical',false,'UseClassDefault',false,...
            'Default','randi(65535,1)');
        groups.PropertyList = setdiff(groups.PropertyList,...
            {'SeedSource','Seed'},'stable');
        groups.PropertyList = [groups.PropertyList,...
            {dSeedSource,dSeed}];
   end
end

methods (Access = protected)
    function varargout = getInputNamesImpl(obj)  %#ok<MANU>
        varargout = {''};
    end

    function varargout = getOutputNamesImpl(obj) 
        if obj.InUseOutputPort
            if ~obj.CoherentOnTransmit && obj.PhaseNoiseOutputPort
                varargout = {'Y','TR','Ph'};
            else
                varargout = {'Y','TR'};
            end
        else
            if ~obj.CoherentOnTransmit && obj.PhaseNoiseOutputPort
                varargout = {'Y','Ph'};
            else
                varargout = {'Y'};
            end
        end
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Transmitter');
    end
    function varargout = getOutputSizeImpl(obj)
        inSz = propagatedInputSize(obj,1);
        varargout = {inSz, inSz, inSz}; 

    end
    function varargout = isOutputFixedSizeImpl(obj)
        inFSz = propagatedInputFixedSize(obj, 1);
        varargout = {inFSz, inFSz, inFSz};
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj,1);
        if obj.InUseOutputPort
            varargout{2} = 'logical'; 
            if ~obj.CoherentOnTransmit && obj.PhaseNoiseOutputPort
                varargout{3} = 'double';
            end
        else
            if ~obj.CoherentOnTransmit && obj.PhaseNoiseOutputPort
                varargout{2} = 'double';
            end
        end
    end
    function varargout = isOutputComplexImpl(obj)
        varargout{1} = true;
        if obj.InUseOutputPort
            varargout{2} = false; 
            if ~obj.CoherentOnTransmit && obj.PhaseNoiseOutputPort
                varargout{3} = true;
            end
        else
            if ~obj.CoherentOnTransmit && obj.PhaseNoiseOutputPort
                varargout{2} = true;
            end
        end
    end
end

end

function yout = addPhaseNoise(y,pnoise)

coder.inline('never');
yout = y.*exp(1i*pnoise);

end