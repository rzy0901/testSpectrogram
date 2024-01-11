classdef AbstractReceiverPreamp < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%AbstractReceiverPreamp   Abstract class for receiver preamp

%   Copyright 2016-2017 The MathWorks, Inc.



%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Gain     Gain (dB)
        %   A scalar containing the gain (in dB) of the pulse receiver. The
        %   default value of this property is 20.
        Gain = 20
        %LossFactor   Loss factor (dB)
        %   A scalar containing the loss factor (in dB) of the pulse
        %   receiver. The default value of this property is 0.
        LossFactor = 0
        %NoiseMethod     Noise specification method
        %   Specify how to compute noise power using one of 'Noise power' |
        %   'Noise temperature', where the default is 'Noise temperature'.
        %
        %   If you set this property to 'Noise temperature', a complex
        %   baseband noise is added to the input signal with noise power
        %   computed from the ReferenceTemperature, NoiseFigure, and
        %   SampleRate properties. The SampleRate property specifies the
        %   sample rate of the input signal and is the same as the noise
        %   bandwidth.
        %
        %   If you set this property to 'Noise power', a noise, whose power
        %   is specified in the NoisePower property, is added to the input
        %   signal. 
        NoiseMethod = 'Noise temperature'
        %NoiseFigure   Noise figure (dB)
        %   A scalar containing the noise figure (in dB) of the pulse
        %   receiver. If the receiver has multiple channels/sensors, the
        %   noise figure applies to each channel/sensor. The default value
        %   of this property is 0. This property is only applicable when
        %   you set the NoiseMethod property to 'Noise temperature'.
        NoiseFigure = 0     
        %ReferenceTemperature   Reference temperature (K)
        %   A scalar containing the reference temperature of the receiver
        %   (in kelvin). If the receiver has multiple channels/sensors, the
        %   reference temperature applies to each channel/sensor. The
        %   default value of this property is 290. This property is only
        %   applicable when you set the NoiseMethod property to 'Noise
        %   temperature'.
        ReferenceTemperature = 290
    end
    
    properties (Nontunable)
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a positive scalar. The
        %   default value of this property is 1e6 (1 MHz). This property is
        %   only applicable when you set the NoiseMethod property to 'Noise
        %   temperature'.
        SampleRate = 1e6
        %NoisePower     Noise power (W)
        %   Specify the noise power (in W) as a positive scalar. The
        %   default value of this property is 1. This property is only
        %   applicable when you set the NoiseMethod property to 'Noise
        %   power'.
        NoisePower = 0
        %NoiseComplexity    Noise complexity
        %   Specify the noise complexity as one of 'Complex' | 'Real',
        %   where the default is 'Complex'. When you set this property to
        %   'Complex', the noise power is evenly divided between real and
        %   imaginary channels.
        NoiseComplexity = 'Complex'
    end 
    
    properties (Hidden, Nontunable)
        %NoiseBandwidth   Noise bandwidth (Hz)
        %   A scalar containing the bandwidth of noise spectrum (in Hz) at
        %   the pulse receiver. If the receiver has multiple
        %   channels/sensors, the noise bandwidth applies to each
        %   channel/sensor. The default value of this property is 1e6.
        NoiseBandwidth = 1e6
    end
     
    properties (Nontunable, Logical) 
        %EnableInputPort    Enable enabling signal input
        %   Set this property to true to add input to specify the receiver
        %   enabling signal. Set this property to false to not add input to
        %   specify the receiver enabling signal. The default value of this
        %   property is false.
        EnableInputPort = false
        %PhaseNoiseInputPort    Enable phase noise input
        %   Set this property to true to add input to specify the phase
        %   noise for each incoming sample. This information can be used to
        %   emulate coherent on receive systems. Set this property to false
        %   to not add input to specify the phase noise. The default value
        %   of this property is false.
        PhaseNoiseInputPort = false
    end
    
    properties (Nontunable)
        %SeedSource   Source of seed for random number generator
        %   Specify how the random numbers are generated as one of 'Auto' |
        %   'Property', where the default is 'Auto'. When you set this
        %   property to 'Auto', the random numbers are generated using the
        %   default MATLAB random number generator. When you set this
        %   property to 'Property', a private random number generator is
        %   used with a seed specified by the value of the Seed property.
        %
        %   To use this object with Parallel Computing Toolbox software,
        %   set this property to 'Auto'.
        SeedSource = 'Auto'
        %Seed     Seed for random number generator
        %   Specify the seed for the random number generator as a
        %   non-negative integer. The integer must be less than 2^32. This
        %   property applies when you set the SeedSource property to
        %   'Property'. The default value of this property is 0.
        Seed = 0
    end    
    
    properties(Constant, Hidden)
        SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        NoiseMethodSet = matlab.system.StringSet(...
            {'Noise power','Noise temperature'});
        NoiseComplexitySet = matlab.system.StringSet(...
            {'Complex','Real'});
    end
    
    properties(Access = protected, Nontunable)
        cNoiseSource
    end
    
    properties(Access = protected)
        pNoiseSamplePower = -1
    end
    
    properties (Access = protected)
        pXSize
    end
        
    methods (Access = protected)

        function obj = AbstractReceiverPreamp(varargin)
            %ReceiverPreamp   Construct the ReceiverPreamp class.
            setProperties(obj, nargin, varargin{:});
        end
        
    end
    
    methods 
        
        function set.Gain(obj, val)
            validateattributes(val,{'double'},{'nonempty','finite',...
                'scalar','real'},'phased.ReceiverPreamp','Gain');
            obj.Gain = val;                  
        end
       
        function set.LossFactor(obj, val)
            validateattributes(val,{'double'},{'nonempty','finite',...
                'scalar','nonnegative'},'phased.ReceiverPreamp',...
                'LossFactor');
            obj.LossFactor = val;                  
        end
        
        function set.NoiseBandwidth(obj, val)
            if isempty(coder.target)	 
                warning(message('phased:system:System:NoiseBandwidthWarning',...	 
                    'NoiseBandwidth'));	 
            end	 
            sigdatatypes.validateFrequency(val,...
                'phased.ReceiverPreamp','NoiseBandwidth',{'scalar'});
            obj.NoiseBandwidth = val;             
        end
        
        function val = get.NoiseBandwidth(obj)
            if isempty(coder.target)	 
                warning(message('phased:system:System:NoiseBandwidthWarning',...	 
                    'NoiseBandwidth'));	 
            end	 
            val = obj.NoiseBandwidth;
        end
        
        function set.NoiseFigure(obj, val)
            validateattributes(val,{'double'},...
                {'scalar','finite','nonnegative'},...
                'phased.ReceiverPreamp','NoiseFigure');
            obj.NoiseFigure = val;             
        end
        
        function set.ReferenceTemperature(obj, val)
            validateattributes(val,{'double'},...
                {'scalar','finite','positive'},...
                'phased.ReceiverPreamp','ReferenceTemperature');
            obj.ReferenceTemperature = val;             
        end
        
        function set.SampleRate(obj, value)
            validateattributes(value,{'double'}, {'scalar',...
                'positive','finite'},...
                'phased.ReceiverPreamp','SampleRate');
            obj.SampleRate = value;
        end
        
        function set.Seed(obj,val)
            validateattributes(val,{'double'},{'scalar','nonnegative',...
                'finite','nonnan','nonempty'},'phased.ReceiverPreamp',...
                'Seed');
            obj.Seed = val;
        end
        
    end
    
    methods(Access = protected, Abstract)
        fs = computeSampleRate(obj)
    end
    
    methods(Access = protected)
        
        function processInputSizeChangeImpl(obj,x,~,~)
            obj.pXSize = size(x);
        end
    
        function setupImpl(obj,x,enrx,pnoise) %#ok<INUSD>
            obj.pXSize = [-1 -1]; % Initialize to -1
            obj.cNoiseSource = phased.internal.NoiseSource(...
                'Distribution','Gaussian',...
                'SeedSource',obj.SeedSource);
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
            if isComplexityPropagated(obj)
                if ~isreal(x)
                    obj.cNoiseSource.OutputComplex = true;
                else
                    obj.cNoiseSource.OutputComplex = false;
                end
            else
                if strcmp(obj.NoiseComplexity,'Complex')
                    obj.cNoiseSource.OutputComplex = true;
                else
                    obj.cNoiseSource.OutputComplex = false;
                end
            end
            
            if (obj.SeedSource(1) == 'P')
                obj.cNoiseSource.Seed = obj.Seed;
            end
            if strcmp(obj.NoiseMethod,'Noise power')
                obj.pNoiseSamplePower = obj.NoisePower;
            else
                Fs = computeSampleRate(obj);
                obj.pNoiseSamplePower = noisepow(Fs,obj.NoiseFigure,...
                    obj.ReferenceTemperature);
            end
            
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = false;  % index == 1
            if obj.EnableInputPort && (index == 2)
                flag = true;
            end
            if obj.PhaseNoiseInputPort
                if obj.EnableInputPort && (index == 3)
                    flag = true;
                else  
                    if index == 2
                        flag = true;
                    end
                end
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
            release(obj.cNoiseSource);
        end

        function resetImpl(obj)
            reset(obj.cNoiseSource);
        end
        
        function flag = isInputSizeLockedImpl(~,~)
            flag = false;
        end
        
        function s = saveObjectImpl(obj)
            ws = warning('off','phased:system:System:NoiseBandwidthWarning');
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.cNoiseSource = saveobj(obj.cNoiseSource);
                s.pNoiseSamplePower = obj.pNoiseSamplePower;
                s.pXSize = obj.pXSize;
            end
            warning(ws);
        end
        
        function s = loadSubObjects(obj,s)
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cNoiseSource = eval(...
                        sprintf('%s.loadobj(s.cNoiseSource)',s.cNoiseSource.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cNoiseSource');
                end
                s = rmfield(s,'isLocked');
            end
        end
        
        function y = stepImpl(obj,x,enrx,pnoise)

            sz_x = [size(x,1) obj.pValidatedNumInputChannels];
            
            if obj.EnableInputPort
                enrx = logical(enrx);
                x(~enrx,:) = 0;
            elseif obj.PhaseNoiseInputPort
                pnoise = enrx;
            end
            
            if obj.PhaseNoiseInputPort
                xc = bsxfun(@times,exp(-1i*pnoise),x);
            else
                xc = x;
            end
            
            y = sqrt(db2pow(obj.Gain-obj.LossFactor))*xc + ...
                step(obj.cNoiseSource,obj.pNoiseSamplePower,sz_x);         
        end
        
        function validateInputsImpl(obj,x,enrx,pnoise)
            
            coder.extrinsic('mat2str');
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
            
            validateNumChannels(obj,x);
            
            sz_x = size(x);
            if obj.EnableInputPort
                cond =  ~islogical(enrx) && ~isa(enrx,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:system:System:incorrectType');
                end
                sz_enrx = size(enrx);
                cond =  ~isequal(sz_enrx, [sz_x(1) 1]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions', 'TR', ...
                        coder.const(mat2str([sz_x(1) 1])), ...
                        coder.const(mat2str(sz_enrx)));
                end
                cond =  ~isreal(enrx);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:step:NeedReal', 'TR');
                end
            end
            
            if obj.PhaseNoiseInputPort
                if obj.EnableInputPort
                    phnoise = pnoise;
                else
                    phnoise = enrx;
                end
                cond =  ~isa(phnoise,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','Ph','double');
                end
                sz_phnoise = size(phnoise);
                cond =  ~isequal(sz_phnoise, [sz_x(1) 1]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions', 'Ph', ...
                        coder.const(mat2str([sz_x(1) 1])), ...
                        coder.const(mat2str(sz_phnoise)));
                end
                cond =  ~isreal(phnoise);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:step:NeedReal', 'Ph');
                end
            end
        end
        
        function num = getNumInputsImpl(obj)
            num = 1;
            if obj.EnableInputPort
                num = num+1;
            end
            if obj.PhaseNoiseInputPort
                num = num+1;
            end
        end
    end
    
    methods (Access = protected)
        function varargout = getOutputSizeImpl(obj)
           varargout{1} = propagatedInputSize(obj,1);
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj, 1);
        end
        function varargout = getOutputDataTypeImpl(obj)
            varargout{1} = propagatedInputDataType(obj,1);
        end
        function varargout = isOutputComplexImpl(obj)  
            varargout{1} = ...
                obj.PhaseNoiseInputPort || ...
                propagatedInputComplexity(obj,1);
        end       
    end
    
end

