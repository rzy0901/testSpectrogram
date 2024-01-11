classdef (Abstract, Hidden) AbstractMultipathChannel < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%   Copyright 2016-2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %OperatingFrequency Signal carrier frequency (Hz)
        %   Specify the carrier frequency of the signal as a positive
        %   scalar in Hz. The default value of this property is 20e3.
        OperatingFrequency = 20e3;  
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate as a positive scalar in Hz. The default
        %   value of this property is 1e3.
        SampleRate = 1e3; 
        %MaximumDelaySource  Source of maximum propagation time delay
        %   Specify how the maximum propagation time delay is defined as
        %   one of 'Auto' | 'Property', where the default is 'Auto'. When
        %   you set this property to 'Auto', the object automatically
        %   allocates the memory to simulate the propagation delay. When
        %   you set this property to 'Property', the maximum propagation
        %   delay time is specified via the MaximumDelay property and
        %   any samples which propagate beyond this time are ignored.
        MaximumDelaySource = 'Auto'
        %MaximumDelay  Maximum propagation time delay (s)
        %   Specify the maximum propagation delay in seconds as a
        %   positive scalar. This property applies when you set the
        %   MaximumDelaySource property to 'Property'. The default
        %   value of this property is 1.
        MaximumDelay = 1
        %InterpolationMethod  Interpolation Method
        %   Specify the method used by the channel to implement signal
        %   fractional delay and doppler time dilation/compression. This
        %   property can be set to one of [{'Linear'} | 'Oversample']. When
        %   this property is set to 'Linear', the input signal is linearly
        %   interpolated directly onto a uniform grid to propagate the
        %   signal. When this property is set to 'Oversample', the input
        %   signal is resampled to a higher rate before linear
        %   interpolation to preserve spectral shape. The default value of
        %   this property is 'Linear'.
        InterpolationMethod = 'Linear'
    end
    
    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from Simulink time engine. Set SampleRateFromInputCheckbox to
        %   false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true;
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties (Access = protected, Nontunable)
        %Sample rate, in MATLAB, specified by property but in Simulink,
        %specified by engine
        pSampleRate
        % Buffer
        cBuffer
    end
    
    properties (Access = protected)
        % logical indices of paths to be used
        pUsedPaths
    end
    
    properties(Constant, Hidden)
        MaximumDelaySourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        InterpolationMethodSet = matlab.system.StringSet({'Linear', ...
      'Oversample'});
    end
    
    methods
        function set.OperatingFrequency(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'positive'},'MultipathChannel','OperatingFrequency');
            obj.OperatingFrequency = value;
        end
        function set.SampleRate(obj,value)
            validateattributes(value,{'double'},{'scalar','positive',...
                'finite'},'MultipathChannel','SampleRate');
            obj.SampleRate = value;
        end
        function set.MaximumDelay(obj,value)
            validateattributes(value,{'double'},{'scalar','positive',...
                'finite'},'MultipathChannel','MaximumDelay');
            obj.MaximumDelay = value;
        end
    end
    
    methods (Access = protected)
        function obj = AbstractMultipathChannel(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access = protected)

        function flag = isInactivePropertyImpl(obj, prop)
            if (obj.MaximumDelaySource(1) == 'A') && ...
                    strcmp(prop,'MaximumDelay')
                flag = true;
            else
                flag = false;
            end
        end
        
        function validateInputSignal(obj,x) 
            validateNumChannels(obj,x);
        end
        
        function validateInputsImpl(obj,x,paths,dop,aloss)
            coder.extrinsic('mat2str');
            coder.extrinsic('num2str');            
            validateInputSignal(obj,x);
            
            % check signal
            cond =  ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','X','double');
            end
            
            % check paths is a real doubles and
            % non-nan.
            cond =  ~isa(paths,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','PATHS','double');              
            end
            
            cond = ~isreal(paths);
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:NeedReal','PATHS');
            end
            
            % Check dop is a real double 
            cond = ~isa(dop,'double');
            if cond
                coder.internal.errorIf(cond,...
                    'MATLAB:system:invalidInputDataType','DOP','double');
            end
           
            cond = ~isreal(dop);
            if cond
                coder.internal.errorIf(cond,...  
                    'phased:step:NeedReal','DOP');
            end
            
            % Check aloss is a double matrix 
            cond =  ~isa(aloss,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','ALOSS','double');              
            end
            
            cond = ~isreal(aloss);
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:NeedReal','ALOSS');
            end
            
            % Check the size of paths
            pathsSize = size(paths);
            xSize = size(x);
            cond = ~isequal(xSize(2),pathsSize(2)) || ~isequal(pathsSize(1),3);
            if cond
                coder.internal.errorIf(cond,...
                    'MATLAB:system:invalidInputDimensions','PATHS',...
                    coder.const(mat2str([3 xSize(2)])),coder.const(mat2str(pathsSize)));
            end

            % Check the size of aloss
            alossSize = size(aloss);
            xSize = size(x);
            cond = ~isequal(xSize(2),alossSize(2)-1);
            if cond
                coder.internal.errorIf(cond,...
                    'MATLAB:system:invalidInputDimensions','ALOSS',...
                    coder.const(['[M ' num2str(xSize(2)+1) ']']),coder.const(mat2str(alossSize)));
            end
            
            % Check the size of dop
            dopSize = size(dop);
            xSize = size(x);
            cond = ~isequal(xSize(2),dopSize(2));
            if cond
                coder.internal.errorIf(cond,...
                  'MATLAB:system:invalidInputDimensions','DOP',...
                  coder.const(mat2str([1 xSize(2)])),coder.const(mat2str(dopSize)));
            end
        end
        
        function setupImpl(obj,x,~,~,~,~)
            sz_x = size(x,1);  
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            obj.pSampleRate = getSampleRate(obj,sz_x,1,obj.SampleRate);
            fs = obj.SampleRate; % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end

            if strcmp(obj.MaximumDelaySource,'Auto')
                obj.cBuffer = phased.internal.DopplerCircularBuffer(...
                    'BufferLength',1);
            else
                % Set the fixed buffer length equal to the maximum delay,
                % in samples. 
                buflen = ceil(...
                    obj.MaximumDelay*obj.pSampleRate);
                obj.cBuffer = phased.internal.DopplerCircularBuffer(...
                        'FixedLengthBuffer',true,'BufferLength',buflen,...
                        'BufferWidthSource','Auto');  
            end
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
            release(obj.cBuffer);
        end
        
        function resetImpl(obj)
            reset(obj.cBuffer);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            if isLocked(obj)
                s.pSampleRate = obj.pSampleRate;
                s.pUsedPaths = obj.pUsedPaths;
                s.cBuffer = saveobj(obj.cBuffer);
            end
        end

        function s = loadSubObjects(obj,s,wasLocked)
            if isfield(s,'isLocked')                                        
                if s.isLocked  
                    obj.cBuffer = phased.internal.CircularBuffer.loadobj(s.cBuffer); 
                    s = rmfield(s,'cBuffer');
                    % recover locked sample rate information
                    obj.pSampleRate = s.SampleRate;
                end
                s = rmfield(s,'isLocked');                                      
            elseif wasLocked
                obj.cBuffer = phased.internal.CircularBuffer.loadobj(s.cBuffer);
                s = rmfield(s,'cBuffer');
                % recover locked sample rate information
                if isfield(s,'pSampleRate')
                    obj.pSampleRate = s.pSampleRate;
                    s = rmfield(s,'pSampleRate');
                else
                    obj.pSampleRate = s.SampleRate;
                end
            end
        end

        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)  %#ok<INUSL>
            if index == 1
                flag = false;
            else % (index == 2,3)
                flag = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end

        function num = getNumInputsImpl(obj)   %#ok<MANU>
            num = 5;
        end
        
        function num = getNumOutputsImpl(obj)   %#ok<MANU>
            num = 1;
        end

    end
    methods (Hidden, Static)
        function flag = isAllowedInSystemBlock(obj) %#ok<INUSD>
            flag = false;
        end
    end
end
