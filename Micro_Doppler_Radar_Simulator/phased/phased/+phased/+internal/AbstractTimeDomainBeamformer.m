classdef (Hidden) AbstractTimeDomainBeamformer < phased.internal.AbstractBeamformer & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%AbstractTimeDomainBeamformer   Abstract class for time domain beamformers

%   Copyright 2010-2017 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from Simulink time engine. Set SampleRateFromInputCheckbox to
        %   false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true
    end
    
    properties (Nontunable)
        %SampleRate     Sample rate (Hz)
        %   Specify the signal sampling rate (in Hz) as a positive scalar.
        %   The default value of this property is 1e6.
        SampleRate = 1e6;
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties(Access = protected, Nontunable)
        cElementDelay;
        pSampleRate
    end

    methods (Access = protected)

        function obj = AbstractTimeDomainBeamformer(varargin)
            obj@phased.internal.AbstractBeamformer(varargin{:});

        end

    end

    methods

        function set.SampleRate(obj,val)
            sigdatatypes.validateFrequency(val,'phased.internal',...
                'SampleRate',{'double','single'},{'scalar'});
            obj.SampleRate = val;
        end

    end
    
    methods (Access = protected)
        function x_steered = steer(obj,x,ang)
            delay = step(obj.cElementDelay,cast(ang,'double'));
            x_steered = delayseq(x,-delay,obj.pSampleRate);
        end
        
        function setupImpl(obj,x)
            setupImpl@phased.internal.AbstractBeamformer(obj,x);
            obj.cElementDelay = phased.ElementDelay(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',cast(obj.PropagationSpeed,'double'));
            %obj.pSampleRate = getSampleRate(obj,size(x,1),1,obj.SampleRate);
            fs = cast(obj.SampleRate,class(x)); % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            obj.pSampleRate = fs;
        end
        
        function releaseImpl(obj)
            release(obj.cElementDelay);
        end
        
        function resetImpl(obj)       
            reset(obj.cElementDelay);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractBeamformer(obj);
            if isLocked(obj)
                s.cElementDelay = saveobj(obj.cElementDelay);
                s.pSampleRate = obj.pSampleRate;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractBeamformer(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cElementDelay = phased.ElementDelay.loadobj(s.cElementDelay);
                    s = rmfield(s,'cElementDelay');
                    % recover locked sample rate information
                    if isfield(s,'pSampleRate')
                        obj.pSampleRate = s.pSampleRate;
                        s = rmfield(s,'pSampleRate');
                    else
                        obj.pSampleRate = s.SampleRate;
                    end
                end
            end
        end
    end

end


