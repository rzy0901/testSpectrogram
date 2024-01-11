classdef (Sealed,StrictDefaults) ReceiverPreamp < phased.internal.AbstractReceiverPreamp 
%ReceiverPreamp   Receiver preamp
%   H = phased.ReceiverPreamp creates a receiver preamp System object, H,
%   that receives the incoming pulses.
%
%   H = phased.ReceiverPreamp(Name,Value) creates a receiver preamp object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) applies the receiver preamp gain and noise to the input
%   signal, X, and returns the resulting output signal, Y. Y has the same
%   dimensions as X.
%
%   Y = step(H,X,TR) uses input TR as the enabling signal when the
%   EnableInputPort property is set to true. TR is a column vector whose
%   number of elements is the same as the number of rows in X. TR can be
%   either of type logical or of type double.
%
%   If TR is of type double, every element of TR that is equal to 0
%   indicates that the preamp is turned off, and no input signal passes
%   through the receiver. Every element of TR that is not equal to zero
%   indicates that the preamp is turned on, and the input passes through.
%   If TR is logical, the preamp is enabled whenever an element of TR is
%   true and disabled whenever an element of TR is false.
%
%   Y = step(H,X,PH) uses input PH as the phase noise for each sample in X
%   when the PhaseNoiseInputPort is set to true. PH is a column vector
%   whose number of elements is the same as the number of rows in X.
%
%   The phase noise is the same for all channels in X. The elements in
%   PH represent the random phases, in radians, added by transmitter
%   to the transmitted pulses. The receiver preamp object removes these
%   random phases from all received samples returned within corresponding
%   pulse intervals. Such setup is often referred to as coherent on
%   receive.
%
%   You can combine optional input arguments when their enabling
%   properties are set. Optional inputs must be listed in the same order
%   as the order of the enabling properties. For example,
%
%   Y = step(H,X,TR,PH)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ReceiverPreamp methods:
%
%   step     - Receive the incoming signal
%   release  - Allow property value and input characteristics changes
%   clone    - Create a pulse receiver object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset the random number generator for noise generation
%
%   ReceiverPreamp properties:
%
%   Gain                 - Gain 
%   LossFactor           - Loss factor
%   NoiseMethod          - Noise specification method
%   NoiseFigure          - Noise figure
%   ReferenceTemperature - Reference temperature
%   SampleRate           - Sample rate
%   NoisePower           - Noise power
%   NoiseComplexity      - Noise complexity
%   EnableInputPort      - Enable enabling signal input
%   PhaseNoiseInputPort  - Enable phase noise input
%   SeedSource           - Source of seed for random number generator
%   Seed                 - Seed for random number generator
%
%   % Examples:
%
%   % Example 1:
%   %   Simulate the reception of a sine wave.
%
%   Fs = 100; t = linspace(0,1-1/Fs,100); x = 1e-6*sin(2*pi*5*t);
%   rx = phased.ReceiverPreamp('NoiseFigure',10,'SampleRate',100);
%   y = rx(x);
%   plot(t,x,t,real(y)), xlabel('Time (s)'), ylabel('Amplitude');
%   legend('Original signal','Received signal');
%
%   % Example 2:
%   %   Simulate the reception of a sine wave with real noise. The SNR is
%   %   3 dB.
%
%   Fs = 100; t = linspace(0,1-1/Fs,100); x = sin(2*pi*5*t);
%   snrdb = 3; npow = 0.5/db2pow(snrdb);
%   rx = phased.ReceiverPreamp('NoiseMethod','Noise power',...
%           'NoisePower',npow,'NoiseComplexity','Real');
%   y = rx(x);
%   plot(t,x,t,real(y)), xlabel('Time (s)'), ylabel('Amplitude');
%   legend('Original signal','Received signal');
%
%   See also phased, phased.Transmitter, phased.Collector.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [2] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods

        function obj = ReceiverPreamp(varargin)
            %ReceiverPreamp   Construct the ReceiverPreamp class.
            obj@phased.internal.AbstractReceiverPreamp(varargin{:});
        end
        
    end
    
    methods(Access = protected)
        
        function flag = isInactivePropertyImpl(obj, prop)
            if (obj.SeedSource(1) == 'A') && ...  %Auto
                    strcmp(prop, 'Seed')
                flag = true;
            elseif strcmp(obj.NoiseMethod,'Noise power') && ...
                    (strcmp(prop, 'ReferenceTemperature') || ...
                    strcmp(prop, 'NoiseFigure') || ...
                    strcmp(prop, 'SampleRate') || ...
                    strcmp(prop, 'PhaseNoiseInputPort') || ...
                    strcmp(prop, 'EnableInputPort'))
                flag = true;
            elseif strcmp(obj.NoiseMethod,'Noise temperature') && ...
                    strcmp(prop, 'NoisePower')
                flag = true;
            else
                flag = false;
            end
        end
               
        function loadObjectImpl(obj,s,~)
            ws = warning('off','phased:system:System:NoiseBandwidthWarning');
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
            warning(ws);
        end
        
        function Fs = computeSampleRate(obj)
            Fs = obj.SampleRate;
        end
        
    end
    
    methods (Access = protected, Static, Hidden)
        function groups = getPropertyGroupsImpl
            groups = matlab.system.display.Section(...
                'phased.ReceiverPreamp');
            % dSampleRate = matlab.system.display.internal.Property(...
            %     'SampleRate','IsObjectDisplayOnly',true);
            dNoiseComplexity = matlab.system.display.internal.Property(...
                'NoiseComplexity','IsObjectDisplayOnly',true);
            dSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false,'UseClassDefault',false,...
                'Default','Property');
            dSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false,'UseClassDefault',false,...
                'Default','randi(65535,1)');
            for m = 1:numel(groups.PropertyList)
                if strcmp(groups.PropertyList{m},'NoiseComplexity')
                    groups.PropertyList{m} = dNoiseComplexity;
                % elseif strcmp(groups.PropertyList{m},'SampleRate')
                %     groups.PropertyList{m} = dSampleRate;
                elseif strcmp(groups.PropertyList{m},'SeedSource')
                    groups.PropertyList{m} = dSeedSource;
                elseif strcmp(groups.PropertyList{m},'Seed')
                    groups.PropertyList{m} = dSeed;
                end
            end
        end
    end
        
    methods (Static,Hidden)
        function a = getAlternateBlock
            a = 'phasedtxrxlib/Receiver Preamp';
        end
    end
end

