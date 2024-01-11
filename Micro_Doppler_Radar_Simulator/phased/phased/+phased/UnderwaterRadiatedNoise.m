classdef (Sealed,StrictDefaults) UnderwaterRadiatedNoise < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.internal.SampleTime
%UnderwaterRadiatedNoise  Underwater radiated noise
%   H = phased.UnderwaterRadiatedNoise creates an acoustic shipping noise
%   System object, H. This object generates a noise signal produced by a
%   ship or subsurface vessel.
%
%   H = phased.UnderwaterRadiatedNoise(Name,Value) returns an shipping
%   noise object, H, with the specified property Name set to the specified
%   Value. You can specify additional name-value pair arguments in any
%   order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,ANG) returns samples of a noise waveform in the
%   columns of vector Y. 
%
%   ANG is a 2-row matrix representing the signal's transmission direction.
%   Each column of ANG specifies the direction of the corresponding noise
%   signal in the form of an [AzimuthAngle; ElevationAngle] pair (in
%   degrees). When ANG represents multiple signals, the DirectionalPattern
%   property contains either one pattern or M patterns where M is the
%   number of columns in ANG. If there is only one pattern, then the
%   multiple noise signals originate from the same source. If there are M
%   patterns, then multiple distinct signals are generated from each
%   corresponding pattern.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   UnderwaterRadiatedNoise methods:
%
%   step     - Generate noise signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create underwater channel object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the object
%
%   UnderwaterRadiatedNoise properties:
%
%   NumSamples               - Number of samples
%   SampleRate               - Sample rate
%   OperatingFrequency       - Operating frequency
%   TonalFrequencies         - Tonal frequencies
%   TonalLevels              - Tonal source levels
%   BroadbandLevels          - Spectrum levels of broadband noise
%   AzimuthAngles            - Azimuth angles
%   ElevationAngles          - Elevation angles
%   DirectionalPattern       - Directional pattern
%   FrequencyVector          - Directional pattern frequency vector
%   SeedSource               - Source of seed for random number generator
%   Seed                     - Seed for random number generator
%
%   % Example:
%   %   Generate and plot a noise signal transmitted at an angle  
%   %   of 10 degrees azimuth. 
%   noise = phased.UnderwaterRadiatedNoise('NumSamples',1024,...
%     'BroadbandLevels',[90 90 90 150 120 100 90],...
%     'SeedSource','Property','Seed',181);
%   y = noise([10; 0]);
%   pwelch(y,[],[],[],noise.SampleRate,'psd','centered');
%
%   See also phased.IsoSpeedUnderwaterPaths, phased.MultipathChannel, 
%   phased.BackscatterSonarTarget

%   Copyright 2017 The MathWorks, Inc.

%   Reference
%   [1] Urick, Principles of Underwater Sound, Peninsula Publishing, 1983

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods (Access = public)
        function obj = UnderwaterRadiatedNoise(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    properties (Nontunable, PositiveInteger) 
        %NumSamples   Number of samples in output
        %   Specify the number of samples in each output as a positive
        %   integer. The default value of this property is 1e3.
        NumSamples = 100;
    end

    properties (Nontunable)
        %SampleRate   Sample rate (Hz)
        %   Specify the sample rate as a positive scalar in Hz. The default
        %   value of this property is 1e3.
        SampleRate = 1e3;
        %OperatingFrequency Signal carrier frequency (Hz)
        %   Specify the carrier frequency of the signal as a positive
        %   scalar in Hz. The default value of this property is 20e3.
        OperatingFrequency = 20e3;  
        %TonalFrequencies    Tonal frequencies (Hz)
        %   Specify line series frequencies as a positive strictly
        %   increasing vector in Hz. The default value of this property is
        %   [19700 20100 20300].
        TonalFrequencies = [19700 20100 20300];
        %TonalLevels    Tonal source levels (dB re 1 uPa)
        %   Specify the source levels for line series as a positive vector
        %   in dB. The default value of this property is [150 150 150].
        TonalLevels = [150 150 150];
        %BroadbandLevel Broadband noise spectrum level (dB re 1 uPa/ 1 Hz)
        %   Specify the spectrum level for broadband noise as a vector in
        %   dB. The noise spectrum is shaped according to these values at
        %   uniformly spaced frequencies in the band [-SampleRate/2
        %   SampleRate/2] + OperatingFrequency. The default value is 130
        %   dB.
        BroadbandLevels = 130;
        %AzimuthAngles    Azimuth angles (deg)
        %   Specify the azimuth angles (in degrees) as a length P vector.
        %   These are the azimuth angles where the custom pattern is
        %   evaluated. P must be greater than 2. The default value of this
        %   property is -180:180.
        AzimuthAngles = -180:180;
        %ElevationAngles   Elevation angles (deg)
        %   Specify the elevation angles (in degrees) as a length Q vector.
        %   These are the elevation angles where the custom pattern is
        %   evaluated. The default value of this property is -90:90.
        ElevationAngles = -90:90;
        %FrequencyVector    Directional pattern frequency vector (Hz)
        %   Specify the frequencies (in Hz) where the frequency responses
        %   of element are measured as a row vector. The elements of the
        %   vector must be increasing. The default value of this property
        %   is [0 100e6].
        FrequencyVector = [0 100e6];
        %DirectionalPattern   Directional pattern (dB)
        %   Specify the directional pattern for the noise source in dB.
        %
        %   For a given source, the pattern is specified as a QxPxK array,
        %   where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and K is the number of
        %   frequencies specified in the FrequencyVector property. Each
        %   page of the array defines a directional pattern across the
        %   region defined by the azimuth and elevation angles at a given
        %   frequency. If K is 1, the same pattern is used across all
        %   frequencies. The default value of this property is a 181x361
        %   matrix with all elements equal to 0.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xPxK
        %   array or a KxP matrix. In this case, each row is a pattern for
        %   a given frequency.
        %
        %   If you want to use the object to represent L sources in the
        %   simulation, put L patterns in a cell array and specify it to
        %   the property. All patterns should use the same format.
        DirectionalPattern = zeros(181,361);
        %SeedSource   Source of seed for random number generator
        %   Specify how the random numbers are generated as one of 'Auto' |
        %   'Property', where the default is 'Auto'. The random numbers are
        %   used to model random TS values. When you set this property to
        %   'Auto', the random numbers are generated using the default
        %   MATLAB random number generator. When you set this property to
        %   'Property', a private random number generator is used with a
        %   seed specified by the value of the Seed property. 
        %
        %   To use this object with Parallel Computing Toolbox software,
        %   set this property to 'Auto'.
        SeedSource = 'Auto'
        %Seed Seed for random number generator
        %   Specify the seed for the random number generator as a
        %   non-negative integer. The integer must be less than 2^32. This
        %   property applies when the SeedSource property is 'Property'.
        %   The default value is 0.
        Seed = 0
    end

    properties(Access = private, Nontunable)
        % Private random stream
        cNoiseSource;
        cNoiseSourcePhase;
        cSubbandDivider;
        cSubbandCombiner;
        pNumSubands
        pSubbandFreqs
        pSampleRate
        pBPFilterLength
        pNumTonals
        pNumSources
        pFilt
    end
    
    properties(Access = private)
        pTime = 0;
        pPhase
    end
    
    properties(Constant, Hidden)
        SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    methods
        function set.OperatingFrequency(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'positive'},'UnderwaterRadiatedNoise','OperatingFrequency');
            obj.OperatingFrequency = value;
        end
        function set.SampleRate(obj,value)
            validateattributes(value,{'double'},{'scalar','positive',...
                'finite'},'UnderwaterRadiatedNoise','SampleRate');
            obj.SampleRate = value;
        end
        function set.TonalFrequencies(obj,value)
            validateattributes(value,{'double'},{'vector','positive',...
                'finite'},'UnderwaterRadiatedNoise','TonalFrequencies');
            obj.TonalFrequencies = value;
        end
        function set.TonalLevels(obj,value)
            validateattributes(value,{'double'},{'vector','real'...
                'finite'},'UnderwaterRadiatedNoise','TonalLevels');
            obj.TonalLevels = value;
        end
        function set.BroadbandLevels(obj,value)
            validateattributes(value,{'double'},{'vector','real'...
                'finite'},'UnderwaterRadiatedNoise','BroadbandLevels');
            obj.BroadbandLevels = value;
        end
        function set.AzimuthAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.UnderwaterRadiatedNoise','AzimuthAngles',...
                {'vector','>=',-180,'<=',180});
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'AzimuthAngles');
            end
            obj.AzimuthAngles = value;
        end
        function set.ElevationAngles(obj,value)
            % allow single elevation cut, i.e., azimuth only pattern
            sigdatatypes.validateAngle(value,...
                'phased.UnderwaterRadiatedNoise','ElevationAngles',...
                {'vector','>=',-90,'<=',90});
            obj.ElevationAngles = value;
        end
        function set.FrequencyVector(obj,value)
            validateattributes( value, { 'double' }, ...
                {'nonempty','finite','row','nonnegative','nondecreasing'}, '', 'FrequencyVector');
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'FrequencyVector');
            end
            obj.FrequencyVector = value;
        end
        function set.DirectionalPattern(obj,value)
            validateattributes(value,{'double','cell'},{'nonempty'},...
                'phased.UnderwaterRadiatedNoise','DirectionalPattern');
            if iscell(value)
                for m = 1:numel(value)
                    validateattributes(value{m},{'double'},...
                        {'real','finite','nonempty','3d'},...
                        'phased.UnderwaterRadiatedNoise',...
                        feval( 'sprintf' ,'DirectionalPattern{%d}',m));
                end
            else
                validateattributes(value,{'double'},...
                    {'real','finite','nonempty','3d'},...
                    'phased.UnderwaterRadiatedNoise',...
                    'DirectionalPattern');
            end
            obj.DirectionalPattern = value;
        end
        function set.Seed(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative',...
                'integer','<',2^32},'phased.UnderwaterRadiatedNoise','Seed');
            obj.Seed = value;
        end
    end

    methods (Access = protected)
        function flag = isInactivePropertyImpl(obj, prop)
            if (obj.SeedSource(1) == 'A') && ...
                    strcmp(prop,'Seed')
                flag = true;
            else
                flag = false;
            end
        end
        
        function validatePropertiesImpl(obj)
            % Check that tonal frequencies and levels are row vectors of
            % the same size.
            validateattributes(obj.TonalFrequencies,{'double'},...
              {'nrows',1},'phased.UnderwaterRadiatedNoise','TonalFrequencies');
            validateattributes(obj.TonalLevels,{'double'},...
              {'nrows',1,'ncols',numel(obj.TonalFrequencies)},...
              'phased.UnderwaterRadiatedNoise','TonalLevels');
            
            % Validate DirectionalPattern
            pattern_size = [numel(obj.ElevationAngles)...
              numel(obj.AzimuthAngles)];
            num_freq = numel(obj.FrequencyVector);
            
            if numel(obj.ElevationAngles) == 1
              % azimuth pattern only
              if iscell(obj.DirectionalPattern)
                  % multiple target
                  if numel(size(obj.DirectionalPattern{1}))==3
                      for m = 1:numel(obj.DirectionalPattern)
                          % multiple 3D azimuth wideband pattern
                          validateattributes(obj.DirectionalPattern{m},{'double'},...
                              {'size',[pattern_size num_freq]},...
                              'phased.UnderwaterRadiatedNoise',...
                              feval( 'sprintf' ,'DirectionalPattern(%d)',m));
                      end
                  elseif size(obj.DirectionalPattern{1},1)==1
                      for m = 1:numel(obj.DirectionalPattern)
                          % multiple 3D azimuth wideband pattern,
                          % share across frequencies
                          validateattributes(obj.DirectionalPattern{m},{'double'},...
                              {'size',[1 pattern_size(2)]},...
                              'phased.UnderwaterRadiatedNoise',...
                              feval( 'sprintf' ,'DirectionalPattern(%d)',m));
                      end
                  else
                      for m = 1:numel(obj.DirectionalPattern)
                          % multiple azimuth wideband pattern
                          validateattributes(obj.DirectionalPattern{m},{'double'},...
                              {'size',[num_freq pattern_size(2)]},...
                              'phased.UnderwaterRadiatedNoise',...
                              feval( 'sprintf' ,'DirectionalPattern(%d)',m));
                      end
                  end
              else
                  % single target
                  if numel(size(obj.DirectionalPattern))==3
                      % multiple 3D azimuth wideband pattern
                      validateattributes(obj.DirectionalPattern,{'double'},...
                          {'size',[pattern_size num_freq]},...
                          'phased.UnderwaterRadiatedNoise',...
                          'DirectionalPattern');
                  elseif size(obj.DirectionalPattern,1)==1
                      % multiple 3D azimuth wideband pattern,
                      % share across frequencies
                      validateattributes(obj.DirectionalPattern,{'double'},...
                          {'size',[1 pattern_size(2)]},...
                          'phased.UnderwaterRadiatedNoise',...
                          'DirectionalPattern(%d)');
                  else
                      % multiple azimuth wideband pattern
                      validateattributes(obj.DirectionalPattern,{'double'},...
                          {'size',[num_freq pattern_size(2)]},...
                          'phased.UnderwaterRadiatedNoise',...
                          'DirectionalPattern');
                  end
              end
            else
              % 3D pattern
              if iscell(obj.DirectionalPattern)
                  if numel(size(obj.DirectionalPattern{1}))==3
                      % multiple target
                      for m = 1:numel(obj.DirectionalPattern)
                          % multiple 3D wideband pattern
                          validateattributes(obj.DirectionalPattern{m},{'double'},...
                              {'size',[pattern_size num_freq]},...
                              'phased.UnderwaterRadiatedNoise',...
                              feval( 'sprintf' ,'DirectionalPattern(%d)',m));
                      end
                  else
                       % multiple target, share pattern across
                       % frequencies
                      for m = 1:numel(obj.DirectionalPattern)
                          % multiple 3D wideband pattern
                          validateattributes(obj.DirectionalPattern{m},{'double'},...
                              {'size',pattern_size},...
                              'phased.UnderwaterRadiatedNoise',...
                              feval( 'sprintf' ,'DirectionalPattern(%d)',m));
                      end
                 end
              else
                  if numel(size(obj.DirectionalPattern))==3
                      % single 3D wideband pattern
                      validateattributes(obj.DirectionalPattern,{'double'},...
                          {'size',[pattern_size num_freq]},...
                          'phased.UnderwaterRadiatedNoise',...
                          'DirectionalPattern');
                  else
                      % single 3D wideband pattern across frequencies
                      validateattributes(obj.DirectionalPattern,{'double'},...
                          {'size',pattern_size},...
                          'phased.UnderwaterRadiatedNoise',...
                          'DirectionalPattern');
                  end
              end
            end
        end
        
        function validateInputsImpl(obj,angle)
            sz_angle = size(angle);
            cond =  ~ismatrix(angle) || isempty(angle);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeMatrix','ANGLE');
            end

            cond =  sz_angle(1) > 2;
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:system:element:NeedTwoRows');
            end

            cond =  ~isreal(angle);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:system:element:InvalidAngle');
            end

            cond =  ~isa(angle,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','ANGLE','double');
            end
            
            if iscell(obj.DirectionalPattern)
                npat = numel(obj.DirectionalPattern);
            else
                npat = 1;
            end

            cond = (npat ~= 1) && (sz_angle(2) ~= npat);
            if cond
                coder.internal.errorIf(cond,'phased:UnderwaterRadiatedNoise:anglePatternMismatch');
            end
            
            sigdatatypes.validateAzElAngle(angle,'','Ang');
        end
        
        function setupImpl(obj,~)
            obj.pSampleRate = obj.SampleRate;
            obj.pNumTonals = length(obj.TonalLevels);
            obj.pTime = (0:obj.NumSamples-1)'/obj.SampleRate;
            if iscell(obj.DirectionalPattern)
                obj.pNumSources = numel(obj.DirectionalPattern);
            else
                obj.pNumSources = 1;
            end
            obj.cNoiseSource = phased.internal.NoiseSource(...
                    'SeedSource',obj.SeedSource,'Distribution','Gaussian');
            if (obj.SeedSource(1) == 'P') %Property
                obj.cNoiseSource.Seed = obj.Seed;
            end
            obj.cNoiseSourcePhase = phased.internal.NoiseSource(...
                    'SeedSource',obj.SeedSource);
            if (obj.SeedSource(1) == 'P') %Property
                obj.cNoiseSourcePhase.Seed = obj.Seed;
            end
            obj.pNumSubands = 64;
            obj.cSubbandDivider = phased.internal.SubbandDivider(...
                'OperatingFrequency',obj.OperatingFrequency,...
                'SampleRate',obj.pSampleRate,...
                'NumSubbands',obj.pNumSubands,'EnableWarning',false);
            obj.cSubbandCombiner = phased.internal.SubbandCombiner(...
                'NumSubbands',obj.pNumSubands,'TimeSignalLengthSource','Inherit',...
                'EnableWarning',false);
            obj.pSubbandFreqs = phased.internal.subbandCenterFrequency(...
                obj.OperatingFrequency,obj.pSampleRate,obj.pNumSubands).';
            obj.pBPFilterLength = 200;
            obj.pPhase = step(obj.cNoiseSourcePhase,1,[obj.pNumSources obj.pNumTonals])*2*pi;
            
            cond = obj.pSampleRate >= 2*obj.OperatingFrequency;
            if cond
                coder.internal.errorIf(cond,'phased:SubbandBeamformer:SampleRateMismatch');
            end
            if numel(obj.BroadbandLevels) > 1
                obj.pFilt = fir2(obj.pBPFilterLength,linspace(0,1,length(obj.BroadbandLevels)),db2mag(obj.BroadbandLevels)/db2mag(max(obj.BroadbandLevels)));
            end
              
        end
        
        function yout = stepImpl(obj,ang)
            % Generate white gaussian phase and noise signal
            yTonals = zeros(obj.NumSamples,obj.pNumSources,'like',1i);
            noiseBroadband = step(obj.cNoiseSource,1,[2*size(obj.pTime,1) obj.pNumSources]);
            maxLevel = max(obj.BroadbandLevels);
            if numel(obj.BroadbandLevels) > 1
                noiseBroadband = filter(obj.pFilt,1,noiseBroadband);
            end
            
            % Generate a signal with baseband spectrum (-fs/2 to fs/2)
            % equal to the specified shape of the frequency sampling
            % bandpass filter. To do this, we demodulate and decimate the
            % signal, maintaining power. The original white random noise
            % signal has twice the length (sample rate) of the output
            % signal.
            noiseBroadband = hilbert(noiseBroadband)/sqrt(2);
            noiseBroadband = bsxfun(@times,noiseBroadband,exp(-1*1i*2*pi*(0:length(noiseBroadband)-1)'/4));
            yBroadband = noiseBroadband(1:2:end,:);
            
            % Scale the signal based on the requested broadband levels.
            yBroadband = db2mag(maxLevel + 10*log10(obj.SampleRate/2))*yBroadband/1e6;
            
            % Generate tonal components
            TonalFrequenciesIdx = (obj.TonalFrequencies > (obj.OperatingFrequency-obj.SampleRate/2)) & ...
              (obj.TonalFrequencies < (obj.OperatingFrequency+obj.SampleRate/2));
            if any(TonalFrequenciesIdx)
                for i = 1:obj.pNumSources
                    argw = 2*pi*bsxfun(@times,obj.TonalFrequencies(TonalFrequenciesIdx)-obj.OperatingFrequency,obj.pTime);
                    arg = bsxfun(@plus,argw,obj.pPhase(i,TonalFrequenciesIdx));
                    yTon = 1/2*exp(1i*arg); 
                    yTonals(:,i) = sqrt(2)*sum(bsxfun(@times,db2mag(obj.TonalLevels(TonalFrequenciesIdx)),yTon),2)/1e6; 
                end
            end
            
            y = yTonals+yBroadband; 
    
            % Scale each channel based on directional factor
            if iscell(obj.DirectionalPattern)
                yout = zeros(obj.NumSamples,obj.pNumSources,'like',1i);
                for i = 1:numel(obj.DirectionalPattern)
                  yout(:,i) = applyDirectionalPattern(obj,y(:,i),obj.DirectionalPattern{i},ang(:,i)); 
                end
            else
                yout = applyDirectionalPattern(obj,y,obj.DirectionalPattern,ang);  
            end
            
            % Update internal time vector
            obj.pTime = obj.pTime + (obj.pTime(end)-obj.pTime(1))+1/obj.SampleRate;
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj); 
            release(obj.cNoiseSource);
            release(obj.cNoiseSourcePhase);
            release(obj.cSubbandDivider);
            release(obj.cSubbandCombiner);
        end
        
        function resetImpl(obj)
            reset(obj.cNoiseSource);
            reset(obj.cNoiseSourcePhase);
            reset(obj.cSubbandDivider);
            reset(obj.cSubbandCombiner);
            obj.pTime = (0:obj.NumSamples-1)'/obj.SampleRate;
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            if isLocked(obj)
                s.pNumSubands = obj.pNumSubands;
                s.pSubbandFreqs = obj.pSubbandFreqs;
                s.pSampleRate = obj.pSampleRate;
                s.pPhase = obj.pPhase;
                s.pBPFilterLength = obj.pBPFilterLength;
                s.pTime = obj.pTime;
                s.pNumTonals = obj.pNumTonals;
                s.pNumSources = obj.pNumSources;
                s.pFilt = obj.pFilt;
                s.cNoiseSource = saveobj(obj.cNoiseSource);
                s.cNoiseSourcePhase = saveobj(obj.cNoiseSourcePhase);
                s.cSubbandDivider = saveobj(obj.cSubbandDivider);
                s.cSubbandCombiner = saveobj(obj.cSubbandCombiner);
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked 
                obj.cNoiseSource = ...
                    phased.internal.NoiseSource.loadobj(s.cNoiseSource);
                s = rmfield(s,'cNoiseSource');
                
                obj.cNoiseSourcePhase = ...
                    phased.internal.NoiseSource.loadobj(s.cNoiseSourcePhase);
                s = rmfield(s,'cNoiseSourcePhase');
      
                obj.cSubbandDivider = ...
                  phased.internal.SubbandDivider.loadobj(s.cSubbandDivider);
                s = rmfield(s,'cSubbandDivider');
                
                obj.cSubbandCombiner = ...
                  phased.internal.SubbandCombiner.loadobj(s.cSubbandCombiner);
                s = rmfield(s,'cSubbandCombiner');
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function num = getNumInputsImpl(obj)   %#ok<MANU>
            num = 1;
        end
        
        function num = getNumOutputsImpl(obj)   %#ok<MANU>
            num = 1;
        end
        
        function flag = isInputComplexityLockedImpl(obj,~) %#ok<INUSD>
            flag = true;
        end

        function flag = isOutputComplexityLockedImpl(obj,~) %#ok<INUSD>
            flag = true;  
        end
    end
    
    methods (Access = protected)
        function yout = applyDirectionalPattern(obj,y,DirectionalPattern,ang)
            az = obj.AzimuthAngles;
            el = obj.ElevationAngles;
            isSingleEl = isequal(numel(el),1); 
            if size(DirectionalPattern,3) > 1 
                g = zeros(1,size(ang,2),obj.pNumSubands);
                ysubband = step(obj.cSubbandDivider,y);
                for m = 1:size(ang,2)
                    for n = 1:obj.pNumSubands
                      % Locate the element of the frequency vector closest to
                      % subband center frequency. Interpolate the pattern
                      % corresponding to that frequency.
                      [~,nidx] = min(abs(obj.pSubbandFreqs(n)-obj.FrequencyVector));
                      g(1,m,n) = interpolatePattern(...
                                        az,el,db2mag(DirectionalPattern(:,:,nidx)),ang(1,m),ang(2,m),isSingleEl);
                    end
                end
                ysubband = bsxfun(@times, ysubband, g);
                yout = step(obj.cSubbandCombiner,ysubband,complex(y));
            elseif isSingleEl && size(DirectionalPattern,1)>1
                % One elevation angle and multiple patterns for frequency
                g = zeros(1,size(ang,2),obj.pNumSubands);
                ysubband = step(obj.cSubbandDivider,y);
                for m = 1:size(ang,2)
                    for n = 1:obj.pNumSubands
                      % Locate the element of the frequency vector closest to
                      % subband center frequency. Interpolate the pattern
                      % corresponding to that frequency.
                      [~,nidx] = min(abs(obj.pSubbandFreqs(n)-obj.FrequencyVector));
                      g(1,m,n) = interpolatePattern(...
                                        az,el,db2mag(DirectionalPattern(nidx,:)),ang(1,m),ang(2,m),isSingleEl);
                    end
                end
                ysubband = bsxfun(@times, ysubband, g);
                yout = step(obj.cSubbandCombiner,ysubband,complex(y));
            else
                if ~coder.target('MATLAB')
                    ysubband = step(obj.cSubbandDivider,1); %codegen
                    step(obj.cSubbandCombiner,ysubband,complex(y)); %codegen
                end
                g = zeros(1,size(ang,2));
                for m = 1:size(ang,2)
                    g(1,m) = interpolatePattern(...
                                        az,el,db2mag(DirectionalPattern),ang(1,m),ang(2,m),isSingleEl);
                end
                yout = bsxfun(@times,y,g);
            end
        end
    end

    methods (Hidden, Static)
        function flag = isAllowedInSystemBlock(obj) %#ok<INUSD>
            flag = false;
        end
    end
    
end

function pat = interpolatePattern(az,el,pattern,az_q,el_q,interpIn1D)
    if interpIn1D
        pat = interpolatePatternRadians1D(phased.internal.deg2rad(az), ...
            pattern, phased.internal.deg2rad(az_q(:)));
    else
        pat = interpolatePatternRadians2D(phased.internal.deg2rad(az), phased.internal.deg2rad(el), ...
            pattern, phased.internal.deg2rad(az_q(:)), phased.internal.deg2rad(el_q(:)));
    end
end

function pat = interpolatePatternRadians2D(azr, elr, pattern, az_qr, el_qr)
 pat = interp2(azr,elr,pattern,...
        az_qr, el_qr,'nearest',0);
end

function pat = interpolatePatternRadians1D(azr, pattern, az_qr)
 pat = interp1(azr,pattern,...
        az_qr,'nearest',0);
end