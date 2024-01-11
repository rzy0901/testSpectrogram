classdef (Hidden, Sealed, StrictDefaults) NoiseSource < matlab.System
%This class is for internal use only. It may be removed in the future.

%NoiseSource   Noise source
%   H = phased.internal.NoiseSource creates a noise source System object,
%   H. The object generates white noise with specified distribution and
%   dimensions.
%
%   H = phased.internal.NoiseSource(Name,Value) creates a noise source
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   X = step(H,K,SZ) generates samples of white noise in X. K is a scalar
%   indicating the scaling factor of the random number. If the distribution
%   is Gaussian, K can be considered as the power of the noise. SZ is a
%   vector specifying the dimensions of X.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   NoiseSource methods:
%
%   step     - Calculate SMI weights (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create SMI weights estimator object with same property
%              values
%   isLocked - Locked status (logical)
%
%   NoiseSource properties:
%
%   SeedSource - Source of seed for random number generator
%   Seed       - Seed for random number generator
%
%   % Example:
%   %   Generates a 100x2 matrix of complex white Gaussian noise samples
%   %   with a noise power of 0.1 watts.
%
%   hns = phased.internal.NoiseSource('Distribution','Gaussian',...
%           'OutputComplex',true);
%   x = step(hns,0.1,[100 2]);
%
%   See also phased.

%   Copyright 2012-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %Distribution   Distribution of noise
        %   Specify distribution of noise as one of 'Uniform' | 'Gaussian',
        %   where the default is 'Uniform'.
        Distribution = 'Uniform'
        %SeedSource   Source of seed for random number generator
        %   Specify how the random numbers are generated as one of 'Auto' |
        %   'Property', where the default is 'Auto'. When you set this
        %   property to 'Auto', the random numbers are generated using the
        %   default MATLAB random number generator. When you set this
        %   property to 'Property', a private random number generator is
        %   used with a seed specified by the value of the Seed property.
        %
        %   To use this object with Parallel Computing Toolbox software or
        %   a GPU, set this property to 'Auto'.
        SeedSource = 'Auto';
        %Seed     Seed for random number generator
        %   Specify the seed for the random number generator as a
        %   non-negative integer. The integer must be less than 2^32. This
        %   property applies when you set the SeedSource property to
        %   'Property'. The default value of this property is 0.
        Seed = 0;
    end
    
    properties(Logical, Nontunable)
        %OutputComplex  Output complex noise
        %   Set this property to true to output complex noise. Set this
        %   property to false to output real noise. The default value of
        %   this property is false.
        %
        %   The complex noise is formed using identical distribution in
        %   both real and imaginary part. The scaling factor is equally
        %   split between the two, i.e., sqrt(K/2).
        OutputComplex = false
    end
    
    properties(Access = protected, Nontunable)
        cRandStream;
    end
    
    properties(Access = private, Logical, Nontunable)
        pUseGlobalStream = true
    end
    
    properties(Constant, Hidden)
        SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        DistributionSet = matlab.system.StringSet({'Uniform',...
            'Gaussian'});
    end
    
    methods 

        function obj = NoiseSource(varargin)

            setProperties(obj, nargin, varargin{:});

        end

    end
    
    methods
        function set.Seed(obj,val)
            validateattributes(val,{'double'},{'scalar','nonnegative',...
                'finite','nonnan','nonempty'},'phasedRandomSource',...
                'Seed');
            obj.Seed = val;
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj,~,~)
            if obj.SeedSource(1) == 'A' %Auto
                obj.pUseGlobalStream = true;
            else
                cond = ~isempty(coder.target);
                if cond
                    rng(obj.Seed);
                else
                    obj.cRandStream = RandStream('mt19937ar','Seed',obj.Seed);
                end
                obj.pUseGlobalStream = false;
            end
        end
            
        function flag = isInputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function resetImpl(obj)
            if ~obj.pUseGlobalStream
                if coder.target('MATLAB')
                    reset(obj.cRandStream,obj.Seed);
                else
                    rng(obj.Seed);
                end
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            if (obj.SeedSource(1) == 'A') && strcmp(prop,'Seed') %Auto and Seed
                flag = true;
            else
                flag = false;
            end
        end
        
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
               
        function validateInputsImpl(~,pow,sz)
            cond =  ~isscalar(pow);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeScalar','K');
            end
            cond =  ~isa(pow,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','K','double');
            end
            cond =  ~isreal(pow);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:ComplexWGNSource:MustBeReal', 'K');
            end
            cond =  ~isrow(sz) || isempty(sz);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeRowVector','SZ');
            end
            cond =  ~isa(sz,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','SZ','double');
            end
            cond =  ~isreal(sz);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:ComplexWGNSource:MustBeReal', 'SZ');
            end
        end
            
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.pUseGlobalStream = obj.pUseGlobalStream;
                if ~obj.pUseGlobalStream
                    if coder.target('MATLAB')
                        s.pRandState = obj.cRandStream.State;
                    else
                        s.pRandState = rng;
                    end
                end
            end
        end

        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                if ~s.pUseGlobalStream
                    if coder.target('MATLAB')
                        obj.cRandStream = RandStream('mt19937ar');
                        obj.cRandStream.State = s.pRandState;
                    else
                        rng(s.pRandState);
                    end
                    s = rmfield(s,'pRandState');
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
        
        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end

        function y = stepImpl(obj,pow,sz)
            % pow must be positive and sz must be integers. 
            if obj.pUseGlobalStream
                if (obj.Distribution(1) == 'U') %Uniform
                    if obj.OutputComplex
                        yr = rand(sz);
                        yi = rand(sz);
                        y = sqrt(pow/2)*complex(yr,yi);
                    else
                        y = sqrt(pow)*rand(sz);    
                    end
                else % Gaussian
                    if obj.OutputComplex
                        yr = randn(sz);
                        yi = randn(sz);
                        y = sqrt(pow/2)*complex(yr,yi);
                    else
                        y = sqrt(pow)*randn(sz);
                    end
                end
            else
                if coder.target('MATLAB')
                    rs = obj.cRandStream;
                    if (obj.Distribution(1) == 'U') %Uniform
                        if obj.OutputComplex
                            yr = rand(rs,sz);
                            yi = rand(rs,sz);
                            y = sqrt(pow/2)*complex(yr,yi);
                        else
                            y = sqrt(pow)*rand(rs,sz);
                        end
                    else % Gaussian
                        if obj.OutputComplex
                            yr = randn(rs,sz);
                            yi = randn(rs,sz);
                            y = sqrt(pow/2)*complex(yr,yi);
                        else
                            y = sqrt(pow)*randn(rs,sz);
                        end
                    end
                else
                    if (obj.Distribution(1) == 'U') %Uniform
                        if obj.OutputComplex
                            yr = rand(sz);
                            yi = rand(sz);
                            y = sqrt(pow/2)*complex(yr,yi);
                        else
                            y = sqrt(pow)*rand(sz);
                        end
                    else % Gaussian
                        if obj.OutputComplex
                            yr = randn(sz);
                            yi = randn(sz);
                            y = sqrt(pow/2)*complex(yr,yi);
                        else
                            y = sqrt(pow)*randn(sz);
                        end
                    end
                end
            end
        end
    end
end


% [EOF]
