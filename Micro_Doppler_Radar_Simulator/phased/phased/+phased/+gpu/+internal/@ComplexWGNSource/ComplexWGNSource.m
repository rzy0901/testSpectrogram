classdef (Hidden, Sealed) ComplexWGNSource < matlab.system.internal.gpu.GPUBase
%This class is for internal use only. It may be removed in the future.

%ComplexWGNSource   Complex white Gaussian noise source
%   H = phased.gpu.internal.ComplexWGNSource creates a complex white
%   Gaussian noise source System object, H. The object generates complex
%   white Gaussian noise with specified dimensions. The object uses the
%   GPU for random number generation and computation.
%
%   Use of this object requires a Parallel Computing Toolbox license.
%
%   H = phased.gpu.internal.ComplexWGNSource(Name,Value) creates a complex
%   white Gaussian noise source object, H, with the specified property Name
%   set to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   X = step(H,POW,SZ) generates samples of complex white Gaussian noise in
%   X whose power is specified in scalar POW (in watts). SZ is a vector
%   specifying the dimensions of X. Note that the return class of X is a
%   gpuArray.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ComplexWGNSource methods:
%
%   step     - Calculate SMI weights (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create SMI weights estimator object with same property
%              values
%   isLocked - Locked status (logical)
%
%   ComplexWGNSource properties:
%
%   SeedSource - Source of seed for random number generator
%   Seed       - Seed for random number generator
%
%   % Example:
%   %   Generates a 100x2 matrix of complex white Gaussian noise samples
%   %   with a noise power of 0.1 watts.
%
%   hns = phased.gpu.internal.ComplexWGNSource;
%   x = step(hns,0.1,[100 2]);
%
%   See also phased.internal.ComplexWGNSource

%   Copyright 2012-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>


%  NOTE: This class inherits from GPUBase not GPUSystem because the output
%  is always a gpuArray. The inputs are (typically) never gpuArrays. This
%  is a helper object for which doubles in -> gpuArrays out is the best
%  protocol. Since this isa GPUBase there is no setupGPUImpl or
%  stepGPUImpl. To be explicit, this class does not follow the protocol
%  doubles in ->doubles out, gpuArrays in -> gpuArrays out like all the
%  other GPU System objects.


    properties (Nontunable)

        %SeedSource   Source of seed for random number generator
        %   Specify how the random numbers are generated as one of 'Auto' |
        %   'Property', where the default is 'Auto'. When you set this
        %   property to 'Auto', the random numbers are generated using the
        %   default MATLAB GPU random number generator. When you set this
        %   property to 'Property', a private random number generator is
        %   used with a seed specified by the value of the Seed property.
        SeedSource = 'Auto';
        %Seed     Seed for random number generator
        %   Specify the seed for the random number generator as a
        %   non-negative integer. This property applies when you set the
        %   SeedSource property to 'Property'. The default value of this
        %   property is 0.
        Seed = 0;
    end
    
    properties(Access = protected, Nontunable)
        cRandStream;
    end

    properties(Access = private, Logical, Nontunable)
        pUseGlobalStream = false
    end

    properties(Constant, Hidden)
        SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    methods 
        function obj = ComplexWGNSource(varargin)
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
              obj.cRandStream = makeLocalRandStream(obj, obj.Seed); 
            end
        end
            
        function flag = isInputComplexityLockedImpl(obj,~) %#ok<MANU>
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~) %#ok<MANU>
            flag = false;
        end
        
        function resetImpl(obj)
            if ~obj.pUseGlobalStream
                reset(obj.cRandStream,obj.Seed);
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            if strncmp(obj.SeedSource,'Auto',1) && strcmp(prop,'Seed')
                flag = true;
            else
                flag = false;
            end
        end
        
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function validateInputsImpl(~,pow,sz)
          if ~isscalar(pow)
            matlab.system.internal.error(...
              'MATLAB:system:inputMustBeScalar','POW');
          end
          if ~isa(pow,'double')
            matlab.system.internal.error(...
              'MATLAB:system:invalidInputDataType','POW','double');
          end          
          if ~isreal(pow)
            matlab.system.internal.error(...
              'phased:ComplexWGNSource:MustBeReal', 'POW');
          end
          
          if ~isrow(sz) || isempty(sz)
            matlab.system.internal.error(...
              'MATLAB:system:inputMustBeRowVector','SZ');
          end
          if ~isa(sz,'double')
            matlab.system.internal.error(...
              'MATLAB:system:invalidInputDataType','SZ','double');
          end          
          if ~isreal(sz)
            matlab.system.internal.error(...
              'phased:ComplexWGNSource:MustBeReal', 'SZ');
          end
        end
            
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.pUseGlobalStream = obj.pUseGlobalStream;
                if ~obj.pUseGlobalStream
                    s.pRandState = obj.cRandStream.State;
                end
            end
        end

        function s = loadSubObjects(obj,s, wasLocked)
            if wasLocked 
              if ~s.pUseGlobalStream,
                    obj.cRandStream = makeLocalRandStream(obj); 
                    obj.cRandStream.State = s.pRandState;
                    s = rmfield(s,'pRandState');
              end
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s, wasLocked);
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
            if obj.pUseGlobalStream,
              rs = parallel.gpu.RandStream.getGlobalStream();
            else
              rs = obj.cRandStream;
            end 
            %The code below is faster than calling randn once for twice the
            %size and making the result complex
            y = sqrt(pow/2)*(  gpuArray.randn(rs,sz, class(pow)) + ...
                               1i*gpuArray.randn(rs,sz, class(pow)));

        end
    end
   
    methods (Access = private)
      function s = makeLocalRandStream(obj, seed)
        %In case the default non-global stream (for Auto Mode) needs
        %changing, just do it here:
        if nargin > 1,
          s = parallel.gpu.RandStream('CombRecursive', ...
              'NormalTransform','Inversion', 'Seed', seed); 
        else
          s = parallel.gpu.RandStream('CombRecursive', ...
              'NormalTransform','Inversion'); 
        end
      end
    end

end


