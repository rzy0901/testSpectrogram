classdef (Hidden) AbstractBarrageJammer < matlab.System & ...
     matlab.system.mixin.Propagates 
%This class is for internal use only. It may be removed in the future.

%ABSTRACTBARRAGEJAMMER Abstract class for barrage jammer

%   Copyright 2016 The MathWorks, Inc.

% References:
% [1] James Ward, "Space-Time Adaptive Processing for Airborne
%     Radar". MIT Lincoln Lab tech report 1015, 1994.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %ERP     Effective radiated power (W)
        %   Specify the effective radiated power (ERP) (in Watts) of the
        %   jamming signal as a positive scalar. The default value of this
        %   property is 5000.
        ERP = 5000;
    end
    
    properties (Nontunable, PositiveInteger) 
        %SamplesPerFrame    Number of samples per frame
        %   Specify the number of samples in the output jamming signal as a
        %   positive integer. This property applies when you set the
        %   SamplesPerFrameSource property to 'Property'. The default value
        %   of this property is 100.
        SamplesPerFrame = 100;
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
        SeedSource = 'Auto';
        %Seed     Seed for random number generator
        %   Specify the seed for the random number generator as a
        %   non-negative integer. The integer must be less than 2^32. This
        %   property applies when you set the SeedSource property to
        %   'Property'. The default value of this property is 0.
        Seed = 0;
    end
            
    properties(Constant, Hidden)
        SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    properties(Access = protected, Nontunable)
        cNoiseSource;
        pSamplesPerFrameViaProp;
    end
    
    properties(Abstract)
        SamplesPerFrameSource
    end
    
    methods
        function set.ERP(obj, val)
            sigdatatypes.validatePower(val,'phased.BarrageJammer',...
                'ERP',{'scalar'});
            obj.ERP = val;             
        end
        
        function set.Seed(obj,val)
            validateattributes(val,{'double'},{'scalar','nonnegative',...
                'finite','nonnan','nonempty'},'phased.BarrageJammer',...
                'Seed');
            obj.Seed = val;
        end
    end
    
    methods (Access = protected)

        function obj = AbstractBarrageJammer(varargin)
            
            setProperties(obj, nargin, varargin{:}, 'ERP');
        end
        
    end
    
    methods (Access = protected)
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if (obj.SamplesPerFrameSource(1) ~= 'P') ... % Not property
                    && strcmp(prop,'SamplesPerFrame')
                flag = true;
            end
            if (obj.SeedSource(1) == 'A') ... % Auto
                    && strcmp(prop, 'Seed')
                flag = true;
            end
        end
        
        function num = getNumInputsImpl(obj)
            if (obj.SamplesPerFrameSource(1) == 'P')
                num = 0;
            else
                num = 1;
            end
        end
        
        function setupImpl(obj,~) 
            obj.cNoiseSource = phased.internal.NoiseSource(...
                'Distribution','Gaussian','OutputComplex',true,...
                'SeedSource',obj.SeedSource);
            if (obj.SeedSource(1) == 'P')
                obj.cNoiseSource.Seed = obj.Seed;
            end
            obj.pSamplesPerFrameViaProp = ...
                (obj.SamplesPerFrameSource(1) == 'P'); %Property
                
        end

        function flag = isInputComplexityLockedImpl(obj,index) 
            flag = false;
            if ~strncmpi(obj.SamplesPerFrameSource,'Property',1) && ...
                    (index == 1)
                flag = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~) %#ok
            flag = false;
        end
        
        function resetImpl(obj)
            reset(obj.cNoiseSource);
        end
        
        function releaseImpl(obj)
            release(obj.cNoiseSource);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.pSamplesPerFrameViaProp = obj.pSamplesPerFrameViaProp;
                s.cNoiseSource = saveobj(obj.cNoiseSource);
            end
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

    end

    methods (Access = protected) %For Simulink propagation and mask
        function varargout = getOutputDataTypeImpl(~)
            varargout{1} = 'double';
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = true;
        end
    end
end




