classdef (Hidden) AbstractTarget < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%ABSTRACTRADARTARGET Define the ABSTRACTRADARTARGET class
% This is an abstract class in support of radar target functionality.

%   Copyright 2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Model Fluctuation model
        %   Specify the statistical model of the target as one of
        %   'Nonfluctuating' | 'Swerling1' | 'Swerling2' | 'Swerling3' |
        %   'Swerling4', where the default is 'Nonfluctuating'.
        %
        %   When you set this property to anything other than
        %   'Nonfluctuating', you need to invoke step method with an extra
        %   UPDATE input.
        Model = 'Nonfluctuating'
        %SeedSource   Source of seed for random number generator
        %   Specify how the random numbers are generated as one of 'Auto' |
        %   'Property', where the default is 'Auto'. The random numbers are
        %   used to model random pattern values. When you set this property
        %   to 'Auto', the random numbers are generated using the default
        %   MATLAB random number generator. When you set this property to
        %   'Property', a private random number generator is used with a
        %   seed specified by the value of the Seed property. This property
        %   applies when you set the Model property to 'Swerling1',
        %   'Swerling2', 'Swerling3', or 'Swerling4'.
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
    
    properties(Access = protected, Nontunable)
        % Private random stream
        cNoiseSource;
        % Private flag about whether the target is fluctuating
        pFluctuate;
        % Private function handle for Swerling models
        pFluctuateFunc;
        %Codegen mode
        pIsCodeGen = false
    end

    properties(Constant, Hidden)
        ModelSet = matlab.system.StringSet({ 'Nonfluctuating', 'Swerling1',...
            'Swerling2', 'Swerling3', 'Swerling4'});
        SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end

    methods (Access = protected)
        function obj = AbstractTarget(varargin)
            setProperties(obj, nargin, varargin{:});
            if isempty(coder.target())
                obj.pIsCodeGen = false;
            else
                obj.pIsCodeGen = true;
            end
        end
    end

    methods
        function set.Seed(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative',...
                'integer','<',2^32},'phased.AbstractTarget','Seed');
            obj.Seed = value;
        end  
    end
    
    methods (Access = protected)
        
        function setupImpl(obj,x) 
            setupImpl@phased.internal.AbstractSampleRateEngine(obj);
            SwerlingModel = obj.Model;
            obj.pFluctuate = (SwerlingModel(1) == 'S');%Swerling
            if obj.pFluctuate
                SwerlingIdx = SwerlingModel(end);
                switch SwerlingIdx
                    case {'1','2'}
                        obj.pFluctuateFunc = @privExponential;
                    case {'3','4'}
                        obj.pFluctuateFunc = @privChiSquare;
                end
            end
            %Generate a random seed
            if obj.pFluctuate
                obj.cNoiseSource = phased.internal.NoiseSource(...
                    'SeedSource',obj.SeedSource);
                if (obj.SeedSource(1) == 'P') %Property
                    obj.cNoiseSource.Seed = obj.Seed;
                end
            end
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractSampleRateEngine(obj);
            if obj.pFluctuate
                reset(obj.cNoiseSource);
            end
        end

        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
            if obj.pFluctuate
                release(obj.cNoiseSource);
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            if isLocked(obj)
                s.pFluctuate = obj.pFluctuate;
                if obj.pFluctuate
                    s.cNoiseSource = saveobj(obj.cNoiseSource);
                    s.pFluctuateFunc = obj.pFluctuateFunc;
                end
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                obj.pFluctuate = s.pFluctuate;
                if s.pFluctuate
                    obj.cNoiseSource = phased.internal.NoiseSource.loadobj(s.cNoiseSource);
                    s = rmfield(s,'cNoiseSource');
                    obj.pFluctuateFunc = s.pFluctuateFunc;
                    s = rmfield(s,'pFluctuateFunc');
                end
                s = rmfield(s,'pFluctuate');
            end
        end
        
        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end

        function fsz_out = isOutputFixedSizeImpl(obj) 
            fsz_out = propagatedInputFixedSize(obj, 1);
        end
        
        function flag = isInputComplexityLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end

        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;  % index == 1
        end

        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if (obj.Model(1) == 'N') %NonFluctuating
                if strcmp(prop,'SeedSource') || strcmp(prop,'Seed')
                    flag = true;
                end
            elseif (obj.SeedSource(1) == 'A') && ... %Auto
                    strcmp(prop,'Seed')
                flag = true;
            end
        end
    end
    
end

function sigma = privExponential(rs,meanval,monostaticpolflag)
% generate exponential random numbers

if nargin < 3
    polflag = false;
else
    polflag = true;
end

if polflag
    if monostaticpolflag
        % polarization with monostatic, meanval symmetric
        temp = step(rs,1,[1 3]);
        u = [temp(1) temp(2);temp(2) temp(3)];
    else
        u = step(rs,1,size(meanval));
    end
    sigma = meanval.*sqrt(-log(u));  % scattering matrix, take sqrt
else %nonpolarized 
    u = step(rs,1,size(meanval));
    sigma = -meanval.*log(u);
end
end

function sigma = privChiSquare(rs,meanval,monostaticpolflag)
% generate chi-square (degree 4) random numbers

    if nargin < 3
        polflag = false;
    else
        polflag = true;
    end
    % first generate standard chi-square y
    % Eq 5.6.5 in Casella and Berger
    if polflag
        if monostaticpolflag
            % polarization with monostatic, meanval symmetric
            temp = -2*sum(log(step(rs,1,[2 3])));
            y = [temp(1);temp(2);temp(2);temp(3)];
        else
            y = -2*sum(log(step(rs,1,[2 numel(meanval(:))])));
        end
        % transform.
        sigma = sqrt(y(:)/4).*meanval(:);   % scattering matrix, take sqrt
        sigma = reshape(sigma,size(meanval));
    else
        %meanval = 1xN
        y = -2*sum(log(step(rs,1,[2 numel(meanval)])));
        sigma = reshape(y,size(meanval)).*meanval/4;
    end
end
