classdef (Hidden) AbstractFreeSpace < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%   Copyright 2015-2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %PropagationSpeed Propagation speed (m/s)
        %   Specify the wave propagation speed (in m/s) in free space as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('LightSpeed')
        %OperatingFrequency Signal carrier frequency (Hz)
        %   Specify the carrier frequency (in Hz) of the narrowband signal
        %   as a scalar. The default value of this property is 3e8 (300
        %   MHz).
        OperatingFrequency = 3e8       
    end

    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from Simulink time engine. Set SampleRateFromInputCheckbox to
        %   false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true
    end
    
    properties (Nontunable)
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a scalar. The default value
        %   of this property is 1e6 (1 MHz).
        SampleRate = 1e6    
    end
    
    properties (Nontunable)
        %MaximumDistanceSource  Source of maximum one-way propagation 
        %                       distance
        %   Specify how the maximum one-way propagation distance is
        %   specified as one of 'Auto' | 'Property', where the default is
        %   'Auto'. When you set this property to 'Auto', FreeSpace
        %   automatically allocates the memory to simulate the propagation
        %   delay. When you set this property to 'Property', the maximum
        %   one-way propagation distance is specified via MaximumDistance
        %   property and any signal that needs to propagation more than
        %   MaximumDistance one way is ignored. 
        %
        %   To use FreeSpace in MATLAB Function Block in Simulink, set this
        %   property to 'Property'.
        MaximumDistanceSource = 'Auto'
        %MaximumDistance    Maximum one-way propagation distance (m)
        %   Specify the maximum one-way propagation distance (in meters) as
        %   a positive scalar. This property indicates the maximum distance
        %   the signal can propagate. If the destination is beyond this
        %   distance from the source, the output at the destination is 0.
        %   This property applies when you set the MaximumDistanceSource
        %   property to 'Property'. The default value of this property is
        %   10e3.
        MaximumDistance = 10e3
        %MaximumNumInputSamplesSource  Source of maximum number of samples
        %                       of the input signal
        %   Specify how the maximum number of samples of the input signal
        %   is specified as one of 'Auto' | 'Property', where the default
        %   is 'Auto'. When you set this property to 'Auto', FreeSpace
        %   automatically allocates the memory to buffer the propagated
        %   signal. When you set this property to 'Property', the maximum
        %   number of samples in the input signal is specified via
        %   MaximumNumInputSamples property and any input signal longer
        %   than that value is truncated. This property applies when you
        %   set the MaximumDistanceSource property to 'Property'. The
        %   default value of this property is 'Auto'.
        %
        %   To use the object in MATLAB Function Block in Simulink with
        %   variable-size signals, set this property to 'Property' and set
        %   the MaximumNumInputSamples property.
        MaximumNumInputSamplesSource = 'Auto'
    end
    
    properties (Nontunable, PositiveInteger)
        %MaximumNumInputSamples Maximum number of samples in input signal
        %   Specify the maximum number of samples in the input signal as a
        %   positive scalar. The input signal is the first input, X, and
        %   the number of samples is number of rows in X. This property
        %   applies when you set the MaximumNumInputSamplesSource property
        %   to 'Property'. The default value of this property is 100.
        MaximumNumInputSamples = 100;
    end

    properties(Constant, Hidden)
        MaximumDistanceSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        MaximumNumInputSamplesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties (Access = private, Nontunable)
        %Internal buffer
        cBuffer
    end
    
    properties (Access = protected, Logical, Nontunable)
        %Whether input is a struct
        pIsInputStruct
        %Codegen mode
        pIsCodeGen = false
    end
    
    properties (Access = protected, Nontunable)
        %Sample rate, in MATLAB, specified by property but in Simulink,
        %specified by engine
        pSampleRate
        %Wavelength
        pLambda
        %Propagation range factor, one-way or round trip
        pRangeFactor
        %Valid field names
        pValidFields 
    end
    
    methods (Access = protected, Abstract)
        [xbuf_in,ndelay] = computeDelayedSignal(obj,x,delay)
        [y,propdelay] = computeMultiplePropagatedSignal(obj,x,...
                ncol_per_path,numOfPropPaths,startLoc,endLoc,baseVel,targetVel,Fs) 
        validateNumberOfPositionPairs(obj,x,pos1size,pos2size)
        setRangeFactor(obj)
    end
    
    methods
        function set.OperatingFrequency(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'positive'},'FreeSpace','OperatingFrequency');
            obj.OperatingFrequency = value;
        end
        function set.SampleRate(obj,value)
            validateattributes(value,{'double'},{'scalar','positive',...
                'finite'},'FreeSpace','SampleRate');
            obj.SampleRate = value;
        end
        function set.PropagationSpeed(obj,value)
            sigdatatypes.validateSpeed(value,'FreeSpace',...
                'PropagationSpeed',{'scalar','positive'});
            obj.PropagationSpeed = value;
        end
        function set.MaximumDistance(obj,value)
            sigdatatypes.validateDistance(value,'FreeSpace',...
                'MaximumDistance',{'scalar','positive'});
            obj.MaximumDistance = value;
        end
    end

    methods (Access = protected)
        function obj = AbstractFreeSpace(varargin)
            setProperties(obj, nargin, varargin{:});
            if isempty(coder.target())
                obj.pIsCodeGen = false;
            else
                obj.pIsCodeGen = true;
            end
        end
    end

    methods (Access = protected)

        function flag = isInactivePropertyImpl(obj, prop)
            if (obj.MaximumDistanceSource(1) == 'A') && ...  %Auto
                    (strcmp(prop, 'MaximumDistance') || ...
                    strcmp(prop, 'MaximumNumInputSamplesSource') || ...
                    strcmp(prop, 'MaximumNumInputSamples'))
                flag = true;
            elseif (obj.MaximumDistanceSource(1) == 'P') && ... %Property
                    (obj.MaximumNumInputSamplesSource(1) == 'A') && ...
                    strcmp(prop,'MaximumNumInputSamples')
                flag = true;
            else
                flag = false;
            end
        end
        
        function validateInputSignal(obj,x) 
            if isstruct(x)
                flag_hasXYZ = isfield(x(1),'X') && isfield(x(1),'Y') && isfield(x(1),'Z');
                flag_hasHV = isfield(x(1),'H') && isfield(x(1),'V');
                cond = ~flag_hasXYZ && ~flag_hasHV;
                if cond
                    coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationStruct');
                end
                % cond = ~isscalar(x);
                % if cond
                %     coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationArrayStruct','X');
                % end
                if flag_hasXYZ
                    xsize = size(x);
                    for m = 1:xsize(2)
                        x_x = x(m).X;
                        x_y = x(m).Y;
                        x_z = x(m).Z;
                        cond =  ~isa(x_x,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).X',m),'double');
                        end
                        cond =  ~iscolumn(x_x) || isempty(x_x);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).X',m));
                        end
                        cond =  ~isa(x_y,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).Y',m),'double');
                        end
                        cond =  ~iscolumn(x_y) || isempty(x_y);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Y',m));
                        end
                        cond =  ~isa(x_z,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).Z',m),'double');
                        end
                        cond =  ~iscolumn(x_z) || isempty(x_z);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Z',m));
                        end
                        cond = numel(x_x)~=numel(x_y) || numel(x_x)~=numel(x_z);
                        if cond
                            coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                                'X,Y,Z',sprintf('X(%d)',m));
                        end
                    end

                end
                if flag_hasHV
                    xsize = size(x);
                    for m = 1:xsize(2)
                        x_h = x(m).H;
                        x_v = x(m).V;
                        cond =  ~isa(x_h,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).H',m),'double');
                        end
                        cond =  ~iscolumn(x_h) || isempty(x_h);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).H',m));
                        end
                        cond =  ~isa(x_v,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).V',m),'double');
                        end
                        cond =  ~iscolumn(x_v) || isempty(x_v);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).V',m));
                        end
                        cond = numel(x_h)~=numel(x_v);
                        if cond
                            coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                                'H,V',sprintf('X(%d)',m));
                        end
                    end
                end
                if flag_hasXYZ && flag_hasHV
                    cond = numel(x_x)~=numel(x_h) ;
                    if cond
                        coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                            'X,Y,Z,H,V','X');
                    end
                end
            else
                cond =  ~isa(x,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','X','double');
                end
            end
            
            validateNumChannels(obj,x);
        end
        
        function validateInputsImpl(obj,x,startLoc,endLoc,baseVel,targetVel)
            coder.extrinsic('mat2str');
            coder.extrinsic('num2str');            
            validateInputSignal(obj,x);
            
            cond =  ~isa(startLoc,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','Pos1','double');
            end
            cond =  ~isreal(startLoc);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:NeedReal', 'Pos1');
            end
            cond = any(any(isnan(startLoc)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Pos1');
            end           
            
            cond =  ~isa(endLoc,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','Pos2','double');
            end
            cond =  ~isreal(endLoc);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:NeedReal', 'Pos2');
            end
            cond = any(any(isnan(endLoc)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Pos2');
            end
            
            startLocSize = size(startLoc);
            endLocSize = size(endLoc);
            validateNumberOfPositionPairs(obj,x,startLocSize,endLocSize);
                       
            cond =   ~isa(baseVel,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','Vel1','double');
            end
            baseVelSize = size(baseVel);
            cond =   ~isequal(baseVelSize,startLocSize);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDimensions','Vel1',...
                     coder.const(mat2str(startLocSize)),coder.const(mat2str(baseVelSize)));
            end
            cond =   ~isreal(baseVel);
            if cond
                coder.internal.errorIf(cond, ...
                      'phased:step:NeedReal', 'Vel1');
            end
            cond = any(any(isnan(baseVel)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Vel1');
            end           

            cond =   ~isa(targetVel,'double');
            if cond
                coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDataType','Vel2','double');
            end
            targetVelSize = size(targetVel);
            cond =   ~isequal(targetVelSize,endLocSize);
            if cond
                coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDimensions','Vel2',...
                       coder.const(mat2str(endLocSize)),coder.const(mat2str(targetVelSize)));
            end
            cond =   ~isreal(targetVel);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:NeedReal', 'Vel2');
            end
            cond = any(any(isnan(targetVel)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Vel2');
            end           
        end

        function setupImpl(obj,x,~,~,~,~)
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
            obj.pIsInputStruct = isstruct(x);
            if obj.pIsInputStruct
                 flag_hasXYZ = isfield(x(1),'X') && isfield(x(1),'Y') && isfield(x(1),'Z');
                 flag_hasHV = isfield(x(1),'H') && isfield(x(1),'V');
                 if flag_hasXYZ
                    if flag_hasHV
                        obj.pValidFields = 'XYZHV';                                
                    else
                        obj.pValidFields = 'XYZ';        
                    end
                else
                    obj.pValidFields = 'HV';  
                end
            end
            
            %obj.pSampleRate = getSampleRate(obj,sz_x,1,obj.SampleRate);
            fs = obj.SampleRate; % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            obj.pSampleRate = fs;
            
            setRangeFactor(obj);
            if strcmp(obj.MaximumDistanceSource,'Auto')
                obj.cBuffer = phased.internal.CircularBuffer(...
                    'BufferLength',1);
            else
                buflen = ceil(obj.pRangeFactor*...
                    obj.MaximumDistance/obj.PropagationSpeed*obj.pSampleRate);
                if isInputDataSizePropagated(obj) % Simulink block
                    obj.cBuffer = phased.internal.CircularBuffer(...
                        'FixedLengthBuffer',true,'BufferLength',buflen,...
                        'MaxNumInputSamplesSource','Property',...
                        'MaxNumInputSamples',getPropagatedNumInputSamples(obj,x),...
                        'BufferWidthSource','Property',...
                        'BufferWidth',getSignalBufferWidth(obj));
                else
                    if strcmp(obj.MaximumNumInputSamplesSource,'Auto')
                        obj.cBuffer = phased.internal.CircularBuffer(...
                            'FixedLengthBuffer',true,'BufferLength',buflen,...
                            'MaxNumInputSamplesSource','Property',...
                            'MaxNumInputSamples',getPropagatedNumInputSamples(obj,x),...
                            'BufferWidthSource','Auto');
                    else
                        obj.cBuffer = phased.internal.CircularBuffer(...
                            'FixedLengthBuffer',true,'BufferLength',buflen,...
                            'MaxNumInputSamplesSource','Property',...
                            'MaxNumInputSamples',obj.MaximumNumInputSamples,...
                            'BufferWidthSource','Auto');
                    end
                end
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)  %#ok<INUSL>
            if index == 1
                flag = false;
            else % (index == 2,3,4,5)
                flag = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
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
                s.pLambda = obj.pLambda;
                s.pRangeFactor = obj.pRangeFactor;
                s.cBuffer = saveobj(obj.cBuffer);
                s.pIsInputStruct = obj.pIsInputStruct;
                s.pValidFields = obj.pValidFields;
                s.pSampleRate = obj.pSampleRate;
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

        function [pos1,pos2,vel1,vel2] = preparePositionVelocity(obj,pos1,pos2,vel1,vel2) %#ok<INUSL>
            % pass through by default, use positions and velocities
            % as passed in. subclass can overwrite the behavior, e.g.
            % two-ray model
        end
        
        function [propdelay,propdistance,rspeed] = computePropagationDelayVelocity(obj,startLoc,endLoc,baseVel,targetVel)
            % propagation distance
            propdistance = sqrt(sum(abs(bsxfun(@minus,startLoc,endLoc)).^2));
            rspeed = calcRadialSpeed(endLoc,targetVel,startLoc,baseVel);
            cond =  any(propdistance == 0);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:phased:FreeSpace:InvalidPropagationDistance');
            end
            propdelay = obj.pRangeFactor*propdistance/obj.PropagationSpeed;
        end
        
        function y = stepImpl(obj,x_inc,startPos,endPos,startVel,endVel)

            Fs = obj.pSampleRate;
            % Prepare the location and velocity
            [startLoc,endLoc,baseVel,targetVel] = preparePositionVelocity(...
                obj,startPos,endPos,startVel,endVel);
            x_in = preProcess(obj,x_inc);
            numOfPropPaths = size(x_in,2);
            
            if obj.pIsInputStruct
                % can be moved to setup
                fn = obj.pValidFields;
                sz_field = size(x_in(1).(fn(1)),1);
                num_fn = numel(fn);
                x = complex(zeros([sz_field num_fn*numel(x_in)])); 
                for n = 1:numel(x_in)
                    for m = coder.unroll(1:numel(fn))
                        x(:,(n-1)*num_fn+m) = x_in(n).(fn(m));
                    end
                end
                ncol_per_path = num_fn;
            else
                x = x_in;
                ncol_per_path = 1;
            end
            
            % Apply propagation loss and phase variation
            
            % only consider the time relative to the current pulse start
            % because all previous pulse time is absorbed into the phase
            % term given by the propagation distance as propagation
            % distance is updated by each pulse.
            
            % Calculate propagation delay in samples and prepare signal for
            % fractional delay
            [tempx,propdelay] = computeMultiplePropagatedSignal(obj,x,ncol_per_path,...
                numOfPropPaths,startLoc,endLoc,baseVel,targetVel,Fs);
            nDelay = propdelay*Fs;
            if  numOfPropPaths ~= 1
                isDelayIntSamples = (rem(propdelay,1/Fs) == 0);
                % For exact samples, round nDelay in case numerical error
                % exists
                nDelay(isDelayIntSamples) = round(nDelay(isDelayIntSamples));
                xbuf_in = tempx;
                
                % cFractionalDelay requires fixed input size, make integer
                % delay trivial
                tempDelay = nDelay;
                tempDelay(isDelayIntSamples) = 0;
                
                [tempxbuf,tempDelay_vec] = computeDelayedSignal(obj,...
                    tempx,reshape(repmat(tempDelay,ncol_per_path,1),1,[]));
                
                tempXColIdx = reshape(repmat(~isDelayIntSamples,ncol_per_path,1),1,[]);
                xbuf_in(:,tempXColIdx) = tempxbuf(:,tempXColIdx);
                tempDelay = tempDelay_vec(1:ncol_per_path:end);
                nDelay(~isDelayIntSamples) = tempDelay(~isDelayIntSamples);
                
                y_out = step(obj.cBuffer,xbuf_in,reshape(repmat(nDelay,ncol_per_path,1),1,[]));
                
            else  
                % single delay, done separately to vectorize struct case
                if ~rem(propdelay,1/Fs)
                    % For exact samples, round nDelay in case numerical
                    % error exists
                    nDelay = round(nDelay);
                    xbuf_in = tempx;
                else
                    [xbuf_in,nDelay] = ...
                        computeDelayedSignal(obj,complex(tempx),nDelay);
                end
                y_out = step(obj.cBuffer,xbuf_in,nDelay);
            end
            
            y_out = postProcess(obj,y_out);
            
            if obj.pIsInputStruct
                odimratio = size(y_out,2)/num_fn;
                if ~obj.pIsCodeGen
                    y = repmat(x_inc(1),1,odimratio);
                else
                    for m = coder.unroll(1:num_fn)
                        ys.(fn(m)) = complex(zeros(sz_field,1));
                    end
                    y = repmat(ys,1,odimratio);
                end
                for chanidx = 1:numel(y)
                    for m = coder.unroll(1:num_fn)
                        y(chanidx).(fn(m)) = y_out(1:sz_field,(chanidx-1)*num_fn+m);
                    end
                end
            else
                y = y_out;
            end
            
        end

        function num = getNumInputsImpl(obj)   %#ok<MANU>
            num = 5;
        end
        
        function x = preProcess(obj,x) %#ok<INUSL>
            % no op, for subclass to overload, e.g., TwoRayChannel
        end
        
        function y = postProcess(obj,y) %#ok<INUSL>
            % no op, for subclass to overload, e.g., TwoRayChannel
        end
        
        function num = getSignalBufferWidth(obj)
            % Used in Simulink, nonpolarized
            sz_x = propagatedInputSize(obj,1);
            num = sz_x(2);
        end
    end
    
    methods (Access = protected)
        
        function sz_out = getOutputSizeImpl(obj)
            sz_out = propagatedInputSize(obj,1);
        end
        
        function dt_out = getOutputDataTypeImpl(obj)
            dt_out = propagatedInputDataType(obj,1);
        end
        
        function cp_out = isOutputComplexImpl(obj) %#ok<MANU>
            cp_out = true;
        end
        
        function fsz_out = isOutputFixedSizeImpl(obj) 
            fsz_out = propagatedInputFixedSize(obj, 1);
        end
        
        function varargout = getInputNamesImpl(obj)   %#ok<MANU>
            varargout = {'X','Pos1','Pos2','Vel1','Vel2'};
        end

        function varargout = getOutputNamesImpl(obj)  %#ok<MANU>
            varargout = {''};
        end
        
    end
    
end

function rspeed = calcRadialSpeed(tgtpos,tgtvel,refpos,refvel)
%calcRadialSpeed    Compute radial speed
%   RSPEED = calcRadialSpeed(POS,VEL,REFPOS,REFVEL) compute the relative
%   speed RSPEED (in m/s) for a target at position POS (in meters) with a
%   velocity VEL (in m/s) relative to the reference position REFPOS (in
%   meters) and reference velocity REFVEL (in m/s).

%   This is the same functionality as radialspeed function. However,
%   because we already done the input validation here, we want to skip the
%   validation to improve the performance. In addition, here all position
%   and velocity are always column vectors and the target and reference can
%   never be colocated, so it simplifies the computation too.

tgtdirec = bsxfun(@minus,tgtpos,refpos);
veldirec = bsxfun(@minus,tgtvel,refvel);

%Take the 2-norm of veldirec and tgtdirec
rn = sqrt(sum(tgtdirec.^2));

% negative sign to ensure that incoming relative speed is positive
rspeed = -(sum(veldirec.*tgtdirec)./rn);

end
