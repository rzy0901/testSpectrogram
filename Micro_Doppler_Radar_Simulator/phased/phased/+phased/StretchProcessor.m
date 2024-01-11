classdef (Sealed,StrictDefaults) StretchProcessor < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%StretchProcessor   Stretch processor for linear FM waveform
%   H = phased.StretchProcessor creates a stretch processor System object,
%   H. This object performs stretch processing on the input data.
%
%   H = phased.StretchProcessor(Name,Value) creates a stretch processor
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) applies the stretch processing to the input X and returns
%   the processed result in Y. The processing is applied along the first
%   dimension where each column in X represents one receiving pulse. Y and
%   X have the same dimensions.
%
%   Y = step(H,X,PRF) uses PRF as the pulse repetition frequency (in
%   Hz), when you set the PRFSource property to 'Input port'. PRF must be a
%   scalar.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   StretchProcessor methods:
%
%   step     - Perform stretch processing (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a stretch processor object with same property values
%   isLocked - Locked status (logical)
%
%   StretchProcessor properties:
%
%   SampleRate       - Sample rate 
%   PulseWidth       - Pulse width 
%   PRFSource        - Source of PRF
%   PRF              - Pulse repetition frequency
%   SweepSlope       - FM sweep slope 
%   SweepInterval    - FM sweep interval
%   PropagationSpeed - Propagation speed
%   ReferenceRange   - Reference range 
%   RangeSpan        - Range span
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Use stretch processing to locate a target at a range of 4950 
%   %   meters. 
%
%   % simulate signal
%   waveform = phased.LinearFMWaveform; x = waveform();
%   c = 3e8; r = 4950; num_sample = r/(c/(2*waveform.SampleRate));
%   x = circshift(x,num_sample);
%
%   % perform stretch processing
%   stretch = getStretchProcessor(waveform,5000,200,c);
%   y = stretch(x);
%
%   % plot the spectrum of the stretch processed signal
%   periodogram(y,[],2048,stretch.SampleRate,'centered');
%   
%   % detect the range
%   [Pxx, F] = periodogram(y,[],2048,stretch.SampleRate,'centered');
%   [~,rngidx] = findpeaks(pow2db(Pxx/max(Pxx)),...
%                   'MinPeakHeight',-5);
%   rngfreq = F(rngidx);
%   re = stretchfreq2rng(rngfreq,stretch.SweepSlope,...
%               stretch.ReferenceRange,c)
%
%   See also phased.LinearFMWaveform, phased.MatchedFilter,
%   stretchfreq2rng.


%   Copyright 2011-2017 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [2] Bassem R. Mahafza, Radar Signal Analysis And Processing Using
%       MATLAB, CRC Press, 2009

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)

    %SampleRate Sample rate (Hz)
    %   Specify the sample rate (in Hz) as a positive scalar. The default
    %   value of this property is 1e6 (1 MHz).
    SampleRate = 1e6;
    %PulseWidth Pulse width (s)
    %   Specify the length of each pulse (in seconds) as a positive scalar.
    %   The default value of this property is 50e-6.
    PulseWidth = 50e-6;
    %PRFSource    Source of PRF
    %   Specify how to determine the PRF for the Stretch processor as one
    %   of 'Auto'| 'Property' | 'Input port', where the default is
    %   'Property'. When you set this  property to 'Auto', the PRF is
    %   inferred from the number of rows in the input signal . When you set
    %   this property to 'Property', the PRF is determined by the value of
    %   the PRF property. When you set this property to 'Input port', the
    %   PRF is determined by the input argument.
    PRFSource = 'Property' 
    %PRF    Pulse repetition frequency (Hz)
    %   Specify the pulse repetition frequency (in Hz) as a positive
    %   scalar. The default value of this property is 1e4 (10 kHz).
    PRF = 1e4
    %SweepSlope FM sweep slope (Hz/s)
    %   Specify the slope of the linear FM sweeping (in Hz/s) as a scalar.
    %   The default value is 2e9.
    SweepSlope = 2e9
    %SweepInterval  FM sweep interval
    %   Specify the linear FM sweeping interval using one of 'Positive' |
    %   'Symmetric', where the default is 'Positive'. If SweepInterval is
    %   'Positive', the waveform sweeps in the interval between 0 and B
    %   where B is the sweeping bandwidth. If SweepInterval is 'Symmetric',
    %   the waveform sweeps in the interval between -B/2 and B/2.
    SweepInterval = 'Positive'
    %PropagationSpeed   Propagation speed (m/s)
    %   Specify the propagation speed (in m/s) of the signal as a scalar.
    %   The default value of this property is the speed of light.
    PropagationSpeed = physconst('LightSpeed');
end

properties
    %ReferenceRange     Reference range (m)
    %   Specify the reference range, i.e., the center of ranges of interest
    %   (in meters) as a positive scalar. The ReferenceRange must be within
    %   the unambiguous range of one pulse. The default value of this
    %   property is 5000. This property is tunable.
    ReferenceRange = 5000
end

properties (Nontunable)
    %RangeSpan  Range span (m)
    %   Specify the length of the interval for ranges of interest (in
    %   meters) as a positive scalar. The range span is centered at the
    %   range value specified in the ReferenceRange property. The default
    %   value of this property is 500.
    RangeSpan = 500
end

properties(Constant, Hidden)
    PRFSourceSet = matlab.system.StringSet({'Auto','Property','Input port'});
    SweepIntervalSet = matlab.system.StringSet({'Positive','Symmetric'});
end

properties (Access=private)
    pMixingPulse
end

properties (Access = private)
    pPRF 
end

properties (Access=private, Nontunable)
    pSampleRate 
end

properties (Access = private, Logical)
    pSizeInitialized = false
end

methods
    function set.PropagationSpeed(obj,val)
        sigdatatypes.validateSpeed(val,...
            '','PropagationSpeed',...
            {'double','single'},{'scalar','positive'});
        obj.PropagationSpeed = val;
    end
    function set.SampleRate(obj, value)
        validateattributes(value,{'double','single'}, {'scalar',...
            'positive','finite'},...
            '','SampleRate');
        obj.SampleRate = value;
    end
    function set.PRF(obj,value)
        validateattributes( value, { 'double','single' }, { 'scalar',...
            'positive', 'finite' }, '', 'PRF');
        obj.PRF = value;
    end
    function set.PulseWidth(obj, value)
        validateattributes(value,{'double','single'},...
            {'scalar','positive','finite'},...
            '','PulseWidth');
        obj.PulseWidth = value;
    end
    function set.SweepSlope(obj, value)
        validateattributes(value,{'double','single'},...
            {'scalar','real','finite'},...
            '','SweepSlope');
        obj.SweepSlope = value;
    end
    function set.ReferenceRange(obj,value)
        sigdatatypes.validateDistance(value,...
            '','ReferenceRange',{'double','single'},{'scalar','positive'});
        obj.ReferenceRange = value;
    end
    function set.RangeSpan(obj,value)
        sigdatatypes.validateDistance(value,...
            '','RangeSpan',{'double','single'},{'scalar','positive'});
        obj.RangeSpan = value;
    end
end

methods
    function obj = StretchProcessor(varargin)
        %StretchProcessor   Construct the StretchProcessor class.
        setProperties(obj,nargin,varargin{:});
    end
end

methods (Hidden, Static)
    function validateStretchRangeSpan(c,prf,pw,refrng,rngspan)
        Rmin = 0;
        Rmax = c*(1/(2*prf)-pw/2);
        Rinterval = refrng+[-1 1]*rngspan/2;
        cond = Rinterval(1)<Rmin || Rinterval(2)>Rmax;
        if cond
            coder.internal.errorIf(cond, ...
                'phased:phased:stretch:InvalidStretchROI',...
                'ReferenceRange','RangeSpan',...
                feval('sprintf','%5.4f',Rmin),...
                feval('sprintf','%5.4f',Rmax),...
                feval('sprintf','%5.4f',Rinterval(1)),...
                feval('sprintf','%5.4f',Rinterval(2)));
        end
    end
end

methods (Access = protected)

    function validatePropertiesImpl(obj)
        if (obj.PRFSource(1) == 'P')
            validatePulseWidthPRF(obj,obj.PRF);
            
            phased.StretchProcessor.validateStretchRangeSpan(...
                obj.PropagationSpeed,obj.PRF,obj.PulseWidth,...
                obj.ReferenceRange,obj.RangeSpan);
        end
    end
    
    function setupImpl(obj,x,prf)
        % setupImpl   Add any one-time calculations for the object
        sz_x = size(x);
        obj.pNumInputChannels = sz_x(2);
        obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        obj.pSampleRate = getSampleRate(obj,size(x,1),1,obj.SampleRate);
        
        switch obj.PRFSource
            case 'Auto'
                prf = obj.pSampleRate/sz_x(1);
            case 'Property'
                prf = obj.PRF;
        end
        obj.pPRF = prf;
        
        if (obj.PRFSource(1) == 'P')
            validateSignalLength(obj,prf,sz_x(1));
            processTunedPropertiesImpl(obj);
        end
    end
    
    function validateSignalLength(obj,prf,len)
        cond =  any(rem(obj.pSampleRate,prf));
        if cond
            coder.internal.errorIf(cond, ...
                'phased:Waveform:NeedRatioInteger', 'SampleRate', 'PRF');
        end

        nrows = round(obj.pSampleRate/prf);
        cond =  len ~= nrows;
        if cond
            coder.internal.errorIf(cond, ...
                'phased:phased:invalidRowNumbers','X',nrows);
        end
    end
    
    function validatePulseWidthPRF(obj,prf)
        cond =  any(obj.PulseWidth.*prf > 1);
        if cond
            coder.internal.errorIf(cond, ...
                'phased:Waveform:NotLessThanOrEqualTo', 'PulseWidth',...
                feval('sprintf', '%5.2e', 1/max(prf)));
        end
            
    end
    
    function s = computeMixingPulse(obj,prf)
        Nsamp = round(obj.pSampleRate/prf); % SampleRate/PRF is integer
        t = (0:Nsamp-1)/obj.pSampleRate;
        c = obj.PropagationSpeed;
        beta = obj.SweepSlope;
        pw = obj.PulseWidth;
        refrng = obj.ReferenceRange;
        rngspan = obj.RangeSpan;
        RangeSpan_max = refrng+rngspan/2;
        RangeSpan_min = refrng-rngspan/2;
        tau_center = 2*obj.ReferenceRange/c;
        tau_max = 2*RangeSpan_max/c+pw;
        tau_min = 2*RangeSpan_min/c;
        
        t_foundmin_idx = find(t>=tau_min,1,'first');
        t_min_idx = t_foundmin_idx(1);
        
        t_foundmax_idx = find(t<=tau_max,1,'last');
        t_max_idx = t_foundmax_idx(1);
        
        s = complex(zeros(size(t(:))));
        t_mix = t(t_min_idx:t_max_idx)-tau_center;
        if (obj.SweepInterval(1) == 'P') %Positive
            if beta > 0
                s(t_min_idx:t_max_idx) = exp(1i*pi*beta*t_mix.^2);
            else
                s(t_min_idx:t_max_idx) = exp(1i*pi*beta*(t_mix.^2-t_mix*2*pw));
            end
        else
            s(t_min_idx:t_max_idx) = exp(1i*pi*beta*(t_mix.^2-t_mix*pw));
        end
    end        

    function processTunedPropertiesImpl(obj)
        if (obj.PRFSource(1) == 'P')
            s = computeMixingPulse(obj,obj.PRF);
            obj.pMixingPulse = s;
        end
    end
    
    function processInputSizeChangeImpl(obj,x,~)
        sz_x = size(x);
        if(obj.PRFSource(1) == 'A')
            prf = obj.pSampleRate/sz_x(1);
            validatePulseWidthPRF(obj,prf);
            
            phased.StretchProcessor.validateStretchRangeSpan(...
                obj.PropagationSpeed,prf,obj.PulseWidth,...
                obj.ReferenceRange,obj.RangeSpan);
            obj.pPRF = prf;
        end
    end
        
    function y = stepImpl(obj,x,prfinput)
        
        classtoUse=class(x);
        if ~obj.pSizeInitialized
            processInputSizeChangeImpl(obj,x);
            obj.pSizeInitialized = true;
        end
        if obj.PRFSource(1) == 'I'
            prf = validateInputPRF(obj,prfinput);
            
            validateSignalLength(obj,prf,size(x,1));
            validatePulseWidthPRF(obj,prf);
            
            phased.StretchProcessor.validateStretchRangeSpan(...
                obj.PropagationSpeed,prf,obj.PulseWidth,...
                obj.ReferenceRange,obj.RangeSpan);
            prf=cast(prf,classtoUse);
        else
            prf = cast(obj.pPRF,classtoUse);
        end
        
        if ~(obj.PRFSource(1) == 'P')
            s = computeMixingPulse(obj,prf);
        else
            s = obj.pMixingPulse;
        end
        
        y = complex(zeros(size(x),classtoUse));
        for m = size(x,2):-1:1
            y(:,m) = x(:,m).*conj(s);
        end
    end
    
    function releaseImpl(obj)
        releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
        obj.pSizeInitialized = false;
    end
    
    function validateInputsImpl(obj,x,prfinput)
        if (obj.PRFSource(1) == 'I')
            validateInputPRFSpec(obj,prfinput);
        end

        cond =   ~isa(x,'float');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','float');
        end
        cond =   ~ismatrix(x) || isempty(x);
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeMatrix','X');
        end
        
        validateNumChannels(obj,x)
    end
    
    function flag = isInputComplexityLockedImpl(obj,index)
        flag = false;
        if ((obj.PRFSource(1) == 'I') && (index == 2))
            flag = true;
        end
    end
    
    function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
        flag = false;
    end
    
    function flag = isInputSizeLockedImpl(obj,index) 
        if index == 1
            if strcmp(obj.PRFSource,'Property')
                flag = true;
            else
                flag = false;
            end
        else
            flag = true;
        end
    end
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        if (strcmp(prop,'PRF') && (~(obj.PRFSource(1) == 'P')))
            flag = true;
        end
    end
    
    function num = getNumOutputsImpl(obj) %#ok<MANU>
        num = 1;
    end
    
    function num = getNumInputsImpl(obj) 
        num = 1;
        if (obj.PRFSource(1) == 'I')
            num = num+1;
        end
    end
       
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
        s.isLocked = isLocked(obj);        
        if isLocked(obj)
            s.pSizeInitialized = obj.pSizeInitialized;
            s.pPRF = obj.pPRF;
            s.pMixingPulse = obj.pMixingPulse;
            s.pSampleRate = obj.pSampleRate;          
        end
    end
    
    function s = loadSubObjects(obj,s) 
        if s.isLocked
            % recover locked sample rate information
            if isfield(s,'pSampleRate')
                obj.pSampleRate = s.pSampleRate;
                s = rmfield(s,'pSampleRate');
            else
                obj.pSampleRate = s.SampleRate;
            end
            if isfield(s,'pPRF')
                obj.pPRF = s.pPRF;
                s = rmfield(s,'pPRF');
            else
                obj.pPRF = s.PRF;
            end
        end
        s = rmfield(s,'isLocked');
    end

    function loadObjectImpl(obj,s,~)
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end
end

methods (Access=protected)
    function validateInputPRFSpec(~,prfArg)
        cond = ~isa(prfArg,'float');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','PRF','float');
        end
        cond = ~isscalar(prfArg);
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeScalar','PRF');
        end
    end
    
    function prf = validateInputPRF(~,prfArg)
        cond = ~isfinite(prfArg);
        if cond
            coder.internal.errorIf(cond, ...
                'phased:step:expectedFinite','PRF');
        end
        cond = ~isreal(prfArg);
        if cond
            coder.internal.errorIf(cond,'phased:stap:NeedReal', 'PRF');
        end
        cond = prfArg <= 0;
        if cond
            coder.internal.errorIf(cond,...
                'phased:step:expectedPositive','PRF');
        end
        prf = prfArg;
    end
end

methods (Static,Hidden,Access=protected)  
    function groups = getPropertyGroupsImpl
        groups = matlab.system.display.Section(...
            'phased.StretchProcessor');
        pPRFSource = matlab.system.display.internal.Property('PRFSource', ...
            'Description', 'Specify PRF as');
        dSampleRate = matlab.system.display.internal.Property(...
            'SampleRate','IsObjectDisplayOnly',true);
        for m = 1:numel(groups.PropertyList)
            if strcmp(groups.PropertyList{m},'SampleRate')
                groups.PropertyList{m} = dSampleRate;
            end
            if strcmp(groups.PropertyList{m},'PRFSource')
                groups.PropertyList{m} = pPRFSource;
            end
        end
    end
    function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:StretchProcessorTitle')),...
          'Text',getString(message('phased:library:block:StretchProcessorDesc')));
    end
end

methods (Access=protected)
    function varargout = getOutputNamesImpl(~)
        varargout = {''};
    end
    
    function varargout = getInputNamesImpl(obj)
        varargout = {'X'};
        if obj.PRFSource(1) == 'I'
            varargout{end+1} = 'PRF';
        end
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Stretch\nProcessor');
    end
    function varargout = getOutputSizeImpl(obj)
        varargout{1} = propagatedInputSize(obj,1);
    end
    function varargout = isOutputFixedSizeImpl(obj)
        varargout{1} = propagatedInputFixedSize(obj, 1);
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj,1);
    end
    function varargout = isOutputComplexImpl(~)
        varargout{1} = true;
    end        
    
end

end

% [EOF]
