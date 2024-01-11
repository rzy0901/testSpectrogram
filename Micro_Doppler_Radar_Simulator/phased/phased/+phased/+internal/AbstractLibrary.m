classdef(Hidden) AbstractLibrary < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime 
%This class is for internal use only. It may be removed in the future.

%   Copyright 2018 The MathWorks, Inc.

%   Reference
%   [1] Cochran et.al., Waveform Libraries, Measures of effectiveness for
%   radar scheduling,IEEE Signal Processing Magazine, 2009.
%
%   [2]Richards, M. A., Fundamentals of Radar Signal Processing. New York:
%   McGraw-Hill, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %SampleRate   Sample rate (Hz)
    %   Specify the sample rate (in Hz) for the waveforms as a positive
    %   scalar. The default value of this property is 1e6 (1 MHz).
    SampleRate = 1e6
end

properties (Nontunable)
    %WaveformSpecification  Specification of each waveform in the library
    %   Specify waveform parameters as a cell array. Each cell defines a
    %   waveform in the library as a cell array and can be specified as 
    %   {WAVEFORMNAME,Name,Value}, where WAVEFORMNAME can be
    %       1.String 
    %       2.Custom function handles
    %
    %  1.String 
    %
    %   If WAVEFORMNAME is string use the following format:
    %   {PULSETYPE,Name,Value}
    %
    %   PULSETYPE, must be one of the following string or character array:
    %   'Rectangular'|'LinearFM'|'SteppedFM'|'PhaseCoded'.
    %
    %   If PULSETYPE is set to 'Rectangular', you can specify additional
    %   name-value pair arguments as:
    %       {...'PRF',PRF} specifies the pulse repetition frequency, PRF,
    %       (in Hz) as a positive scalar. The default value is 1e4.
    %       {...'PulseWidth',PW} or {...,'DutyCycle',DC} specifies the
    %       length of each pulse, PW, (in seconds) as a positive scalar,
    %       default value is 50e-6 or, the duty cycle, DC, of the pulse as
    %       a positive scalar between 0 and 1, default value is 0.5.
    %       {...'FrequencyOffset',FO} specifies a frequency offset, FO, (in
    %       Hz). The default value is 0.
    %
    %   The FrequencyOffset can be used to provide a shift in the frequency
    %   of the generated pulse. This frequency offset lets you to hop the
    %   signal and mitigate effect of jamming.
    %
    %   If PULSETYPE is set to 'LinearFM', you can specify additional
    %   name-value pair arguments as:
    %       {...'PRF',PRF} specifies the pulse repetition frequency, PRF,
    %       (in Hz) as a positive scalar. The default value is 1e4.
    %       {...'PulseWidth',PW} or {...,'DutyCycle',DC} specifies the
    %       length of each pulse, PW, (in seconds) as a positive scalar,
    %       default value is 50e-6 or, the duty cycle, DC, of the pulse as
    %       a positive scalar between 0 and 1, default value is 0.5.
    %       {...'SweepBandwidth',SBW} specifies the bandwidth of the sweep,
    %       SBW, (in Hz) as a positive scalar. The default value is 1e5.
    %       {...'SweepDirection',SDIR} specifies the direction of the
    %       sweep, SDIR, as one of 'Up' | 'Down'. The default value is
    %       'Up'.
    %       {...'SweepInterval',SI} specifies the sweeping interval, SI,
    %       using one of 'Positive' | 'Symmetric'. The default value is
    %       'Positive'.
    %       {...'Envelope',ENV} specifies the envelope function, ENV, as
    %       one of 'Rectangular' | 'Gaussian'. The default value is
    %       'Rectangular'.
    %       {...'FrequencyOffset',FO} specifies a frequency offset, FO, (in
    %       Hz). The default value is 0.
    %
    %   If PULSETYPE is set to 'SteppedFM', you can specify additional
    %   name-value pair arguments as:
    %       {...'PRF',PRF} specifies the pulse repetition frequency, PRF,
    %       (in Hz) as a positive scalar. The default value is 1e4.
    %       {...'PulseWidth',PW} or  {...,'DutyCycle',DC} specifies the
    %       length of each pulse, PW, (in seconds) as a positive scalar,
    %       default value is 50e-6 or, the duty cycle, DC, of the pulse as
    %       a positive scalar between 0 and 1, default value is 0.5.
    %       {...'NumSteps',NSTEPS} specifies the number of frequency steps,
    %       NSTEPS, as a positive integer. The default value is 5.
    %       {...'FrequencyStep',FSTEPS} specifies the linear frequency step
    %       size, FSTEPS, (in Hz) as a positive scalar. The default value
    %       is 2e4.
    %       {...'FrequencyOffset',FO} specifies a frequency offset, FO, (in
    %       Hz). The default value is 0.
    %
    %   If PULSETYPE is set to 'PhaseCoded', you can specify additional
    %   name-value pair arguments as:
    %       {...'PRF',PRF} specifies the pulse repetition frequency, PRF,
    %       (in Hz) as a positive scalar. The default value is 1e4.
    %       {...'Code',C} specifies the type of code, C, used in phase
    %       modulation as one of 'Frank' | 'P1' | 'P2' | 'Px' |
    %       'Zadoff-Chu' | 'P3' | 'P4' | 'Barker'. The default value is
    %       'Frank'.
    %       {...'SequenceIndex',SIDX} specifies the sequence index, SIDX,
    %       used in Zadoff-Chu code as a positive integer. This property is
    %       applicable when you set the Code property to 'Zadoff-Chu'. The
    %       default value is 1.
    %       {...'ChipWidth',CW} specifies the duration of each chip, CW,
    %       (in seconds) as a positive scalar. The default value is 1e-5.
    %       {...'NumChips',NC} specifies the number of chips, NC, as a
    %       positive integer. The default value is 4.
    %       {...'FrequencyOffset',FO} specifies a frequency offset, FO, (in
    %       Hz). The default value is 0.
    %
    %   2. Custom function handles
    %
    %   If WAVEFORMNAME is a CustomWaveform use the following
    %   format: {<function handle>,ARGS} 
    %      
    %   The first input argument of the custom waveform function that
    %   generates samples must be sample rate (in Hz). Input arguments
    %   of the function other than sample rate can be specified in ARGS.
    %   The function must have at least one output argument to return the
    %   samples of each pulse in a column vector.
    %
    %   The default value of WaveformSpecification is a 2-element library
    %   {{'Rectangular','PRF',1e4,'PulseWidth', 50e-6},...
    %    {'LinearFM','PRF',1e4,'PulseWidth',50e-6,'SweepBandwidth',1e5,...
    %    'SweepDirection','Up','SweepInterval','Positive'}}
    WaveformSpecification = {{'Rectangular',...
                    'PRF',1e4,'PulseWidth', 50e-6},...
                    {'LinearFM','PRF',1e4,'PulseWidth',50e-6,....
                    'SweepBandwidth',1e5,'SweepDirection','Up',...
                    'SweepInterval','Positive'}};
end


properties (Access = protected)
    cWaveformSpecification
    pFrequencyOffset        % Property to hold the hop frequency
end

properties (Access = protected, Nontunable)
    pIsFunctionHandle
end

properties (Access = protected)
    pPreviousIndex          % To keep track of index
end

methods (Access = protected)
    % Constructor
    function obj = AbstractLibrary(varargin)
        setProperties(obj, nargin,varargin{:});
    end
    
    function validateDuplicateSetting (~,~)
    % For validation
    end
end

methods
    
    function set.WaveformSpecification(obj,value)
        
        validateattributes(value,{'cell'},{'row'},'','WaveformSpecification');
        validateWaveformProperties(obj,value); % validate WaveformSpecification entries
        validateDuplicateSetting (obj,value);
        
        obj.WaveformSpecification = value;
    end
    
    function set.SampleRate(obj, value)
        validateattributes(value,{'double'}, {'scalar',...
            'positive','finite'},...
            '','SampleRate');
        obj.SampleRate = value;
    end
end

methods (Access = protected)
    function setupImpl(obj,~)
        % Perform one-time calculations, such as computing constants
        coder.extrinsic('phased.internal.AbstractLibrary.createWaveformGenerator');
        
        inp = obj.WaveformSpecification;
        samplerate = obj.SampleRate;
        num = numel(inp);
        isfh = false(1,num);
        
        freqoffset = zeros(1,num);
        wavgenparam = cell(1,num);
        for m = 1:num
            if isa(inp{m}{1},'function_handle')
                isfh(m) = true;
                wavgenparam{m} = inp{m}{1};
            else
                [wavgenparam{m},freqoffset(m)] = coder.const(...
                    @phased.PulseWaveformLibrary.createWaveformGenerator,...
                    inp{m},samplerate);
            end
        end
        
        obj.pIsFunctionHandle = isfh;
        
        obj.cWaveformSpecification = cell(1,num);
        for m = coder.unroll(1:num)
            if isfh(m)
                obj.cWaveformSpecification{m} = wavgenparam{m};
                % dummy, makes cWaveformSpecification homogeneous for codegen
            else
                wavetype = validatestring(inp{m}{1},{'Rectangular',...
                    'LinearFM','SteppedFM','PhaseCoded'},...
                    'PulseWaveformLibrary');
                switch wavetype
                    case 'Rectangular'
                        obj.cWaveformSpecification{m} = ...
                            phased.RectangularWaveform(wavgenparam{m}{:});
                    case 'LinearFM'
                        obj.cWaveformSpecification{m} = ...
                            phased.LinearFMWaveform(wavgenparam{m}{:});
                    case 'SteppedFM'
                        obj.cWaveformSpecification{m} = ...
                            phased.SteppedFMWaveform(wavgenparam{m}{:});
                    case 'PhaseCoded'
                        obj.cWaveformSpecification{m} = ...
                            phased.PhaseCodedWaveform(wavgenparam{m}{:});
                end
            end
        end
        obj.pFrequencyOffset = freqoffset;
        
        for m = coder.unroll(1:num)
            % Error check to make sure SampleRate > Bandwidth of the selected
            % waveform
            if ~isfh(m)
                bw = bandwidth(obj,m);
            else
                bw = samplerate/2;
            end
            cond = ~(samplerate>=bw);
            if cond
                coder.internal.errorIf(cond,...
                    'phased:PulseWaveformLibrary:UnderSample',...
                    'SampleRate','Bandwidth',m,'WaveformSpecification');
            end
        end
    end
    
    function validateDuplicateWarningCond(~,value)
        coder.unroll();
        for idx = 1:numel(value)
            wavetype = value{idx}{1};
            
            if ~isa(wavetype,'function_handle') && ...
                    ~strcmpi(wavetype,'PhaseCoded')
                pw = find(strcmpi(value{idx},'PulseWidth'),1,'last');
                dc = find(strcmpi(value{idx},'DutyCycle'),1,'last');
                
                if ~isempty(pw)&&~isempty(dc) % If both PulseWidth and DutyCycle is specified
                    coder.internal.warning(...
                        'phased:PulseWaveformLibrary:SpecifyEitherPulseWidthOrDutyCycle',...
                        'PulseWidth','DutyCycle',idx,...
                        'WaveformSpecification');
                end
                
            elseif ~isa(wavetype,'function_handle') && ...
                    strcmpi(wavetype,'PhaseCoded')
                
                si_idx = find(strcmpi(value{idx},'SequenceIndex'),1);
                [~,~,~,PhaseCode] = parseInput_phasecoded(value{idx});
                Code = validatestring(PhaseCode,{'Barker','Frank',...
                    'P1','P2','P3','P4','Px','Zadoff-Chu'},...
                    'PulseWaveformLibrary');
                if ~strcmpi(Code,'Zadoff-Chu') && ~isempty(si_idx)
                    coder.internal.warning(....
                        'phased:PulseWaveformLibrary:NonRelevantParameter',....
                        'SequenceIndex',idx,'WaveformSpecification');
                end
            end
        end
    end
    
    function [PRF,ChipWidth,NumChips] = validatePhaseCodedProperties(obj,idx)
        inp =  obj.WaveformSpecification;
        
        [PRF,ChipWidth,NumChips,Code,SequenceIndex] =...
            parseInput_phasecoded(inp{idx});
        
        temp = obj.SampleRate*ChipWidth;
        cond =  abs(temp-round(temp)) > eps(temp);
        if cond
            coder.internal.errorIf(cond, ...
                'phased:PulseWaveformLibrary:NeedProductInteger',...
                'SampleRate', 'ChipWidth',idx,'WaveformSpecification');
        end
        
        cond = (strcmp(Code,'Frank') || strcmp(Code,'P1') || ...
            strcmp(Code,'P2') || strcmp(Code,'Px')) && ...
            rem(sqrt(NumChips),1);
        if cond
            coder.internal.errorIf( cond, ...
                'phased:PulseWaveformLibrary:NeedSquare','NumChips',...
                idx,'WaveformSpecification');
        end
        cond =  strcmp(Code,'P2') && rem(NumChips,2) ;
        if cond
            coder.internal.errorIf(cond, ...
                'phased:PulseWaveformLibrary:NeedEven','NumChips',...
                idx,'WaveformSpecification');
        end
        cond = strcmp(Code,'Zadoff-Chu') && gcd(NumChips,SequenceIndex)~= 1;
        if cond
            coder.internal.errorIf(cond, ...
                'phased:PulseWaveformLibrary:NeedPrime','NumChips',...
                'SequenceIndex',idx,'WaveformSpecification');
        end
        cond =  strcmp(Code,'Barker') && all(NumChips ~= [2 3 4 5 7 11 13]);
        if cond
            coder.internal.errorIf(cond, ...
                'phased:PulseWaveformLibrary:InvalidBarkerLength',...
                'NumChips',idx,'WaveformSpecification');
        end
        
    end
    
    function validatePropertiesImpl(obj)
        coder.extrinsic('sprintf');
        inp =  obj.WaveformSpecification;
        num = numel(inp);
        
        prf = populatePRF(obj);
        coder.unroll;
        for idx = 1:num
            wavetype = inp{idx}{1};
            if ~isa(wavetype,'function_handle')
                
                wavetype = validatestring(wavetype,{'Rectangular',...
                    'LinearFM','SteppedFM','PhaseCoded'},...
                    'PulseWaveformLibrary');
                
                if ~strcmpi(wavetype,'PhaseCoded')
                    pw = find(strcmpi(inp{idx},'PulseWidth'),1,'last');
                    dc = find(strcmpi(inp{idx},'DutyCycle'),1,'last');
                    
                    switch wavetype
                        case 'Rectangular'
                            [~,PulseWidth,DutyCycle] =...
                                parseInput_rect(inp{idx});
                        case 'LinearFM'
                            [~,PulseWidth,DutyCycle] =...
                                parseInput_lfm(inp{idx});
                        case 'SteppedFM'
                            [~,PulseWidth,DutyCycle] =...
                                parseInput_steppedfm(inp{idx});
                        otherwise
                            error('phased:PulseWaveformLibrary:InvalidWaveform',...
                                idx,'WaveformSpecification');
                    end
                    if ~isempty(dc)&&isempty(pw)
                        duration = DutyCycle/prf(idx);
                    else
                        duration = PulseWidth;
                    end
                    cond = any(duration*prf(idx)> 1);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:PulseWaveformLibrary:NotLessThanOrEqualTo',...
                            'PulseWidth', sprintf('%5.2e', 1/prf(idx)),...
                            idx, 'WaveformSpecification');
                    end
                else
                    [~,ChipWidth,NumChips] = ...
                        validatePhaseCodedProperties(obj,idx);
                    duration = ChipWidth*NumChips;
                    cond = any(duration*prf(idx)> 1);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:PulseWaveformLibrary:NotLessThanOrEqualTo', ...
                            'ChipWidth', sprintf( '%5.2e', 1/(prf(idx)*NumChips)),...
                            idx,'WaveformSpecification');
                    end
                end
                % SampleRate and PRF check
                
                val = obj.SampleRate./prf(idx);
                cond =  any(abs(val-round(val))> eps(val));
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:PulseWaveformLibrary:NeedRatioInteger',...
                        'SampleRate', 'PRF', idx, 'WaveformSpecification');
                end
            end
            
        end
    end
    
    function PRF = populatePRF(obj) % For Simulink, to resolve propagator issues.
        
        inp = obj.WaveformSpecification;
        n = numel(inp);
        prf = zeros(1,n);
        
        for idx = coder.unroll(1:numel(inp))
            y = inp{idx};
            wavetype = y{1};
            if isa(wavetype,'function_handle')
                y = getCustomWaveformOutput(obj,idx); % Output of function handle
                wavprf = obj.SampleRate/size(y,1); % PRF from sample rate and size of y
            else
                wavetype = validatestring(wavetype,{'Rectangular',...
                    'LinearFM','SteppedFM','PhaseCoded'},...
                    'PulseWaveformLibrary');
                switch wavetype
                    case 'Rectangular'
                        wavprf = parseInput_rect(inp{idx});
                    case 'LinearFM'
                        wavprf = parseInput_lfm(inp{idx});
                    case 'SteppedFM'
                        wavprf = parseInput_steppedfm(inp{idx});
                    case 'PhaseCoded'
                        wavprf = parseInput_phasecoded(inp{idx});
                    otherwise
                        error('phased:PulseWaveformLibrary:InvalidWaveform',...
                            idx,'WaveformSpecification');
                end
            end
            prf(idx) = wavprf;
        end
        PRF = prf;
    end
    
    function wavestring = validateWaveformString(~,inp,idx)
        
        wavetype = inp{idx}{1};
        
        wav_string = validatestring(wavetype,....
            {'LinearFM','Rectangular','SteppedFM','PhaseCoded'},...
            'PulseWaveformLibrary');
        
        param = {'OutputFormat','NumSamples','PRFSelectionInputPort',...
            'PRFOutputPort','NumPulses','DurationSpecification'};  % List of parameters not valid for WaveformLibrary Configuration
        
        y = inp{idx};
        
        for i = 1:numel(y)
            
            for j = 1:numel(param)
                cond =  any(strcmpi(y{i},param{j}));
                if cond
                    coder.internal.errorIf (cond,...
                        'phased:PulseWaveformLibrary:InvalidConfiguration',...
                        param{j},idx,'WaveformSpecification');
                end
            end
            
            cond = any(strcmpi(y{i},'SampleRate'));
            if cond  % SampleRate shouldn't be specified for individual waveforms to avoid multirate processing
                coder.internal.errorIf(cond,...
                    'phased:PulseWaveformLibrary:NonRelevantSampleRate',...
                    'SampleRate',idx,'WaveformSpecification');
            end
            
        end
        
        wavestring = wav_string;
    end
    
    function validateWaveformProperties(obj,inp)
        
        num = numel(inp);
        for idx = coder.unroll(1:num)
            cond = ~iscell(inp{idx}) ;
            if cond
                coder.internal.errorIf(cond,...
                    'phased:PulseWaveformLibrary:InvalidDatatype',...
                    idx,'WaveformSpecification')
            end
            
            cond = ~(isa(inp{idx}{1},'function_handle')|| ...
                ischar(inp{idx}{1}));
            
            if cond
                coder.internal.errorIf(cond,...
                    'phased:PulseWaveformLibrary:InvalidWaveform',...
                    idx,'WaveformSpecification');
            end
            
            wavetype = inp{idx}{1}; % String/FunctionHandle
            
            if ~isa(wavetype,'function_handle')
                wav_string = validateWaveformString(obj,inp,idx);
                
                y = inp{idx};
                
                switch wav_string
                    
                    case 'Rectangular'
                        
                        % List of parameters not valid for RectangularWaveform
                        param = {'SweepBandwidth','SweepInterval',.....
                            'SweepDirection','Envelope',....
                            'FrequencyStep','NumSteps','NumChips',....
                            'ChipWidth','Code','SequenceIndex'};
                        for i = 1:numel(y)
                            for j = 1:numel(param)
                                cond =  any(strcmpi(y{i},param{j}));
                                if cond
                                    coder.internal.errorIf(cond,...
                                        'phased:PulseWaveformLibrary:InvalidParameter',....
                                        param{j},idx,'WaveformSpecification');
                                end
                            end
                        end
                        
                        [PRF,PulseWidth,DutyCycle,FrequencyOffset] =...
                            parseInput_rect(y);
                        validateattributes(PRF,{'double'},...
                            {'scalar', 'positive', 'finite' },'',...
                            'PRF',idx);
                        
                        validateattributes(PulseWidth,{'double'},...
                            {'scalar','positive','finite'},'',...
                            'PulseWidth',idx);
                        
                        validateattributes(DutyCycle,{'double'},...
                            {'scalar','positive','>',0,'<',1},'',...
                            'DutyCycle',idx);
                        
                        validateattributes(FrequencyOffset,{'double'},...
                            {'scalar','finite'},'',...
                            'FrequencyOffset',idx);
                        
                    case 'LinearFM'
                        
                        % List of parameters not valid for LinearFMWaveform
                        param = {'FrequencyStep','NumSteps','NumChips',....
                            'ChipWidth','Code','SequenceIndex'};
                        
                        for i = 1:numel(y)
                            for j = 1:numel(param)
                                cond =  any(strcmpi(y{i},param{j}));
                                if cond
                                    coder.internal.errorIf (cond,....
                                        'phased:PulseWaveformLibrary:InvalidParameter',...
                                        param{j},idx,'WaveformSpecification');
                                end
                            end
                        end
                        
                        [PRF,PulseWidth,DutyCycle,SweepBandwidth,...
                            SweepInterval,Envelope,SweepDirection,...
                            FrequencyOffset] = parseInput_lfm(y);
                        
                        validateattributes(PRF, {'double'},...
                            {'scalar', 'positive', 'finite' },'',...
                            'PRF',idx);
                        
                        validateattributes(PulseWidth,{'double'},...
                            {'scalar','positive','finite'},'',...
                            'PulseWidth',idx);
                        
                        validateattributes(DutyCycle,{'double'},...
                            {'scalar','positive','>',0,'<',1},'',...
                            'DutyCycle',idx);
                        
                        validateattributes(SweepBandwidth,{'double'},...
                            {'scalar','positive','finite'},'',...
                            'SweepBandwidth',idx);
                        
                        validatestring(SweepInterval,...
                            {'Positive','Symmetric'},'',...
                            'SweepInterval',idx);
                        
                        validatestring(SweepDirection,...
                            {'Up','Down'},'','SweepDirection',idx);
                        
                        validatestring(Envelope,....
                            {'Rectangular','Gaussian'},'','Envelope',idx);
                        
                        validateattributes(FrequencyOffset,{'double'},...
                            {'scalar','finite'},'',...
                            'FrequencyOffset',idx);
                        
                    case 'SteppedFM'
                        
                        % List of parameters not valid for SteppedFMWaveform
                        param = {'SweepBandwidth','SweepInterval',...
                            'SweepDirection','Envelope',....
                            'NumChips','ChipWidth','Code','SequenceIndex'};
                        
                        for i = 1:numel(y)
                            for j = 1:numel(param)
                                cond =  any(strcmpi(y{i},param{j}));
                                if cond
                                    coder.internal.errorIf (cond,...
                                        'phased:PulseWaveformLibrary:InvalidParameter',...
                                        param{j},idx,'WaveformSpecification');
                                end
                            end
                        end
                        
                        [PRF,PulseWidth,DutyCycle,FrequencyStep,...
                            NumSteps,FrequencyOffset] = ...
                            parseInput_steppedfm(y);
                        
                        validateattributes(PRF, {'double'},...
                            {'scalar', 'positive', 'finite' },'',...
                            'PRF',idx);
                        
                        validateattributes(PulseWidth,{'double'},...
                            {'scalar','positive','finite'},'',...
                            'PulseWidth',idx);
                        
                        validateattributes(DutyCycle,{'double'},...
                            {'scalar','positive','>',0,'<',1},'',...
                            'DutyCycle',idx);
                        
                        validateattributes(NumSteps,{'double'},...
                            {'scalar','positive','finite'},'',...
                            'NumSteps',idx);
                        
                        validateattributes(FrequencyStep,{'double'},...
                            {'scalar', 'positive', 'finite' },'',...
                            'FrequencyStep',idx);
                        
                        validateattributes(FrequencyOffset,{'double'},...
                            {'scalar','finite'},'',...
                            'FrequencyOffset',idx);
                        
                    case 'PhaseCoded'
                        
                        % List of parameters not valid for PhaseCodedWaveform
                        param = {'SweepBandwidth','SweepInterval',...
                            'SweepDirection','Envelope',....
                            'FrequencyStep','NumSteps'};
                        
                        for i = 1:numel(y)
                            for j = 1:numel(param)
                                cond =  any(strcmpi(y{i},param{j}));
                                if cond
                                    coder.internal.errorIf (cond,....
                                        'phased:PulseWaveformLibrary:InvalidParameter',...
                                        param{j},idx,'WaveformSpecification');
                                end
                            end
                        end
                        
                        [PRF,ChipWidth,NumChips,Code,SequenceIndex,...
                            FrequencyOffset] = ...
                            parseInput_phasecoded(y);
                        
                        validateattributes(PRF, {'double' },...
                            {'scalar', 'positive', 'finite'},'',...
                            'PRF',idx);
                        
                        validateattributes(ChipWidth,{'double'},...
                            {'scalar','finite','nonnan','nonempty',...
                            'positive'},'','ChipWidth',idx);
                        
                        validateattributes(NumChips,{'double'},...
                            {'scalar', 'positive', 'finite'},'',...
                            'NumChips',idx);
                        
                        Code = validatestring(Code,...
                            {'Barker','Frank','P1','P2','P3','P4','Px',...
                            'Zadoff-Chu'},'','Code',idx);
                        
                        validateattributes(FrequencyOffset,{'double'},...
                            {'scalar','finite'},'','FrequencyOffset',idx);
                        
                        if strcmpi(Code,'Zadoff-Chu')
                            validateattributes(SequenceIndex,...
                                {'double'},{'scalar', 'positive',...
                                'finite'},'','SequenceIndex',idx);
                        end
                end
            end
        end
    end
    
    function y = getCustomWaveformOutput(obj,idx)   % To get custom waveform output
        
        inp = obj.WaveformSpecification{idx};
        fh = inp{1};
        yout = fh(obj.SampleRate,inp{2:end});  % If function has SampleRate as required input and 'n' optional inputs (SampleRate, varargin)
        
        if isvector(yout)
            y = yout(:);
        else
            coder.internal.assert(false, ...
                'phased:PulseWaveformLibrary:MustBeVector',idx,...
                'WaveformSpecification');
        end
    end
    
    function resetImpl(obj)
        for idx = 1:length(obj.cWaveformSpecification)
            if ~isa(obj.WaveformSpecification{idx}{1},'function_handle')
                reset(obj.cWaveformSpecification{idx});
            end
        end
        obj.pPreviousIndex = 0;
    end
    
    function releaseImpl(obj)
        for idx = 1:length(obj.cWaveformSpecification)
            if ~isa(obj.WaveformSpecification{idx}{1},'function_handle')
                release(obj.cWaveformSpecification{idx});
            end
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        s.isLocked = isLocked(obj);
        if isLocked(obj)
            s.pFrequencyOffset = obj.pFrequencyOffset;
            s.pPreviousIndex = obj.pPreviousIndex;
            s.pIsFunctionHandle = obj.pIsFunctionHandle;
            for idx=1:length(obj.WaveformSpecification)
                if ~isa(obj.WaveformSpecification{idx}{1},'function_handle')
                    s.cWaveformSpecification{idx} = ...
                        saveobj(obj.cWaveformSpecification{idx});
                else
                    s.cWaveformSpecification{idx} = ...
                        obj.cWaveformSpecification{idx};
                end
            end
        end
    end
    
    function s = loadSubObjects(obj,s)
        if isfield(s,'isLocked')
            if s.isLocked
                for idx= 1:length(s.WaveformSpecification)
                    if ~isa(s.WaveformSpecification{idx}{1},...
                            'function_handle')
                        obj.cWaveformSpecification{idx} = ...
                            phased.internal.AbstractPulseWaveform.loadobj...
                            (s.cWaveformSpecification{idx});
                    else
                        obj.cWaveformSpecification{idx}= ...
                            s.cWaveformSpecification{idx};
                    end
                end
                s = rmfield(s,'cWaveformSpecification');
            end
            s = rmfield(s,'isLocked');
        end
    end
    
    function loadObjectImpl(obj,s,~)
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end
    
end

methods(Hidden,Access = protected)
    
    function bw = bandwidth(obj,idx)
        narginchk(2,2)
        
        validateattributes(idx,{'double'},{'scalar', 'positive','finite'},...
            '','IDX');
        
        if ~obj.pIsFunctionHandle(idx)
            bw = bandwidth(obj.cWaveformSpecification{idx});
        else
            bw = 0;
        end
    end
    
end

methods (Hidden, Static)
    function [tempWaveformParam,offset] = createWaveformGenerator(inp,fs)
        wavetype = inp{1};
        
        pwidx = find(strcmpi(inp,'PulseWidth'),1);
        dcidx = find(strcmpi(inp,'DutyCycle'),1);
        
        wav_string = validatestring(wavetype,{'Rectangular',...
            'LinearFM','SteppedFM','PhaseCoded'},...
            'PulseWaveformLibrary');
        switch wav_string
            case 'Rectangular'
                [PRF,PulseWidth,DutyCycle,FrequencyOffset] = ...
                    parseInput_rect(inp);
                
                offset = FrequencyOffset;
                
                if isempty(pwidx) && ~isempty(dcidx)            % Case DurationSpecification - 'DutyCycle'
                    tempWaveformParam = {'PRF',PRF,...
                        'DurationSpecification','Duty cycle',...
                        'DutyCycle',DutyCycle,'SampleRate',fs};
                else                                            % Case DurationSpecification - 'PulseWidth'
                    tempWaveformParam = {'PRF',PRF,...
                        'PulseWidth',PulseWidth,'SampleRate',fs};
                end
            case 'LinearFM'
                
                [PRF,PulseWidth,DutyCycle,SweepBandwidth,...
                    SweepInterval,Envelope,SweepDirection,...
                    FrequencyOffset] = ...
                    parseInput_lfm(inp);
                
                offset = FrequencyOffset;
                
                SweepDirection = validatestring(SweepDirection,{'Up',...
                    'Down'},'PulseWaveformLibrary');
                SweepInterval = validatestring(SweepInterval,{'Positive',...
                    'Symmetric'},'PulseWaveformLibrary');
                Envelope = validatestring(Envelope,{'Rectangular',...
                    'Gaussian'},'PulseWaveformLibrary');
                
                if isempty(pwidx) && ~isempty(dcidx)
                    tempWaveformParam = {...
                        'PRF',PRF,...
                        'DurationSpecification','Duty cycle',...
                        'DutyCycle',DutyCycle,...
                        'SampleRate',fs,...
                        'SweepBandwidth',SweepBandwidth,...
                        'SweepDirection', SweepDirection,...
                        'Envelope',Envelope,...
                        'SweepInterval',SweepInterval};
                else
                    tempWaveformParam = {...
                        'PRF',PRF,...
                        'PulseWidth',PulseWidth,...
                        'SampleRate',fs,...
                        'SweepBandwidth',SweepBandwidth,...
                        'SweepDirection', SweepDirection,...
                        'Envelope',Envelope,...
                        'SweepInterval',SweepInterval};
                end
            case 'SteppedFM'
                
                [PRF,PulseWidth,DutyCycle,FrequencyStep,...
                    NumSteps,FrequencyOffset] = ...
                    parseInput_steppedfm(inp);
                
                offset = FrequencyOffset;
                
                if isempty(pwidx) && ~isempty(dcidx)
                    tempWaveformParam = {...
                        'PRF',PRF,...
                        'DurationSpecification','Duty cycle',...
                        'DutyCycle',DutyCycle,...
                        'SampleRate',fs,...
                        'FrequencyStep',FrequencyStep,...
                        'NumSteps',NumSteps};
                else
                    tempWaveformParam = {...
                        'PRF',PRF,...
                        'DurationSpecification','Pulse width',...
                        'PulseWidth',PulseWidth,...
                        'SampleRate',fs,...
                        'FrequencyStep',FrequencyStep,...
                        'NumSteps',NumSteps};
                end
            case 'PhaseCoded'
                
                [PRF,ChipWidth,NumChips,Code,SequenceIndex,...
                    FrequencyOffset] = ...
                    parseInput_phasecoded(inp);
                
                offset = FrequencyOffset;
                
                Code = validatestring(Code,....
                    {'Barker','Frank','P1','P2','P3','P4','Px',...
                    'Zadoff-Chu'},'PulseWaveformLibrary');
                
                if strcmpi(Code,'Zadoff-Chu')
                    tempWaveformParam = {...
                        'PRF',PRF,...
                        'ChipWidth',ChipWidth,...
                        'Code',Code,...
                        'SequenceIndex',SequenceIndex,...
                        'NumChips',NumChips,...
                        'SampleRate',fs};
                else
                    tempWaveformParam =  {...
                        'PRF',PRF,...
                        'ChipWidth',ChipWidth,...
                        'Code',Code,...
                        'NumChips',NumChips,...
                        'SampleRate',fs};
                end
        end
    end
end

end

function [PRF,PulseWidth,DutyCycle,FrequencyOffset] = ...
    parseInput_rect(varargin)

    defaultPRF = 1e4;
    defaultPulseWidth = 50e-6;
    defaultDutyCycle = 0.5;
    defaultFrequencyOffset = 0;

    if isempty(coder.target)
        p = inputParser;
        p.FunctionName = 'PulseWaveformLibrary';
        p.addParameter('PRF',defaultPRF);
        p.addParameter('PulseWidth', defaultPulseWidth);
        p.addParameter('DutyCycle', defaultDutyCycle);
        p.addParameter('FrequencyOffset', defaultFrequencyOffset);
        p.parse(varargin{1}{2:end});
        PRF = p.Results.PRF;
        PulseWidth = p.Results.PulseWidth;
        DutyCycle = p.Results.DutyCycle;
        FrequencyOffset = p.Results.FrequencyOffset;
    else
        parms = struct( ...
            'PRF',uint32(0), ...
            'PulseWidth',uint32(0), ...
            'DutyCycle',uint32(0),...
            'FrequencyOffset',uint32(0));
        % Select parsing options.
        poptions = struct( ...
            'CaseSensitivity',false, ...
            'PartialMatching','unique', ...
            'StructExpand',false, ...
            'IgnoreNulls',false);
        pstruct = coder.internal.parseParameterInputs(...
            parms,poptions,varargin{1}{2:end});
        PRF = coder.internal.getParameterValue(...
            pstruct.PRF,defaultPRF,varargin{1}{2:end});
        PulseWidth = coder.internal.getParameterValue(...
            pstruct.PulseWidth,defaultPulseWidth,varargin{1}{2:end});
        DutyCycle = coder.internal.getParameterValue(...
            pstruct.DutyCycle,defaultDutyCycle,varargin{1}{2:end});
        FrequencyOffset = coder.internal.getParameterValue(...
            pstruct.FrequencyOffset,defaultFrequencyOffset,...
            varargin{1}{2:end});
    end

end

function [PRF,PulseWidth,DutyCycle,SweepBandwidth,SweepInterval,...
    Envelope,SweepDirection,...
    FrequencyOffset] = parseInput_lfm(varargin)

    defaultPRF = 1e4;
    defaultPulseWidth = 50e-6;
    defaultDutyCycle = 0.5;
    defaultFrequencyOffset = 0;
    defaultSweepBandwidth = 1e5;
    defaultSweepDirection = 'Up';
    defaultSweepInterval = 'Positive';
    defaultEnvelope = 'Rectangular';

    if isempty(coder.target)
        p = inputParser;
        p.FunctionName = 'PulseWaveformLibrary';
        p.addParameter('PRF',defaultPRF);
        p.addParameter('PulseWidth', defaultPulseWidth);
        p.addParameter('DutyCycle', defaultDutyCycle);
        p.addParameter('FrequencyOffset', defaultFrequencyOffset);
        p.addParameter('SweepBandwidth',defaultSweepBandwidth);
        p.addParameter('SweepDirection', defaultSweepDirection);
        p.addParameter('SweepInterval', defaultSweepInterval);
        p.addParameter('Envelope', defaultEnvelope);
        p.parse(varargin{1}{2:end});
        PRF = p.Results.PRF;
        PulseWidth = p.Results.PulseWidth;
        DutyCycle = p.Results.DutyCycle;
        FrequencyOffset = p.Results.FrequencyOffset;
        SweepBandwidth = p.Results.SweepBandwidth;
        SweepDirection = p.Results.SweepDirection;
        SweepInterval = p.Results.SweepInterval;
        Envelope = p.Results.Envelope;
    else
        parms = struct( ...
            'PRF',uint32(0), ...
            'PulseWidth',uint32(0), ...
            'DutyCycle',uint32(0),...
            'FrequencyOffset',uint32(0),...
            'SweepBandwidth',uint32(0),...
            'SweepDirection', uint32(0),...
            'SweepInterval', uint32(0),...
            'Envelope', uint32(0));
        % Select parsing options.
        poptions = struct( ...
            'CaseSensitivity',false, ...
            'PartialMatching','unique', ...
            'StructExpand',false, ...
            'IgnoreNulls',false);
        pstruct = coder.internal.parseParameterInputs(...
            parms,poptions,varargin{1}{2:end});
        PRF = coder.internal.getParameterValue(...
            pstruct.PRF,defaultPRF,varargin{1}{2:end});
        PulseWidth = coder.internal.getParameterValue(...
            pstruct.PulseWidth,defaultPulseWidth,varargin{1}{2:end});
        DutyCycle = coder.internal.getParameterValue(...
            pstruct.DutyCycle,defaultDutyCycle,varargin{1}{2:end});
        FrequencyOffset = coder.internal.getParameterValue(...
            pstruct.FrequencyOffset,defaultFrequencyOffset,...
            varargin{1}{2:end});
        SweepBandwidth = coder.internal.getParameterValue(...
            pstruct.SweepBandwidth,defaultSweepBandwidth,...
            varargin{1}{2:end});
        SweepDirection = coder.internal.getParameterValue(...
            pstruct.SweepDirection,defaultSweepDirection,...
            varargin{1}{2:end});
        SweepInterval  = coder.internal.getParameterValue(...
            pstruct.SweepInterval,defaultSweepInterval,...
            varargin{1}{2:end});
        Envelope = coder.internal.getParameterValue(...
            pstruct.Envelope,defaultEnvelope,varargin{1}{2:end});
    end

end

function [PRF,PulseWidth,DutyCycle,FrequencyStep,NumSteps,...
    FrequencyOffset] = parseInput_steppedfm(varargin)

    defaultPRF = 1e4;
    defaultPulseWidth = 50e-6;
    defaultDutyCycle = 0.5;
    defaultFrequencyOffset = 0;
    defaultFrequencyStep = 2e4;
    defaultNumSteps = 5;

    if isempty(coder.target)
        p = inputParser;
        p.FunctionName = 'PulseWaveformLibrary';
        p.addParameter('PRF',defaultPRF);
        p.addParameter('PulseWidth', defaultPulseWidth);
        p.addParameter('DutyCycle', defaultDutyCycle);
        p.addParameter('FrequencyOffset', defaultFrequencyOffset);
        p.addParameter('FrequencyStep',defaultFrequencyStep);
        p.addParameter('NumSteps', defaultNumSteps);
        p.parse(varargin{1}{2:end});
        PRF = p.Results.PRF;
        PulseWidth = p.Results.PulseWidth;
        DutyCycle = p.Results.DutyCycle;
        FrequencyOffset = p.Results.FrequencyOffset;
        FrequencyStep = p.Results.FrequencyStep;
        NumSteps = p.Results.NumSteps;
    else
        parms = struct( ...
            'PRF',uint32(0), ...
            'PulseWidth',uint32(0), ...
            'DutyCycle',uint32(0),...
            'FrequencyOffset',uint32(0),...
            'FrequencyStep',uint32(0),...
            'NumSteps', uint32(0));
        % Select parsing options.
        poptions = struct( ...
            'CaseSensitivity',false, ...
            'PartialMatching','unique', ...
            'StructExpand',false, ...
            'IgnoreNulls',false);
        pstruct = coder.internal.parseParameterInputs(...
            parms,poptions,varargin{1}{2:end});
        PRF = coder.internal.getParameterValue(...
            pstruct.PRF,defaultPRF,varargin{1}{2:end});
        PulseWidth = coder.internal.getParameterValue(...
            pstruct.PulseWidth,defaultPulseWidth,varargin{1}{2:end});
        DutyCycle = coder.internal.getParameterValue(...
            pstruct.DutyCycle,defaultDutyCycle,varargin{1}{2:end});
        FrequencyOffset = coder.internal.getParameterValue(...
            pstruct.FrequencyOffset,defaultFrequencyOffset,...
            varargin{1}{2:end});
        FrequencyStep = coder.internal.getParameterValue(...
            pstruct.FrequencyStep,defaultFrequencyStep,...
            varargin{1}{2:end});
        NumSteps = coder.internal.getParameterValue(...
            pstruct.NumSteps,defaultNumSteps,varargin{1}{2:end});
    end

end

function [PRF,ChipWidth,NumChips,Code,SequenceIndex,...
    FrequencyOffset] = parseInput_phasecoded(varargin)

    defaultPRF = 1e4;
    defaultChipWidth = 1e-5;
    defaultCode = 'Frank';
    defaultFrequencyOffset = 0;
    defaultSequenceIndex = 1;
    defaultNumChips = 4;

    if isempty(coder.target)
        p = inputParser;
        p.FunctionName = 'PulseWaveformLibrary';
        p.addParameter('PRF',defaultPRF);
        p.addParameter('ChipWidth', defaultChipWidth);
        p.addParameter('Code', defaultCode);
        p.addParameter('FrequencyOffset', defaultFrequencyOffset);
        p.addParameter('SequenceIndex',defaultSequenceIndex);
        p.addParameter('NumChips',defaultNumChips);
        p.parse(varargin{1}{2:end});
        PRF = p.Results.PRF;
        ChipWidth = p.Results.ChipWidth;
        NumChips = p.Results.NumChips;
        FrequencyOffset = p.Results.FrequencyOffset;
        Code = p.Results.Code;
        Code = validatestring(Code,{'Barker','Frank','P1','P2','P3',...
            'P4','Px','Zadoff-Chu'},'PulseWaveformLibrary');
        if strcmpi(Code,'Zadoff-Chu')
            SequenceIndex = p.Results.SequenceIndex;
        else
            SequenceIndex = [];
        end

    else
        parms = struct( ...
            'PRF',uint32(0), ...
            'ChipWidth',uint32(0), ...
            'Code',uint32(0),...
            'FrequencyOffset',uint32(0),...
            'SequenceIndex',uint32(0),...
            'NumChips', uint32(0));
        % Select parsing options.
        poptions = struct( ...
            'CaseSensitivity',false, ...
            'PartialMatching','unique', ...
            'StructExpand',false, ...
            'IgnoreNulls',false);
        pstruct = coder.internal.parseParameterInputs(...
            parms,poptions,varargin{1}{2:end});
        PRF = coder.internal.getParameterValue(...
            pstruct.PRF,defaultPRF,varargin{1}{2:end});
        ChipWidth = coder.internal.getParameterValue(...
            pstruct.ChipWidth,defaultChipWidth,varargin{1}{2:end});
        NumChips = coder.internal.getParameterValue(...
            pstruct.NumChips,defaultNumChips,varargin{1}{2:end});
        FrequencyOffset = coder.internal.getParameterValue(...
            pstruct.FrequencyOffset,defaultFrequencyOffset,...
            varargin{1}{2:end});
        Code = coder.internal.getParameterValue(...
            pstruct.Code,defaultCode,varargin{1}{2:end});
        Code = validatestring(Code,{'Barker','Frank','P1','P2','P3',...
            'P4','Px','Zadoff-Chu'},'PulseWaveformLibrary');
        if strcmpi(Code,'Zadoff-Chu')
            SequenceIndex  = coder.internal.getParameterValue(...
                pstruct.SequenceIndex,defaultSequenceIndex,...
                varargin{1}{2:end});
        else
            SequenceIndex  = [];
        end
    end

end
