classdef (Sealed,StrictDefaults)PulseCompressionLibrary < ...
        matlab.system.mixin.CustomIcon &...
        phased.internal.AbstractLibrary 
    
%PulseCompressionLibrary   Library of range processing techniques.
%   H = phased.PulseCompressionLibrary creates a pulse compression library
%   System object, H. This object performs match filtering or stretch
%   processing on the input data based on a predefined list and returns
%   the range processed data.
%
%   H = phased.PulseCompressionLibrary(Name,Value) creates a pulse
%   compression library object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   [Y,RNG] = step(H,X,IDX) performs pulse compression on the input X using
%   the IDX-th entry of the list and returns the range processed result in
%   Y. IDX must be a positive integer.
%
%   X is an input signal. X must be a KxL matrix, a KxN matrix, or
%   a KxNxL array where K denotes the number of fast time samples, L is
%   the number of pulses, and N is the number of channels (antenna elements
%   or beams).  
%
%   Y is either an MxL matrix, an MxN matrix, or an MxNxL array containing
%   the complex range response of the input X, where the number of
%   dimensions in Y will match the number of dimensions in X. M is
%   determined by either the number of rows in X when matched filtering is
%   performed and RangeFFTLength parameter of StretchProcessor is not
%   specified, or value specified in the RangeFFTLength parameter of
%   StretchProcessor.
%
%   RNG is a length-M column vector containing the range samples at which
%   the range response is evaluated.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   PulseCompressionLibrary methods:
%
%   step         - Returns range processed data using IDX-th pulse 
%                  compression
%   release      - Allows property value and input characteristics changes
%   clone        - Creates a pulse compression library object with same
%                  property values
%   isLocked     - Locked status (logical)
%   plotResponse - Plots the range response
%
%   PulseCompressionLibrary properties:
%
%   SampleRate                  - Sample rate 
%   PropagationSpeed            - Propagation speed
%   WaveformSpecification       - Specify the type of pulsed waveform and
%                                 its parameters                         
%   ProcessingSpecification     - Specify the type of processing and its
%                                 parameters
%
%   % Example:
%   %   Use PulseCompressionLibrary to range process the waveforms in 
%   %   PulseWaveformLibrary using different range processing methods.
%
%   % PulseWaveformLibrary
%   waveform1 = {'Rectangular','PRF',1e4, 'PulseWidth', 50e-6};
%   waveform2 = {'LinearFM','PRF',1e4,'PulseWidth',50e-6,...
%       'SweepBandwidth',1e5,'SweepDirection','Up',...
%       'SweepInterval', 'Positive'};
%   wav = phased.PulseWaveformLibrary('WaveformSpecification',...
%       {waveform1,waveform2},'SampleRate',1e6);
%   rect = wav(1);
%   lfm = wav(2);
% 
%   % perform range processing
%   process1 = {'MatchedFilter','Coefficients',getMatchedFilter(wav,1)};
%   process2 = {'StretchProcessor','ReferenceRange',5000,...
%        'RangeSpan',200,'RangeWindow','Hamming'};
% 
%   ProcessLib = phased.PulseCompressionLibrary('WaveformSpecification',...
%       {waveform1,waveform2},'ProcessingSpecification',...
%       {process1,process2},'SampleRate',1e6);
% 
%   out_rect = ProcessLib(rect,1);
%   out_lfm = ProcessLib(lfm,2);
%   
%   See also phased,phased.RangeResponse, phased.RangeDopplerResponse,
%   phased.MatchedFilter,phased.StretchProcessor.

%   Copyright 2018 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %ProcessingSpecification Pulse compression specifications
    %   Specify pulse processing parameters as a cell array. Each cell
    %   defines a pulse processing technique as a cell array and can be
    %   specified as {PROCESSTYPE,Name,Value}.
    %
    %   PROCESSTYPE, must be one of the following string or character
    %   array: 'MatchedFilter'|'StretchProcessor'.
    %
    %   If PROCESSTYPE is set to 'MatchedFilter', you can specify
    %   additional name-value pair arguments as:
    %       {...'Coefficients',COEFF} specifies matched filter coefficient,
    %       COEFF, as a column vector. If not specified COEFF is calculated
    %       using the WaveformSpecification property. In case of Stepped FM
    %       waveform which contains multiple pulses, COEFF corresponds to
    %       each pulse until IDX changes, where IDX is the index input.
    %       {...'SpectrumWindow',SW} specifies SW, the window used for
    %       spectrum weighting as one of 'None'| 'Hamming' | 'Chebyshev' |
    %       'Hann' | 'Kaiser' | 'Taylor'. The default value is 'None'.
    %       {...'SidelobeAttenuation',SLB} specifies the sidelobe
    %       attenuation level, SLB, (in dB) of a Chebyshev or Taylor window
    %       as a positive scalar. The default value is 30. This parameter
    %       applies when you set the SpectrumWindow to 'Chebyshev' or
    %       'Taylor'.
    %       {...'Beta',Beta} specifies the parameter, Beta, that affects
    %       the Kaiser window sidelobe attenuation as a nonnegative scalar.
    %       The default value is 0.5. This parameter applies when you set
    %       the SpectrumWindow to 'Kaiser'.
    %       {...'Nbar',Nbar} specifies the number of nearly constant level
    %       sidelobes, Nbar, adjacent to the mainlobe in a Taylor window as
    %       a positive integer. The default value is 4. This parameter
    %       applies when you set the SpectrumWindow to 'Taylor'.
    %       {...'SpectrumRange',SR} specifies the spectrum region, SR, on
    %       which the spectrum window is applied as a 1x2 vector in the
    %       form of [StartFrequency EndFrequency] (in Hz). The default
    %       value is [0 1e5]. This parameter applies when you set the
    %       SpectrumWindow to a value other than 'None'.
    %   
    %       Note that both StartFrequency and EndFrequency are measured in
    %       baseband, i.e., within [-Fs/2 Fs/2] where Fs is the sample rate
    %       that you specify in the SampleRate property. StartFrequency
    %       cannot be larger than EndFrequency.
    %
    %   If PROCESSTYPE is set to 'StretchProcessor', you can specify
    %   additional name-value pair arguments as:
    %       {...'ReferenceRange',REFRNG} specifies the center of ranges of
    %       interest, REFRNG, (in meters) as a positive scalar. The
    %       ReferenceRange must be within the unambiguous range of one
    %       pulse. The default value is 5000.
    %       {...'RangeSpan',RNGSPAN} specifies the length of the interval
    %       for ranges of interest, RNGSPAN, (in meters) as a positive
    %       scalar. The range span is centered at the range value specified
    %       in the ReferenceRange parameter. The default value is 500.
    %       {...'RangeFFTLength',LEN} specifies the FFT length, LEN, in the
    %       range domain as a positive integer. If not specified default
    %       value is same as length of input data.
    %       {...'RangeWindow',RW} specifies window used for range
    %       processing, RW, as one of 'None'| 'Hamming' | 'Chebyshev' |
    %       'Hann' | 'Kaiser' | 'Taylor'. The default value is 'None'.
    %       {...'SidelobeAttenuation',SLB} specifies the sidelobe
    %       attenuation level, SLB, in dB) of a Chebyshev or Taylor window
    %       as a positive scalar. The default value is 30. This parameter
    %       applies when you set the RangeWindow to 'Chebyshev' or
    %       'Taylor'.
    %       {...'Beta',Beta} specifies the parameter, Beta, that affects
    %       the Kaiser window sidelobe attenuation as a nonnegative scalar.
    %       The default value is 0.5. This parameter applies when you set
    %       the RangeWindow to 'Kaiser'.
    %       {...'Nbar',Nbar} specifies the number of nearly constant level
    %       sidelobes, Nbar, adjacent to the mainlobe in a Taylor window as
    %       a positive integer. The default value is 4. This parameter
    %       applies when you set the rangeWindow to 'Taylor'.
    %
    %   The default value of ProcessingSpecification is a 2-element
    %   library {{'MatchedFilter','SpectrumWindow','None'},...
    %   {'StretchProcessor','RangeSpan',200,'ReferenceRange',5000,...
    %    'RangeWindow','None'}}
    ProcessingSpecification = {{'MatchedFilter',...
        'SpectrumWindow','None'},{'StretchProcessor',...
        'RangeSpan',200,'ReferenceRange',5000,...
        'RangeWindow','None'}};
    
    %PropagationSpeed   Propagation speed (m/s)
    %   Specify the propagation speed (in m/s) of the signal as a scalar.
    %   The default value of this property is the speed of light.
    PropagationSpeed = physconst('LightSpeed');
  
end

properties (Access = protected)
    pFreqOffset                 % frequency offset
    pCubeDim = [-1 -1 -1]       % Cube dimension
    pRngWinCoeff                % RangeWinCoeff
    pRngFFTLength               % Range FFT length 
end

properties(Access = private)     
     cProcessingSpecification   % ProcessingSpecification    
     pRngGrid
     pRngOffset     
     pIndex
     pPulseIdx
end

properties (Nontunable, Logical)
    %SampleRateFromInputCheckbox Inherit sample rate
    %   Set SampleRateFromInputCheckbox to true to derive sample rate
    %   from Simulink time engine. Set SampleRateFromInputCheckbox to
    %   false to specify the sample rate. This property applies when
    %   used in Simulink.
    SampleRateFromInputCheckbox = true
end

properties (Constant, Hidden)
    SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
        'SystemBlock', 'SampleRateFromInputCheckbox',...
        'getSampleRateInSimulation',false})
end

properties (Access = protected, Logical, Nontunable)
    pInput3DFlag    % indicate whether the input is a cube or not
end

methods
    function set.ProcessingSpecification(obj,value)
        validateattributes(value,{'cell'},{'row'},'','PulseCompressionSpecification');
        
        validateProcessingparams(obj,value); 
        
        obj.ProcessingSpecification = value;
    end
    
    function set.PropagationSpeed(obj,value)
        sigdatatypes.validateSpeed(value,...
            '','PropagationSpeed',...
            {'scalar','positive'});
        obj.PropagationSpeed = value;
    end

end

methods
    % Constructor
    
    function obj = PulseCompressionLibrary(varargin)
        obj@phased.internal.AbstractLibrary(varargin{:});
    end
    
    function varargout = plotResponse(obj,x,idx,varargin)
    %plotResponse   Plot range response
    %   plotResponse(H,X,IDX) plots the range response of the input signal,
    %   X, for the IDX-th processing in dB scale.
    %
    %   X must be a matrix where each column contains the signal from one
    %   pulse. IDX must be a positive scalar.
    %
    %   plotResponse(...,PIDX) specifies the index of the pulse to plot.
    %   PIDX must be a scalar and its default value is 1.
    %
    %   plotResponse(...,Name,Value) plots the range response with the
    %   specified parameter Name set to the specified value.
    %   The parameter Names are
    %                   Unit: The unit of the plot, using one of
    %                         | 'db' | 'mag' | 'pow' |. The default value
    %                         is 'db'.
    %
    %   % Example:
    %   %   Plot the range response of an LFM signal. The signal contains
    %   %   the return from three targets two are approximately 2000 m  
    %   %   away and third is approximately 3500 m away. 
    %
    %   % Load example data
    %   load('RangeResponseExampleData','lfmdata');
    %   fs = lfmdata.fs;
    %   propspeed = lfmdata.propspeed;
    %   rxdata = lfmdata.rxdata;
    %   
    %   % Create range response 
    %   w1 = {'LinearFM','PulseWidth',40/fs,'PRF',fs/320,...
    %        'SweepBandwidth',(20*fs)/40};
    %   p1 = { 'MatchedFilter','SpectrumWindow','None'};
    %   idx = 1;
    %   rngresp = phased.PulseCompressionLibrary(...
    %       'WaveformSpecification',{w1},...
    %       'ProcessingSpecification',{p1},...
    %       'SampleRate',fs,...
    %       'PropagationSpeed',propspeed);
    %   
    %   % Plot range response of processed data
    %   plotResponse(rngresp,rxdata,idx,'Unit','db');
    %
    %   See also phased.RangeResponse, phased.RangeDopplerResponse,
    %   phased.RangeResponse/plotResponse,
    %   phased.RangeDopplerResponse/plotResponse,
    %   phased.AngleDopplerResponse/plotResponse.

        if ~isempty(coder.target)
            coder.internal.assert(false, ...
                'phased:Waveform:CodegenNotSupported','plotResponse');
        end

        narginchk(3,inf);
        if (mod(nargin,2) ~= 0)
            pulse_idx = 1;            
        else
            pulse_idx = varargin{1};
            varargin(1) = [];
        end
        
        if ~ismatrix(x) || isempty(x)
            error(message(...
                'phased:step:MatrixNonEmpty','X'));
        end
        
        rngresp = clone(obj);
        release(rngresp);
        
        for i = 1:pulse_idx
            [resp,rng_grid] = step(rngresp,x,idx);
        end
        
        unit = 'db';
        sigutils.pvparse(varargin{:});
        unit = validatestring(unit,{'db','mag','pow'},...
            'plotResponse','Unit');
        
        hresp = phased.internal.RangePattern(...
            'Range',rng_grid,...
            'Pattern',abs(resp));
        if nargout
            varargout{1} = ...
                plot(hresp,'Units',unit);
        else
            plot(hresp,'Units',unit);
        end
    end
end



methods(Access = protected)
    
    function setupImpl(obj,x,~)
        
        setupImpl@phased.internal.AbstractLibrary(obj);
        
        coder.extrinsic('phased.PulseCompressionLibrary.processparams');
        coder.extrinsic('phased.PulseCompressionLibrary.getWinCoeff');
        coder.extrinsic('phased.PulseCompressionLibrary.MatchedFilterparams');
        
        obj.pNumInputChannels = size(x,2);
        obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        
        if ndims(x) == 3
            obj.pInput3DFlag = true;
            obj.pCubeDim = size(x);
        else
            obj.pInput3DFlag = false;
        end
        
        CompressInfo = obj.ProcessingSpecification;
        WavSpec = obj.WaveformSpecification;
        n = numel(CompressInfo);
        m = numel(WavSpec);
        
        cond = m ~= n;
        if cond
            coder.internal.errorIf(cond,...
                'phased:PulseCompressionLibrary:InvalidSpecifications',...
                'WaveformSpecification','ProcessingSpecification');
        end
        
        FreqOffset = zeros(1,n);
        FFT_len = zeros(1,n);
        obj.pRngFFTLength = zeros(1,n);
        temp = cell(1,n);
        params = cell(1,n);
        
        rng_winCoeff = cell(1,n);
        grid = cell(1,n);
        
        fs = obj.SampleRate; % property/method duality
        cond = ~isscalar(fs) || (fs<=0);
        if cond
            coder.internal.errorIf(cond,...
                'phased:phased:invalidSampleTime');
        end
        
        for m = coder.unroll(1:n)
            
            if strcmpi(CompressInfo{m}{1},'StretchProcessor')
                params{m} =  coder.const(...
                    @phased.PulseCompressionLibrary.processparams,...
                    CompressInfo{m},WavSpec{m},fs);                
            else
                if ~isa(WavSpec{m}{1},'function_handle')
                    params{m} = coder.const(...
                        @phased.PulseCompressionLibrary.processparams,...
                        CompressInfo{m},WavSpec{m},fs);
                else
                    params{m} = coder.const(@phased.PulseCompressionLibrary.MatchedFilterparams,...
                        CompressInfo{m},fs);
                end
            end
        end
        wav_prf = populatePRF(obj);
        
        for idx = coder.unroll(1:n)
            
            samples = obj.SampleRate/wav_prf(idx);
            if ~strcmpi(CompressInfo{idx}{1},'StretchProcessor')
                temp{idx} = ...
                    phased.MatchedFilter(params{idx}{:},...
                    'CoefficientsSource','Input port',...
                    'MaximumNumInputSamplesSource','Property',...
                    'MaximumNumInputSamples',samples);
                
                grid{idx} = zeros(samples,1);
                grid{idx} = obj.PropagationSpeed*...
                    (0:samples-1).'/obj.SampleRate/2;
                rng_winCoeff{idx} = 0;
                
            else
                pw = find(strcmpi(obj.WaveformSpecification{idx},...
                    'PulseWidth'),1);
                dc = find(strcmpi(obj.WaveformSpecification{idx},...
                    'DutyCycle'),1);
                
                [PRF,PulseWidth,DutyCycle,SweepBandwidth,~,...
                    ~,~,FrequencyOffset] = parseInput_lfm(...
                    obj.WaveformSpecification{idx});
                
                FreqOffset(idx) = FrequencyOffset;
                
                if isempty(pw) && ~isempty(dc)
                    PulseWidth = DutyCycle/PRF;
                end
                
                [~,ReferenceRange,RangeFFTLength,...
                    RangeWindow,Beta,Nbar,SidelobeAttenuation] = ...
                    parseInput_Stretch(CompressInfo{idx});
                
                SweepSlope = SweepBandwidth/PulseWidth;
                
                RangeWindow = validatestring(RangeWindow,{'None',...
                    'Hamming','Chebyshev','Hann','Kaiser','Taylor'},'PulseCompressionLibrary');
                
                temp{idx} = ...
                    phased.StretchProcessor(params{idx}{:},...
                    'PropagationSpeed',obj.PropagationSpeed);
                
                num_rng_samples = obj.SampleRate/PRF;
                rng_winCoeff{idx} = ...
                    coder.const(phased.PulseCompressionLibrary.getWinCoeff(RangeWindow,num_rng_samples,...
                    SidelobeAttenuation, ...
                    Nbar,Beta));
                
                if isempty(RangeFFTLength)
                    FFT_len(idx) = num_rng_samples;
                    grid{idx} = zeros(num_rng_samples,1);
                else
                    FFT_len(idx) = coder.const(RangeFFTLength);
                    grid{idx} = zeros(RangeFFTLength,1);
                end
                
                grid{idx} = stretchfreq2rng(...
                    obj.fftshiftfreqgrid(FFT_len(idx),...
                    obj.SampleRate),SweepSlope,ReferenceRange,...
                    obj.PropagationSpeed);
            end
        end
        
        obj.pRngWinCoeff = rng_winCoeff;
        obj.pRngFFTLength = FFT_len;
        obj.pFreqOffset = FreqOffset;
        obj.cProcessingSpecification = temp;
        obj.pRngGrid = grid;
    end

    function num = getNumInputsImpl(~)
        num = 2; % Idx & Data
    end
    
    function num = getNumOutputsImpl(~)
        num = 2; % Data out & range samples
    end
    
   
    function validatePropertiesImpl(obj)
        
      processInfo = obj.ProcessingSpecification;
      n = numel(processInfo);
      inp = obj.WaveformSpecification;
            
      validateDuplicateWarningCond(obj,inp); 
      for idx = coder.unroll(1:n)
          if strcmpi(processInfo{idx}{1},'MatchedFilter')
              
              [SpectrumRange,SpectrumWindow] = ...
                  parseInput_MatchedFilter(processInfo{idx});

              if ~strcmpi(SpectrumWindow,'None')
                  Fs = obj.SampleRate;
                  
                  cond = any(SpectrumRange < -Fs/2) || ...
                      any(SpectrumRange > Fs/2);
                  if cond
                      coder.internal.errorIf(cond,'phased:PulseCompressionLibrary:InvalidBaseBand',...
                          idx,'ProcessingSpecification',num2str(-Fs/2), num2str(Fs/2),...
                          num2str(SpectrumRange(1)), num2str(SpectrumRange(2)));
                  end
              end
          else
              
              wave = validatestring(obj.WaveformSpecification{idx}{1},...
                  {'LinearFM','Rectangular','SteppedFM','PhaseCoded'},...
                  'WaveformSpecification');
              
              if ~strcmpi(wave,'LinearFM')
                  coder.internal.assert(false,...
                      'phased:PulseCompressionLibrary:StretchProcessingNotSupported');
              end
              
              
              pw = find(strcmpi(inp{idx},'PulseWidth'),1);
              dc = find(strcmpi(inp{idx},'DutyCycle'),1);
              
              [PRF,PulseWidth,DutyCycle] = parseInput_lfm(...
                  inp{idx});
              
              if isempty(pw) && ~isempty(dc)
                  PulseWidth = DutyCycle/PRF;
                  
              end
              
              [RangeSpan,ReferenceRange] = ...
                  parseInput_Stretch(processInfo{idx});
              
              phased.PulseCompressionLibrary.validateStretchRangeSpan(...
                  obj.PropagationSpeed,PRF,PulseWidth,...
                  ReferenceRange,RangeSpan,idx)
          end
      end  
    end

    function validateInputsImpl(obj,x,idx)
        
        cond = ndims(x) > 3;
        if cond
            coder.internal.errorIf(cond,'phased:PulseCompressionLibrary:InvalidDimensions');
        end
        
        cond = ~isa(x,'double');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','double');
        end
        
        validateNumChannels(obj,x);
        if ndims(x) == 3 && obj.pCubeDim(3) ~= -1
            validateNumPages(obj,x,obj.pCubeDim(3));
        end
        
        cond = ~isa(idx,'double');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','IDX','double');
        end
        
        cond = ~isscalar(idx);
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeScalar','IDX');
        end
      
    end
    
    function idx = validateInputIdx(obj,value)
        idx = sigdatatypes.validateIndex(value,'','IDX',{'<=',...
            numel(obj.ProcessingSpecification)});
    end
    
    function validateProcessingparams(~,value)
        
        n = numel(value);
        coder.unroll();
        for idx = 1:n 
            
            cond = ~iscell(value{idx}) ;
            if cond
                coder.internal.errorIf(cond,...
                    'phased:PulseWaveformLibrary:InvalidDatatype',...
                    idx,'ProcessingSpecification')
            end
            
            cond = ~(strcmpi(value{idx}{1},'StretchProcessor') || ...
                strcmpi(value{idx}{1},'MatchedFilter'));
            
            if cond
                coder.internal.errorIf(cond,...
                    'phased:PulseCompressionLibrary:InvalidProcessing',...
                    idx,'ProcessingSpecification','MatchedFilter','StretchProcessor');
            end
            
            if ~strcmpi(value{idx}{1},'StretchProcessor')
                
                param = {'CoefficientsSource','SampleRate','GainOutputPort',...
                    'CustomSpectrumWindow','MaximumNumInputSamplesSource',...
                    'SweepSlope','SampleRate','PRFSource','PRF',.....
                    'PulseWidth','SweepInterval','PropagationSpeed',...
                    'ReferenceRange','RangeSpan','RangeFFTLength',...
                    'RangeWindow','MaximumNumInputSamples'};
                
                y = value{idx};
                
                for i = 1:numel(y)
                    for j = 1:numel(param)
                        cond =  any(strcmpi(y{i},param{j}));
                        if cond
                            coder.internal.errorIf(cond,...
                                'phased:PulseWaveformLibrary:InvalidParameter',....
                                param{j},idx,'ProcessingSpecification');
                        end
                    end
                end
                
                [SpectrumRange,SpectrumWindow,Beta,Nbar,...
                    SidelobeAttenuation,Coefficients] = ...
                    parseInput_MatchedFilter(value{idx});
                
                validateattributes( SpectrumRange, { 'double' }, ...
                    { 'real', 'finite', 'size', [ 1, 2 ] }, '',...
                    'SpectrumRange',idx);
                
                cond = SpectrumRange(1) > SpectrumRange(2);
                if cond
                    coder.internal.errorIf(cond,...
                        'phased:PulseCompressionLibrary:InvalidSpectrumRange',...
                        idx,'ProcessingSpecification');
                end
                
                validatestring(SpectrumWindow,{'None',...
                    'Hamming','Chebyshev','Hann','Kaiser','Taylor'},...
                    '','SpectrumWindow',idx);
                
                validateattributes(Beta, { 'double' },...
                    { 'scalar', 'nonnegative', 'finite' }, '', 'Beta',idx);
                
                validateattributes(SidelobeAttenuation, { 'double' },...
                    { 'scalar','positive', 'finite' }, '',...
                    'SidelobeAttenuation',idx);
                
                validateattributes(Nbar, { 'double' },...
                    { 'scalar', 'positive', 'finite' }, '', 'Nbar',idx);
                
                if ~isempty(Coefficients)
                    validateattributes(Coefficients, { 'double' }, ...
                        {'vector', 'finite' }, '', 'Coefficients',idx);
                    
                    cond = ~iscolumn(Coefficients) || isempty(Coefficients);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:inputMustBeColVector','Coeff');
                    end
                end
            elseif strcmpi(value{idx}{1},'StretchProcessor')
                
                param = {'SweepSlope','SampleRate','PRFSource','PRF',.....
                    'PulseWidth','SweepInterval','PropagationSpeed',...
                    'CoefficientsSource','SampleRate','GainOutputPort',...
                    'CustomSpectrumWindow','MaximumNumInputSamplesSource',...
                    'Coefficients','SpectrumWindow','SpectrumRange',...
                    'MaximumNumInputSamples'};
                
                y = value{idx};
                
                for i = 1:numel(y)
                    for j = 1:numel(param)
                        cond =  any(strcmpi(y{i},param{j}));
                        if cond
                            coder.internal.errorIf(cond,...
                                'phased:PulseWaveformLibrary:InvalidParameter',....
                                param{j},idx,'ProcessingSpecification');
                        end
                    end
                end
                            
                [RangeSpan,ReferenceRange,RangeFFTLength,...
                    RangeWindow,Beta,Nbar,SidelobeAttenuation] = ...
                    parseInput_Stretch(value{idx});

                validateattributes(RangeSpan,{'double'},{'scalar',...
                    'positive'},'','RangeSpan',idx);
                
                validateattributes(ReferenceRange,{'double'},{'scalar',...
                    'positive'},'','ReferenceRange',idx);
                
                if ~isempty(RangeFFTLength)
                    validateattributes(RangeFFTLength,{'double'},{'scalar',...
                        'positive','finite'},'','RangeFFTLength',idx);
                end
                validatestring(RangeWindow,{'None',...
                    'Hamming','Chebyshev','Hann','Kaiser','Taylor'},...
                    '','RangeWindow',idx);
                
                validateattributes(Beta, { 'double' },...
                    { 'scalar', 'nonnegative', 'finite' }, '', 'Beta',idx);
                
                validateattributes(SidelobeAttenuation, { 'double' },...
                    { 'scalar','positive', 'finite' }, '',...
                    'SidelobeAttenuation',idx);
                
                validateattributes(Nbar, { 'double' },...
                    { 'scalar', 'positive', 'finite' }, '', 'Beta',idx);
            end
        end
        
    end
    
    function [rng,rng_grid] = stepImpl(obj,x_in,idx)
        
        currentIdx = idx;
        cPulseCompressSpec = obj.cProcessingSpecification;
        rng_winCoeff = obj.pRngWinCoeff;
        
        inp = obj.WaveformSpecification;
        CompressInfo = obj.ProcessingSpecification;
        N = numel(CompressInfo);
        rng = complex(0);
        rng_grid = 0;
        
        sz_x = size(x_in);
        idx = validateInputIdx(obj,idx);
        assert(idx<=N);
        for m = coder.unroll(1:N)
            if m == idx
                if strcmpi(CompressInfo{m}{1},'StretchProcessor')
                    x = x_in;
                    % Reshape Cube inputs before Stretch Processing
                    if numel(sz_x) == 3
                        x_inp = reshape(x,sz_x(1),[]);
                    else
                        x_inp = x;
                    end
                    y = step(cPulseCompressSpec{m},x_inp); % StretchProcessing
                    
                    if numel(size(x)) == 3
                        yout = reshape(y,sz_x);
                    else
                        yout = y;
                    end
                    
                    r_fft_input = bsxfun(@times,rng_winCoeff{m},yout); % Range Windowing
                    
                    rng = fftshift(fft(complex(r_fft_input),obj.pRngFFTLength(m),1)); % Range Processing
                    rng_grid = obj.pRngGrid{m};
                else
                    
                    if  m == obj.pIndex
                        pulseIdx = obj.pPulseIdx+1;
                    else
                        pulseIdx = 1;
                    end
                    
                    if strcmpi(inp{m}{1},'SteppedFM')
                        [~,~,~,~,NumSteps] = parseInput_steppedfm(inp{m});
                        
                        n = mod(pulseIdx-1,NumSteps)+1;
                    else
                        n = 1;
                    end
                    
                    [~,~,~,~,~,Coefficients...
                        ] = parseInput_MatchedFilter(CompressInfo{m});
                    
                    if ~isempty(Coefficients)
                        Coeffs = Coefficients;
                    else
                        wav = phased.PulseWaveformLibrary(...
                            'SampleRate',obj.SampleRate,...
                            'WaveformSpecification',{inp{m}});
                        if ~isa(inp{m}{1},'function_handle')
                            Coeffs = getMatchedFilter(wav,1,n);
                        else
                            Coeffs = getMatchedFilter(wav,1);
                        end
                    end
                    Coeff = complex(Coeffs);
                    RangeOffset = numel(Coeff)-1;
                    x_rng = step(cPulseCompressSpec{m},x_in,Coeff); % MatchedFiltering
                    
                    rng = complex(zeros(sz_x));
                    rng(1:size(x_in,1)- RangeOffset,:,:) = ...
                        x_rng(RangeOffset+1:sz_x(1),:,:);
                    
                    grid = obj.pRngGrid{m};
                    rng_grid = grid(1:sz_x(1));
                    obj.pPulseIdx = pulseIdx;
                    obj.pIndex = currentIdx;
                    break;
                end
            end
        end
    end
    
    function resetImpl(obj)
        for idx = 1:length(obj.cProcessingSpecification)
            reset(obj.cProcessingSpecification{idx});
        end
        obj.pIndex = 0;
        obj.pPulseIdx = 1;
    end
    
    function releaseImpl(obj)
        for idx = 1:length(obj.cProcessingSpecification)
            release(obj.cProcessingSpecification{idx});
            
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractLibrary(obj); 
        s.isLocked = isLocked(obj);
        if isLocked(obj)
            s.pFreqOffset = obj.pFreqOffset;
            s.pRngGrid = obj.pRngGrid;
            s.pRngWinCoeff = obj.pRngWinCoeff;
            s.pInput3DFlag = obj.pInput3DFlag;
            s.pCubeDim = obj.pCubeDim;
            s.pRngFFTLength = obj.pRngFFTLength;
            s.pIndex = obj.pIndex;
            s.pPulseIdx = obj.pPulseIdx;
            for idx=1:length(obj.ProcessingSpecification)
                s.cProcessingSpecification{idx} = ...
                    saveobj(obj.cProcessingSpecification{idx});
            end
        end
    end   
    function s = loadSubObjects(obj,s)
        if isfield(s,'isLocked')
            if s.isLocked
                for idx = 1:length(s.ProcessingSpecification)
                    if strcmpi(s.ProcessingSpecification{idx}{1},'MatchedFilter')
                        obj.cProcessingSpecification{idx} = ...
                            phased.MatchedFilter.loadobj(...
                            s.cProcessingSpecification{idx});
                    else
                        obj.cProcessingSpecification{idx} = ...
                            phased.StretchProcessor.loadobj(...
                            s.cProcessingSpecification{idx});
                    end
                end
                s = rmfield(s,'cProcessingSpecification');
                s = rmfield(s,{'pFrequencyOffset','pPreviousIndex','pIsFunctionHandle',...
                    'cWaveformSpecification'});
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

methods (Access = protected)
    
    function flag = isInputComplexityLockedImpl(~,index) 
        flag = true;
        if index == 1
            flag = false;
        end
    end
    
    function flag = isInputSizeLockedImpl(~,index)  
        % Return true if input size is not allowed to change while
        % system is running        
        if index == 1
            flag = false;
        else
            flag = true;
        end
    end
    
    function varargout = getOutputSizeImpl(obj)
        szX = propagatedInputSize(obj,1);
        n = numel(obj.ProcessingSpecification);
        len = zeros(1,n);
        
        for idx = 1:n
            if strcmpi(obj.ProcessingSpecification{idx}{1},...
                    'StretchProcessor')
                [~,~,rng_FFT] = parseInput_Stretch(obj.ProcessingSpecification{idx});
                if ~isempty(rng_FFT)
                    len(idx) = rng_FFT;
                else
                    len(idx) = szX(1);
                end
            else
                len(idx) = szX(1);
            end
        end
        
        P = max(len);
        
        if numel(szX)==3
            respSz = [P szX(2:3)];
        else
            respSz = [P szX(2)];
        end
        
        varargout = {respSz,[P 1]};
    end
    
    function varargout = getOutputNamesImpl(~)
        varargout = {'Y','Range'};
    end
    
    function varargout = getInputNamesImpl(~)
        varargout = {'X','Idx'};
    end
    
    function varargout = isOutputFixedSizeImpl(~)
        varargout{1} = false;
        varargout{2} = false;
    end
    
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj,1);
        varargout{2} = 'double';
    end
    
    function varargout = isOutputComplexImpl(~)
        varargout{1} = true;
        varargout{2} = false;
    end

end


methods (Static, Hidden, Access = protected)
    
    function groups = getPropertyGroupsImpl
        props = {'PropagationSpeed',...
            'WaveformSpecification',....
            'ProcessingSpecification'...
            'SampleRateFromInputCheckbox',...
            'SampleRate'...            
            };
        groups = matlab.system.display.Section('Title',...
            'Parameters', ...
            'PropertyList', props);
        
    end
    
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:PulseCompressionLibraryTitle')),...
            'Text',getString(message('phased:library:block:PulseCompressionLibraryDesc')));
    end
end

methods (Access=protected)
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Pulse Compression\nLibrary');
    end
end

methods(Hidden,Static)
    
    function validateStretchRangeSpan(c,prf,pw,refrng,rngspan,idx)
        Rmin = 0;
        Rmax = c*(1/(2*prf)-pw/2);
        Rinterval = refrng+[-1 1]*rngspan/2;
        cond = Rinterval(1)<Rmin || Rinterval(2)>Rmax;
        if cond
            coder.internal.errorIf(cond, ...
                'phased:PulseCompressionLibrary:InvalidStretchROI',...
                'ReferenceRange','RangeSpan',idx,'ProcessingSpecification',...
                feval('sprintf','%5.4f',Rmin),...
                feval('sprintf','%5.4f',Rmax),...
                feval('sprintf','%5.4f',Rinterval(1)),...
                feval('sprintf','%5.4f',Rinterval(2)));
        end
    end
    
    function params = processparams(inp,wavinp,fs)
        if ~strcmpi(inp{1},'MatchedFilter')
                      
            pw = find(strcmpi(wavinp,'PulseWidth'),1);
            dc = find(strcmpi(wavinp,'DutyCycle'),1);
            
            [PRF,PulseWidth,DutyCycle,SweepBandwidth,SweepInterval...
               ] = parseInput_lfm(wavinp);
           
            if isempty(pw) && ~isempty(dc)
                PulseWidth = DutyCycle/PRF;
            end
            
            [RangeSpan,ReferenceRange] = ...
                parseInput_Stretch(inp);
            
            SweepSlope = SweepBandwidth/PulseWidth;
            
            params = {...
                'SampleRate',fs,...
                'PulseWidth',PulseWidth,...
                'PRFSource','Auto',...
                'SweepSlope',SweepSlope,....
                'SweepInterval',SweepInterval,...
                'RangeSpan',RangeSpan,...
                'ReferenceRange',ReferenceRange};
        else
                params = phased.PulseCompressionLibrary.MatchedFilterparams(inp,fs);

        end
    end
    
    function params = MatchedFilterparams(inp,fs)

        
        [SpectrumRange,SpectrumWindow,Beta,Nbar,...
            SidelobeAttenuation...
            ] = parseInput_MatchedFilter(inp);
        SpectrumWindow = validatestring(SpectrumWindow,{'None',...
            'Hamming','Chebyshev','Hann','Kaiser','Taylor'},...
            'PulseProcessingLibrary');
        switch SpectrumWindow
            case 'None'
                params = {'SpectrumWindow','None'};
            case 'Hamming'
                params = {'SampleRate',fs,...
                    'SpectrumRange',SpectrumRange,...
                    'SpectrumWindow','Hamming'};
            case 'Hann'
                params = {'SampleRate',fs,...
                    'SpectrumRange',SpectrumRange,...
                    'SpectrumWindow','Hann'};
            case 'Taylor'
                params = {'SampleRate',fs,...
                    'SpectrumRange',SpectrumRange,...
                    'SpectrumWindow','Taylor',...
                    'Nbar',Nbar,...
                    'SidelobeAttenuation',SidelobeAttenuation};
            case 'Chebyshev'
                params = {'SampleRate',fs,...
                    'SpectrumRange',SpectrumRange,...
                    'SpectrumWindow','Chebyshev',...
                    'SidelobeAttenuation',SidelobeAttenuation};
            case 'Kaiser'
                params = {'SampleRate',fs,...
                    'SpectrumRange',SpectrumRange,...
                    'SpectrumWindow','Kaiser',...
                    'Beta',Beta};
        end
    end

    function winCoeff = getWinCoeff(rangeWindow,num_rng_samples,...
            rangeSidelobeAttenuation,Nbar,Beta)
        switch rangeWindow
            case 'None'
                winCoeff = ones(num_rng_samples,1);
            case 'Hamming'
                winCoeff = hamming(num_rng_samples);
            case 'Hann'
                winCoeff = hann(num_rng_samples);
            case 'Kaiser'
                winCoeff = kaiser(num_rng_samples,Beta);
            case 'Chebyshev'
                winCoeff = chebwin(num_rng_samples,...
                    rangeSidelobeAttenuation);
            case 'Taylor'
                winCoeff = taylorwin(num_rng_samples,...
                    Nbar,-rangeSidelobeAttenuation);
        end
    end

    function freq_grid = fftshiftfreqgrid(N,Fs)
        %fftshiftfreqgrid   Generate frequency grid
        %   freq_grid = fftshiftfreqgrid(N,Fs) generate an N point
        %   frequency grid according to sample rate Fs. This grid matches
        %   the operation used in fftshift.
        %
        %   % Example:
        %   %   Create a 16 point frequency grid for a sample rate of 10 
        %   %   Hz.
        %   
        %   freq_grid = fftshiftfreqgrid(16,10)
            
        % set 'CenterDC' in psdfreqvec to true preserves Nyquist point,
        % which does not match our processing to the data since we use
        % fftshift.
        
        % adopted from psdfreqvec 
        
            freq_res = Fs/N;
            freq_grid = (0:N-1).'*freq_res;
            Nyq = Fs/2;
            half_res = freq_res/2;
            if rem(N,2) % odd
                idx = 1:(N-1)/2;
                halfpts = (N+1)/2;
                freq_grid(halfpts) = Nyq-half_res;
                freq_grid(halfpts+1) = Nyq+half_res;
            else
                idx = 1:N/2;
                hafpts = N/2+1;
                freq_grid(hafpts) = Nyq;
            end
            freq_grid(N) = Fs-freq_res;
            freq_grid = fftshift(freq_grid);
            freq_grid(idx) = freq_grid(idx)-Fs;
    end
end
end
function [RangeSpan,ReferenceRange,...
    RangeFFTLength,RangeWindow,Beta,Nbar,SidelobeAttenuation] = ...
    parseInput_Stretch(varargin)
    
    fft_len = find(strcmpi(varargin,...
                        'RangeFFTLength'),1);
    if ~isempty(fft_len)
        defaultRangeFFTLength = 1024;
    else
        defaultRangeFFTLength = [];
    end
    defaultRangeSpan = 500;
    defaultReferenceRange = 5000;
    defaultBeta = 0.5;
    defaultNbar = 4;
    defaultSidelobeAttenuation = 30;
    defaultRangeWindow = 'None';
       
    if isempty(coder.target)
        p = inputParser;
        p.FunctionName = 'PulseCompressionLibrary';
        p.addParameter('RangeSpan',defaultRangeSpan);
        p.addParameter('ReferenceRange',defaultReferenceRange);
        p.addParameter('Beta',defaultBeta);
        p.addParameter('Nbar',defaultNbar);
        p.addParameter('SidelobeAttenuation',defaultSidelobeAttenuation);
        p.addParameter('RangeFFTLength',defaultRangeFFTLength);
        p.addParameter('RangeWindow',defaultRangeWindow);
        p.parse(varargin{1}{2:end});

        RangeSpan = p.Results.RangeSpan;
        ReferenceRange = p.Results.ReferenceRange;
        Beta = p.Results.Beta;
        Nbar = p.Results.Nbar;
        SidelobeAttenuation = p.Results.SidelobeAttenuation;
        RangeFFTLength = p.Results.RangeFFTLength;
        RangeWindow = p.Results.RangeWindow;
    else
        parms = struct(...
            'RangeSpan',uint32(0),...
            'ReferenceRange',uint32(0),...
            'RangeFFTLength',uint32(0),...
            'RangeWindow',uint32(0),...
            'Beta',uint32(0),...
            'Nbar',uint32(0),...
            'SidelobeAttenuation',uint32(0));
        poptions = struct( ...
            'CaseSensitivity',false, ...
            'PartialMatching','unique', ...
            'StructExpand',false, ...
            'IgnoreNulls',false);   
        
        pstruct = coder.internal.parseParameterInputs(...
            parms,poptions,varargin{1}{2:end});

        RangeSpan = coder.internal.getParameterValue(pstruct.RangeSpan,...
            defaultRangeSpan,varargin{1}{2:end});
        ReferenceRange = coder.internal.getParameterValue(...
            pstruct.ReferenceRange,defaultReferenceRange,varargin{1}{2:end});
        Beta = coder.internal.getParameterValue(...
            pstruct.Beta,defaultBeta,...
            varargin{1}{2:end});
        Nbar = coder.internal.getParameterValue(...
            pstruct.Nbar,defaultNbar,...
            varargin{1}{2:end});
        SidelobeAttenuation = coder.internal.getParameterValue(...
            pstruct.SidelobeAttenuation,defaultSidelobeAttenuation,...
            varargin{1}{2:end});
        RangeFFTLength = coder.internal.getParameterValue(pstruct.RangeFFTLength,...
            defaultRangeFFTLength,varargin{1}{2:end});
        RangeWindow = coder.internal.getParameterValue(...
            pstruct.RangeWindow,defaultRangeWindow,varargin{1}{2:end});
    end
end

function [SpectrumRange,SpectrumWindow,Beta,Nbar,...
    SidelobeAttenuation,Coefficients] = parseInput_MatchedFilter(varargin)

    coeff_idx = find(strcmpi(varargin,...
        'Coefficients'),1);

    if ~isempty(coeff_idx)
        defaultCoefficients = [1;1];
    else
        defaultCoefficients = [];
    end

    defaultSpectrumRange = [0 1e5];
    defaultSpectrumWindow = 'None';
    defaultBeta = 0.5;
    defaultNbar = 4;
    defaultSidelobeAttenuation = 30;

    if isempty(coder.target)
        p = inputParser;
        p.FunctionName = 'PulseCompressionLibrary';
        p.addParameter('SpectrumRange',defaultSpectrumRange);
        p.addParameter('SpectrumWindow',defaultSpectrumWindow);
        p.addParameter('Beta',defaultBeta);
        p.addParameter('Nbar',defaultNbar);
        p.addParameter('SidelobeAttenuation',defaultSidelobeAttenuation);
        p.addParameter('Coefficients',defaultCoefficients);
        p.parse(varargin{1}{2:end});
       
        SpectrumRange = p.Results.SpectrumRange;
        SpectrumWindow = p.Results.SpectrumWindow;
        Beta = p.Results.Beta;
        Nbar = p.Results.Nbar;
        SidelobeAttenuation = p.Results.SidelobeAttenuation;
        Coefficients = p.Results.Coefficients;
    else
        parms = struct( ...
            'SpectrumRange',uint32(0), ...
            'SpectrumWindow',uint32(0), ...        
            'Beta',uint32(0),...
            'Nbar',uint32(0),...
            'SidelobeAttenuation', uint32(0),...
            'Coefficients', uint32(0));
        % Select parsing options.
        poptions = struct( ...
            'CaseSensitivity',false, ...
            'PartialMatching','unique', ...
            'StructExpand',false, ...
            'IgnoreNulls',false);
        pstruct = coder.internal.parseParameterInputs(...
            parms,poptions,varargin{1}{2:end});
        SpectrumRange = coder.internal.getParameterValue(...
            pstruct.SpectrumRange,defaultSpectrumRange,varargin{1}{2:end});
        SpectrumWindow = coder.internal.getParameterValue(...
            pstruct.SpectrumWindow,defaultSpectrumWindow,varargin{1}{2:end});
        Beta = coder.internal.getParameterValue(...
            pstruct.Beta,defaultBeta,...
            varargin{1}{2:end});
        Nbar = coder.internal.getParameterValue(...
            pstruct.Nbar,defaultNbar,...
            varargin{1}{2:end});
        SidelobeAttenuation = coder.internal.getParameterValue(...
            pstruct.SidelobeAttenuation,defaultSidelobeAttenuation,...
            varargin{1}{2:end});
        Coefficients  = coder.internal.getParameterValue(...
            pstruct.Coefficients,defaultCoefficients,...
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
        p.FunctionName = 'PulseCompressionLibrary';
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

