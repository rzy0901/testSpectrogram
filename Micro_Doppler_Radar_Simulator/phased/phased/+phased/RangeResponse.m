classdef (Sealed,StrictDefaults) RangeResponse < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%RangeResponse   Range response
%   H = phased.RangeResponse creates a range response System object, H.
%   This object calculates the range response of the input data.
%
%   H = phased.RangeResponse(Name,Value) creates a range response object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   The RangeResponse object generates the response by processing the input
%   signal in the range domain using either a matched filter or dechirp
%   operation.
%
%   Step method syntax:
%   
%   [RESP,RANGE] = step(H,X) calculates the range response of the input
%   data X. This syntax applies when you set the RangeMethod property to
%   'FFT' and then the DechirpInput property to false. This syntax is most
%   commonly used with FMCW signals.
%
%   X is a dechirped input signal. X must be a KxL matrix, a KxN matrix, or
%   a KxNxL array where K denotes the number of fast time samples, L is the
%   number of dechirped frequency sweeps, and N is the number of channels
%   (antenna elements or beams). Each column and page of X is processed as
%   an independent Kx1 signal.
%
%   RESP is either an MxL matrix, an MxN matrix, or an MxNxL array
%   containing the complex range response of the input X, where the number
%   of dimensions in RESP will match the number of dimensions in X. M is
%   determined by either the number of rows in X when you set the
%   RangeFFTLengthSource property to 'Auto', or the value specified in the
%   RangeFFTLength property when you set the RangeFFTLengthSource property
%   to 'Property'.
%
%   RANGE is a length M column vector containing the range samples at which
%   the range response is evaluated.
%
%   [RESP,RANGE] = step(H,X,XREF) uses input XREF as the reference signal
%   to dechirp the input signal X. This syntax applies when you set the
%   RangeMethod property to 'FFT' and then the DechirpInput property to
%   true. This syntax is most commonly used with FMCW signals and the
%   reference signal is, in general, the transmitted signal.
%
%   X is an input signal to be dechirped by the RangeResponse object. X
%   must be a KxL matrix, a KxN matrix, or a KxNxL array where K denotes
%   the number of fast time samples, L is the number of frequency sweeps,
%   and N is the number of channels (antenna elements or beams). XREF must
%   be a column vector whose number of rows is the same as the number of
%   rows of X. XREF is the reference signal used to dechirp X.
%
%   In this syntax, the number of rows in RESP is the quotient of the
%   number of rows in X and the value specified in the DecimationFactor
%   property.
%
%   [RESP,RANGE] = step(H,X,COEFF) uses COEFF as the matched filter
%   coefficients. This method applies when you set the RangeMethod property
%   to 'Matched filter'. This syntax is most commonly used with pulsed
%   signals.
%
%   X is an input signal to be match filtered by the RangeResponse object.
%   X must be a KxL matrix, a KxN matrix, or a KxNxL array where K denotes
%   the number of fast time samples, L is the number of pulses, and N is
%   the number of channels (antenna elements or beams). COEFF must be a
%   column vector containing the matched filter coefficients.
%
%   In this syntax, the number of rows of RESP equals the number of rows of
%   X.
%
%   RANGE is a length K column vector containing the range samples at which
%   the range response is evaluated.
%   
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   RangeResponse methods:
%
%   step         - Calculate range response
%   release      - Allow property value and input characteristics changes
%   clone        - Create range response object with same property values
%   isLocked     - Locked status (logical)
%   plotResponse - Plot range response
%
%   RangeResponse properties:
%       
%   RangeMethod                  - Range processing method
%   PropagationSpeed             - Propagation speed
%   SampleRate                   - Sample rate
%   SweepSlope                   - FM sweep slope
%   DechirpInput                 - Dechirp input signal
%   DecimationFactor             - Decimation factor for dechirped signal
%   RangeFFTLengthSource         - Source of FFT length in range processing
%   RangeFFTLength               - FFT length in range processing
%   RangeWindow                  - Range processing window
%   RangeSidelobeAttenuation     - Range sidelobe attenuation level 
%   CustomRangeWindow            - Custom range processing window
%   ReferenceRangeCentered       - Set reference range at center
%   ReferenceRange               - Reference range
%   MaximumNumInputSamplesSource - Source of maximum number of samples
%                                  of the input signal
%   MaximumNumInputSamples       - Maximum number of samples in input 
%                                  signal
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Calculate the range response from a pulsed radar transmitting an
%   %   LFM waveform using the matched filter approach. The signal includes
%   %   three target returns. Two are approximately 2000 m away and the
%   %   third is approximately 3500 m away.
%
%   % Load example data
%   load('RangeResponseExampleData','lfmdata');
%   fs = lfmdata.fs;
%   propspeed = lfmdata.propspeed;
%   rxdata = lfmdata.rxdata;
%   mfcoeffs = lfmdata.mfcoeffs;
%   noisepower = lfmdata.noisepower;
%   
%   % Create range response for matched filter processing
%   rngresp = phased.RangeResponse(...
%       'SampleRate',fs,...
%       'PropagationSpeed',propspeed);
%   
%   % Perform range processing
%   [resp,rng_grid] = rngresp(rxdata,mfcoeffs);
%   
%   % Calculate noise power after range processing
%   mfgain = mfcoeffs'*mfcoeffs;
%   noisepower_proc = mfgain*noisepower;
%   
%   subplot(2,1,1);
%   plot(rng_grid,pow2db(abs(rxdata).^2./noisepower));
%   ylabel('SNR (dB)'); title('Before Range Processing');
%   xlim(rng_grid([1 end]));
%   xlabel('Range (m)');
%   
%   subplot(2,1,2);
%   plot(rng_grid,pow2db(abs(resp).^2./noisepower_proc));
%   ylabel('SNR (dB)'); title('After Range Processing');
%   xlim(rng_grid([1 end]));
%   xlabel('Range (m)');
%
%   See also phased, phased.RangeDopplerResponse,
%   phased.AngleDopplerResponse, phased.MatchedFilter, dechirp.

%   Copyright 2016-2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)

        %RangeMethod    Range processing method
        %   Specify the method of range processing as one of 'Matched
        %   filter' | 'FFT', where the default is 'Matched filter'. When
        %   you set this property to 'Matched filter', the range processing
        %   is achieved by applying a matched filter to the incoming
        %   signal. When you set this property to 'FFT', the range
        %   processing is achieved by applying FFT to the input signal.
        %
        %   The matched filter approach is often used with pulsed signals,
        %   where the matched filter is the time reverse of the transmitted
        %   signal.  The FFT approach are often used with FMCW signals.
        RangeMethod = 'Matched filter';
        %PropagationSpeed   Propagation speed (m/s)
        %   Specify the propagation speed (in m/s) of the signal as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('lightspeed');
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
        %   Specify the sample rate (in Hz) as a positive scalar. The
        %   default value of this property is 1e6 (1 MHz).
        SampleRate = 1e6;
        %SweepSlope     FM sweep slope (Hz/s)
        %   Specify the slope of the linear FM sweeping (in Hz/s) as a
        %   scalar. This property applies when you set the RangeMethod
        %   property to 'FFT'. The default value is 1e9.
        SweepSlope = 1e9;
        
    end
    
    properties (Nontunable,Logical)
        %DechirpInput     Dechirp input signal
        %   Set this property to true to dechirp the input signal first
        %   before range processing. Set this property to false to indicate
        %   that the input signal is already dechirped and no dechirp
        %   operation is necessary. This property applies when you set the
        %   RangeMethod property to 'FFT'. The default value of this
        %   property is false.
        DechirpInput = false;
    end
    
    properties (Nontunable,PositiveInteger)
        %DecimationFactor   Decimation factor for dechirped signal
        %   Specify the decimation factor for the dechirped signal as a
        %   positive integer. This property applies when you set the
        %   RangeMethod property to 'FFT' and then the DechirpInput
        %   property to true. The default value of this property is 1,
        %   indicating no decimation.
        %
        %   When processing FMCW signals, it is often possible to decimate
        %   the dechirped signal to alleviate the requirement on A/D
        %   converter. The decimation algorithm implemented here uses a
        %   30th order FIR filter generated by FIR1(30,1/R) where R is the
        %   value of the property.
        DecimationFactor = 1;
    end
    
    properties (Nontunable)
        %RangeFFTLengthSource   Source of FFT length in range processing
        %   Specify how to determine the FFT length in range processing as
        %   one of 'Auto' | 'Property', where the default is 'Auto'. When
        %   you set this property to 'Auto', the FFT length is equal to the
        %   number of rows of the input signal. When you set this property
        %   to 'Property', the FFT length is specified in RangeFFTLength
        %   property. This property applies when you set the RangeMethod
        %   property to 'FFT'.
        RangeFFTLengthSource = 'Auto'
    end
    
    properties (Nontunable, PositiveInteger)
        %RangeFFTLength     FFT length in range processing
        %   Specify the FFT length in range domain as a positive integer.
        %   This property applies when you set the RangeMethod property to
        %   'FFT' and then the RangeFFTLengthSource property to 'Property'.
        %   The default value of this property is 1024.
        RangeFFTLength = 1024
    end
    
    properties (Nontunable)
        %RangeWindow    Range processing window
        %   Specify the window used for range processing using one of
        %   'None' | 'Hamming' | 'Chebyshev' | 'Hann' | 'Kaiser' | 'Taylor'
        %   | 'Custom', where the default is 'None'. This property applies
        %   when you set the RangeMethod property to 'FFT'.
        %
        %   Note that if you set the RangeWindow property to 'Taylor', the
        %   generated Taylor window has 4 nearly constant sidelobes
        %   adjacent to the mainlobe.
        RangeWindow = 'None'
        %RangeSidelobeAttenuation  Range sidelobe attenuation level 
        %   Specify the sidelobe attenuation level (in dB) of a Kaiser,
        %   Chebyshev or Taylor window in range processing as a positive
        %   scalar. This property applies when you set the 'RangeMethod' to
        %   'FFT' and then the RangeWindow property to 'Kaiser',
        %   'Chebyshev', or 'Taylor'. The default value of this property is
        %   30.
        RangeSidelobeAttenuation = 30
        %CustomRangeWindow   Custom range processing window
        %   Specify the user-defined window for range processing using a
        %   function handle or a cell array. The default value of this
        %   property is @hamming. This property applies when you set the
        %   RangeMethod property to 'FFT' and then the RangeWindow property
        %   to 'Custom'.
        %
        %   If CustomRangeWindow is a function handle, the specified
        %   function takes the window length as the input and generates
        %   appropriate window coefficients.
        %
        %   If CustomRangeWindow is a cell array, then the first cell must
        %   be a function handle. The specified function takes the window
        %   length as the first input argument, with other additional input
        %   arguments if necessary, and generates appropriate window
        %   coefficients. The remaining entries in the cell array are the
        %   additional input arguments to the function, if any.
        CustomRangeWindow = @hamming
    end
    
    properties (Nontunable, Logical)
        %ReferenceRangeCentered     Set reference range at center
        %   Set this property to true to set the reference range to the
        %   center of the range grid. Set this property to false to
        %   set the reference range to the beginning of the range grid. The
        %   default value of this property is true. This property only
        %   applies when you set the RangeMethod to 'FFT'.
        ReferenceRangeCentered = true
    end
    
    properties 
        %ReferenceRange     Reference range (m)
        %   Specify the reference range of the range grid as a nonnegative
        %   scalar. The default value is 0. If you set the RangeMethod
        %   property to 'Matched filter', the reference range marks the
        %   start of the range grid. If you set the RangeMethod property to
        %   'FFT', the position of the reference range is determined by the
        %   ReferenceRangeCentered property. If you set the
        %   ReferenceRangeCentered property to true, the reference range
        %   marks the center of the range grid. If you set the
        %   ReferenceRangeCentered property to false, the reference range
        %   marks the start of the range grid. This property is tunable.
        ReferenceRange = 0
    end
    
    properties (Nontunable)
        %MaximumNumInputSamplesSource  Source of maximum number of samples
        %                       of the input signal
        %   Specify how the maximum number of samples of the input signal
        %   is specified as one of 'Auto' | 'Property', where the default
        %   is 'Auto'. When you set this property to 'Auto', the object
        %   automatically allocates the memory to buffer the input signal.
        %   When you set this property to 'Property', the maximum number of
        %   samples in the input signal is specified via
        %   MaximumNumInputSamples property and any input signal longer
        %   than that value is truncated. This property applies when you
        %   set the RangeMethod property to 'Matched Filter'. The default
        %   value of this property is 'Auto'.
        %
        %   To use the object in MATLAB Function Block in Simulink with
        %   variable-size signal, set this property to 'Property' and set
        %   the MaximumNumInputSamples property.
        MaximumNumInputSamplesSource = 'Auto'
    end

    properties (Nontunable, PositiveInteger)
        %MaximumNumInputSamples Maximum number of samples in input signal
        %   Specify the maximum number of samples in the input signal as a
        %   positive scalar. The input signal is the first input, X, and
        %   the number of samples is number of rows in X. This property
        %   applies when you set the RangeMethod property to 'Matched
        %   Filter' and the MaximumNumInputSamplesSource property to
        %   'Property'. The default value of this property is 100.
        MaximumNumInputSamples = 100;
    end
    
    properties (Constant, Hidden)
        RangeMethodSet = matlab.system.internal.StringSetGF(...
            {'Matched filter','FFT'},{'Dechirp'},{'FFT'});
        RangeWindowSet = matlab.system.StringSet({'None','Hamming',...
            'Chebyshev','Hann','Kaiser','Taylor','Custom'});
        RangeFFTLengthSourceSet = dsp.CommonSets.getSet(...
            'AutoOrProperty');
        MaximumNumInputSamplesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties (Access = private, Nontunable, Logical)
        pUseMatchedFilter = false
    end
    
    properties (Access = private)
        cMatchedFilter
        cRangeFFT
    end
    
    properties (Access = private, Nontunable)
        pRangeWinCoeff
        pSampleRate
        pRangeGrid
        pRangeOffset = 0
        pRangeFFTLength
        pDecimationFilterCoefficients
        pValidatedNumInputPages
    end
    
    properties (Access = private)
        pNumPages = -1
    end
    
    properties (Access = private, Logical)
        pSizeInitialized = false
    end
    
    methods
        function set.SampleRate(obj, value)
            validateattributes(value,{'double','single'}, {'scalar',...
                'positive','finite'},...
                '','SampleRate');
            obj.SampleRate = value;
        end
        function set.RangeSidelobeAttenuation(obj,value)
            validateattributes( value, { 'double','single' },...
                { 'scalar', 'positive', 'finite' },...
                '', 'RangeSidelobeAttenuation');
            obj.RangeSidelobeAttenuation = value;
        end
        function set.CustomRangeWindow(obj,value)
            cond = ~isa(value,'function_handle') && ~isa(value,'cell');
            if cond
                coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidCustomWindow','CustomRangeWindow');
            end
            cond = isa(value,'cell') && ~isa(value{1},'function_handle');
            if cond
                coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidCustomWindowCell','CustomRangeWindow');
            end
            obj.CustomRangeWindow = value;
        end
        function set.PropagationSpeed(obj,val)
            sigdatatypes.validateSpeed(val,...
                '','PropagationSpeed',...
                {'double','single'},{'scalar','positive'});
            obj.PropagationSpeed = val;
        end
        function set.SweepSlope(obj, value)
            validateattributes(value,{'double','single'},...
                {'scalar','real','finite'},...
                '','SweepSlope');
            obj.SweepSlope = value;
        end
        function set.ReferenceRange(obj,val)
            sigdatatypes.validateDistance(val,'','ReferenceRange',...
                {'double','single'},{'scalar'});
            obj.ReferenceRange = val;
        end
        
    end

    methods

        function obj = RangeResponse(varargin)
            setProperties(obj, nargin, varargin{:});
        end
        
        function varargout = plotResponse(obj,x,varargin)
        %plotResponse   Plot range response
        %   plotResponse(H,X) plots the range response of the input signal,
        %   X, in dB scale. This syntax applies when you set the
        %   RangeMethod property to 'FFT' and then the DechirpInput
        %   property to false.
        %
        %   X must be a matrix where each column contains the dechirped
        %   signal from one frequency sweep.
        %
        %   plotResponse(H,X,XREF) plots the range response after
        %   performing a dechirp operation on X using the reference signal
        %   specified in XREF. This syntax applies when you set the
        %   RangeMethod property to 'FFT' and then the DechirpInput
        %   property to true.
        %
        %   X must be a matrix where each column contains the signal from
        %   one frequency sweep, where these signals have not yet been
        %   dechirped. XREF must be a column vector whose number of rows is
        %   the same as the number of rows of X.
        %
        %   plotResponse(H,X,COEFF) plots the range response after
        %   performing a matched filter operation on X using the matched
        %   filter coefficients specified in COEFF. This syntax applies
        %   when you set the RangeMethod property to 'Matched filter'.
        %
        %   X must be a matrix where each column contains the signal from
        %   one pulse. COEFF must be a column vector containing the matched
        %   filter coefficients.
        %
        %   plotResponse(...,Name,Value) plots the range response with the
        %   specified parameter Name set to the specified value.
        %   The parameter Names are
        %                   Unit: The unit of the plot, using one of 
        %                         | 'db' | 'mag' | 'pow' |. The default 
        %                         value is 'db'.
        %
        %   % Example:
        %   %   Plot the range response of an FMCW signal. The signal is
        %   %   not dechirped. The signal contains the return from one
        %   %   target which is approximately 2200 m away.
        %   
        %   % Load example data
        %   load('RangeResponseExampleData','fmcwdata');
        %   fs = fmcwdata.fs;
        %   propspeed = fmcwdata.propspeed;
        %   rxdata = fmcwdata.rxdata;
        %   refsig = fmcwdata.refsig;
        %   sweepslope = fmcwdata.sweepslope;
        %   
        %   % Create range response for dechirp and FFT processing
        %   rngresp = phased.RangeResponse(...
        %       'RangeMethod','FFT',...
        %       'SweepSlope',sweepslope,...
        %       'DechirpInput',true,...
        %       'SampleRate',fs,...
        %       'PropagationSpeed',propspeed);
        %   
        %   % Plot range response of processed data
        %   plotResponse(rngresp,rxdata,refsig,'Unit','db');
        %
        %   See also phased.RangeResponse, phased.RangeDopplerResponse,
        %   phased.RangeDopplerResponse/plotResponse,
        %   phased.AngleDopplerResponse/plotResponse.

            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','plotResponse');
            end

            narginchk(2,inf);
            
            if ~ismatrix(x) || isempty(x)
                error(message(...
                    'phased:step:MatrixNonEmpty','X'));
            end
            
            rngresp = clone(obj);
            release(rngresp);
            
            if (rngresp.RangeMethod(1) == 'F') && ... %Dechirp
                    ~rngresp.DechirpInput
                [resp,rng_grid] = step(rngresp,x);
            else
                if isempty(varargin)
                    if (rngresp.RangeMethod(1) == 'M') %Matched filter
                        error(message(...
                            'phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                            'Coeff','RangeMethod','Matched filter'));
                    else
                        error(message(...
                            'phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                            'XRef','DechirpInput','true'));
                    end
                end
                xref = varargin{1};
                varargin(1) = [];
                [resp,rng_grid] = step(rngresp,x,xref);
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
    
    methods (Access = protected)
        function resetImpl(obj)
            if (obj.RangeMethod(1) == 'M') %Matched filter
                reset(obj.cMatchedFilter);
            else
                reset(obj.cRangeFFT);
            end
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
            obj.pSizeInitialized = false;
            obj.pNumPages = -1;
            
            if (obj.RangeMethod(1) == 'M') %Matched filter
                release(obj.cMatchedFilter);
            else
                release(obj.cRangeFFT);
            end
        end
                
        function validateInputsImpl(obj,x,xref)
            cond = ~isa(x,'float');
            if cond
                coder.internal.errorIf(cond, ...
                  'MATLAB:system:invalidInputDataType','X','float');
            end
            cond = ndims(x)>3 || isempty(x);
            if cond               
                 coder.internal.errorIf(cond, ...
                     'phased:step:MatrixOr3DNonEmpty','X');
            end
            cond = size(x,1)<2;
            if cond               
                 coder.internal.errorIf(cond, ...
                     'phased:step:LowRowNum','X',2);
            end
            if (obj.RangeMethod(1) == 'M') %Matched filter
                cond = ~isa(xref,'float');
                if cond
                    coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDataType','Coeff','float');
                end
                cond = ~iscolumn(xref) || isempty(xref);
                if cond
                    coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeColVector','Coeff');
                end
            elseif obj.DechirpInput
                cond = ~isa(xref,'float');
                if cond
                    coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDataType','XRef','float');
                end
                cond = ~iscolumn(xref) || isempty(xref);
                if cond
                    coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeColVector','XRef');
                end
                cond = size(x,1)~=size(xref,1);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:NumRowsMismatch','X','XRef');
                end
            end
            
            validateNumChannels(obj,x);
            if ndims(x) == 3 && obj.pNumPages ~= -1
                validateNumPages(obj,x,obj.pNumPages);
            end             
        end
        
        function setupImpl(obj,x,xref)
            classtouse = class(x);
            coder.extrinsic('phased.RangeResponse.getWinCoeff');
            coder.extrinsic('fir1');
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            obj.pNumPages = size(x,3);
            obj.pValidatedNumInputPages = size(x,3);
            
            sz_x = size(x);
            
            if isempty(coder.target)
                cRW = obj.CustomRangeWindow;
            else
                %'function_handle' as property not supported in codegen
                cRW = [];
                cond = (obj.RangeWindow(2) == 'u'); %'Custom'
                if cond
                    coder.internal.errorIf(cond, ...
                     'phased:phased:MatchedFilter:NoCodegenCustom','RangeWindow','Custom');
                end
            end
            
            obj.pUseMatchedFilter = (obj.RangeMethod(1) == 'M'); %Matched filter
            
            fs = obj.SampleRate; % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            if obj.pUseMatchedFilter
                %obj.pSampleRate = getSampleRate(obj,sz_x(1),1,obj.SampleRate);
                obj.pSampleRate = fs;
                num_rng_samples = getPropagatedNumInputSamples(obj,x);
                
                if strcmp(obj.MaximumNumInputSamplesSource,'Auto')
                    maxNumInputSamples = num_rng_samples;
                else
                    maxNumInputSamples = obj.MaximumNumInputSamples;
                end
                
                obj.cMatchedFilter = phased.MatchedFilter(...
                    'CoefficientsSource','Input port',...
                    'MaximumNumInputSamplesSource','Property',...
                    'MaximumNumInputSamples',maxNumInputSamples);
                obj.pRangeGrid = cast(obj.PropagationSpeed*...
                    (0:maxNumInputSamples-1).'/obj.pSampleRate/2,classtouse);
                obj.pRangeOffset = numel(xref)-1;
            else
                if ~obj.DechirpInput
                    %obj.pSampleRate = getSampleRate(obj,sz_x(1),1,obj.SampleRate);
                    obj.pSampleRate = fs;
                    num_rng_samples = sz_x(1);
                else
                    %obj.pSampleRate = getSampleRate(obj,sz_x(1),obj.DecimationFactor,obj.SampleRate);
                    obj.pSampleRate = fs/obj.DecimationFactor;
                    num_rng_samples = ceil(sz_x(1)/obj.DecimationFactor);
                    
                    if obj.DecimationFactor > 1
                        obj.pDecimationFilterCoefficients = ...
                            coder.internal.const(fir1(30,1/obj.DecimationFactor));
                    end
                end
                

                obj.pRangeWinCoeff = coder.internal.const(...
                   obj.getWinCoeff(obj.RangeWindow,num_rng_samples,...
                                                 obj.RangeSidelobeAttenuation, ...
                                                 cRW));
                
                if (obj.RangeFFTLengthSource(1) == 'A')%Auto
                    obj.pRangeFFTLength = num_rng_samples;
                else
                    obj.pRangeFFTLength = obj.RangeFFTLength;
                end
                
                if obj.ReferenceRangeCentered
                    obj.pRangeGrid = beat2range(...
                        coder.internal.const(...
                            cast(obj.fftshiftfreqgrid(obj.pRangeFFTLength,obj.pSampleRate),classtouse)),...
                        obj.SweepSlope,obj.PropagationSpeed);
                else
                    if obj.SweepSlope > 0
                        obj.pRangeGrid = beat2range(...
                            cast(((0:obj.pRangeFFTLength-1).'/obj.pRangeFFTLength*obj.pSampleRate),classtouse),...
                            obj.SweepSlope,obj.PropagationSpeed);
                    else
                        obj.pRangeGrid = beat2range(...
                            cast((obj.pRangeFFTLength-1:-1:0).'/obj.pRangeFFTLength*obj.pSampleRate,classtouse),...
                            -obj.SweepSlope,obj.PropagationSpeed);
                    end
                end
                
                obj.cRangeFFT = dsp.FFT('FFTLengthSource','Property',...
                    'FFTLength',obj.pRangeFFTLength);
                
            end
        end
        
        function [rngresp,rnggrid] = stepImpl(obj,x,xref)
            
            classtouse = class(x);
            if ~obj.pSizeInitialized
                processInputSizeChangeImpl(obj,x);
                obj.pSizeInitialized = true;
            end
            
            sz_x = size(x);
            
            if obj.pUseMatchedFilter
                x_rng = complex(zeros(sz_x));
                for p = 1:obj.pNumPages
                    for m = 1:sz_x(2)
                        x_rng(:,m,p) = step(obj.cMatchedFilter,x(:,m,p),xref);
                    end
                end
                rngresp = complex(zeros(sz_x,classtouse));
                rngresp(1:sz_x(1)-obj.pRangeOffset,:) = x_rng(obj.pRangeOffset+1:sz_x(1),:);
                rnggrid = obj.pRangeGrid(1:sz_x(1));
            else
                if obj.DechirpInput
                    x_dechirp = reshape(dechirp(reshape(x,sz_x(1),[]),xref),sz_x);
                    if obj.DecimationFactor>1
                        x_rng = decfilt(obj.pDecimationFilterCoefficients,...
                                obj.DecimationFactor,x_dechirp);
                    else
                        x_rng = x_dechirp;
                    end
                else
                    x_rng = x;
                end
                r_fft_input = bsxfun(@times,obj.pRangeWinCoeff,x_rng);
                if obj.ReferenceRangeCentered
                    rngresp = fftshift(step(obj.cRangeFFT,complex(r_fft_input)),1);
                else
                    rngresp = step(obj.cRangeFFT,complex(r_fft_input));
                end
                rnggrid = obj.pRangeGrid;
            end
            rnggrid = rnggrid+cast(obj.ReferenceRange,classtouse);
        end
        
        function num = getNumInputsImpl(obj)
            num = 2;
            if strcmp(obj.RangeMethod,'FFT') && ...
                ~obj.DechirpInput
                num = num-1;
            end
        end
        
        function num = getNumOutputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if strcmp(obj.RangeMethod,'Matched filter')
                if strcmp(prop,'DechirpInput') || ...
                        strcmp(prop,'SweepSlope') || ...
                        strcmp(prop,'DecimationFactor') || ...
                        strcmp(prop,'RangeFFTLengthSource') || ...
                        strcmp(prop,'RangeFFTLength') || ...
                        strcmp(prop,'RangeWindow') || ...
                        strcmp(prop,'RangeSidelobeAttenuation') || ...
                        strcmp(prop,'CustomRangeWindow') || ...
                        strcmp(prop,'ReferenceRangeCentered')
                    flag = true;
                end
                if strcmp(prop,'MaximumNumInputSamples') && ...
                        strcmp(obj.MaximumNumInputSamplesSource,'Auto')
                    flag = true;
                end
            else
                if strcmp(prop,'DecimationFactor') && ...
                        ~obj.DechirpInput
                    flag = true;
                end
                
                if strcmp(prop,'RangeFFTLength') && ...
                        strcmp(obj.RangeFFTLengthSource,'Auto')
                    flag = true;
                end
                if strcmp(prop,'CustomRangeWindow') && ...
                        ~strcmp(obj.RangeWindow,'Custom')
                    flag = true;
                end
                if strcmp(prop,'RangeSidelobeAttenuation') && ...
                        ~(strcmp(obj.RangeWindow,'Chebyshev') || ...
                        strcmp(obj.RangeWindow,'Taylor') || ...
                        strcmp(obj.RangeWindow,'Kaiser'))
                    flag = true;
                end
                if strcmp(prop,'MaximumNumInputSamples') || ...
                        strcmp(prop,'MaximumNumInputSamplesSource')
                    flag = true;
                end
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            if isLocked(obj)
                s.pSizeInitialized = obj.pSizeInitialized;
                s.pNumPages = obj.pNumPages;
                s.pValidatedNumInputPages = obj.pValidatedNumInputPages;
                s.pUseMatchedFilter = obj.pUseMatchedFilter;
                s.pRangeWinCoeff = obj.pRangeWinCoeff;
                s.pRangeGrid = obj.pRangeGrid;
                s.pSampleRate = obj.pSampleRate;
                s.pRangeOffset = obj.pRangeOffset;
                s.pRangeFFTLength = obj.pRangeFFTLength;
                s.pDecimationFilterCoefficients = obj.pDecimationFilterCoefficients; 
                s.cMatchedFilter = saveobj(obj.cMatchedFilter);
                s.cRangeFFT = saveobj(obj.cRangeFFT);
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                obj.cMatchedFilter = phased.MatchedFilter.loadobj(s.cMatchedFilter);
                s = rmfield(s,'cMatchedFilter');
                if isfield(s,'cRangeFFT')
                    obj.cRangeFFT = dsp.FFT.loadobj(s.cRangeFFT);
                    s = rmfield(s,'cRangeFFT');
                end
                % recover locked sample rate information
                if isfield(s,'pSampleRate')
                    obj.pSampleRate = s.pSampleRate;
                    s = rmfield(s,'pSampleRate');
                else
                    obj.pSampleRate = s.SampleRate;
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

    end
    
    methods (Static,Hidden,Access=protected)  
        function groups = getPropertyGroupsImpl
            groups = matlab.system.display.Section(...
                'PropertyList',{'RangeMethod','PropagationSpeed',...
                'SampleRateFromInputCheckbox','SampleRate','SweepSlope','DechirpInput',...
                'DecimationFactor','RangeFFTLengthSource',...
                'RangeFFTLength','RangeWindow','RangeSidelobeAttenuation',...
                'CustomRangeWindow','ReferenceRangeCentered','ReferenceRange','MaximumNumInputSamplesSource',...
                'MaximumNumInputSamples'});
            dDecimationFactor = matlab.system.display.internal.Property(...
                'DecimationFactor','IsObjectDisplayOnly',true); % need multirate support
            dRangeWindow = matlab.system.display.internal.Property(...
                'RangeWindow', 'StringSetValues', {'None','Hamming',...
                'Chebyshev','Hann','Kaiser','Taylor'});
            dMaximumNumInputSamplesSource = matlab.system.display.internal.Property(...
                'MaximumNumInputSamplesSource','IsGraphical',false);
            dMaximumNumInputSamples = matlab.system.display.internal.Property(...
                'MaximumNumInputSamples','IsGraphical',false);
            for m = 1:numel(groups.PropertyList)
                if strcmp(groups.PropertyList{m},'DecimationFactor')
                    groups.PropertyList{m} = dDecimationFactor;
                elseif strcmp(groups.PropertyList{m},'RangeWindow')
                    groups.PropertyList{m} = dRangeWindow;
                elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamplesSource')
                    groups.PropertyList{m} = dMaximumNumInputSamplesSource;
                elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamples')
                    groups.PropertyList{m} = dMaximumNumInputSamples;
                end
            end
        end
        
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:RangeResponseTitle')),...
                'Text',getString(message('phased:library:block:RangeResponseDesc')));
        end
    end

    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)
            if strcmp(obj.RangeMethod,'Matched filter')
                varargout = {'X','Coeff'};
            else
                if obj.DechirpInput
                    varargout = {'X','XRef'};
                else
                    varargout = {'X'};
                end
            end
        end
        
        function varargout = getOutputNamesImpl(~)
                    varargout = {'Resp','Range'};
        end
               
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Range Response');
        end
        function flag = isInputSizeLockedImpl(obj,ind)
            if strcmp(obj.RangeMethod,'Matched filter') && ind==1
                flag = false;
            else
                flag = true;
            end
        end
        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            if strcmp(obj.RangeMethod,'Matched filter')
                P = szX(1);
            else
                if obj.RangeFFTLengthSource(1) == 'A' %Auto
                    P = szX(1);
                else
                    P = obj.RangeFFTLength;
                end
            end
            if numel(szX)==3
                respSz = [P szX(2:3)];
            else
                respSz = [P szX(2)];
            end
            varargout = {respSz,[P 1]};
        end
        function varargout = isOutputFixedSizeImpl(obj) 
            varargout = {propagatedInputFixedSize(obj, 1), propagatedInputFixedSize(obj, 1)};
        end
        function varargout = getOutputDataTypeImpl(obj)  %#ok<MANU>
            varargout = {propagatedInputDataType(obj,1),propagatedInputDataType(obj,1)};
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout = {true,false};
        end            
    end
    
    methods (Static, Hidden)
        
        function winCoeff = getWinCoeff(rangeWindow,num_rng_samples,...
                        rangeSidelobeAttenuation, customRangeWindow)
                switch rangeWindow
                    case 'None'
                        winCoeff = ones(num_rng_samples,1);
                    case 'Hamming'
                        winCoeff = hamming(num_rng_samples);
                    case 'Hann'
                        winCoeff = hann(num_rng_samples);
                    case 'Kaiser'
                        winCoeff = kaiser(num_rng_samples,...
                            signal.internal.kaiserBeta(rangeSidelobeAttenuation));
                    case 'Chebyshev'
                        winCoeff = chebwin(num_rng_samples,...
                            rangeSidelobeAttenuation);
                    case 'Taylor'
                        winCoeff = taylorwin(num_rng_samples,...
                            4,-rangeSidelobeAttenuation);
                    case 'Custom'
                        if isa(customRangeWindow,'function_handle')
                            winCoeff = ...
                                customRangeWindow(num_rng_samples);
                        else
                            winCoeff = ...
                                customRangeWindow{1}(num_rng_samples,...
                                customRangeWindow{2:end});
                        end
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
        
            % freq_grid = fftshift(psdfreqvec(...
            %     'Npts',N,'Fs',Fs,'Range','whole'));
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

function odata = decfilt(b,r,idata)
%decfilt    Decimation filter
%   y = decfilt(b,r,x) performs decimation filter on input x. b is the
%   FIR decimation filter coefficients and r is the decimation factor.
%
%   % Example:
%   %   Perform a decimation filter with factor of 4.
%
%   t = 0:1/100:1; x = sin(2*pi*10*t);
%   r = 4; b = fir1(10,1/r); 
%   y = decfilt(b,r,x);
%   plot(t,x,t(1:r:end),y);

nfilt = numel(b)+1;

szin = size(idata);
nd = szin(1);
nc = prod(szin(2:end));
idata = idata(:,:);

nout = ceil(nd/r);
itemp = bsxfun(@minus,2*idata(1,:),idata((nfilt+1):-1:2,:));
[~,zi]=filter(b,1,itemp);
[odata_temp,zf] = filter(b,1,idata,zi);
itemp2 = complex(zeros(2*nfilt,nc));
itemp2(:) = bsxfun(@minus,2*idata(nd,:),idata((nd-1):-1:(nd-2*nfilt),:));
itemp_out = filter(b,1,itemp2,zf);
% finally, select only every r'th point from the interior of the lowpass
% filtered sequence
gd = (numel(b)-1)/2;
list = round(gd(1)+1.25):r:nd;
odata_temp1 = odata_temp(list,:);
lod = size(odata_temp1,1);
nlen = nout - lod;
nbeg = r - (nd - list(numel(list)));
odata = [odata_temp1; itemp_out(nbeg:r:nbeg+nlen*r-1,:)];
szout = [nout szin(2:end)];
odata = reshape(odata,szout);
end

% [EOF]
