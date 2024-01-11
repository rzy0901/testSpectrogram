classdef (Sealed,StrictDefaults) RangeDopplerResponse < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%RangeDopplerResponse   Range-Doppler response
%   H = phased.RangeDopplerResponse creates a range-Doppler response
%   System object, H. This object calculates the range-Doppler response of
%   the input data.
%
%   H = phased.RangeDopplerResponse(Name,Value) creates a range-Doppler
%   response object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The RangeDopplerResponse object generates the response by first
%   processing the input signal in the range domain using either a matched
%   filter or dechirp operation, and then processing in the Doppler domain
%   using an FFT.
%
%   Step method syntax:
%   
%   [RESP,RANGE,DOP] = step(H,X) calculates the range-Doppler response of
%   the input data X. This syntax applies when you set the RangeMethod
%   property to 'FFT' and then the DechirpInput property to false. This
%   syntax is most commonly used with FMCW signals.
%
%   X is a dechirped input signal. X must be a KxL matrix or a KxNxL array
%   where K denotes the number of fast time samples, L is the number of
%   dechirped frequency sweeps, and N is the number of channels (antenna
%   elements or beams). When X is a 3-dimensional array, each column of X
%   is processed as an independent KxL signal.
%
%   By default, all sweeps in X are assumed to be consecutive. If the
%   sweeps are not consecutive, the PRF of the input data must be provided
%   by setting the PRFSource property to either 'Property' or 'Input port'.
% 
%   RESP is either an MxP matrix or an MxNxP array containing the complex
%   range-Doppler response of the input X. The number of dimensions in RESP
%   will match the number of dimensions in X. When X and RESP are
%   3-dimensional arrays, N is the number of channels. M is determined by
%   either the number of rows in X when you set the RangeFFTLengthSource
%   property to 'Auto', or the value specified in the RangeFFTLength
%   property when you set the RangeFFTLengthSource property to 'Property'.
%   P is determined by either the number of frequency sweeps in X when you
%   set the DopplerFFTLengthSource property to 'Auto', or the value
%   specified in the DopplerFFTLength property when you set the
%   DopplerFFTLengthSource property to 'Property'.
%
%   RANGE is a length M column vector containing the range samples at which
%   the range-Doppler response is evaluated. DOP is a length P column
%   vector containing either Doppler or speed samples at which the
%   range-Doppler response is evaluated. Whether DOP contains Doppler or
%   speed samples depends on the DopplerOutput property.
%
%   [RESP,RANGE,DOP] = step(H,X,XREF) uses input XREF as the reference
%   signal to dechirp the input signal X. This syntax applies when you set
%   the RangeMethod property to 'FFT' and then the DechirpInput property to
%   true. This syntax is most commonly used with FMCW signals and the
%   reference signal is, in general, the transmitted signal.
%
%   X is an input signal to be dechirped by the RangeDopplerResponse
%   object. X must be a KxL matrix or an KxNxL array where K denotes the
%   number of fast time samples, L is the number of frequency sweeps, and N
%   is the number of channels (antenna elements or beams). XREF must be a
%   column vector whose number of rows is the same as the number of rows of
%   X. XREF is the reference signal used to dechirp X.
%
%   By default, all sweeps in X are assumed to be consecutive. If the
%   sweeps are not consecutive, the PRF of the input data must be provided
%   by setting the PRFSource property to either 'Property' or 'Input port'.
%
%   In this syntax, the number of rows in RESP is the quotient of the
%   number of rows in X and the value specified in the DecimationFactor
%   property.
%
%   [RESP,RANGE,DOP] = step(H,X,COEFF) uses COEFF as the matched filter
%   coefficients. This method applies when you set the RangeMethod property
%   to 'Matched filter'. This syntax is most commonly used with pulsed
%   signals.
%
%   X is an input signal to be match filtered by the RangeDopplerResponse
%   object. X must be a KxL matrix or a KxNxL array where K denotes the
%   number of fast time samples, L is the number of pulses, and N is the
%   number of channels (antenna elements or beams). COEFF must be a column
%   vector containing the matched filter coefficients.
%
%   By default, all pulses in X are assumed to be consecutive. If the
%   pulses are not consecutive, the PRF of the input data must be provided
%   by setting the PRFSource property to either 'Property' or 'Input port'.
%
%   In this syntax, RESP is either a KxP matrix or a KxNxP array, where the
%   number of rows of RESP equals the number of rows in X. When X is a
%   3-dimensional array, N is the number of channels. P is determined by
%   either the number of pulses in X when you set the
%   DopplerFFTLengthSource property to 'Auto', or the value specified in
%   the DopplerFFTLength property when you set the DopplerFFTLengthSource
%   property to 'Property'.
%
%   RANGE is a length K column vector containing the range samples at which
%   the range response is evaluated. DOP is a length P column vector
%   containing either Doppler or speed samples at which the range-Doppler
%   response is evaluated. Whether DOP contains Doppler or speed samples
%   depends on the DopplerOutput property.
%   
%   [RESP,RANGE,DOP] = step(H,X,...,PRF) uses PRF as the pulse repetition
%   frequency for the signal in X. This method applies when you set the
%   PRFSource property to 'Input port'.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   RangeDopplerResponse methods:
%
%   step         - Calculate range-Doppler response
%   release      - Allow property value and input characteristics changes
%   clone        - Create range-Doppler response object with same property
%                  values
%   isLocked     - Locked status (logical)
%   plotResponse - Plot range-Doppler response
%
%   RangeDopplerResponse properties:
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
%   PRFSource                    - Source of pulse repetition frequency
%   PRF                          - Pulse repetition frequency used in
%                                  Doppler processing
%   DopplerFFTLengthSource       - Source of FFT length in Doppler 
%                                  processing
%   DopplerFFTLength             - FFT length in Doppler processing
%   DopplerWindow                - Doppler processing window
%   DopplerSidelobeAttenuation   - Doppler sidelobe attenuation level
%   CustomDopplerWindow          - Custom Doppler processing window
%   DopplerOutput                - Doppler output
%   OperatingFrequency           - Signal carrier frequency
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
%   %   Calculate the range-Doppler response from a pulsed radar
%   %   transmitting a rectangular waveform using the matched filter
%   %   approach. The signal includes three target returns. Two are
%   %   approximately 2000 m away and the third is approximately 3500 m
%   %   away. In addition, two targets are stationary relative to the radar
%   %   while the third is moving away from the radar at approximately
%   %   100 m/s.
%   
%   % Load example data
%   load('RangeDopplerResponseExampleData','rectdata');
%   fs = rectdata.fs;
%   propspeed = rectdata.propspeed;
%   fc = rectdata.fc;
%   rxdata = rectdata.rxdata;
%   mfcoeffs = rectdata.mfcoeffs;
%   noisepower = rectdata.noisepower;
%   
%   % Create range-Doppler response for matched filter processing.
%   % Interpolate to 1024 Doppler bins.
%   rngdopresp = phased.RangeDopplerResponse(...
%       'DopplerFFTLengthSource','Property',...
%       'DopplerFFTLength',1024,...
%       'DopplerOutput','Speed',...
%       'OperatingFrequency',fc,...
%       'SampleRate',fs,...
%       'PropagationSpeed',propspeed);
%   
%   % Perform range-Doppler processing
%   [resp,rng_grid,dop_grid] = rngdopresp(rxdata,mfcoeffs);
%   
%   % Calculate noise power after range-Doppler processing
%   mfgain = mfcoeffs'*mfcoeffs;
%   noisepower_proc = mfgain*noisepower;
%   
%   imagesc(dop_grid,rng_grid,pow2db(abs(resp).^2/noisepower_proc));
%   xlabel('Speed (m/s)'); ylabel('Range (m)'); title('Range Doppler Map');
%   set(get(colorbar,'YLabel'),'String','SNR (dB)'); caxis([0 60]);
%
%   See also phased, phased.RangeResponse, phased.AngleDopplerResponse,
%   phased.MatchedFilter, dechirp.

%   Copyright 2012-2017 The MathWorks, Inc.

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
    
    properties (Nontunable)
        %PRFSource  Source of pulse repetition frequency
        %   Specify how to determine the pulse repetition frequency (PRF)
        %   of the processed input signal as one of 'Auto' | 'Property' |
        %   'Input port', where the default is 'Auto'. When you set this
        %   property to 'Auto', the PRF is inferred from the number of rows
        %   in the input signal and the value of the SampleRate property.
        %   When you set this property to 'Property', the PRF is specified
        %   in the PRF property. When you set this property to 'Input
        %   port', the PRF is specified as an input argument.
        PRFSource = 'Auto'
    end
    
    properties (Nontunable)
        %PRF  Pulse repetition frequency of the input signal (Hz)
        %   Specify the pulse repetition frequency (PRF) of the input
        %   signal (in Hz) as a positive scalar. The default value of this
        %   property is 10e3 (10 kHz). This property only applies when you
        %   set the PRFSource property to 'Property'.
        PRF = 10e3
    end
    
    properties (Nontunable)
        %DopplerFFTLengthSource  Source of FFT length in Doppler processing
        %   Specify how to determine the FFT length in Doppler processing
        %   as one of 'Auto' | 'Property', where the default is 'Auto'.
        %   When you set this property to 'Auto', the FFT length is equal
        %   to the number of rows of the input signal. When you set this
        %   property to 'Property', the FFT length is specified in
        %   DopplerFFTLength property.
        DopplerFFTLengthSource = 'Auto'
    end
    
    properties (Nontunable, PositiveInteger)
        %DopplerFFTLength     FFT length in Doppler processing
        %   Specify the FFT length in Doppler processing as a positive
        %   integer. This property applies when you set the
        %   DopplerFFTLengthSource property to 'Property'. The default
        %   value of this property is 1024.
        DopplerFFTLength = 1024
    end
    
    properties (Nontunable)
        %DopplerWindow    Doppler processing window
        %   Specify the window used for Doppler processing using one of
        %   'None' | 'Hamming' | 'Chebyshev' | 'Hann' | 'Kaiser' | 'Taylor'
        %   | 'Custom', where the default is 'None'.
        %
        %   Note that if you set the DopplerWindow property to 'Taylor',
        %   the generated Taylor window has 4 nearly constant sidelobes
        %   adjacent to the mainlobe.
        DopplerWindow = 'None'
        %DopplerSidelobeAttenuation  Doppler sidelobe attenuation level 
        %   Specify the sidelobe attenuation level (in dB) of a Kaiser,
        %   Chebyshev or Taylor window in Doppler processing as a positive
        %   scalar. This property applies when you set the DopplerWindow
        %   property to 'Kaiser', 'Chebyshev', or 'Taylor'. The default
        %   value of this property is 30.
        DopplerSidelobeAttenuation = 30
        %CustomDopplerWindow   Custom Doppler processing window
        %   Specify the user-defined window for Doppler processing using a
        %   function handle or a cell array. The default value of this
        %   property is @hamming. This property applies when you set the
        %   DopplerWindow property to 'Custom'.
        %
        %   If CustomDopplerWindow is a function handle, the specified
        %   function takes the window length as the input and generates
        %   appropriate window coefficients.
        %
        %   If CustomDopplerWindow is a cell array, then the first cell
        %   must be a function handle. The specified function takes the
        %   window length as the first input argument, with other
        %   additional input arguments if necessary, and generates
        %   appropriate window coefficients. The remaining entries in the
        %   cell array are the additional input arguments to the function,
        %   if any.
        CustomDopplerWindow = @hamming
        %DopplerOutput  Doppler output
        %   Specify the Doppler domain output as one of 'Frequency' |
        %   'Speed', where the default is 'Frequency'. If you set this
        %   property to 'Frequency', the Doppler domain output of step,
        %   DOP, is the Doppler shift (in Hz). If you set this
        %   property to 'Speed', the Doppler domain output is the
        %   corresponding radial speed (in m/s).
        DopplerOutput = 'Frequency'
    end
    
    properties (Nontunable)
        %OperatingFrequency     Signal carrier frequency (Hz)
        %   Specify the signal carrier frequency (in Hz) as a scalar. This
        %   property applies when you set the DopplerOutput property to
        %   'Speed'.The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 3e8;
    end
    
    properties (Constant, Hidden)
        RangeMethodSet = matlab.system.internal.StringSetGF(...
            {'Matched filter','FFT'},{'Dechirp'},{'FFT'});
        DopplerOutputSet = matlab.system.StringSet({'Frequency','Speed'});
        RangeWindowSet = matlab.system.StringSet({'None','Hamming',...
            'Chebyshev','Hann','Kaiser','Taylor','Custom'});
        DopplerWindowSet = matlab.system.StringSet({'None','Hamming',...
            'Chebyshev','Hann','Kaiser','Taylor','Custom'});
        RangeFFTLengthSourceSet = dsp.CommonSets.getSet(...
            'AutoOrProperty');
        DopplerFFTLengthSourceSet = dsp.CommonSets.getSet(...
            'AutoOrProperty');
        MaximumNumInputSamplesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        PRFSourceSet = matlab.system.StringSet({'Auto','Property','Input port'});
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties (Access = private, Nontunable, Logical)
        pUseMatchedFilter = false
        pOutputSpeed = false;
    end
    
    properties (Access = private)
        cRangeResponse
        cDopplerFFT
    end
    
    properties (Access = private, Nontunable)
        pDopplerWinCoeff
        pSampleRate
        pDopplerFFTLength
        pValidatedNumInputPages
    end
    
    properties (Access = private)
        pDopplerGrid
        pDopplerSampleRate
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
        function set.OperatingFrequency(obj,val)
            
            sigdatatypes.validateFrequency(val,'',...
                'OperatingFrequency',{'double','single'},{'scalar'});
            obj.OperatingFrequency = val;
        end
        function set.PRF(obj, value)
            validateattributes(value,{'double','single'}, {'scalar',...
                'positive','finite'},...
                '','PRF');
            obj.PRF = value;
        end        
        function set.RangeSidelobeAttenuation(obj,value)
            validateattributes( value, { 'double','single' },...
                { 'scalar', 'positive', 'finite' },...
                '', 'RangeSidelobeAttenuation');
            obj.RangeSidelobeAttenuation = value;
        end
        function set.DopplerSidelobeAttenuation(obj,value)
            validateattributes( value, { 'double','single' },...
                { 'scalar', 'positive', 'finite' },...
                '', 'DopplerSidelobeAttenuation');
            obj.DopplerSidelobeAttenuation = value;
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
        function set.CustomDopplerWindow(obj,value)
            cond = ~isa(value,'function_handle') && ~isa(value,'cell');
            if cond
                coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidCustomWindow','CustomDopplerWindow');
            end
            cond = isa(value,'cell') && ~isa(value{1},'function_handle');
            if cond
                coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidCustomWindowCell','CustomDopplerWindow');
            end
            obj.CustomDopplerWindow = value;
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

        function obj = RangeDopplerResponse(varargin)
            setProperties(obj, nargin, varargin{:});
        end
        
        function varargout = plotResponse(obj,x,varargin)
        %plotResponse   Plot range-Doppler response
        %   plotResponse(H,X) plots the range-Doppler response of the input
        %   signal, X, in dB scale. This syntax applies when you set the
        %   RangeMethod property to 'FFT' and then the DechirpInput
        %   property to false.
        %
        %   X must be a matrix where each column contains the dechirped
        %   signal from one frequency sweep.
        %
        %   By default, all sweeps in X are assumed to be consecutive. If
        %   the sweeps are not consecutive, the PRF of the input data must
        %   be provided by setting the PRFSource property to either
        %   'Property' or 'Input port'.
        %
        %   plotResponse(H,X,XREF) plots the range-Doppler response after
        %   performing a dechirp operation on X using the reference signal
        %   specified in XREF. This syntax applies when you set the
        %   RangeMethod property to 'FFT' and then the DechirpInput
        %   property to true.
        %
        %   X must be a matrix where each column contains the signal from
        %   one frequency sweep. These signals are not yet dechirped. XREF
        %   must be a column vector whose number of rows is the same as the
        %   number of rows of X.
        %
        %   By default, all sweeps in X are assumed to be consecutive. If
        %   the sweeps are not consecutive, the PRF of the input data must
        %   be provided by setting the PRFSource property to either
        %   'Property' or 'Input port'.
        %
        %   plotResponse(H,X,COEFF) plots the range-Doppler response after
        %   performing a matched filter operation on X using the matched
        %   filter coefficients specified in COEFF. This syntax applies
        %   when you set the RangeMethod property to 'Matched filter'.
        %
        %   X must be a matrix where each column contains the signal from
        %   one pulse. COEFF must be a column vector containing the matched
        %   filter coefficients.
        %
        %   By default, all pulses in X are assumed to be consecutive. If
        %   the pulses are not consecutive, the PRF of the input data must
        %   be provided by setting the PRFSource property to either
        %   'Property' or 'Input port'.
        %
        %   plotResponse(H,X,...,PRF) uses PRF as the pulse repetition
        %   frequency for the signal in X. This method applies when you set
        %   the PRFSource property to 'Input port'.
        %
        %   plotResponse(...,Name,Value) plots the range-Doppler response
        %   with the specified parameter Name set to the specified value.
        %   The parameter Names are
        %                   Unit: The unit of the plot, using one of 
        %                         | 'db' | 'mag' | 'pow' |. The default 
        %                         value is 'db'.
        %       NormalizeDoppler: Set this to true to normalize the Doppler
        %                         frequencies. Set this to false to plot
        %                         the range-Doppler response without 
        %                         normalizing the Doppler frequency. The
        %                         default value is false. This parameter
        %                         applies when you set the DopplerOutput
        %                         property to 'Frequency'.
        %
        %   % Example:
        %   %   Plot the range-Doppler response of an FMCW signal. The
        %   %   signal is not dechirped. The signal contains the return
        %   %   from one target which is approximately 2200 m away and has
        %   %   a normalized Doppler frequency of about -0.36 relative to
        %   %   the radar.
        %   
        %   % Load example data
        %   load('RangeDopplerResponseExampleData','fmcwdata');
        %   fs = fmcwdata.fs;
        %   propspeed = fmcwdata.propspeed;
        %   rxdata = fmcwdata.rxdata;
        %   refsig = fmcwdata.refsig;
        %   sweepslope = fmcwdata.sweepslope;
        %   
        %   % Create range-Doppler response for dechirp and FFT processing
        %   rngdopresp = phased.RangeDopplerResponse(...
        %       'RangeMethod','FFT',...
        %       'SweepSlope',sweepslope,...
        %       'DechirpInput',true,...
        %       'SampleRate',fs,...
        %       'PropagationSpeed',propspeed);
        %   
        %   % Plot range-Doppler response of processed data
        %   plotResponse(rngdopresp,rxdata,refsig, ...
        %       'Unit','db','NormalizeDoppler',true);
        %
        %   See also phased.RangeDopplerResponse,
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
            
            hrdresp = clone(obj);
            release(hrdresp);
            
            if (hrdresp.RangeMethod(1) == 'F') && ... %Dechirp
                    ~hrdresp.DechirpInput
                if (hrdresp.PRFSource(1) == 'I') && ... % Input port
                        isempty(varargin)
                    error(message(...
                        'phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                        'PRF','PRFSource','Input port'));
                end
                if (hrdresp.PRFSource(1) == 'I')
                    prf = varargin{1};
                    varargin(1) = [];
                    [resp,rng_grid,dop_grid] = step(hrdresp,x,prf);
                else
                    [resp,rng_grid,dop_grid] = step(hrdresp,x);
                end
            else
                if isempty(varargin)
                    if (hrdresp.RangeMethod(1) == 'M') %Matched filter
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
                if (hrdresp.PRFSource(1) == 'I') && ... % Input port
                        isempty(varargin)
                    error(message(...
                        'phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                        'PRF','PRFSource','Input port'));
                end
                if (hrdresp.PRFSource(1) == 'I')
                    prf = varargin{1};
                    varargin(1) = [];
                    [resp,rng_grid,dop_grid] = step(hrdresp,x,xref,prf);
                else
                    [resp,rng_grid,dop_grid] = step(hrdresp,x,xref);
                end
            end
            
            unit = 'db';
            normalizedoppler = [];
            sigutils.pvparse(varargin{:});
            unit = validatestring(unit,{'db','mag','pow'},...
                'plotResponse','Unit');
            
            if ~hrdresp.pOutputSpeed
                if isempty(normalizedoppler)
                    normalizedoppler = false;
                    validateattributes(normalizedoppler,{'logical'},{'scalar'},...
                        'plotResponse','NormalizeDoppler');
                end
                if (hrdresp.RangeMethod(1) == 'M') %Matched filter
                    fs_dop = hrdresp.SampleRate/numel(rng_grid);
                else
                    if ~hrdresp.DechirpInput
                        fs_dop = hrdresp.SampleRate/numel(rng_grid);
                    else
                        fs_dop = hrdresp.SampleRate/numel(rng_grid)/...
                            hrdresp.DecimationFactor;
                    end
                end
                hresp = phased.internal.RangeDopplerPattern(...
                    'Doppler',dop_grid(:)',...
                    'Range',rng_grid(:)',...
                    'Pattern',abs(resp),...
                    'PRF',fs_dop);
                if nargout
                    varargout{1} = ...
                        plot(hresp,'Units',unit,'NormalizeDoppler',normalizedoppler);
                else
                    plot(hresp,'Units',unit,'NormalizeDoppler',normalizedoppler);
                end
            else
                if ~isempty(normalizedoppler)
                    error(message(...
                        'phased:AngleDopplerResponse:plotResponse:InvalidParameter',...
                        'NormalizeDoppler','DopplerOutput','Speed'));
                end
                hresp = phased.internal.RangeSpeedPattern(...
                    'Speed',dop_grid(:)',...
                    'Range',rng_grid(:)',...
                    'Pattern',abs(resp));
                if nargout
                    varargout{1} = ...
                        plot(hresp,'Units',unit);
                else
                    plot(hresp,'Units',unit);
                end
            end
             
        end

    end
    
    methods (Access = protected)
        function resetImpl(obj)
            reset(obj.cRangeResponse);
            reset(obj.cDopplerFFT);
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
            obj.pSizeInitialized = false;
            obj.pNumPages = -1;
            
            release(obj.cRangeResponse);
            release(obj.cDopplerFFT);
        end

        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object
            obj.cRangeResponse.ReferenceRange = obj.ReferenceRange;
        end
                
        function validateInputsImpl(obj,x,opt1,opt2)
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
            cond = size(x,2)<2;
            if cond               
                 coder.internal.errorIf(cond, ...
                     'phased:step:LowColNum','X',2);
            end
            if (obj.RangeMethod(1) == 'M') %Matched filter
                coeff = opt1;
                cond = ~isa(coeff,'float');
                if cond
                    coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDataType','Coeff','float');
                end
                cond = ~iscolumn(coeff) || isempty(coeff);
                if cond
                    coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeColVector','Coeff');
                end
                if strcmp(obj.PRFSource,'Input port')
                    prf = opt2;
                    validateattributes(prf,{'double','single'},{'real','scalar','positive','finite'},'','PRF');
                end
            else
                if obj.DechirpInput
                    xref = opt1;
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
                    if strcmp(obj.PRFSource,'Input port')
                        prf = opt2;
                        validateattributes(prf,{'double','single'},{'real','scalar','positive','finite'},'','PRF');
                    end
                else
                    if strcmp(obj.PRFSource,'Input port')
                        prf = opt1;
                        validateattributes(prf,{'double','single'},{'real','scalar','positive','finite'},'','PRF');
                    end
                end
            end
            
            validateNumChannels(obj,x);
            if ndims(x) == 3 && obj.pNumPages ~= -1
                validateNumPages(obj,x,obj.pNumPages);
            end             
        end
        
        function setupImpl(obj,x,~,prf)
            classtouse = class(x);
            coder.extrinsic('phased.RangeDopplerResponse.getWinCoeff');
            
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputPages = size(x,3);
            
            obj.pNumInputChannels = obj.pValidatedNumInputChannels;
            obj.pNumPages = obj.pValidatedNumInputPages;
            
            if obj.pValidatedNumInputPages>1 % is a cube
                num_dop_samples = obj.pValidatedNumInputPages;
            else
                num_dop_samples = obj.pValidatedNumInputChannels;
            end
            
            fs = obj.SampleRate; % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            
            if strcmp(obj.RangeMethod,'Matched filter')
                
                if strcmp(obj.MaximumNumInputSamplesSource,'Auto')
                    maxNumInputSamples = getPropagatedNumInputSamples(obj,x);
                else
                    maxNumInputSamples = obj.MaximumNumInputSamples;
                end
                
                obj.cRangeResponse = phased.RangeResponse(...
                    'RangeMethod',obj.RangeMethod,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'SampleRate',fs,...
                    'MaximumNumInputSamplesSource','Property',...
                    'MaximumNumInputSamples',maxNumInputSamples,...
                    'ReferenceRange',obj.ReferenceRange);
            else
                obj.cRangeResponse = phased.RangeResponse(...
                    'RangeMethod',obj.RangeMethod,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'SampleRate',fs,...
                    'SweepSlope',obj.SweepSlope,...
                    'DechirpInput',obj.DechirpInput,...
                    'RangeFFTLengthSource',obj.RangeFFTLengthSource,...
                    'RangeWindow',obj.RangeWindow,...
                    'ReferenceRangeCentered',obj.ReferenceRangeCentered,...
                    'ReferenceRange',obj.ReferenceRange);
                if obj.DechirpInput
                    obj.cRangeResponse.DecimationFactor = obj.DecimationFactor;
                end
                if ~strcmp(obj.RangeFFTLengthSource,'Auto')
                    obj.cRangeResponse.RangeFFTLength = obj.RangeFFTLength;
                end
                if strcmp(obj.RangeWindow,'Custom')
                    obj.cRangeResponse.CustomRangeWindow = obj.CustomRangeWindow;
                end
                if (strcmp(obj.RangeWindow,'Chebyshev') || ...
                        strcmp(obj.RangeWindow,'Taylor') || ...
                        strcmp(obj.RangeWindow,'Kaiser'))
                    obj.cRangeResponse.RangeSidelobeAttenuation = obj.RangeSidelobeAttenuation;
                end
            end
            
            if isempty(coder.target)
                cDW = obj.CustomDopplerWindow;
            else
                %'function_handle' as property not supported in codegen
                cDW = [];
                cond = (obj.DopplerWindow(2) == 'u'); %'Custom'
                if cond
                    coder.internal.errorIf(cond, ...
                     'phased:phased:MatchedFilter:NoCodegenCustom','DopplerWindow','Custom');
                end

            end
            
            obj.pUseMatchedFilter = (obj.RangeMethod(1) == 'M'); %Matched filter
            obj.pOutputSpeed = (obj.DopplerOutput(1) == 'S'); %Speed
            
            sz_x = size(x);
            if obj.pUseMatchedFilter
                obj.pSampleRate = fs;
                num_rng_samples = getPropagatedNumInputSamples(obj,x);
            else
                if ~obj.DechirpInput
                    obj.pSampleRate = fs;
                    num_rng_samples = sz_x(1);
                else
                    obj.pSampleRate = fs/obj.DecimationFactor;
                    num_rng_samples = ceil(sz_x(1)/obj.DecimationFactor);
                end
            end
            
            switch obj.PRFSource
                case 'Auto'
                    prf = obj.pSampleRate/num_rng_samples;
                case 'Property'
                    prf = obj.PRF;
                    checkPRF(obj,x,prf);
                case 'Input port'
                    prf = 1; % Normalized frequency
            end
            obj.pDopplerSampleRate = prf;
            obj.pDopplerWinCoeff = coder.internal.const(...
                   obj.getWinCoeff(obj.DopplerWindow,num_dop_samples,...
                             obj.DopplerSidelobeAttenuation, cDW));

            
            if (obj.DopplerFFTLengthSource(1) == 'A') %Auto
                obj.pDopplerFFTLength = num_dop_samples;
            else
                obj.pDopplerFFTLength = obj.DopplerFFTLength;
            end
            
            obj.cDopplerFFT = dsp.FFT('FFTLengthSource','Property',...
                'FFTLength',obj.pDopplerFFTLength);

            dopplerGridtemp = cast(obj.fftshiftfreqgrid(obj.pDopplerFFTLength,...
                obj.pDopplerSampleRate),classtouse);
            
            if ~obj.pUseMatchedFilter  %dechirp
                dopplerGrid = -dopplerGridtemp; % due to order of dechirp mixing
            else
                dopplerGrid = dopplerGridtemp;
            end
            
            if obj.pOutputSpeed
                obj.pDopplerGrid = dop2speed(dopplerGrid,...
                    obj.PropagationSpeed/obj.OperatingFrequency)/2;
            else
                obj.pDopplerGrid = dopplerGrid;
            end
        end
        
        function processInputSizeChangeImpl(obj,x,~,~)
            classtouse = class(x);
            if obj.pUseMatchedFilter && strcmp(obj.PRFSource,'Auto')
                num_rng_samples = size(x,1);
                obj.pDopplerSampleRate = obj.pSampleRate/num_rng_samples;
                dopplerGrid = cast(obj.fftshiftfreqgrid(obj.pDopplerFFTLength,...
                    obj.pDopplerSampleRate),classtouse);

                if obj.pOutputSpeed
                    obj.pDopplerGrid = dop2speed(dopplerGrid,...
                        obj.PropagationSpeed/obj.OperatingFrequency)/2;
                else
                    obj.pDopplerGrid = dopplerGrid;
                end
            end
        end

        function [rdresp_out,rng_grid,dop_grid] = stepImpl(obj,x,varargin)
            
            if ~obj.pSizeInitialized
                processInputSizeChangeImpl(obj,x);
                obj.pSizeInitialized = true;
            end
            
            % Range processing
            if obj.pUseMatchedFilter || ...
                    (strcmp(obj.RangeMethod,'FFT') && obj.DechirpInput)
                xref = varargin{1};
                [rresp,rng_grid] = step(obj.cRangeResponse,x(:,1:obj.pValidatedNumInputChannels,1:obj.pValidatedNumInputPages),xref);
                
                if strcmp(obj.PRFSource,'Input port')
                    prf = varargin{2};
                else
                    prf = obj.pDopplerSampleRate;
                end
            else
                [rresp,rng_grid] = step(obj.cRangeResponse,x(:,1:obj.pValidatedNumInputChannels,1:obj.pValidatedNumInputPages));
                
                if strcmp(obj.PRFSource,'Input port')
                    prf = varargin{1};
                else
                    prf = obj.pDopplerSampleRate;
                end
            end
            
            % Doppler processing
            if ~strcmp(obj.PRFSource,'Auto')
                checkPRF(obj,x,prf);
            end
            
            isCube = obj.pValidatedNumInputPages>1;
            if isCube
                rresp_rpc = permute(rresp,[1 3 2]); % (rng,ch,pulse) -> (rng,pulse,ch)
            else
                rresp_rpc = rresp; % (rng,pulse) -> (rng,pulse)
            end
            
            szResp = size(rresp_rpc);
            if obj.pUseMatchedFilter
                d_fft_input = bsxfun(@times,obj.pDopplerWinCoeff.',rresp_rpc);
                rdresp = fftshift(fft(complex(d_fft_input),obj.pDopplerFFTLength,2),2);
            else
                if isCube
                    rdresp = zeros(szResp(1),obj.pDopplerFFTLength,obj.pValidatedNumInputChannels,'like',rresp_rpc);
                else
                    rdresp = zeros(szResp(1),obj.pDopplerFFTLength,'like',rresp_rpc);
                end
                
                d_fft_input = shiftdim(rresp_rpc,1);
                d_fft_input = bsxfun(@times,obj.pDopplerWinCoeff,d_fft_input);
                rdtemp = fftshift(step(obj.cDopplerFFT,complex(d_fft_input)),1);
                if isCube
                    rdtemp = shiftdim(rdtemp,2);
                    rdtemp = reshape(rdtemp,szResp(1),obj.pDopplerFFTLength,obj.pValidatedNumInputChannels);
                else
                    rdtemp = shiftdim(rdtemp,1);
                    rdtemp = reshape(rdtemp,szResp(1),[]);
                end
                rdresp(:) = rdtemp;
            end
            
            if isCube
                rdresp_out = permute(rdresp,[1 3 2]); % (rng,dop,ch) -> (rng,ch,dop)
            else
                rdresp_out = rdresp; % (rng,dop) -> (rng,dop)
            end
            
            dop_grid = obj.pDopplerGrid;
            if strcmp(obj.PRFSource,'Input port')
                dop_grid = dop_grid*prf;
            end
        end
        
        function num = getNumInputsImpl(obj)
            num = 3;
            if strcmp(obj.RangeMethod,'FFT') && ...
                ~obj.DechirpInput
                num = num-1;
            end
            if ~strcmp(obj.PRFSource,'Input port')
                num = num-1;
            end
        end
        
        function num = getNumOutputsImpl(obj) %#ok<MANU>
            num = 3;
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
                        strcmp(prop,'RangeGridInterval') || ...
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
                           
            if strcmp(prop,'OperatingFrequency') && ...
                    strcmp(obj.DopplerOutput,'Frequency')
                flag = true;
            end
            
            if strcmp(prop,'CustomDopplerWindow') && ...
                    ~strcmp(obj.DopplerWindow,'Custom')
                flag = true;
            end
            if strcmp(prop,'DopplerSidelobeAttenuation') && ...
                    ~(strcmp(obj.DopplerWindow,'Chebyshev') || ...
                    strcmp(obj.DopplerWindow,'Taylor') || ...
                    strcmp(obj.DopplerWindow,'Kaiser'))
                flag = true;
            end
            if strcmp(prop,'DopplerFFTLength') && ...
                    strcmp(obj.DopplerFFTLengthSource,'Auto')
                flag = true;
            end
            if strcmp(prop,'PRF') && ...
                    ismember(obj.PRFSource,{'Auto' 'Input port'})
                flag = true;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            if isLocked(obj)
                s.pSizeInitialized = obj.pSizeInitialized;
                s.pNumPages = obj.pNumPages;
                s.pValidatedNumInputPages = obj.pValidatedNumInputPages;
                s.pUseMatchedFilter = obj.pUseMatchedFilter;
                s.cRangeResponse = saveobj(obj.cRangeResponse);
                s.pOutputSpeed = obj.pOutputSpeed;
                s.pDopplerWinCoeff = obj.pDopplerWinCoeff;
                s.pDopplerGrid = obj.pDopplerGrid;
                s.pSampleRate = obj.pSampleRate;
                s.pDopplerSampleRate = obj.pDopplerSampleRate;
                s.pDopplerFFTLength = obj.pDopplerFFTLength;
                s.cDopplerFFT = saveobj(obj.cDopplerFFT);
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                if isfield(s,'cRangeResponse')
                    obj.cRangeResponse = phased.RangeResponse.loadobj(s.cRangeResponse);
                    s = rmfield(s,'cRangeResponse');
                else
                    % Range processing moved into RangeResponse in R2017a
                    % Use saved configuration to construct RangeResponse
                    % object
                    
                    % Remove properties on RangeDopplerResponse not
                    % available on RangeResponse object
                    fn = {'PRFSource','PRF','DopplerFFTLengthSource',...
                        'DopplerFFTLength','DopplerWindow','DopplerSidelobeAttenuation',...
                        'CustomDopplerWindow','DopplerOutput','OperatingFrequency',...
                        'pOutputSpeed','pDopplerWinCoeff','pDopplerGrid','pDopplerSampleRate',...
                        'pDopplerFFTLength','cDopplerFFT'};
                    fn = fn(isfield(s,fn));
                    srngresp = rmfield(s,fn);
                    obj.cRangeResponse = phased.RangeResponse.loadobj(srngresp);

                    % Properties removed from RangeDopplerResponse after
                    % refactor to use RangeResponse
                    fn = {'pRangeWinCoeff','pRangeGrid','pRangeOffset',...
                        'pRangeFFTLength','pDecimationFilterCoefficients',...
                        'cMatchedFilter','cRangeFFT'};
                    fn = fn(isfield(s,fn));
                    s = rmfield(s,fn);
                end
                
                if isfield(s,'cDopplerFFT')
                    obj.cDopplerFFT = dsp.FFT.loadobj(s.cDopplerFFT);
                    s = rmfield(s,'cDopplerFFT');
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
            rangegroups = matlab.system.display.Section(...
                'Title','Range Settings',...
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
            for m = 1:numel(rangegroups.PropertyList)
                if strcmp(rangegroups.PropertyList{m},'DecimationFactor')
                    rangegroups.PropertyList{m} = dDecimationFactor;
                elseif strcmp(rangegroups.PropertyList{m},'RangeWindow')
                    rangegroups.PropertyList{m} = dRangeWindow;
                elseif strcmp(rangegroups.PropertyList{m},'MaximumNumInputSamplesSource')
                    rangegroups.PropertyList{m} = dMaximumNumInputSamplesSource;
                elseif strcmp(rangegroups.PropertyList{m},'MaximumNumInputSamples')
                    rangegroups.PropertyList{m} = dMaximumNumInputSamples;
                end
            end
            
            dopplergroups = matlab.system.display.Section(...
                'Title','Doppler Settings',...
                'PropertyList',{'PRFSource','PRF','DopplerFFTLengthSource',...
                'DopplerFFTLength','DopplerWindow',...
                'DopplerSidelobeAttenuation','CustomDopplerWindow',...
                'DopplerOutput','OperatingFrequency'});
            dDopplerWindow = matlab.system.display.internal.Property(...
                'DopplerWindow', 'StringSetValues', {'None','Hamming',...
                'Chebyshev','Hann','Kaiser','Taylor'});
            for m = 1:numel(dopplergroups.PropertyList)
                if strcmp(dopplergroups.PropertyList{m},'DopplerWindow')
                    dopplergroups.PropertyList{m} = dDopplerWindow;
                end
            end
            
            groups = [rangegroups,dopplergroups];
        end
        
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:RangeDopplerResponseTitle')),...
                'Text',getString(message('phased:library:block:RangeDopplerResponseDesc')));
        end
    end

    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)
            if strcmp(obj.RangeMethod,'Matched filter')
                varargout = {'X','Coeff','PRF'};
            else
                if obj.DechirpInput
                    varargout = {'X','XRef','PRF'};
                else
                    varargout = {'X','PRF'};
                end
            end
            if ~strcmp(obj.PRFSource,'Input port')
                varargout = varargout(~ismember(varargout,'PRF'));
            end
        end
        
        function varargout = getOutputNamesImpl(~)
                    varargout = {'Resp','Range','Dop'};
        end
               
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Range-Doppler\nResponse');
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
            if obj.DopplerFFTLengthSource(1) == 'A' %Auto
                Q = szX(end);
            else
                Q = obj.DopplerFFTLength;
            end
            if numel(szX)==3
                respSz = [P szX(2) Q];
            else
                respSz = [P Q];
            end
            varargout = {respSz,[P 1],[Q 1]};
        end
        function varargout = isOutputFixedSizeImpl(obj) 
            varargout = {propagatedInputFixedSize(obj, 1), propagatedInputFixedSize(obj, 1), true};
        end
        function varargout = getOutputDataTypeImpl(obj)  %#ok<MANU>
            varargout = {propagatedInputDataType(obj,1),propagatedInputDataType(obj,1),propagatedInputDataType(obj,1)};
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout = {true,false,false};
        end            
    end
    
    methods(Access=private)
        function checkPRF(obj,x,specifiedPRF)
            
            if obj.pUseMatchedFilter
                numRngSmps = getPropagatedNumInputSamples(obj,x);
            else
                sz_x = size(x);
                if ~obj.DechirpInput
                    numRngSmps = sz_x(1);
                else
                    numRngSmps = ceil(sz_x(1)/obj.DecimationFactor);
                end
            end
            
            % Maximum PRF that makes sense based on size of input data
            % - The length of the input data (number of range samples),
            %   sets the minimum PRI that makes sense for the input signal
            %   (the signal could be longer, if for example, the input data
            %   has been range gated). Setting the minimum PRI is
            %   equivalent to setting the maximum PRF.
            maxPRF = obj.pSampleRate/numRngSmps;
            
            cond = specifiedPRF>maxPRF;
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:phased:prfMustBeLessThanSampleRateOverNumSamples');
            end
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

% [EOF]
