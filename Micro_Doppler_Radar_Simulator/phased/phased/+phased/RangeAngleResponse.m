classdef (Sealed,StrictDefaults) RangeAngleResponse < phased.internal.AbstractNarrowbandArrayProcessing & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%RangeAngleResponse   Range-angle response
%   H = phased.RangeAngleResponse creates a range-angle response
%   System object, H. This object calculates the range-angle response of
%   the input data.
%
%   H = phased.RangeAngleResponse(Name,Value) creates a range-angle
%   response object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The RangeAngleResponse object generates the response by first
%   processing the input signal in the range domain using either a matched
%   filter or dechirp operation, and then processing along azimuth angles.
%
%   Step method syntax:
%   
%   [RESP,RANGE,ANG] = step(H,X) calculates the range-angle response of
%   the input data X. This syntax applies when you set the RangeMethod
%   property to 'FFT' and then the DechirpInput property to false. This
%   syntax is most commonly used with FMCW signals.
%
%   X is a dechirped input signal. X must be a KxN matrix or a KxNxL array
%   where K denotes the number of fast time samples, L is the number of
%   dechirped frequency sweeps, and N is the number of channels (antenna
%   elements or beams). 
%
%   RESP is either an MxP matrix or an MxPxL array containing the complex
%   range-angle response of the input X. The number of dimensions in RESP
%   will match the number of dimensions in X. When X and RESP are
%   3-dimensional arrays, L is the number of frequency sweeps. M is
%   determined by either the number of rows in X when you set the
%   RangeFFTLengthSource property to 'Auto', or the value specified in the
%   RangeFFTLength property when you set the RangeFFTLengthSource property
%   to 'Property'. P is determined by the NumAngleSamples property.
%
%   RANGE is a length M column vector containing the range samples, in
%   meters, at which the range-angle response is evaluated. ANG is a length
%   P column vector containing angles, in degrees, at which the range-angle
%   response is evaluated.
%
%   [RESP,RANGE,ANG] = step(H,X,XREF) uses input XREF as the reference
%   signal to dechirp the input signal X. This syntax applies when you set
%   the RangeMethod property to 'FFT' and then the DechirpInput property to
%   true. This syntax is most commonly used with FMCW signals and the
%   reference signal is, in general, the transmitted signal.
%
%   X is an input signal to be dechirped by the RangeAngleResponse
%   object. X must be a KxN matrix or an KxNxL array where K denotes the
%   number of fast time samples, L is the number of frequency sweeps, and N
%   is the number of channels (antenna elements or beams). XREF must be a
%   column vector whose number of rows is the same as the number of rows of
%   X. XREF is the reference signal used to dechirp X.
%
%   In this syntax, the number of rows in RESP is the quotient of the
%   number of rows in X and the value specified in the DecimationFactor
%   property.
%
%   [RESP,RANGE,ANG] = step(H,X,COEFF) uses COEFF as the matched filter
%   coefficients. This method applies when you set the RangeMethod property
%   to 'Matched filter'. This syntax is most commonly used with pulsed
%   signals.
%
%   X is an input signal to be match filtered by the RangeAngleResponse
%   object. X must be a KxN matrix or a KxNxL array where K denotes the
%   number of fast time samples, L is the number of pulses, and N is the
%   number of channels (antenna elements or beams). COEFF must be a column
%   vector containing the matched filter coefficients.
%
%   In this syntax, RESP is either a KxP matrix or a KxPxL array, where the
%   number of rows of RESP equals the number of rows in X. When X is a
%   3-dimensional array, L is the number of pulses. P is determined by
%   the NumAngleSamples property.
%
%   RANGE is a length K column vector containing the range samples at which
%   the range response is evaluated. ANG is a length P column vector
%   containing the angles at which the range-angle response is evaluated.
%   
%   [RESP,RANGE,ANG] = step(...,EL) uses input EL as the elevation angle
%   (in degrees) to calculate the range-angle response when you set the
%   ElevationAngleSource property to 'Input port'.
%
%   You can combine optional input and output arguments when their enabling
%   properties are set. Optional inputs and outputs must be listed in the
%   same order as the order of the enabling properties. For example,
% 
%   [RESP,RANGE,ANG] = step(H,X,XREF,EL)
%
%   or
%   
%   [RESP,RANGE,ANG] = step(H,X,COEFF,EL)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   RangeAngleResponse methods:
%
%   step         - Calculate range-angle response
%   release      - Allow property value and input characteristics changes
%   clone        - Create range-angle response object with same property
%                  values
%   isLocked     - Locked status (logical)
%   plotResponse - Plot range-angle response
%
%   RangeAngleResponse properties:
%       
%   SensorArray                  - Sensor array
%   RangeMethod                  - Range processing method
%   PropagationSpeed             - Propagation speed
%   OperatingFrequency           - Operating frequency
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
%   ElevationAngleSource         - Source of elevation angle
%   ElevationAngle               - Elevation angle
%   AngleSpan                    - Angle span
%   NumAngleSamples              - Number of angles 
%
%   % Example:
%   %   Calculate the range-angle response from a pulsed radar
%   %   transmitting a rectangular waveform using the matched filter
%   %   approach. The signal includes three target returns. Two are
%   %   approximately 2000 m away and the third is approximately 3500 m
%   %   away. In addition, two targets are stationary relative to the radar
%   %   while the third is moving away from the radar at approximately
%   %   100 m/s.
%   
%   % Load example data
%   load('RangeAngleResponseExampleData','rectdata');
%   fs = rectdata.fs;
%   propspeed = rectdata.propspeed;
%   fc = rectdata.fc;
%   rxdata = rectdata.rxdata;
%   mfcoeffs = rectdata.mfcoeffs;
%   noisepower = rectdata.noisepower;
%   antennaarray = rectdata.antennaarray;
%   
%   % Create range-angle response for matched filter processing.
%   % Interpolate to 1024 Doppler bins.
%   rngangresp = phased.RangeAngleResponse(...
%       'SensorArray',antennaarray,'OperatingFrequency',fc,...
%       'SampleRate',fs,'PropagationSpeed',propspeed);
%   
%   % Perform range-angle processing
%   [resp,rng_grid,ang_grid] = rngangresp(rxdata,mfcoeffs);
%   
%   % Calculate noise power after range-angle processing
%   mfgain = mfcoeffs'*mfcoeffs;
%   noisepower_proc = mfgain*noisepower;
%   
%   imagesc(ang_grid,rng_grid,pow2db(abs(resp).^2/noisepower_proc));
%   xlabel('Angle (degrees)'); ylabel('Range (m)'); 
%   title('Range Angle Map');
%   set(get(colorbar,'YLabel'),'String','SNR (dB)'); caxis([0 60]);
%
%   See also phased, phased.RangeResponse, phased.AngleDopplerResponse,
%   phased.MatchedFilter, dechirp.

%   Copyright 2018 The MathWorks, Inc.

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
        %ElevationAngleSource   Source of elevation angle
        %   Specify how to determine the elevation angle when calculating
        %   the angle-Doppler response as one of 'Property' | 'Input port',
        %   where the default is 'Property'. When you set this property to
        %   'Property', the elevation angle is determined by the
        %   ElevationAngle property. When you set this property to 'Input
        %   port', the elevation angle is determined by the input argument.
        ElevationAngleSource = 'Property';
        %ElevationAngle     Elevation angle (deg)
        %   Specify the elevation angle (in degrees) used to calculate the
        %   range-angle response as a scalar. The angle must be between
        %   -90 and 90. This property applies when you set the
        %   ElevationAngleSource property to 'Property'. The default value
        %   of this property is 0.
        ElevationAngle = 0;
        %AngleSpan      Angle span (deg)
        %   Specify the angle span (in degrees) used to calculate the
        %   range-angle response as a 2-element row vector in the form of
        %   [min_angle max_angle]. The default value of this property is
        %   [-90 90].
        AngleSpan = [-90 90]
        %NumAngleSamples    Number of angle bins
        %   Specify the number of samples in angular domain used to
        %   calculate the angle-Doppler response as a positive integer.
        %   This value must be greater than 2. The default value of this
        %   property is 256.
        NumAngleSamples = 256;
    end
    
    properties (Constant, Hidden)
        RangeMethodSet = matlab.system.internal.StringSetGF(...
            {'Matched filter','FFT'},{'Dechirp'},{'FFT'});
        RangeWindowSet = matlab.system.StringSet({'None','Hamming',...
            'Chebyshev','Hann','Kaiser','Taylor','Custom'});
        RangeFFTLengthSourceSet = dsp.CommonSets.getSet(...
            'AutoOrProperty');
        MaximumNumInputSamplesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        ElevationAngleSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
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
        cRangeResponse
        cSteeringVector
    end
    
    properties (Access = private, Nontunable)
        pSampleRate
        pValidatedNumInputPages
    end
    
    properties (Access = private)
        pAngleGrid
        pNumPages = -1
end
    
    properties (Access = private, Logical)
        pSizeInitialized = false
    end
    
    methods
        function set.SampleRate(obj, value)
            validateattributes(value,{'double'}, {'scalar',...
                'positive','finite'},...
                '','SampleRate');
            obj.SampleRate = value;
        end
        function set.RangeSidelobeAttenuation(obj,value)
            validateattributes( value, { 'double' },...
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
        function set.SweepSlope(obj, value)
            validateattributes(value,{'double'},...
                {'scalar','real','finite'},...
                '','SweepSlope');
            obj.SweepSlope = value;
        end
        function set.ElevationAngle(obj,val)
            sigdatatypes.validateAngle(val,'phased.AngleDopplerResponse',...
                'ElevationAngle',{'scalar','>=',-90,'<=',90});
            obj.ElevationAngle = val;
        end
        function set.AngleSpan(obj,val)
            sigdatatypes.validateAngle(val,'phased.AngleDopplerResponse',...
                'AngleSpan',{'size',[1 2],'increasing','>=',-90,'<=',90});
            obj.AngleSpan = val;
        end
        function set.NumAngleSamples(obj,val)
            sigdatatypes.validateIndex(val,'phased.AngleDopplerResponse',...
                'NumAngleSamples',{'scalar','>=',2});
            obj.NumAngleSamples = val;
        end
        function set.ReferenceRange(obj,val)
            sigdatatypes.validateDistance(val,'','ReferenceRange',...
                {'scalar'});
            obj.ReferenceRange = val;
        end
        
    end

    methods

        function obj = RangeAngleResponse(varargin)
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
        end
        
        function varargout = plotResponse(obj,x,varargin)
        %plotResponse   Plot range-angle response
        %   plotResponse(H,X) plots the range-angle response of the input
        %   signal, X, in dB scale. This syntax applies when you set the
        %   RangeMethod property to 'FFT' and then the DechirpInput
        %   property to false.
        %
        %   X must be a matrix where each column contains the dechirped
        %   signal from one frequency sweep.
        %
        %   plotResponse(H,X,XREF) plots the range-angle response after
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
        %   plotResponse(H,X,COEFF) plots the range-angle response after
        %   performing a matched filter operation on X using the matched
        %   filter coefficients specified in COEFF. This syntax applies
        %   when you set the RangeMethod property to 'Matched filter'.
        %
        %   X must be a matrix where each column contains the signal from
        %   one pulse. COEFF must be a column vector containing the matched
        %   filter coefficients.
        %
        %   plotResponse(...,Name,Value) plots the range-angle response
        %   with the specified parameter Name set to the specified value.
        %   The parameter Names are
        %                   Unit: The unit of the plot, using one of 
        %                         | 'db' | 'mag' | 'pow' |. The default 
        %                         value is 'db'.
        %       CoordinateSystem: Set the style of the plot as one of
        %                         'rectangular' | 'polar', where the
        %                         default is 'rectangular'.
        %
        %   % Example:
        %   %   Plot the range-angle response of an FMCW signal. The
        %   %   signal is not dechirped. The signal contains the return
        %   %   from one target which is approximately 2200 m away and has
        %   %   a normalized Doppler frequency of about -0.36 relative to
        %   %   the radar.
        %   
        %   % Load example data
        %   load('RangeAngleResponseExampleData','fmcwdata');
        %   fs = fmcwdata.fs;
        %   propspeed = fmcwdata.propspeed;
        %   rxdata = fmcwdata.rxdata;
        %   refsig = fmcwdata.refsig;
        %   sweepslope = fmcwdata.sweepslope;
        %   antennaarray = fmcwdata.antennaarray;
        %   fc = fmcwdata.fc;
        %   
        %   % Create range-angle response for dechirp and FFT processing
        %   rngangresp = phased.RangeAngleResponse(...
        %       'SensorArray',antennaarray,'RangeMethod','FFT',...
        %       'SweepSlope',sweepslope,...
        %       'DechirpInput',true,'ReferenceRangeCentered',false,...
        %       'SampleRate',fs,'OperatingFrequency',fc,...
        %       'PropagationSpeed',propspeed,'AngleSpan',[-45 45]);
        %   
        %   % Plot range-angle response of processed data
        %   plotResponse(rngangresp,rxdata,refsig,...
        %       'Unit','db','CoordinateSystem','polar');
        %
        %   See also phased.RangeAngleResponse,
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
            
            haresp = clone(obj);
            release(haresp);
            
            if (haresp.RangeMethod(1) == 'F') && ... %Dechirp
                    ~haresp.DechirpInput
                if (haresp.ElevationAngleSource(1) == 'I') && ... % Input port
                        isempty(varargin)
                    error(message(...
                        'phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                        'EL','ElevationAngleSource','Input port'));
                end
                if (haresp.ElevationAngleSource(1) == 'I')
                    elang = varargin{1};
                    varargin(1) = [];
                    [resp,rng_grid,ang_grid] = step(haresp,x,elang);
                else
                    [resp,rng_grid,ang_grid] = step(haresp,x);
                end
            else
                if isempty(varargin)
                    if (haresp.RangeMethod(1) == 'M') %Matched filter
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
                if (haresp.ElevationAngleSource(1) == 'I') && ... % Input port
                        isempty(varargin)
                    error(message(...
                        'phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                        'EL','ElevationAngleSource','Input port'));
                end
                if (haresp.ElevationAngleSource(1) == 'I')
                    elang = varargin{1};
                    varargin(1) = [];
                    [resp,rng_grid,ang_grid] = step(haresp,x,xref,elang);
                else
                    [resp,rng_grid,ang_grid] = step(haresp,x,xref);
                end
            end
            
            unit = 'db';
            coordinatesystem = 'rectangular';
            sigutils.pvparse(varargin{:});
            unit = validatestring(unit,{'db','mag','pow'},...
                'plotResponse','Unit');
            coordinatesystem = validatestring(coordinatesystem,{'rectangular','polar'},...
                'plotResponse','CoordinateSystem');
            
            hresp = phased.internal.RangeAnglePattern(...
                'Angle',ang_grid(:)',...
                'Range',rng_grid(:)',...
                'Pattern',abs(resp));
            if nargout
                varargout{1} = ...
                    plot(hresp,'Units',unit,'Style',coordinatesystem);
            else
                plot(hresp,'Units',unit,'Style',coordinatesystem);
            end
             
        end

    end
    
    methods (Access = protected)
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            reset(obj.cRangeResponse);
            reset(obj.cSteeringVector);
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            obj.pSizeInitialized = false;
            obj.pNumPages = -1;
            
            release(obj.cRangeResponse);
            release(obj.cSteeringVector);
        end
                
        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object
            obj.cRangeResponse.ReferenceRange = obj.ReferenceRange;
        end
                
        function validateInputsImpl(obj,x,opt1,opt2)
            cond = ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                  'MATLAB:system:invalidInputDataType','X','double');
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
                cond = ~isa(coeff,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDataType','Coeff','double');
                end
                cond = ~iscolumn(coeff) || isempty(coeff);
                if cond
                    coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeColVector','Coeff');
                end
                if strcmp(obj.ElevationAngleSource,'Input port')
                    elang = opt2;
                    cond = ~isa(elang,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:invalidInputDataType','El','double');
                    end
                    cond = ~isscalar(elang);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:inputMustBeScalar','El');
                    end
                    cond = ~isreal(elang);
                    if cond
                        coder.internal.errorIf(cond,'phased:AngleDopplerResponse:NeedReal','El');
                    end
                end
            else
                if obj.DechirpInput
                    xref = opt1;
                    cond = ~isa(xref,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:invalidInputDataType','XRef','double');
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
                    if strcmp(obj.ElevationAngleSource,'Input port')
                        elang = opt2;
                        cond = ~isa(elang,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:invalidInputDataType','El','double');
                        end
                        cond = ~isscalar(elang);
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:inputMustBeScalar','El');
                        end
                        cond = ~isreal(elang);
                        if cond
                            coder.internal.errorIf(cond,'phased:AngleDopplerResponse:NeedReal','El');
                        end
                    end
                else
                    if strcmp(obj.ElevationAngleSource,'Input port')
                        elang = opt1;
                        cond = ~isa(elang,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:invalidInputDataType','El','double');
                        end
                        cond = ~isscalar(elang);
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:inputMustBeScalar','El');
                        end
                        cond = ~isreal(elang);
                        if cond
                            coder.internal.errorIf(cond,'phased:AngleDopplerResponse:NeedReal','El');
                        end
                    end
                end
            end
            
            validateNumChannels(obj,x);
            if ndims(x) == 3 && obj.pNumPages ~= -1
                validateNumPages(obj,x,obj.pNumPages);
            end             
        end
        
        function setupImpl(obj,x,~,~)
            coder.extrinsic('phased.RangeAngleResponse.getWinCoeff');
            
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputPages = size(x,3);
            
            obj.pNumInputChannels = obj.pValidatedNumInputChannels;
            obj.pNumPages = obj.pValidatedNumInputPages;
            
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
            
            obj.pUseMatchedFilter = (obj.RangeMethod(1) == 'M'); %Matched filter
            
            if obj.pUseMatchedFilter
                obj.pSampleRate = fs;
            else
                if ~obj.DechirpInput
                    obj.pSampleRate = fs;
                else
                    obj.pSampleRate = fs/obj.DecimationFactor;
                end
            end
            
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed);
            
            obj.pAngleGrid = linspace(obj.AngleSpan(1),obj.AngleSpan(2),obj.NumAngleSamples).';
            
        end
        
        function [raresp,rng_grid,ang_grid] = stepImpl(obj,x,varargin)
            
            if ~obj.pSizeInitialized
                processInputSizeChangeImpl(obj,x);
                obj.pSizeInitialized = true;
            end
            
            % Range processing
            if obj.pUseMatchedFilter || ...
                    (strcmp(obj.RangeMethod,'FFT') && obj.DechirpInput)
                xref = varargin{1};
                [rresp,rng_grid] = step(obj.cRangeResponse,x(:,1:obj.pValidatedNumInputChannels,1:obj.pValidatedNumInputPages),xref);
                if strcmp(obj.ElevationAngleSource,'Property')
                    elang = obj.ElevationAngle;
                else
                    elang = varargin{2};
                end
            else
                [rresp,rng_grid] = step(obj.cRangeResponse,x(:,1:obj.pValidatedNumInputChannels,1:obj.pValidatedNumInputPages));
                if strcmp(obj.ElevationAngleSource,'Property')
                    elang = obj.ElevationAngle;
                else
                    elang = varargin{1};
                end
            end
            cond = (elang > 90) || (elang < -90);
            if cond
                coder.internal.errorIf(cond,'phased:AngleDopplerResponse:Step:OutOfBoundElAng','El');
            end
            
            % Angular processing
            isCube = obj.pValidatedNumInputPages>1;           
            
            ang_grid = obj.pAngleGrid;
            ang_stv = step(obj.cSteeringVector,obj.OperatingFrequency,...
                [ang_grid elang*ones(size(ang_grid))].');
            
            szResp = size(rresp);

            if isCube
                raresp = zeros(szResp(1),obj.NumAngleSamples,obj.pValidatedNumInputPages,'like',rresp);
                for m = 1:obj.pValidatedNumInputPages
                    raresp(:,:,m) = rresp(:,:,m)*conj(ang_stv);
                end
            else
                raresp = rresp*conj(ang_stv);
            end
                        
        end
        
        function num = getNumInputsImpl(obj)
            num = 3;
            if strcmp(obj.RangeMethod,'FFT') && ...
                ~obj.DechirpInput
                num = num-1;
            end
            if strcmp(obj.ElevationAngleSource,'Property')
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
            if strcmp(obj.ElevationAngleSource,'Input port') && ...
                    strcmp(prop, 'ElevationAngle')
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
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.pAngleGrid = obj.pAngleGrid;
                s.pSampleRate = obj.pSampleRate;
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                obj.cRangeResponse = phased.RangeResponse.loadobj(s.cRangeResponse);
                s = rmfield(s,'cRangeResponse');
                obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                s = rmfield(s,'cSteeringVector');
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
            groups = getPropertyGroupsImpl@phased.internal.AbstractNarrowbandArrayProcessing('subarray');
            rangegroups = matlab.system.display.Section(...
                'Title','Range Settings',...
                'PropertyList',{'RangeMethod',...
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
            angleprops = {...
                'ElevationAngleSource',...
                'ElevationAngle',...
                'AngleSpan',...
                'NumAngleSamples'};
                       
            groups(1).PropertyList = [groups(1).PropertyList,rangegroups.PropertyList,angleprops];
        end
        
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:RangeAngleResponseTitle')),...
                'Text',getString(message('phased:library:block:RangeAngleResponseDesc')));
        end
    end

    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)
            if strcmp(obj.RangeMethod,'Matched filter')
                varargout = {'X','Coeff','El'};
            else
                if obj.DechirpInput
                    varargout = {'X','XRef','El'};
                else
                    varargout = {'X','El'};
                end
            end
            if ~strcmp(obj.ElevationAngleSource,'Input port')
                varargout = varargout(~ismember(varargout,'El'));
            end
        end
        
        function varargout = getOutputNamesImpl(~)
                    varargout = {'Resp','Range','Ang'};
        end
               
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Range-Angle\nResponse');
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
            Q = obj.NumAngleSamples;
            if numel(szX)==3
                respSz = [P Q szX(3)];
            else
                respSz = [P Q];
            end
            varargout = {respSz,[P 1],[Q 1]};
        end
        function varargout = isOutputFixedSizeImpl(obj) 
            varargout = {propagatedInputFixedSize(obj, 1), propagatedInputFixedSize(obj, 1), true};
        end
        function varargout = getOutputDataTypeImpl(obj)  %#ok<MANU>
            varargout = {'double','double','double'};
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout = {true,false,false};
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
        
    end
end

% [EOF]
