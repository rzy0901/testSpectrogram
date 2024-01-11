classdef (Sealed, StrictDefaults) RangeDopplerScope < matlab.System   
% RangeDopplerScope Visualize range-Doppler patterns
%   SCOPE = phased.RangeDopplerScope returns a System object, SCOPE, that
%   can display range-Doppler response patterns.
%
%   SCOPE = phased.RangeDopplerScope ('Name', Value,) returns a
%   range-Doppler scope System object, SCOPE, with each specified property
%   name set to the specified value. You can specify name-value pair
%   arguments in any order as (Name 1, Value 1, ..., NameN, ValueN).
%
%   Step method syntax when IQDataInput property is false:
%
%   step(SCOPE,X,RANGE,DOP) displays the range-Doppler response, X. X is
%   KxL matrix  containing the range-Doppler response where K denotes the
%   number of fast time samples, L is the number of pulses in case of
%   pulsed signals or K denotes the number of fast time samples, L is the
%   number of dechirped frequency sweeps in case of FMCW signals.
%
%   RANGE is a length-K column vector containing the range values at which
%   the range-Doppler response is evaluated. DOP is a length-L column
%   vector containing either Doppler or speed values at which the
%   range-Doppler response is evaluated.
%
%   Step method syntax when IQDataInput property is true:
%
%   step(SCOPE,X) calculates and displays the range-Doppler response of X.
%   This syntax applies when you set the RangeMethod property to 'FFT' and
%   the DechirpInput property to false. This syntax is most commonly used
%   with FMCW signals.
%
%   X is a dechirped input signal. X must be a KxL matrix where K denotes
%   the number of fast time samples and L is the number of dechirped
%   frequency sweeps.
%
%   By default, all sweeps in X are assumed to be consecutive. If the
%   sweeps are not consecutive, the PRF of the input data must be provided
%   by setting the PRFSource property to 'Property'.
%
%   step(SCOPE,X,XREF) uses input XREF as the reference signal to dechirp
%   the input signal X to calculate and display the range-Doppler response
%   of X. This syntax applies when you set the RangeMethod property to
%   'FFT' and the DechirpInput property to true. This syntax is most
%   commonly used with FMCW signals and the reference signal is, in
%   general, the transmitted signal.
%
%   X is an input signal to be dechirped by the RangeDopplerScope object. X
%   must be a KxL matrix where K denotes the number of fast time samples
%   and L is the number of frequency sweeps.
%
%   By default, all sweeps in X are assumed to be consecutive. If the
%   sweeps are not consecutive, the PRF of the input data must be provided
%   by setting the PRFSource property to 'Property'.
%
%   step(SCOPE,X,COEFF) uses COEFF as the matched filter coefficients, to
%   matched filer the input X and display the range-Doppler response of X.
%   This method applies when you set the RangeMethod property to 'Matched
%   filter'. This syntax is most commonly used with pulsed signals.
%
%   By default, all pulses in X are assumed to be consecutive. If the
%   pulses are not consecutive, the PRF of the input data must be provided
%   by setting the PRFSource property to 'Property'.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(SCOPE,x) and y = SCOPE(x) are
%   equivalent
%
%
%   RangeDopplerScope methods:
%
%   step         - Display range-Doppler response
%   release      - Allow property value and input characteristics changes
%   clone        - Create range-Doppler scope object with same property
%                  values
%   isLocked     - Locked status (logical)
%   show         - Turn on visibility of the scope
%   hide         - Turn off visibility of the scope
%   isVisible    - Return visibility of the scope (logical)
%
%   RangeDopplerScope properties:
%
%   Name                         - Scope window name
%   Position                     - Scope window position
%   IQDataInput                  - Type of input
%   ResponseUnits                - Output response unit
%   RangeLabel                   - Range axis label
%   DopplerLabel                 - Doppler axis label
%   RangeMethod                  - Range processing method
%   RangeUnits                   - Range unit
%   PropagationSpeed             - Propagation speed
%   SampleRate                   - Sample rate
%   SweepSlope                   - FM sweep slope
%   DechirpInput                 - Dechirp input signal
%   RangeFFTLength               - FFT length in range processing
%   ReferenceRangeCentered       - Set reference range at center
%   ReferenceRange               - Reference range
%   PRFSource                    - Source of pulse repetition frequency
%   PRF                          - Pulse repetition frequency used in
%                                  Doppler processing
%   DopplerFFTLength             - FFT length in Doppler processing
%   DopplerOutput                - Doppler output
%   OperatingFrequency           - Signal carrier frequency
%   NormalizeDoppler             - Normalize Doppler grid
%   SpeedUnits                   - Doppler speed unit
%   FrequencyUnits               - Doppler frequency unit
%
%   % Example:
%   %   Calculate and visualize the range-Doppler response from a pulsed 
%   %   radar transmitting a rectangular waveform using the matched filter
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
%   % Create range-Doppler scope for matched filter processing and 
%   % visualization. Interpolate to 1024 Doppler bins.
%  
%   rngdopScope = phased.RangeDopplerScope(...
%       'IQDataInput',true,'RangeMethod','Matched filter',...
%       'Name','Range-Doppler Scope',...
%       'Position',[560 375 560 420],'ResponseUnits','db',...
%       'RangeUnits','m','DopplerFFTLength',1024,...
%       'DopplerOutput','Speed',...
%       'OperatingFrequency',fc,...
%       'SampleRate',fs,...
%       'PropagationSpeed',propspeed);
%   
%    rngdopScope(rxdata,mfcoeffs);
%
%   See also phased, phased.AngleDopplerScope, phased.RangeAngleScope,
%   phased.RangeDopplerResponse.

%   Copyright 2018 The MathWorks, Inc.

    properties(Nontunable,Logical)
        %IQDataInput Specify the type of input
        %   Specify whether the input is I/Q (raw) data or processed data.
        %   If you set this property true, it implies that raw data is
        %   passed and processing is done along range and Doppler domain
        %   before visualizing. Set this property false if you have
        %   processed data. The default value is true.
        IQDataInput = true;
        
        %NormalizeDoppler Normalize Doppler
        %   Set the NormalizeDoppler property true if you want the Doppler
        %   frequency to be normalized. Set this to false to plot the
        %   range-Doppler response without normalizing the Doppler
        %   frequency. This parameter applies when 'IQDataInput' property
        %   is true and you set the DopplerOutput property to 'Frequency'.
        %   The default value is false.
        NormalizeDoppler = false;
    end

    properties (Nontunable)
        %ResponseUnits Unit of the response
        %   Specify the unit of the plot, using one of | 'db' | 'magnitude'
        %   | 'power'|. The default is 'db'.
        ResponseUnits = 'db';
        
        %RangeMethod    Range processing method
        %   Specify the method of range processing as one of 'Matched
        %   filter' | 'FFT', where the default is 'Matched filter'. This
        %   property applies when you set the 'IQDataInput' property true.
        %   When you set this property to 'Matched filter', the range
        %   processing is achieved by applying a matched filter to the
        %   incoming signal. When you set this property to 'FFT', the range
        %   processing is achieved by applying FFT to the input signal.
        %
        %   The matched filter approach is often used with pulsed signals,
        %   where the matched filter is the time reverse of the transmitted
        %   signal.  The FFT approach are often used with FMCW signals.
        RangeMethod = 'Matched filter';
        
        %PropagationSpeed   Propagation speed (m/s)
        %   Specify the propagation speed (in m/s) of the signal as a
        %   scalar.  This property applies when you set the 'IQDataInput'
        %   property true. The default value of this property is the speed
        %   of light.
        PropagationSpeed = physconst('lightspeed');
    end
       
    properties (Nontunable)
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a positive scalar.  This
        %   property applies when you set the 'IQDataInput' property true.
        %   The default value of this property is 1e6 (1 MHz).
        SampleRate = 1e6;
        
        %SweepSlope     FM sweep slope (Hz/s)
        %   Specify the slope of the linear FM sweeping (in Hz/s) as a
        %   scalar. This property applies when you set the RangeMethod
        %   property to 'FFT' and 'IQDataInput' proeprty true. The default
        %   value is 1e9.
        SweepSlope = 1e9;
    end
    
    properties (Nontunable,Logical)
        %DechirpInput     Dechirp input signal
        %   Set this property to true to dechirp the input signal first
        %   before range processing. Set this property to false to indicate
        %   that the input signal is already dechirped and no dechirp
        %   operation is necessary. This property applies when you set the
        %   RangeMethod property to 'FFT' and 'IQDataInput' property true.
        %   The default value of this property is false.
        DechirpInput = false;
    end
    
    properties (Nontunable, PositiveInteger)
        %RangeFFTLength     FFT length in range processing
        %   Specify the FFT length in range domain as a positive integer.
        %   This property applies when you set the RangeMethod property to
        %   'FFT' and 'IQDataInput' property true. The default value of
        %   this property is 1024.
        RangeFFTLength = 1024
    end
    
    properties (Nontunable, Logical)
        %ReferenceRangeCentered     Set reference range at center
        %   Set this property to true to set the reference range to the
        %   center of the range grid. Set this property to false to set the
        %   reference range to the beginning of the range grid. The default
        %   value of this property is true. This property only applies when
        %   you set the RangeMethod to 'FFT' and 'IQDataInput' property
        %   true.
        ReferenceRangeCentered = true
    end
    
    properties 
        %ReferenceRange     Reference range (m)
        %   Specify the reference range of the range grid as a nonnegative
        %   scalar.  This property applies when you set the 'IQDataInput'
        %   property true. The default value is 0. If you set the
        %   RangeMethod property to 'Matched filter', the reference range
        %   marks the start of the range grid. If you set the RangeMethod
        %   property to 'FFT', the position of the reference range is
        %   determined by the ReferenceRangeCentered property. If you set
        %   the ReferenceRangeCentered property to true, the reference
        %   range marks the center of the range grid. If you set the
        %   ReferenceRangeCentered property to false, the reference range
        %   marks the start of the range grid. This property is tunable.
        ReferenceRange = 0
    end
    
    properties (Nontunable)
        %PRFSource  Source of pulse repetition frequency
        %   Specify how to determine the pulse repetition frequency (PRF)
        %   of the processed input signal as one of 'Auto' | 'Property',
        %   where the default is 'Auto'. This property applies when you set
        %   the 'IQDataInput' property true. When you set this property to
        %   'Auto', the PRF is inferred from the number of rows in the
        %   input signal and the value of the SampleRate property. When you
        %   set this property to 'Property', the PRF is specified in the
        %   PRF property.
        PRFSource = 'Auto'
    end
    
    properties (Nontunable)
        %PRF  Pulse repetition frequency of the input signal (Hz)
        %   Specify the pulse repetition frequency (PRF) of the input
        %   signal (in Hz) as a positive scalar. The default value of this
        %   property is 10e3 (10 kHz). This property only applies when you
        %   set the PRFSource property to 'Property' and 'IQDataInput'
        %   property true.
        PRF = 10e3
    end
    
    properties (Nontunable, PositiveInteger)
        %DopplerFFTLength     FFT length in Doppler processing
        %   Specify the FFT length in Doppler processing as a positive
        %   integer. This property applies when you set the 'IQDataInput'
        %   property true. The default value of this property is 1024.
        DopplerFFTLength = 1024
    end
    
    properties (Nontunable)
        %DopplerOutput  Doppler output
        %   Specify the Doppler domain output as one of 'Frequency' |
        %   'Speed', where the default is 'Frequency'.  This property
        %   applies when you set the 'IQDataInput' property true. If you
        %   set this property to 'Frequency', the Doppler domain output of
        %   step, DOP, is the Doppler shift (in Hz). If you set this
        %   property to 'Speed', the Doppler domain output is the
        %   corresponding radial speed (in m/s).
        DopplerOutput = 'Frequency'
    end
    
    properties (Nontunable)
        %OperatingFrequency     Signal carrier frequency (Hz)
        %   Specify the signal carrier frequency (in Hz) as a scalar. This
        %   property applies when you set the DopplerOutput property to
        %   'Speed'. The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 3e8;
    end

    properties (Nontunable)
        %RangeUnits  Range unit
        %   Specify the range unit, using one of | 'm' | 'km' |
        %   'mi'|'nmi'|. This property applies when 'IQDataInput' property
        %   is true. The default is 'm'.
        RangeUnits = 'm';
        
        %SpeedUnits  Doppler speed unit
        %   Specify the Doppler speed unit, using one of | 'm/s' | 'km/h' |
        %   'mph'|'kt'|. This property applies only when 'IQDataInput'
        %   property is true and the DopplerOutput property is 'Speed'. The
        %   default is 'm/s'.
        SpeedUnits = 'm/s';
        
        %FrequencyUnits  Doppler frequency unit
        %   Specify the Doppler frequency unit, using one of | 'Hz' | 'kHz'
        %   | 'MHz'|.  This property applies only when 'IQDataInput'
        %   property is true and the DopplerOutput property is 'Frequency'
        %   and Normalized Doppler property is false. The default is 'Hz'.
        FrequencyUnits = 'Hz';
    end
    
    properties
        %RangeLabel Range-axis label
        %   Specify the range-axis label as a string.  This property
        %   applies when you set the 'IQDataInput' property false. The
        %   default value is 'Range (m)'. This property is tunable.
        RangeLabel = 'Range (m)';
        
        %DopplerLabel Doppler-axis label
        %   Specify the Doppler-axis label as a string. This property
        %   applies when you set the 'IQDataInput' property false. The
        %   default value is 'Doppler Frequency (Hz)'. This property is
        %   tunable.
        DopplerLabel = 'Doppler Frequency (Hz)';
        
        %Name Caption to display on scope window
        %   Specify the caption to display on the scope window as any
        %   string. The default value is 'Range-Doppler Scope'. This
        %   property is tunable.
        Name = 'Range-Doppler Scope';
        
        %Position Scope window position in pixels
        %   Specify the size and location of the scope window in pixels, as
        %   a four-element double vector of the form: [left bottom width
        %   height]. The default value of this property is dependent on the
        %   screen resolution, and is such that the window is positioned in
        %   the center of the screen, with a width and height of 800 and
        %   450 pixels respectively. This property is tunable.
        Position = [560 375 800 450];
    end
       
    properties(Access = private)
        cScope
        cResponse
        pRangeConversion
        pSpeedConversion
        pFrequencyConversion
    end
    
    properties (Constant, Hidden)
        ResponseUnitsSet = matlab.system.StringSet({'db','power',...
            'magnitude'});
        RangeUnitsSet = matlab.system.StringSet({'m','km','mi','nmi'});
        SpeedUnitsSet = matlab.system.StringSet({'m/s','km/h','mph','kt'});
        FrequencyUnitsSet = matlab.system.StringSet({'Hz','kHz','MHz'});
        RangeMethodSet = matlab.system.StringSet({'Matched filter','FFT'});
        DopplerOutputSet = matlab.system.StringSet({'Frequency','Speed'});
        PRFSourceSet = matlab.system.StringSet({'Auto','Property'});
    end
    
    methods
        % Constructor
        function obj = RangeDopplerScope(varargin)
            % Support name-value pair arguments when constructing the
            % RangeDopplerScope object
            setProperties(obj, nargin, varargin{:});
            
            % matlabshared.scopes.MatrixViewer - Visualization
            obj.cScope = matlabshared.scopes.MatrixViewer(...
                'XDataMode','Custom',...
                'YDataMode','Custom',...
                'AxisOrigin','Lower left corner',...
                'Title', getString(message('phased:scopes:rngdopLabel')),...
                'Colormap','parula'...
                );
        end
        
        function show(obj)
            %show    Show scope window
            %   SHOW(H) turns on the visibility of the scope window
            %   associated with the System object H.
            show(obj.cScope);
        end
        
        function hide(obj)
            %hide    Hide scope window
            %   HIDE(H) turns off the visibility of the scope window
            %   associated with the System object H.
            hide(obj.cScope);
        end

        function value = isVisible(obj)
            %isVisible    Visibility of the Scope window
            %   ISVISIBLE(H) Returns the visibility of the scope window
            %   associated with the System object H.
            value = isVisible(obj.cScope);
        end
    end

    methods
        function set.SampleRate(obj, value)
            validateattributes(value,{'double','single'}, {'scalar',...
                'positive','finite'},...
                '','SampleRate');
            obj.SampleRate = value;
        end
        function set.OperatingFrequency(obj,value)           
            sigdatatypes.validateFrequency(value,'',...
                'OperatingFrequency',{'double','single'},{'scalar'});
            obj.OperatingFrequency = value;
        end
        function set.PRF(obj, value)
            validateattributes(value,{'double','single'}, {'scalar',...
                'positive','finite'},...
                '','PRF');
            obj.PRF = value;
        end
        function set.PropagationSpeed(obj,value)
            sigdatatypes.validateSpeed(value,...
                '','PropagationSpeed',...
                {'double','single'},{'scalar','positive'});
            obj.PropagationSpeed = value;
        end
        function set.SweepSlope(obj, value)
            validateattributes(value,{'double','single'},...
                {'scalar','real','finite'},...
                '','SweepSlope');
            obj.SweepSlope = value;
        end
        function set.ReferenceRange(obj,value)
            sigdatatypes.validateDistance(value,'','ReferenceRange',...
                {'double','single'},{'scalar'});
            obj.ReferenceRange = value;
        end
        function set.Position(obj, value)
            if ~matlabshared.scopes.Validator.Position(value)
                error(message('Spcuilib:scopes:InvalidPosition'));
            end         
            obj.Position = value;        
        end
        function set.Name(obj, value)
             validateattributes(value, {'char', 'string'}, ...
                {'nonsparse'}, '', 'Name');
            obj.Name = value;
        end
        function set.RangeLabel(obj,value)
             validateattributes(value, {'char', 'string'}, ...
                {'nonsparse'}, '', 'RangeLabel');
            obj.RangeLabel = value;
        end
        function set.DopplerLabel(obj,value)
             validateattributes(value, {'char', 'string'}, ...
                {'nonsparse'}, '', 'DopplerLabel');
            obj.DopplerLabel = value;
        end     
    end
    
    methods(Access = protected)
        %% Common functions    
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and System block dialog
            if obj.IQDataInput
                flag = false;
                if strcmp(obj.RangeMethod,'Matched filter')
                    
                    if strcmp(prop,'DechirpInput') || ...
                            strcmp(prop,'SweepSlope') || ...
                            strcmp(prop,'RangeFFTLengthSource') || ...
                            strcmp(prop,'RangeFFTLength') || ...
                            strcmp(prop,'ReferenceRangeCentered') 
                        flag = true;
                    end
                end
                
                if strcmp(prop,'RangeLabel') || ...
                        strcmp(prop,'DopplerLabel')
                    flag = true;
                end
                
                if strcmp(prop,'PRF') && ...
                        ismember(obj.PRFSource,{'Auto'})
                    flag = true;
                end
                
                if  strcmp(obj.DopplerOutput,'Frequency') && ...
                        (strcmp(prop,'OperatingFrequency')|| ...
                        strcmp(prop,'SpeedUnits') || ...
                        (strcmp(prop,'FrequencyUnits') && obj.NormalizeDoppler))
                    flag = true;
                end
                
                if strcmp(obj.DopplerOutput,'Speed') && ...
                        (strcmp(prop,'NormalizeDoppler') || ...
                        strcmp(prop,'FrequencyUnits'))
                    flag = true;
                end
                           
            else
                flag = true;
                if strcmp(prop,'IQDataInput') || ...
                        strcmp(prop,'ResponseUnits') || ...
                        strcmp(prop,'Position') || ...
                        strcmp(prop,'Name') || ...
                        strcmp(prop,'RangeLabel') || ...
                        strcmp(prop,'DopplerLabel')
                    flag = false;
                end        
            end

        end

        function setupImpl(obj,varargin)
                        
            % Perform one-time calculations, such as computing constants
            
            obj.cScope.Name = obj.Name;
            obj.cScope.Position = obj.Position;
            
            % Set the Label to be shown for the response
            switch obj.ResponseUnits
                case 'magnitude'
                    obj.cScope.ColorBarLabel = getString(message('phased:scopes:MagLabel'));
                case 'power'
                    obj.cScope.ColorBarLabel = getString(message('phased:scopes:PowLabel'));
                case 'db'
                    obj.cScope.ColorBarLabel = getString(message('phased:scopes:PowdBLabel'));
            end
             
            if obj.IQDataInput
                % RangeProcessing : 'Matched filter', DopplerOutput : 'Speed'
                if strcmp(obj.RangeMethod,'Matched filter') && ...
                        strcmp(obj.DopplerOutput,'Speed')
                    obj.cResponse = phased.RangeDopplerResponse(...
                        'RangeMethod',obj.RangeMethod,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'SampleRate',obj.SampleRate,...
                        'PRFSource',obj.PRFSource,...
                        'ReferenceRange',obj.ReferenceRange,...
                        'DopplerOutput','Speed',...
                        'DopplerFFTLengthSource','Property',...
                        'DopplerFFTLength', obj.DopplerFFTLength,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'OperatingFrequency',obj.OperatingFrequency);
                    
                % RangeProcessing : 'Matched filter', DopplerOutput : 'Frequency'    
                elseif strcmp(obj.RangeMethod,'Matched filter') && ...
                        ~strcmp(obj.DopplerOutput,'Speed')
                    obj.cResponse = phased.RangeDopplerResponse(...
                        'RangeMethod',obj.RangeMethod,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'SampleRate',obj.SampleRate,...
                        'PRFSource',obj.PRFSource,...
                        'DopplerOutput','Frequency',...
                        'ReferenceRange',obj.ReferenceRange,...
                        'DopplerFFTLengthSource','Property',...
                        'DopplerFFTLength', obj.DopplerFFTLength',...
                        'PropagationSpeed',obj.PropagationSpeed);
                    
                % RangeProcessing : 'FFT', DopplerOutput : 'Speed'
                elseif ~strcmp(obj.RangeMethod,'Matched filter') && ...
                        strcmp(obj.DopplerOutput,'Speed')
                    obj.cResponse = phased.RangeDopplerResponse(...
                        'RangeMethod',obj.RangeMethod,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'DechirpInput',obj.DechirpInput,...
                        'SweepSlope',obj.SweepSlope,...
                        'SampleRate',obj.SampleRate,...
                        'PRFSource',obj.PRFSource,...
                        'DopplerOutput','Speed',...
                        'RangeFFTLengthSource','Property',...
                        'RangeFFTLength', obj.RangeFFTLength,...
                        'ReferenceRangeCentered',obj.ReferenceRangeCentered,...
                        'ReferenceRange',obj.ReferenceRange,...
                        'DopplerFFTLengthSource','Property',...
                        'DopplerFFTLength', obj.DopplerFFTLength,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'OperatingFrequency',obj.OperatingFrequency);
                
                % RangeProcessing : 'FFT', DopplerOutput : 'Frequency'
                elseif ~strcmp(obj.RangeMethod,'Matched filter') && ...
                        ~strcmp(obj.DopplerOutput,'Speed')
                    obj.cResponse = phased.RangeDopplerResponse(...
                        'RangeMethod',obj.RangeMethod,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'SweepSlope',obj.SweepSlope,...
                        'DechirpInput',obj.DechirpInput,...
                        'SampleRate',obj.SampleRate,...
                        'PRFSource',obj.PRFSource,...
                        'DopplerOutput','Frequency',...
                        'RangeFFTLengthSource','Property',...
                        'RangeFFTLength', obj.RangeFFTLength,...
                        'ReferenceRangeCentered',obj.ReferenceRangeCentered,...
                        'ReferenceRange',obj.ReferenceRange,...
                        'DopplerFFTLengthSource','Property',...
                        'DopplerFFTLength', obj.DopplerFFTLength,...
                        'PropagationSpeed',obj.PropagationSpeed);
                end
                
                if ~strcmp(obj.PRFSource,'Auto')
                    obj.cResponse.PRF = obj.PRF;
                end
                
                % Set the XLabel (Doppler axis Label)
                if (obj.NormalizeDoppler && ...
                        strcmp(obj.DopplerOutput,'Frequency'))
                    obj.cScope.XLabel = getString(message('phased:scopes:NormFreqLabel'));
                elseif strcmp(obj.DopplerOutput,'Frequency')
                    
                    switch obj.FrequencyUnits
                        case 'Hz'  % Hz
                            obj.pFrequencyConversion = 1;
                            obj.cScope.XLabel = getString(message('phased:scopes:HzLabel'));
                        case 'kHz' % Hz to kHz
                            obj.pFrequencyConversion = unitsratio('km', 'm');
                            obj.cScope.XLabel = getString(message('phased:scopes:kHzLabel'));
                        case 'MHz' % Hz to MHz
                            obj.pFrequencyConversion = unitsratio('km', 'mm');
                            obj.cScope.XLabel = getString(message('phased:scopes:MHzLabel'));
                    end 
                    
                elseif strcmp(obj.DopplerOutput,'Speed')
                    
                    TimeConversion = 1/3600; % seconds to hours
                    switch obj.SpeedUnits
                        case 'm/s'  % m/s
                            obj.pSpeedConversion = 1;
                            obj.cScope.XLabel = getString(message('phased:scopes:msecLabel'));
                        case 'km/h' % Hz to kHz
                            obj.pSpeedConversion = ...
                                unitsratio('km', 'm')/TimeConversion;
                            obj.cScope.XLabel = getString(message('phased:scopes:kmhLabel'));
                        case 'mph' % Hz to MHz
                            obj.pSpeedConversion = ...
                                unitsratio('mi', 'm')/TimeConversion;
                            obj.cScope.XLabel = getString(message('phased:scopes:mphLabel'));
                        case 'kt'
                            obj.pSpeedConversion = ...
                                unitsratio('nm', 'm')/TimeConversion;
                            obj.cScope.XLabel = getString(message('phased:scopes:ktLabel'));
                    end
                end
                               
                % Set YLabel (Range axis Label)
                switch obj.RangeUnits
                    case 'm' % m
                        obj.pRangeConversion = 1;
                        obj.cScope.YLabel = getString(message('phased:scopes:mLabel'));
                    case 'km' % m to km
                        obj.pRangeConversion = unitsratio('km', 'm');
                        obj.cScope.YLabel = getString(message('phased:scopes:kmLabel'));
                    case 'mi' % m to mi
                        obj.pRangeConversion = unitsratio('mi', 'm');
                        obj.cScope.YLabel = getString(message('phased:scopes:milesLabel'));                    
                    case 'nmi' % m to nmi
                        obj.pRangeConversion = unitsratio('nm', 'm');
                        obj.cScope.YLabel = getString(message('phased:scopes:nmiLabel'));                    
                end
            else
                % Set the XLabel (Doppler axis Label) and YLabel (Range
                % axis Label)
                obj.cScope.XLabel = obj.DopplerLabel;
                obj.cScope.YLabel = obj.RangeLabel;
            end

            % Set data cursor labels
            setCursorDataLabels(obj.cScope,["Doppler","Range","Intensity"]);
        end
        
        function validateInputsImpl(obj,varargin)
            % Validate inputs to the step method at initialization
            x = varargin{1};
            sz_x = size(x);
            
            % Input data check
            cond = ~isa(x,'float');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','float');
            end
            
            cond = ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','X');
            end
            if obj.IQDataInput             % 'IQDataInput' is true     
                % Check input to RangeDopplerResponse
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
                
                if strcmp(obj.RangeMethod,'Matched filter') % Matched filter
                    coeff = varargin{2};
                    cond = ~isa(coeff,'float');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:invalidInputDataType',...
                            'Coeff','float');
                    end
                    cond = ~iscolumn(coeff) || isempty(coeff);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:inputMustBeColVector','Coeff');
                    end
                else                                        % FFT
                    if obj.DechirpInput
                        xref = varargin{2};
                        cond = ~isa(xref,'float');
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:invalidInputDataType',...
                                'XRef','float');
                        end
                        cond = ~iscolumn(xref) || isempty(xref);
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:inputMustBeColVector','XRef');
                        end
                        cond = size(x,1)~=size(xref,1);
                        if cond
                            coder.internal.errorIf(cond,...
                                'phased:phased:NumRowsMismatch','X','XRef');
                        end
                    end
                end
            else                                  % 'IQDataInput' is false    
                rng_grid = varargin{2};
                dop_grid = varargin{3};
                
                cond = ~iscolumn(dop_grid) || isempty(dop_grid);
                if cond
                   coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeColVector','DOP'); 
                end
                
                cond = ~iscolumn(rng_grid) || isempty(rng_grid);
                if cond
                   coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeColVector','RANGE'); 
                end
                                
                cond = (sz_x(1) ~= numel(rng_grid));
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions',...
                        'RANGE',sz_x(1),numel(dop_grid));
                end
                
                cond = (sz_x(2) ~= numel(dop_grid));
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions',...
                        'DOP',sz_x(2),numel(dop_grid));
                end   
            end
        end
        
        function stepImpl(obj,varargin)
            % Type of response 
            switch obj.ResponseUnits
                case 'db'
                    unit = {'db'};
                case 'power'
                    unit = {'pow'};
                case 'magnitude'
                    unit = {'mag'};
            end
            
            if obj.IQDataInput                    % 'IQDataInput' is true            
                if ~strcmp(obj.RangeMethod, 'FFT')
                    % data,xref (MatchedFilter case)
                    [rdresp_out,rng_grid,dop_grid] = ...
                        obj.cResponse(varargin{:});
                else        
                    if ~obj.DechirpInput
                        % data (FFT case, no dechirp)
                        [rdresp_out,rng_grid,dop_grid] = ...
                            obj.cResponse(varargin{:});
                    else
                        % data,xref (FFT case, dechirp)
                        [rdresp_out,rng_grid,dop_grid] = ...
                            obj.cResponse(varargin{:});
                    end
                end
                
                % Range conversion m/km/mi/nmi
                rng_grid = rng_grid.*obj.pRangeConversion;
                
               
                if strcmp(obj.DopplerOutput,'Speed')
                    % If DopplerOutput is 'Speed' Doppler conversion
                    % (m/s)/(km/h)/(mph)/(knot)
                    
                    dop_grid = dop_grid.*obj.pSpeedConversion;
                    
                elseif strcmp(obj.DopplerOutput,'Frequency') && ...
                        ~obj.NormalizeDoppler
                    % If DopplerOutput is 'Frequency' Doppler conversion
                    % Hz/kHz/MHz
                    
                    dop_grid = dop_grid.*obj.pFrequencyConversion;
                end
                
            else                                  % 'IQDataInput' is false   
                
                % Get the response/ range grid/ Doppler grid
                rdresp_out = varargin{1};
                rng_grid = varargin{2};
                dop_grid = varargin{3};
            end
            
            % NormalizeDoppler
            if obj.IQDataInput && (strcmp(obj.DopplerOutput,'Frequency')...
                    && obj.NormalizeDoppler)
                % Doppler Sample rate
                fs_dop = obj.SampleRate/numel(rng_grid); 
                dop_grid = dop_grid./fs_dop;
                
                % Visualize
                plotMatrixdata(obj,rdresp_out,rng_grid,dop_grid,unit{:});
            else
                % Visualize
                plotMatrixdata(obj,rdresp_out,rng_grid,dop_grid,unit{:});
            end
        end
        
        function plotMatrixdata(obj,resp,rng_grid,...
                dop_grid,varargin)
            
            unit = varargin{:};
            
            % Matrix Viewer
            scope = obj.cScope;
                        
            response = phased.internal.computePlotPattern(...
                abs(resp),false,unit);
            
            scope.CustomYData = rng_grid;
            scope.CustomXData = dop_grid;
            
            % Visualization (Launch Scope)
            scope(response);
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            reset(obj.cResponse);
            reset(obj.cScope);
        end
        
        function releaseImpl(obj)
            % Release resources, such as file handles
            if obj.IQDataInput
                release(obj.cResponse);
            end
            release(obj.cScope);
        end
        
        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object
            
            obj.cScope.Name = obj.Name;
            obj.cScope.Position = obj.Position;
            
            if ~obj.IQDataInput
                obj.cScope.YLabel = obj.RangeLabel;
                obj.cScope.XLabel = obj.DopplerLabel;
            else
                obj.cResponse.ReferenceRange = obj.ReferenceRange;
            end
                
        end
               
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.cResponse = saveobj(obj.cResponse);
                s.cScope = saveobj(obj.cScope);
                s.pSpeedConversion = obj.pSpeedConversion;
                s.pFrequencyConversion = obj.pFrequencyConversion;
                s.pRangeConversion = obj.pRangeConversion;
            end  
        end

        function loadObjectImpl(obj,s,wasLocked)       
            s = loadSubObjects(obj,s,wasLocked);            
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
                if wasLocked
                    obj.cResponse = phased.RangeDopplerResponse.loadobj(s.cResponse);
                    obj.cScope = matlabshared.scopes.MatrixViewer.loadobj(s.cScope);
                    s = rmfield(s,'cScope');
                    s = rmfield(s,'cResponse');
                end

        end

        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end
        
        function flag = isInputComplexityLockedImpl(~,index)
            flag = true;  % index == 2 || index == 3
            if index == 1
                flag = false;
            end
        end 
        
        function num = getNumInputsImpl(obj)
            % Define total number of inputs for system with optional inputs
            num = 3;
            if obj.IQDataInput         
                if strcmp(obj.RangeMethod,'FFT') && ...
                        ~obj.DechirpInput
                    num = num-2;
                elseif (strcmp(obj.RangeMethod,'FFT') && ...
                        obj.DechirpInput) || ...
                        strcmp(obj.RangeMethod,'Matched filter')
                    num = num-1;
                end
            end
        end
        
        function num = getNumOutputsImpl(~)
            % Define total number of outputs for system with optional
            % outputs
            num = 0;
        end       
    end

    methods(Static, Access = protected)

        function groups = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
                                   
            groupScope = matlab.system.display.Section(...
                'Title','Scope Settings',...
                'PropertyList',{...
                'Name',...
                'Position',...
                'IQDataInput',...
                'ResponseUnits'...
                'RangeLabel',...
                'DopplerLabel'});
            
            rangegroups = matlab.system.display.Section(...
                'Title','Range Settings',...
                'PropertyList',{'RangeMethod','RangeUnits',...
                'PropagationSpeed','SampleRate','SweepSlope',...
                'DechirpInput','RangeFFTLength',...
                'ReferenceRangeCentered','ReferenceRange'});
            
            
            dopplergroups = matlab.system.display.Section(...
                'Title','Doppler Settings',...
                'PropertyList',{'PRFSource','PRF','DopplerFFTLength',...
                'DopplerOutput','OperatingFrequency',...
                'NormalizeDoppler','SpeedUnits','FrequencyUnits'});
            
            groupProcess = [rangegroups,dopplergroups];
            
            ResponseGroup = matlab.system.display.SectionGroup(...
                'Title','Main', ...
                'Sections',groupScope);
            
            ProcessGroup = matlab.system.display.SectionGroup(...
                'Title','Processing Settings', ...
                'Sections',groupProcess);
            
            groups = [ResponseGroup ProcessGroup];
        end

    end
    
    methods (Access = ?matlab.unittest.TestCase)      
        function this = getMatrixViewer(obj)        
            this = obj.cScope;
        end
    end
    
    methods (Static, Hidden)
        function flag = isAllowedInSystemBlock(~)
            flag = false;
        end
    end
end

% [EOF]