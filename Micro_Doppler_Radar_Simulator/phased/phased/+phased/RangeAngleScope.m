classdef (Sealed, StrictDefaults) RangeAngleScope < ...
        phased.internal.AbstractNarrowbandArrayProcessing 
% RangeAngleScope Visualize range-angle patterns
%   SCOPE = phased.RangeAngleScope returns a System object, SCOPE, that
%   displays range-angle response patterns in rectangular coordinates.
%
%   SCOPE = phased.RangeAngleScope ('Name', Value,) returns a range-angle
%   scope System object, SCOPE, with each specified property name set to
%   the specified value. You can specify name-value pair arguments in any
%   order as (Name 1, Value 1, ..., NameN, ValueN).
%
%   Step method syntax when IQDataInput property is false:
%
%   step(SCOPE,X,RANGE,ANG) displays the range-angle response pattern of X,
%   in the RangeAngleScope figure. X is KxP matrix  containing the
%   range-angle response where K denotes the number of fast time samples, P
%   is the number of angle samples.
%
%   RANGE is a length-M column vector containing the range values, in
%   meters, at which the range-angle response is evaluated. ANG is a
%   length-P column vector containing angles, in degrees, at which the
%   range-angle response is evaluated.
%
%   Step method syntax when IQDataInput property is true:
%
%   step(SCOPE,X) calculates and displays the range angle map of input X.
%   This syntax applies when you set the RangeMethod property to 'FFT' and
%   the DechirpInput property to false. This syntax is most commonly used
%   with FMCW signals.
%
%   X is a dechirped input signal. X must be a KxL matrix, where K denotes
%   the number of fast time samples, and L is the number of subarrays if
%   SensorArray contains subarrays, or the number of elements, otherwise.
%
%   step(SCOPE,X,XREF) uses input XREF as the reference signal to dechirp
%   the input signal X to calculate and display the range-angle response of
%   X. This syntax applies when you set the RangeMethod property to 'FFT'
%   and the DechirpInput property to true. This syntax is most commonly
%   used with FMCW signals and the reference signal is, in general, the
%   transmitted signal.
%
%   X is an input signal to be dechirped by the RangeAngleScope object. X
%   must be a KxL matrix, where K denotes the number of fast time samples,
%   and L is the number of subarrays if SensorArray contains subarrays, or
%   the number of elements, otherwise.
%
%   step(SCOPE,X,COEFF) uses COEFF as the matched filter coefficients, to
%   matched filer the input X and display the range-angle response of X.
%   This method applies when you set the RangeMethod property to 'Matched
%   filter'. This syntax is most commonly used with pulsed signals.
%
%   X must be a KxL matrix, where K denotes the number of fast time
%   samples, and L is the number of channels (antenna elements or
%   beams),that is, the number of subarrays if SensorArray contains
%   subarrays, or the number of elements otherwise.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(SCOPE, x) and y = SCOPE(x) are
%   equivalent
%
%   RangeAngleScope methods:
%
%   step         - Display range-angle response
%   release      - Allow property value and input characteristics changes
%   clone        - Create range-angle scope object with same property
%                  values
%   isLocked     - Locked status (logical)
%   show         - Turn on visibility of the scope
%   hide         - Turn off visibility of the scope
%   isVisible    - Return visibility of the scope (logical)
%
%   RangeAngleScope properties:
%
%   Name                         - Scope window name
%   Position                     - Scope window position 
%   IQDataInput                  - Type of input
%   ResponseUnits                - Output response unit
%   RangeLabel                   - Range axis label
%   AngleLabel                   - Angle axis label
%   SensorArray                  - Sensor array
%   RangeMethod                  - Range processing method
%   PropagationSpeed             - Propagation speed
%   OperatingFrequency           - Operating frequency
%   RangeUnits                   - Range unit
%   SampleRate                   - Sample rate
%   SweepSlope                   - FM sweep slope
%   DechirpInput                 - Dechirp input signal
%   RangeFFTLength               - FFT length in range processing
%   ReferenceRangeCentered       - Set reference range at center
%   ReferenceRange               - Reference range
%   ElevationAngle               - Elevation angle
%   AngleSpan                    - Angle span
%   NumAngleSamples              - Number of angles 
%
%   % Example:
%   %   Calculate and visualize the range-angle response from a pulsed
%   %   radar transmitting a rectangular waveform using the matched filter
%   %   approach. The signal includes three target returns. Two are 
%   %   approximately 2000 m away and the third is approximately 3500 m 
%   %   away. In addition, two targets are stationary relative to the  
%   %   radar while the third is moving away from the radar at 
%   %   approximately 100 m/s. The signals arrive at an 8-element uniform
%   %   linear array.
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
%   % Create range-angle scope for processing and visualization.
%   
%   rngangscope = phased.RangeAngleScope(...
%       'IQDataInput',true,'RangeMethod','Matched filter',...
%       'Name','Range-Angle Scope','ResponseUnits','magnitude',...
%       'Position',[560 375 560 420],'RangeUnits','m',...
%       'SensorArray',antennaarray,'OperatingFrequency',fc,...
%       'SampleRate',fs,'PropagationSpeed',propspeed);
%
%   rngangscope(rxdata,mfcoeffs);
%
%   See also phased, phased.RangeDopplerScope, phased.AngleDopplerScope,
%   phased.RangeAngleResponse.

%   Copyright 2018 The MathWorks, Inc.
  
    properties (Nontunable)

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
    end
   
    properties (Nontunable)
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a positive scalar. This
        %   property applies when you set the 'IQDataInput' property true.
        %   The default value of this property is 1e6 (1 MHz).
        SampleRate = 1e6;
        %SweepSlope     FM sweep slope (Hz/s)
        %   Specify the slope of the linear FM sweeping (in Hz/s) as a
        %   scalar. This property applies when you set the RangeMethod
        %   property to 'FFT' and 'IQDataInput' property true. The default
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
        %   scalar. This property applies when you set the 'IQDataInput'
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
        %ElevationAngle     Elevation angle (deg)
        %   Specify the elevation angle (in degrees) used to calculate the
        %   range-angle response as a scalar. The angle must be between -90
        %   and 90. This property applies when you set the 'IQDataInput'
        %   property true. The default value of this property is 0.
        ElevationAngle = 0;
        %AngleSpan      Angle span (deg)
        %   Specify the angle span (in degrees) used to calculate the
        %   range-angle response as a 2-element row vector in the form of
        %   [min_angle max_angle]. This property applies when you set the
        %   'IQDataInput' property true. The default value of this property
        %   is [-90 90].
        AngleSpan = [-90 90]
        %NumAngleSamples    Number of angle bins
        %   Specify the number of samples in angular domain used to
        %   calculate the angle-Doppler response as a positive integer.
        %   This value must be greater than 2. This property applies when
        %   you set the 'IQDataInput' property true. The default value of
        %   this property is 256.
        NumAngleSamples = 256;
    end
    
    properties(Nontunable,Logical)
        %IQDataInput Specify the type of input
        %   Specify whether the input is I/Q (raw) data or processed data.
        %   If you set this property true, it implies that raw data is
        %   passed and processing is done along range and angle domain
        %   before visualizing. Set this property false if you have
        %   processed data. The default value is true.
        IQDataInput = true;
    end

    properties (Nontunable)
        %ResponseUnits Unit of the response
        %   Specify the unit of the plot, using one of | 'db' | 'magnitude'
        %   | 'power'|. The default is 'db'.
        ResponseUnits = 'db';
    end
    
    properties (Nontunable)
        %RangeUnits  Range unit
        %   Specify the range unit, using one of | 'm' | 'km' | 'mi'|
        %   'nmi'|. This property applies when you set the 'IQDataInput'
        %   property true. The default is 'm'.
        RangeUnits = 'm';
    end

    properties
        %RangeLabel Range-axis label
        %   Specify the range-axis label as a string. The default value is
        %   'Range (m)'. This property applies when you set the
        %   'IQDataInput' property false. This property is tunable.
        RangeLabel = 'Range (m)';
        
        %AngleLabel Angle-axis label
        %   Specify the angle-axis label as a string. This property applies
        %   when you set the 'IQDataInput' property false. The default
        %   value is 'Angle (degrees)'. This property is tunable.
        AngleLabel = 'Angle (degrees)';
        
        %Name Caption to display on scope window
        %   Specify the caption to display on the scope window as any
        %   string. The default value is 'Range-Angle Scope'. This property
        %   is tunable.
        Name = 'Range-Angle Scope';
        
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
    end
  
    properties (Constant, Hidden)
        ResponseUnitsSet = matlab.system.StringSet({'db','power',...
            'magnitude'});
        RangeUnitsSet = matlab.system.StringSet({'m','km','mi','nmi'});
        RangeMethodSet = matlab.system.StringSet({'Matched filter','FFT'});
    end
    
    methods
        % Constructor
        function obj = RangeAngleScope(varargin)
            % Support name-value pair arguments when constructing
            % RangeAngleScope object
 
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
            
            % matlabshared.scopes.MatrixViewer - Visualization 
            obj.cScope = matlabshared.scopes.MatrixViewer(...
                'XDataMode','Custom',...
                'YDataMode','Custom',...
                'AxisOrigin','Lower left corner',...
                'XLabel',getString(message('phased:scopes:AngLabel')),...
                'Colormap','parula',...
                'Title',getString(message('phased:scopes:rngangLabel'))...
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
            validateattributes(value,{'double'}, {'scalar',...
                'positive','finite'},...
                '','SampleRate');
            obj.SampleRate = value;
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
        function set.AngleLabel(obj,value)
             validateattributes(value, {'char', 'string'}, ...
                {'nonsparse'}, '', 'AngleLabel');
            obj.AngleLabel = value;
        end
        function set.RangeLabel(obj,value)
             validateattributes(value, {'char', 'string'}, ...
                {'nonsparse'}, '', 'RangeLabel');
            obj.RangeLabel = value;
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
                        strcmp(prop,'AngleLabel')
                    flag = true;
                end
                                           
            else
                flag = true;
                
                if strcmp(prop,'IQDataInput') || ...
                        strcmp(prop,'ResponseUnits') || ...
                        strcmp(prop,'Name') || ...
                        strcmp(prop,'Position') || ...
                        strcmp(prop,'RangeLabel') || ...
                        strcmp(prop,'AngleLabel')
                    flag = false;
                end
                
            end
        end
   
        function validateInputsImpl(obj,varargin)
            % Validate inputs to the step method at initialization
            x = varargin{1};
            sz_x = size(x);
            
            % Input data check
            cond = ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','double');
            end
            
            cond = ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','X');
            end
            
            if obj.IQDataInput
                % Check input to RangeAngleResponse
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
                
                N = getDOF(obj.SensorArray);
                
                cond = sz_x(2) ~= N;
                if cond
                    coder.internal.errorIf(cond,...
                        'phased:phased:invalidColumnNumbers','X', N);
                end
                               
                if (obj.RangeMethod(1) == 'M') %Matched filter
                    coeff = varargin{2};
                    cond = ~isa(coeff,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:invalidInputDataType',...
                            'Coeff','double');
                    end
                    cond = ~iscolumn(coeff) || isempty(coeff);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:inputMustBeColVector','Coeff');
                    end
                else
                    if obj.DechirpInput
                        xref = varargin{2};
                        cond = ~isa(xref,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:invalidInputDataType',...
                                'XRef','double');
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
                
            else
                rng_grid = varargin{2};
                ang_grid = varargin{3};
                
                cond = ~iscolumn(ang_grid) || isempty(ang_grid);
                if cond
                   coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeColVector','ANG'); 
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
                        'RANGE',sz_x(1),numel(rng_grid));
                end
                
                cond = (sz_x(2) ~= numel(ang_grid));
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions',...
                        'ANG',sz_x(2),numel(ang_grid));
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
                   % RangeProcessing : 'Matched filter'
                if strcmp(obj.RangeMethod,'Matched filter')
                    obj.cResponse = phased.RangeAngleResponse(...
                        'SensorArray',obj.SensorArray,...
                        'RangeMethod',obj.RangeMethod,...
                        'SampleRate',obj.SampleRate,...
                        'ReferenceRange',obj.ReferenceRange,...
                        'ElevationAngle',obj.ElevationAngle,...
                        'NumAngleSamples',obj.NumAngleSamples,...
                        'AngleSpan',obj.AngleSpan,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'OperatingFrequency',obj.OperatingFrequency);
                    
                else
                    % RangeProcessing : 'FFT'
                    obj.cResponse = phased.RangeAngleResponse(...
                        'SensorArray',obj.SensorArray,...
                        'RangeMethod',obj.RangeMethod,...
                        'DechirpInput',obj.DechirpInput,...
                        'SweepSlope',obj.SweepSlope,...
                        'SampleRate',obj.SampleRate,...
                        'RangeFFTLengthSource','Property',...
                        'RangeFFTLength',obj.RangeFFTLength,...
                        'ReferenceRangeCentered',obj.ReferenceRangeCentered,...
                        'ReferenceRange',obj.ReferenceRange,...
                        'ElevationAngle',obj.ElevationAngle,...
                        'NumAngleSamples',obj.NumAngleSamples,...
                        'AngleSpan',obj.AngleSpan,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'OperatingFrequency',obj.OperatingFrequency);
                end
                
                % Set YLabel (Range axis Label)
                switch obj.RangeUnits
                    case 'm'
                        obj.pRangeConversion = 1;
                        obj.cScope.YLabel = getString(message('phased:scopes:mLabel'));
                    case 'km'
                        obj.pRangeConversion = unitsratio('km', 'm');
                        obj.cScope.YLabel = getString(message('phased:scopes:kmLabel'));
                    case 'mi'
                        obj.pRangeConversion = unitsratio('mi', 'm');
                        obj.cScope.YLabel = getString(message('phased:scopes:milesLabel'));
                    case 'nmi'
                        obj.pRangeConversion = unitsratio('nm', 'm');
                        obj.cScope.YLabel = getString(message('phased:scopes:nmiLabel'));
                end
            else          
                obj.cScope.XLabel = obj.AngleLabel;
                obj.cScope.YLabel = obj.RangeLabel;
                
            end

            % Set data cursor labels
            setCursorDataLabels(obj.cScope,["Angle","Range","Intensity"]);
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
            
            % Style of response - 'Rectangular'           
            if obj.IQDataInput                   % 'IQDataInput' is true
                
                % Process raw data using phased.RangeAngleResponse
                if ~strcmp(obj.RangeMethod, 'FFT')
                    % data,xref (MatchedFilter case)
                    [raresp_out,rng_grid,ang_grid] = ...
                        obj.cResponse(varargin{:});
                else
                    if ~obj.DechirpInput
                        % data (FFT case, no dechirp)
                        [raresp_out,rng_grid,ang_grid] = ...
                            obj.cResponse(varargin{:});
                    else
                        % data,xref (FFT case, dechirp)
                        [raresp_out,rng_grid,ang_grid] = ....
                            obj.cResponse(varargin{:});
                    end
                end
                rng_grid = rng_grid.*obj.pRangeConversion;
            else                               % 'IQDataInput' is false
                
                % Get the response/ range grid/ angle grid
                raresp_out = varargin{1};
                rng_grid = varargin{2};
                ang_grid = varargin{3};
            end
            
            % Visualization
            plotMatrixdata(obj,raresp_out,rng_grid,ang_grid,unit{:});
        end
        
        function plotMatrixdata(obj,resp,rng_grid,ang_grid,varargin)
            
            unit = varargin{:};
            
            % Matrix Viewer
            scope = obj.cScope;
            
            % Compute response based on 'ResponseUnits'
            response = phased.internal.computePlotPattern(...
                abs(resp),false,unit);
            
            scope.CustomXData = ang_grid;
            scope.CustomYData = rng_grid;
            
            % Visualization (Launch Scope)
            scope(response)
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            reset(obj.cResponse);
            reset(obj.cScope);
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            if obj.IQDataInput
                releaseImpl(obj.cResponse);
            end
            releaseImpl(obj.cScope);
        end

        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object     
            obj.cScope.Name = obj.Name;
            obj.cScope.Position = obj.Position;
            
            if ~obj.IQDataInput
                obj.cScope.YLabel = obj.RangeLabel;
                obj.cScope.XLabel = obj.AngleLabel;
            else
                obj.cResponse.ReferenceRange = obj.ReferenceRange;
            end            
        end
        
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@...
                phased.internal.AbstractNarrowbandArrayProcessing(obj);
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.cScope = saveobj(obj.cScope);
                s.cResponse = saveobj(obj.cResponse);
                s.pRangeConversion = obj.pRangeConversion;
            end
        end

        function loadObjectImpl(obj,s,~)
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@...
                phased.internal.AbstractNarrowbandArrayProcessing(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cScope = matlabshared.scopes.MatrixViewer.loadobj(s.cScope);
                    obj.cResponse = phased.RangeAngleResponse.loadobj(s.cResponse);
                    s = rmfield(s,'cScope');
                    s = rmfield(s,'cResponse');
                end
                s = rmfield(s,'isLocked');
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
            num = 3 ;
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
        
        function group = getPropertyGroupsLongImpl(obj)
            allGroups = getPropertyGroupsLongImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            group = allGroups;
            % Shuffle to group similar properties
            if obj.IQDataInput
                group.PropertyList = [allGroups.PropertyList(2:5),...
                    allGroups.PropertyList(1),...
                    allGroups.PropertyList(6:end)];
            end
        end
        
    end

    methods(Static, Access = protected)
        
        function groups = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
            groups = ...
                getPropertyGroupsImpl@phased.internal.AbstractNarrowbandArrayProcessing('subarray');
            
            groupScope = matlab.system.display.Section(...
                'Title','Scope Settings',...
                'PropertyList',{...
                'Name',...
                'Position',...
                'IQDataInput',...
                'ResponseUnits',...
                'RangeLabel',...
                'AngleLabel'});
            
            rangeProps = [{'RangeMethod','RangeUnits'},groups(1).PropertyList,...
                {'SampleRate','SweepSlope',...
                'DechirpInput','RangeFFTLength',...
                'ReferenceRangeCentered','ReferenceRange'}];
                
            rangegroups = matlab.system.display.Section(...
                'Title','Range Settings',...
                'PropertyList',rangeProps);
  
            anglegroups = matlab.system.display.Section(...
                'Title','Angle Settings',...
                'PropertyList',{'ElevationAngle',...
                'AngleSpan',...
                'NumAngleSamples'});
            
            groupProcess = [rangegroups,anglegroups];
            
            group = matlab.system.display.SectionGroup(...
                'Title','Main',...
                'Section',groupScope);
            
            ProcessGroup = matlab.system.display.SectionGroup(...
                'Title','Processing Settings', ...
                'Sections',groupProcess);
            
            sensorGroup = groups(2);

            groups = [group ProcessGroup sensorGroup];

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