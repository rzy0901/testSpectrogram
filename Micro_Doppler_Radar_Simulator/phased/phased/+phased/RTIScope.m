classdef (Sealed, StrictDefaults) RTIScope < matlab.System   
%RTIScope Display scrolling range intensity vectors or arrays
%   RangeScope = phased.RTIScope returns a System object, RangeScope, that
%   can show a scrolling plot of range intensity vectors versus time.
%
%   RangeScope = phased.RTIScope('Name', Value, ...) returns a RTIScope
%   System object, RangeScope, with each specified property name set to the
%   specified value. You can specify name-value pair arguments in any order
%   as (Name 1, Value 1, ..., Name N, Value N).
%
%   Step method syntax when IQDataInput is false:
%
%   step(RangeScope, X) displays the real signal, X. X is an N by M matrix
%   where N is the number of intensity bins in a range intensity vector and
%   M is the number of intensity vectors. N should be greater than 1 and M
%   is equal or greater than 1. Each column of the matrix represents an
%   intensity vector from successive times. The time between intensity
%   vectors is specified in the TimeResolution property.
%
%   Step method syntax when IQDataInput is true:
%
%   step(RangeScope,X) calculates the range response of the input data X,
%   and displays the processed signal. This syntax applies when you set the
%   RangeMethod property to 'FFT' and then the DechirpInput property to
%   false. This syntax is most commonly used with FMCW signals.
%
%   X is a dechirped input signal. X is an N by M matrix where N is the
%   number of intensity bins in an intensity vector and M is the number of
%   intensity vectors. N should be greater than 1 and M is equal or greater
%   than 1. Each column of the matrix represents an intensity vector from
%   successive times. The time between intensity vectors are specified in
%   the TimeResolution property.
%   
%   step(RangeScope,X,XREF) uses input XREF as the reference signal to
%   dechirp the input signal X. This syntax applies when you set the
%   RangeMethod property to 'FFT' and then the DechirpInput property to
%   true. This syntax is most commonly used with FMCW signals and the
%   reference signal is, in general, the transmitted signal.
%
%   X is an input signal to be dechirped by the RTIScope object. X is an N
%   by M matrix where N is the number of intensity bins in an intensity
%   vector and M is the number of intensity vectors. N should be greater
%   than 1 and M is equal or greater than 1. Each column of the matrix
%   represents an intensity vector from successive times. The time between
%   intensity vectors are specified in the TimeResolution property. XREF
%   must be a column vector whose number of rows is the same as the number
%   of rows of X. XREF is the reference signal used to dechirp X.
%
%   step(RangeScope,X,COEFF) uses COEFF as the matched filter coefficients.
%   This method applies when you set the RangeMethod property to 'Matched
%   filter'. This syntax is most commonly used with pulsed signals.
%
%   X is an input signal to be match filtered by the RTIScope object. X is
%   an N by M matrix where N is the number of intensity bins in an
%   intensity vector and M is the number of intensity vectors. N should be
%   greater than 1 and M is equal or greater than 1. Each column of the
%   matrix represents an intensity vector from successive times. The time
%   between intensity vectors are specified in the TimeResolution property.
%   COEFF must be a column vector containing the matched filter
%   coefficients.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(RangeScope, x) and y =
%   RangeScope(x) are equivalent.
%
%   RTIScope methods:
%
%   step      - Display signal in the RTIScope figure 
%   release   - Allow property value and input characteristics changes, and
%               release RTIScope resources
%   clone     - Create RTIScope object with same property values
%   isLocked  - Display locked status (logical)
%   show      - Turn on visibility of RTIScope figure
%   hide      - Turn off visibility of RTIScope figure
%   isVisible - Return visibility of the scope (logical)
%
%   RTIScope properties:
%
%   Name                    - Scope window name
%   Position                - Scope window position
%   IQDataInput             - Type of input
%   RangeLabel              - Range Label
%   RangeResolution         - Range sample spacing
%   RangeOffset             - Range display offset
%   TimeResolution          - Time resolution    
%   TimeSpan                - Time span
%   IntensityUnits          - Intensity units
%   RangeMethod             - Range processing method
%   PropagationSpeed        - Propagation speed
%   SampleRate              - Sample rate
%   SweepSlope              - FM sweep slope
%   DechirpInput            - Dechirp input signal
%   RangeFFTLength          - FFT length in range processing
%   ReferenceRangeCentered  - Set reference range at center
%
%   % Example: 
%   %   Plot a range-time intensity for three simulated targets.
%    
%   % Function to simulate received range bins
%   txpow = 200; gain = 2e8; std = 5;
%   rangePow = @(bins,range) ...
%       gain.*exp(-0.5*((bins-range)./std).^2).* ...
%       txpow./(range.^4)./(sqrt(2*pi).*std);
% 
%   % Setup the scope
%   scope = phased.RTIScope( ...
%               'IQDataInput',false,...
%               'Name','Range-Time Intensity Scope',...
%               'Position',[560 375 560 420],...
%               'RangeLabel','Range (m)', ...
%               'RangeResolution',1, ... %Assume range bin is 1 meter
%               'TimeResolution',0.05,'TimeSpan',5, ...
%               'IntensityUnits','magnitude');
%   
%   % Create range bins for three targets, two moving in opposite 
%   % directions
%   x = 0:1000;
%   ranges(:,1) = 250:10:850;
%   ranges(:,2) = 850:-10:250;
%   ranges(:,3)  = 400;
%   for k = 1:size(ranges,1)
%       y = ((rangePow(x,ranges(k,1))+ ...
%             rangePow(x,ranges(k,2)) + ...
%             rangePow(x,ranges(k,3))));
%       scope(y.');
%   end
%
%   See also phased, phased.DTIScope, phased.IntensityScope

%   Copyright 2018 The MathWorks, Inc.

    properties (Logical, Nontunable)
        %IQDataInput Specify the type of input
        %   Specify whether the input is I/Q (raw) data or processed data.
        %   If you set this property true, it implies that raw data is
        %   passed and range processing is done before visualizing. Set
        %   this property false if you have processed data. The default
        %   value is false.
        IQDataInput = false
    end

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
        %   signal.  The FFT approach is often used with FMCW signals.
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
        %   property to 'FFT' and 'IQDataProperty' true. The default value
        %   is 1e9.
        SweepSlope = 1e9;
        %PropagationSpeed   Propagation speed (m/s)
        %   Specify the propagation speed (in m/s) of the signal as a
        %   scalar. This property applies when you set the 'IQDataInput'
        %   property true. The default value of this property is the speed
        %   of light.
        PropagationSpeed = physconst('lightspeed');
    end

    properties (Nontunable,Logical)
        %DechirpInput     Dechirp input signal
        %   Set this property to true to dechirp the input signal first
        %   before range processing. Set this property to false to indicate
        %   that the input signal is already dechirped and no dechirp
        %   operation is necessary. This property applies when you set the
        %   RangeMethod property to 'FFT' and 'IQDataInput' true. The
        %   default value of this property is false.
        DechirpInput = false;
    end

    properties (Nontunable, PositiveInteger)
        %RangeFFTLength     FFT length in range processing
        %   Specify the FFT length in range domain as a positive integer.
        %   This property applies when you set the RangeMethod property to
        %   'FFT' and 'IQDataInput' true. The default value of this
        %   property is 1024.
        RangeFFTLength = 1024
    end

    properties (Nontunable, Logical)
        %ReferenceRangeCentered     Set reference range at center
        %   Set this property to true to set the reference range to the
        %   center of the range grid. Set this property to false to set the
        %   reference range to the beginning of the range grid. The default
        %   value of this property is true. This property only applies when
        %   you set the RangeMethod to 'FFT' and 'IQDataInput' true.
        ReferenceRangeCentered = true
    end

    properties(Nontunable)
        %RangeOffset Range offset
        %   Specify the offset to apply to the range. This
        %   property is a numeric scalar with default value 0 m.
        RangeOffset = 0;      
        %RangeResolution  Range resolution
        %   Specify the spacing between the range samples  as a finite
        %   numeric scalar. The default value of this property is
        %   1 m.
        RangeResolution = 1;
    end

    properties (Nontunable)
        %TimeResolution Time resolution (s)
        %   Specify the time resolution of each intensity line as a
        %   positive scalar in seconds. The default is 1 ms.
        TimeResolution = 1e-3
        %TimeSpan Time span (s)
        %   Specify the time span of the intensity display as a positive
        %   scalar in seconds. The default is 100 ms.
        TimeSpan = 100e-3
        %IntensityUnits Intensity units
        %   Specify the unit of the intensity, using one of | 'db' |
        %   'magnitude' | 'power'|. The default is 'db'.
        IntensityUnits = 'db';  
    end

    properties 
        %RangeLabel Range-axis label
        %   Specify the range-axis label as a string. The default value is
        %   'Range (m)'. This property is tunable.
        RangeLabel = 'Range (m)';
        %Name Specify scope window title
        %   Specify the caption to display on the scope window as any
        %   string. The default value of this property is 'Range Time
        %   Intensity Scope'. This property is tunable.
        Name = 'Range Time Intensity Scope';
        %Position Scope window position in pixels
        %   Specify the size and location of the scope window in pixels,
        %   as a four-element double vector of the form: [left bottom width
        %   height]. The default value of obj property is dependent on the
        %   screen resolution, and is such that the window is positioned in
        %   the center of the screen, with a width and height of 800 and
        %   450 pixels respectively. This property is tunable.
        Position = [560 375 800 450];
    end

    properties (Access = private)
       cResponse
       cScope
    end

    properties (Constant, Hidden)
        RangeMethodSet = matlab.system.StringSet({'Matched filter','FFT'});
        IntensityUnitsSet = matlab.system.StringSet({'db','power',...
                'magnitude'});
    end

    methods   

        function obj = RTIScope(varargin)
            setProperties(obj, nargin, varargin{:});
            
            % matlabshared.scopes.MatrixViewer - Visualization
            obj.cScope = matlabshared.scopes.MatrixViewer(...
                'XDataMode','Offset and resolution',...
                'YDataMode','Span and resolution',...
                'AxisOrigin','Lower left corner',...
                'ColorBarLocation','northoutside',...
                'Colormap','parula',...
                'Title',getString(message('phased:scopes:rtiLabel')),...
                'YLabel',getString(message('phased:scopes:timeLabel')));
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
        function set.TimeResolution(obj, value)
            validateattributes(value,{'double'}, ...
                {'positive','real','scalar','finite','nonnan'},...
                '','TimeResolution');
            obj.TimeResolution = value;
        end

        function set.TimeSpan(obj, value)
            validateattributes(value,{'double'}, ...
                {'positive','real','scalar','finite','nonnan'},...
                '','TimeSpan');
            obj.TimeSpan = value;
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

        function set.SweepSlope(obj, value)
            validateattributes(value,{'double','single'},...
                {'scalar','real','finite'},...
                '','SweepSlope');
            obj.SweepSlope = value;
        end  
    end
    
    methods (Access = protected)
        %% Common functions  
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog

            if ~obj.IQDataInput  
                flag = true;

                if strcmp(prop,'IQDataInput')
                    flag = false;
                end
                
                if strcmp(prop,'Name') || strcmp(prop,'RangeOffset') ||...
                        strcmp(prop,'RangeResolution') || ...
                        strcmp(prop,'RangeOffset') ||...
                        strcmp(prop,'RangeLabel') || ...
                        strcmp(prop,'TimeSpan') ||...
                        strcmp(prop,'TimeResolution') || ...
                        strcmp(prop,'IntensityUnits')||...
                        strcmp(prop,'Position') || strcmp(prop,'Title')
                    flag = false;
                end
            else
                flag = false;
                
                if strcmp(obj.RangeMethod,'Matched filter') && ...
                        (strcmp(prop,'RangeFFTLength') ...
                        || strcmp(prop,'DechirpInput') || ...
                        strcmp(prop,'SweepSlope') ||...
                        strcmp(prop,'ReferenceRangeCentered'))
                    flag = true;
                end                       
            end
        end

        function setupImpl(obj,varargin)

            % MatrixViewer
            scope = obj.cScope;
            
            % Set properties
            scope.XOffset = obj.RangeOffset;
            scope.XResolution = obj.RangeResolution;        
            scope.YSpan = obj.TimeSpan;
            scope.YResolution = obj.TimeResolution;
            switch obj.IntensityUnits
                case 'magnitude'
                    obj.cScope.ColorBarLabel = getString(message('phased:scopes:MagLabel'));
                case 'power'
                    obj.cScope.ColorBarLabel = getString(message('phased:scopes:PowLabel'));
                case 'db'
                    obj.cScope.ColorBarLabel = getString(message('phased:scopes:PowdBLabel'));
            end
            scope.XLabel = obj.RangeLabel; %'Range(meters)'
            scope.Name = obj.Name;
            scope.Position = obj.Position;
                             
            if obj.IQDataInput                 % IQDataInput - true 
                
                % Range processing - Matched Filter
                if strcmp(obj.RangeMethod,'Matched filter')
                    obj.cResponse = phased.RangeResponse(...
                        'RangeMethod',obj.RangeMethod,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'SampleRate',obj.SampleRate);
                else
                    % Range processing - FFT
                    obj.cResponse = phased.RangeResponse(...
                        'RangeMethod',obj.RangeMethod,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'SampleRate',obj.SampleRate,...
                        'SweepSlope',obj.SweepSlope,...
                        'DechirpInput',obj.DechirpInput,...
                        'RangeFFTLengthSource','Property',...
                        'RangeFFTLength',obj.RangeFFTLength,...
                        'ReferenceRangeCentered',obj.ReferenceRangeCentered);
                end
            end

            % Set data cursor labels
            setCursorDataLabels(scope,["Range","Time","Intensity"]);
        end

        function validateInputsImpl(obj,x,xref)
            
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
            
            cond = sz_x(1) < 2;
            if cond
                coder.internal.errorIf(cond,...
                    'phased:scopes:InvalidInputDimension');
            end
            
            if obj.IQDataInput
                if strcmp(obj.RangeMethod,'Matched filter') % Matched filter
                    coeff = xref;
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
                        cond = sz_x(1) ~= size(xref,1);
                        if cond
                            coder.internal.errorIf(cond,...
                                'phased:phased:NumRowsMismatch','X','XRef');
                        end
                    end
                end
            end
        end

        function stepImpl(obj,varargin) 
                            
            data = varargin{1}; 
            
            if obj.IQDataInput  % IQDataInput - true
                
                % Range processing - Matched Filter
                if ~strcmp(obj.RangeMethod,'FFT')
                    xref = varargin{2};
                    % data,xref (MatchedFilter case)
                    [resp,~] = obj.cResponse(data,xref);
                                     
                else  % Range processing - FFT
                    
                    if obj.DechirpInput
                        xref = varargin{2};
                        % data,xref (FFT case, dechirp)
                        [resp,~] = obj.cResponse(data,xref);
                    else
                        % data (FFT case, no dechirp)
                        [resp,~] = obj.cResponse(data);
                    end                    
                end
                resp = abs(resp);
            else               % IQDataInput - false               
                resp = abs(data);
            end
            
            switch obj.IntensityUnits
                case 'magnitude'
                    plotResp = resp;
                case 'power'
                    plotResp = (resp).^2;                   
                case 'db'
                    plotResp = mag2db(resp);
            end
            
            % Visualization (Launch Scope)
            obj.cScope(plotResp.');  
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            reset(obj.cResponse);
            reset(obj.cScope);
        end

        function releaseImpl(obj)
            release(obj.cScope);
            if obj.IQDataInput
                release(obj.cResponse);
            end
        end

        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object
            obj.cScope.Name = obj.Name;
            obj.cScope.Position = obj.Position;
        end

        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.cResponse = obj.cResponse;
                s.cScope = obj.cScope;
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s
            s = loadSubObjects(obj,s,wasLocked);        
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                obj.cScope = matlabshared.scopes.MatrixViewer.loadobj(s.cScope);
                obj.cResponse = phased.RangeResponse.loadobj(s.cResponse);
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
            num = 2;
            if ~obj.IQDataInput
                num = num-1;
            elseif (obj.IQDataInput && strcmp(obj.RangeMethod,'FFT') && ...
                    ~obj.DechirpInput)
                num = num-1;
            end

        end

        function num = getNumOutputsImpl(~)
            % Define total number of outputs for system with optional
            % outputs
            num = 0;
        end
    end

    methods(Access = protected, Static,Hidden)
        
        function groups = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
            
            groupRange =  matlab.system.display.SectionGroup(...
                'Title','Processing Settings',...
                'PropertyList',{'RangeMethod',...
                'PropagationSpeed',...
                'SampleRate',...
                'SweepSlope',...
                'DechirpInput',...
                'RangeFFTLength',...
                'ReferenceRangeCentered'});
            
            scopeSection = matlab.system.display.Section(...
                'Title','Scope Settings',...
                'PropertyList',{...
                'Name',...
                'Position',...
                'IQDataInput',...
                'RangeLabel',...
                'RangeResolution',...
                'RangeOffset',...
                'TimeResolution',...
                'TimeSpan',...
                'IntensityUnits'...
                });
                
           groupDisp = matlab.system.display.SectionGroup(...
                'Title','Main',...
                'Sections',scopeSection);   

           groups = [groupDisp,groupRange];
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