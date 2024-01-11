classdef (Sealed, StrictDefaults) AngleDopplerScope < ...
        phased.internal.AbstractNarrowbandArrayProcessing 
% AngleDopplerScope Visualize angle-Doppler patterns
%   SCOPE = phased.AngleDopplerScope returns a System object, SCOPE, that
%   displays angle-Doppler response patterns.
%
%   SCOPE = phased.AngleDopplerScope ('Name', Value,) returns an
%   angle-Doppler scope System object, SCOPE, with each specified property
%   name set to the specified value. You can specify name-value pair
%   arguments in any order as (Name 1, Value 1, ..., NameN, ValueN).
%
%   Step method syntax when IQDataInput property is false:
%
%   step(SCOPE,X,ANG,DOP) displays the angle-Doppler response pattern, X. X
%   is a PxQ matrix containing the complex angle-Doppler response, where P
%   is the number of Doppler samples and Q is the number of angle samples
%   respectively.
%
%   ANG is a length-Q column vector containing the angle values at which
%   the angle-Doppler response is evaluated. DOP is a length-P column
%   vector containing the Doppler values at which the angle-Doppler
%   response is evaluated.
%
%   Step method syntax when IQDataInput property is true:
%
%   step(SCOPE,X) calculates and displays the angle-Doppler response of
%   input data X. X must be a matrix or a column vector. If X is a matrix,
%   the number of rows must equal the number of subarrays if SensorArray
%   contains subarrays, or the number of elements, otherwise. If X is a
%   vector, the number of rows must be an integer multiple of the number of
%   elements or subarrays of the array specified in the SensorArray
%   property. In addition, the multiplier must be at least 2.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(SCOPE,x) and y = SCOPE(x) are
%   equivalent.
%
%   AngleDopplerScope methods:
%
%   step         - Display angle-Doppler response
%   release      - Allow property value and input characteristics changes
%   clone        - Create angle-Doppler scope object with same property
%                  values
%   isLocked     - Locked status (logical)
%   show         - Turn on visibility of the scope
%   hide         - Turn off visibility of the scope
%   isVisible    - Return visibility of the scope (logical)
%
%   AngleDopplerScope properties:
%
%   Name                 - Scope window name
%   Position             - Scope window position
%   IQDataInput          - Type of input
%   ResponseUnits        - Output response unit
%   AngleLabel           - Angle axis label
%   DopplerLabel         - Doppler axis label
%   SensorArray          - Sensor array
%   PropagationSpeed     - Signal propagation speed
%   OperatingFrequency   - Operating frequency
%   PRF                  - Pulse repetition frequency
%   ElevationAngle       - Elevation angle
%   NumAngleSamples      - Number angular bins
%   NumDopplerSamples    - Number Doppler bins
%   NormalizeDoppler     - Normalize Doppler grid
%   FrequencyUnits       - Doppler frequency unit
%
%   % Example:
%   %   Calculate and visualize the angle Doppler response of the 190th
%   %   cell of a collected data cube.
%
%   load STAPExampleData;
%
%   % 190th cell of data cube
%   x = shiftdim(STAPEx_ReceivePulse(190,:,:));
%
%   % Set up scope to process and visualize the angle-Doppler response
%   response = phased.AngleDopplerScope(...
%               'IQDataInput', true,...
%               'Name','Angle-Doppler Scope',...
%               'Position', [560 375 560 420],...
%               'NormalizeDoppler', false,...
%               'ResponseUnits',  'db',...
%               'SensorArray',STAPEx_HArray,...
%               'OperatingFrequency',STAPEx_OperatingFrequency,...
%               'PropagationSpeed',STAPEx_PropagationSpeed,...
%               'PRF',STAPEx_PRF);
%   response(x);
%   
%   See also phased, phased.RangeDopplerScope, phased.RangeAngleScope,
%   phased.AngleDopplerResponse.

%   Copyright 2018 The MathWorks, Inc. 

     properties (Nontunable)
        %PRF     Pulse repetition frequency (Hz)
        %   Specify the pulse repetition frequency (PRF) (in Hz) of the
        %   input signal as a positive scalar.  This property applies when
        %   you set the 'IQDataInput' property true. The default value of
        %   this property is 1.
        PRF = 1;
        
        %ElevationAngle     Elevation angle (deg)
        %   Specify the elevation angle (in degrees) used to calculate the
        %   angle-Doppler response as a scalar. The angle must be between
        %   -90 and 90.   This property applies when you set the
        %   'IQDataInput' property true. The default value of this property
        %   is 0.
        ElevationAngle = 0;
        
        %NumAngleSamples    Number of angle bins
        %   Specify the number of samples in the angular domain used to
        %   calculate the angle-Doppler response as a positive integer.
        %   This value must be greater than 2.  This property applies when
        %   you set the 'IQDataInput' property true. The default value of
        %   this property is 256.
        NumAngleSamples = 256;
        
        %NumDopplerSamples    Number of Doppler bins
        %   Specify the number of samples in the Doppler domain used to
        %   calculate the angle-Doppler response as a positive integer.
        %   This value must be greater than 2.  This property applies when
        %   you set the 'IQDataInput' property true. The default value of
        %   this property is 256.
        NumDopplerSamples = 256;
    end
    
    properties(Nontunable,Logical)
        %IQDataInput Specify the type of input
        %   Specify whether the input is I/Q (raw) data or processed data.
        %   If you set this property true, it implies that raw data is
        %   passed and processing is done along angle and Doppler domain
        %   before visualizing. Set this property false if you have
        %   processed data. The default value is true.
        IQDataInput = true;
        
        %NormalizeDoppler Normalize Doppler
        %   Set the NormalizeDoppler property to true if you want the
        %   Doppler frequency to be normalized. Set this to false to plot
        %   the angle-Doppler response without normalizing the Doppler
        %   frequency.  This property applies when you set the
        %   'IQDataInput' property true. The default value is false.
        NormalizeDoppler = false;
    end

    properties (Nontunable)
        %ResponseUnits Unit of the response
        %   Specify the unit of the plot, using one of | 'db' | 'magnitude'
        %   | 'power'|. The default is 'db'.
        ResponseUnits = 'db';
        
        %FrequencyUnits  Doppler frequency unit
        %   Specify the Doppler frequency unit, using one of | 'Hz' | 'kHz'
        %   | 'MHz'|. This property applies only when 'IQDataInput'
        %   property is true and NormalizedDoppler property is false. The
        %   default is 'Hz'.
        FrequencyUnits = 'Hz';
    end

    properties
        %AngleLabel Angle-axis label
        %   Specify the angle-axis label as a string. This property applies
        %   when you set the 'IQDataInput' property false. The default
        %   value is 'Angle (degrees)'. This property is tunable.
        AngleLabel = 'Angle (degrees)';
        
        %DopplerLabel Doppler-axis label
        %   Specify the Doppler-axis label as a string. This property
        %   applies when you set the 'IQDataInput' property false. The
        %   default value is 'Doppler Frequency (Hz)'. This property is
        %   tunable.
        DopplerLabel = 'Doppler Frequency (Hz)';
        
        %Name Caption to display on scope window
        %   Specify the caption to display on the scope window as any
        %   string. The default value is 'Angle Doppler Scope'. This
        %   property is tunable.
        Name = 'Angle Doppler Scope';

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
        pFrequencyConversion
    end
    
    properties (Constant, Hidden) 
        ResponseUnitsSet = matlab.system.StringSet({'db','power','magnitude'}); 
        FrequencyUnitsSet = matlab.system.StringSet({'Hz','kHz','MHz'});
    end
    
    methods
        % Constructor
        function obj = AngleDopplerScope(varargin)
            % Support name-value pair arguments when constructing
            % AngleDopplerScope object
                     
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
            
            % matlabshared.scopes.MatrixViewer - Visualization
            obj.cScope = matlabshared.scopes.MatrixViewer(...
                'XDataMode','Custom',...
                'YDataMode','Custom',...
                'AxisOrigin','Lower left corner',...
                'Colormap','parula',...
                'Title',getString(message('phased:scopes:angdopLabel'))...
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
        function set.Position(obj, value)
            if ~matlabshared.scopes.Validator.Position(value)
                error(message('Spcuilib:scopes:InvalidPosition'));
            end
            obj.Position = value;
        end       
        function set.PRF(obj,val)
            sigdatatypes.validateFrequency(val,'phased.AngleDopplerResponse',...
                'PRF',{'scalar'});
            obj.PRF = val;
        end
        function set.ElevationAngle(obj,val)
            sigdatatypes.validateAngle(val,'phased.AngleDopplerResponse',...
                'ElevationAngle',{'scalar','>=',-90,'<=',90});
            obj.ElevationAngle = val;
        end
        function set.NumAngleSamples(obj,val)
            sigdatatypes.validateIndex(val,'phased.AngleDopplerResponse',...
                'NumAngleSamples',{'scalar','>=',2});
            obj.NumAngleSamples = val;
        end
        function set.NumDopplerSamples(obj,val)
            sigdatatypes.validateIndex(val,'phased.AngleDopplerResponse',...
                'NumDopplerSamples',{'scalar','>=',2});
            obj.NumDopplerSamples = val;
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
        function set.DopplerLabel(obj,value)
             validateattributes(value, {'char', 'string'}, ...
                {'nonsparse'}, '', 'DopplerLabel');
            obj.DopplerLabel = value;
        end
    end
    
    
    methods(Access = protected)
        
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and System block dialog
            if obj.IQDataInput
                flag = false;             
                if   strcmp(prop,'AngleLabel') || ...
                        strcmp(prop,'DopplerLabel')
                    flag = true;
                end
                
                if obj.NormalizeDoppler && strcmp(prop,'FrequencyUnits')
                    flag = true;
                end
            else
                flag = true;
                
                if strcmp(prop,'IQDataInput') || ...
                        strcmp(prop,'ResponseUnits') || ...
                        strcmp(prop,'Name') || ...
                        strcmp(prop,'Position') || ...
                        strcmp(prop,'AngleLabel') || ...
                        strcmp(prop,'DopplerLabel')
                    flag = false;
                end
            end
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

            if obj.IQDataInput
                % Check input to AngleDopplerResponse
                sz_x = size(x);
                N = getDOF(obj.SensorArray);
                if sz_x(2) == 1
                    temp = sz_x(1)/N;
                    cond = rem(temp,1) || temp < 2;
                    if cond
                        coder.internal.errorIf(cond,...
                            'phased:AngleDopplerResponse:InvalidWeightDimension', N, 2*N);
                    end
                else
                    cond = sz_x(1) ~= N;
                    if cond
                        coder.internal.errorIf(cond,...
                            'phased:AngleDopplerResponse:InvalidDataDimension', N);
                    end
                end
            else
                % Grid values check
                ang_grid = varargin{2};
                dop_grid = varargin{3};
                
                cond = ~iscolumn(ang_grid) || isempty(ang_grid);
                if cond
                   coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeColVector','ANG'); 
                end
                
                cond = ~iscolumn(dop_grid) || isempty(dop_grid);
                if cond
                   coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeColVector','DOP'); 
                end
                                
                cond = (sz_x(1) ~= numel(ang_grid));
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions',...
                        'ANG',sz_x(1),numel(ang_grid));
                end
                
                cond = (sz_x(2) ~= numel(dop_grid));
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDimensions',...
                        'DOP',sz_x(2),numel(dop_grid));
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
                
                obj.cResponse = phased.AngleDopplerResponse(...
                    'SensorArray',obj.SensorArray,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'OperatingFrequency',obj.OperatingFrequency,...
                    'PRFSource','Property',...
                    'PRF',obj.PRF,...
                    'ElevationAngleSource','Property',...
                    'ElevationAngle',obj.ElevationAngle,...
                    'NumAngleSamples',obj.NumAngleSamples,...
                    'NumDopplerSamples',obj.NumDopplerSamples);
                              
                % Set the YLabel
                if ~obj.NormalizeDoppler
                    switch obj.FrequencyUnits
                        case 'Hz'  % Hz
                            obj.pFrequencyConversion = 1;
                            obj.cScope.YLabel = getString(message('phased:scopes:HzLabel'));
                        case 'kHz' % Hz to kHz
                            obj.pFrequencyConversion = unitsratio('km','m');
                            obj.cScope.YLabel = getString(message('phased:scopes:kHzLabel'));
                        case 'MHz' % Hz to MHz
                            obj.pFrequencyConversion = unitsratio('km','mm');
                            obj.cScope.YLabel = getString(message('phased:scopes:MHzLabel'));
                    end
                else
                    obj.cScope.YLabel = getString(message('phased:scopes:NormFreqLabel'));
                end
                
                obj.cScope.XLabel = getString(message('phased:scopes:AngLabel'));
            else
                obj.cScope.XLabel = obj.AngleLabel;
                obj.cScope.YLabel = obj.DopplerLabel;
            end

            % Set data cursor labels
            setCursorDataLabels(obj.cScope,["Angle","Doppler","Intensity"]);
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

            if obj.IQDataInput             % 'IQDataInput' is true
                
                % Process raw data using phased.AngleDopplerResponse
                [adresp_out,ang_grid,dop_grid] = ...
                    obj.cResponse(varargin{:});
                
                % NormalizeDoppler
                if obj.NormalizeDoppler
                    dop_grid = dop_grid/obj.PRF;
                end
                
                if ~obj.NormalizeDoppler
                    dop_grid = dop_grid.*obj.pFrequencyConversion;
                end
                
            else                            % 'IQDataInput' is false
                
                % Get the response/ angle grid/ Doppler grid
                adresp_out = varargin{1};
                ang_grid = varargin{2};
                dop_grid = varargin{3};
            end
            
            % Visualize
            plotMatrixData(obj,adresp_out,ang_grid,dop_grid,unit{:});
        end
        
        function plotMatrixData(obj,resp,ang_grid,...
                dop_grid,varargin)

            unit = varargin{:};
            
            % Matrix Viewer
            scope = obj.cScope;             
            
            % Compute response based on 'ResponseUnits'
            response = phased.internal.computePlotPattern(...
                abs(resp),false,unit);
            
            scope.CustomXData = ang_grid;
            scope.CustomYData = dop_grid;
            
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
                obj.cScope.XLabel = obj.AngleLabel;
                obj.cScope.YLabel = obj.DopplerLabel;
            end         
        end
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            if isLocked(obj)
                s.cScope = saveobj(obj.cScope);
                s.cResponse = saveobj(obj.cResponse);
                s.pFrequencyConversion = obj.pFrequencyConversion;
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
            s = loadSubObjects@phased.internal.AbstractNarrowbandArrayProcessing(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cResponse = phased.AngleDopplerResponse.loadobj(s.cResponse);
                    obj.cScope = matlabshared.scopes.MatrixViewer.loadobj(s.cScope);
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
            num = 1;
            if ~obj.IQDataInput
                num = num+2;
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
            groups = ...
                getPropertyGroupsImpl@phased.internal.AbstractNarrowbandArrayProcessing('subarray');
            
            groupScope = matlab.system.display.Section(...
                'Title','Scope Settings',...
                'PropertyList',{...
                'Name',...
                'Position',...
                'IQDataInput',... 
                'ResponseUnits'...
                'AngleLabel',...
                'DopplerLabel'...
                });
            
            props = [groups(1).PropertyList,...
                {...
                'PRF',...
                'ElevationAngle',...
                'NumAngleSamples',...      
                'NumDopplerSamples',...
                'NormalizeDoppler',...
                'FrequencyUnits'}];
                       
            groupProcess = matlab.system.display.Section(...         
                'PropertyList',props);
            
            processGroup = matlab.system.display.SectionGroup(...
                'Title','Processing Settings', ...
                'Sections',groupProcess);
          
            responseGroup = matlab.system.display.SectionGroup(...
                'Title','Main', ...
                'Sections',groupScope);
            
            sensorGroup = groups(2);
            
            groups = [responseGroup processGroup sensorGroup];        
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