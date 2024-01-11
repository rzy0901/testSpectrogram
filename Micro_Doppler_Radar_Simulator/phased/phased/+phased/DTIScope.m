classdef (Sealed, StrictDefaults) DTIScope < matlab.System   
%DTIScope Display scrolling Doppler intensity vectors or arrays
%   DopplerScope = phased.DTIScope returns a System object, DopplerScope,
%   that can show a scrolling plot of Doppler intensity vectors versus
%   time.
%
%   DopplerScope = phased.DTIScope('Name', Value, ...) returns a DTIScope
%   System object, DopplerScope, with each specified property name set to
%   the specified value. You can specify name-value pair arguments in any
%   order as (Name 1, Value 1, ..., Name N, Value N).
%
%   Step method syntax when IQDataInput is false:
%
%   step(DopplerScope, X) displays the real signal, X. X is an N by M
%   matrix where N is the number of Doppler intensity bins in an intensity
%   vector and M is the number of intensity vectors. N should be greater
%   than 1 and M is equal or greater than 1. Each column of the matrix
%   represents a Doppler intensity vector from successive times. The time
%   between intensity vectors are specified in the TimeResolution property.
%   
%   Step method syntax when IQDataInput is true:
%
%   step(DopplerScope,X) Doppler processes the input raw data X to extract
%   the Doppler information and displays the processed signal.
%
%   X is an N by M matrix where N is the number of intensity bins in an
%   intensity vector and M is the number of intensity vectors. N should be
%   greater than 1 and M is equal or greater than 1. Each column of the
%   matrix represents an intensity vector from successive times. The time
%   between intensity vectors is specified in the TimeResolution property.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(DopplerScope, x) and y =
%   DopplerScope(x) are equivalent.
%
%   DTIScope methods:
%
%   step      - Display signal in the DTIScope figure 
%   release   - Allow property value and input characteristics changes, and
%               release DTIScope resources
%   clone     - Create DTIScope object with same property values
%   isLocked  - Display locked status (logical)
%   show      - Turn on visibility of DTIScope figure
%   hide      - Turn off visibility of DTIScope figure
%   isVisible - Return visibility of the scope (logical)
%
%   DTIScope properties:
%
%   Name                - Scope window name 
%   Position            - Scope window position
%   IQDataInput         - Type of input
%   DopplerResolution   - Doppler sample spacing 
%   DopplerOffset       - Doppler display offset 
%   TimeResolution      - Time resolution    
%   TimeSpan            - Time span 
%   IntensityUnits      - Intensity units
%   DopplerOutput       - Doppler output
%   PropagationSpeed    - Signal propagation speed
%   OperatingFrequency  - Operating frequency       
%   DopplerFFTLength    - FFT length in Doppler processing
%
%   % Example:
%   %   Use phased.DTIScope to visualize the Doppler information of the  
%   %   detection output of a complete radar system simulation. The radar
%   %   scenario contains a stationary single-element monostatic radar and 
%   %   three targets. In addition, one target is stationary relative to 
%   %   the radar and two targets are moving one approaching the radar at
%   %   approxiamtely 150 m/s and the other receding from the radar at 
%   %   approximately 150 m/s.
%   
%   % Load data
%   load('RTIDTIExampleData.mat')
%   rx_pulses = zeros(numel(fast_time),num_pulse_int);
% 
%   % DTIScope
%   dtiscope = phased.DTIScope('IQDataInput',false,...
%     'DopplerOutput','Speed',...
%     'PropagationSpeed',c,...
%     'OperatingFrequency',fc,...
%     'Name','Doppler-Time Display',...
%     'DopplerResolution',DopplerRes, ...
%     'DopplerOffset',-prf/2,...
%     'TimeResolution',0.05,...
%     'TimeSpan',5,...
%     'IntensityUnits','magnitude',...
%     'Position',[560 375 560 420]);
% 
%   pri = 1/prf;
%   nsteps = 200;
% 
%   % Transmit 200 pulses. Coherently process groups of 10 pulses at a
%   % time.
%   for k = 1:nsteps
%     
%       for m = 1:num_pulse_int
%         	[ant_pos,ant_vel] = radarplatform(pri);
%           [tgt_pos,tgt_vel] = targetplatforms(pri);
%        	sig = waveform();
%           [s,tx_status] = transmitter(sig);
%           [~,tgt_ang] = rangeangle(tgt_pos,ant_pos);
%           tsig = radiator(s,tgt_ang);
%           tsig = channels(tsig,ant_pos,tgt_pos,ant_vel,tgt_vel);
%           rsig = targets(tsig);
%           rsig = collector(rsig,tgt_ang);
%           rx_pulses(:,m) = preamplifier(rsig,~(tx_status>0));
%       end
%     
%       rx_pulses = gain(rx_pulses);
%       dshift = fft(rx_pulses.');
%       dshift = fftshift(abs(dshift),1);
%       dtiscope(mean(dshift,2));
%     
%       radarplatform(.05);
%       targetplatforms(.05);
%   end
%
%   See also phased, phased.RTIScope, phased.IntensityScope

%   Copyright 2018 The MathWorks, Inc.

    properties (Logical, Nontunable)
        %IQDataInput Specify the type of input
        %   Specify whether the input is I/Q (raw) data or processed data.
        %   If you set this property true, it implies that raw data is
        %   passed and Doppler processing is done before visualizing. Set
        %   this property false if you have processed data. The default
        %   value is false.
        IQDataInput = false
    end

    properties(Nontunable, PositiveInteger)
        %DopplerFFTLength     FFT length in Doppler processing
        %   Specify the FFT length in Doppler processing as a positive
        %   integer.  This Property applies when you set the 'IQDataInput'
        %   property true. The default value of this property is 1024.
        DopplerFFTLength = 1024
    end

    properties (Nontunable)
        %DopplerOutput  Doppler output
        %   Specify the Doppler domain as one of 'Frequency' | 'Speed',
        %   where the default is 'Frequency'. If you set this property to
        %   'Frequency', the Doppler domain is the Doppler shift (in Hz).
        %   If you set this property to 'Speed', the Doppler domain is the
        %   corresponding radial speed (in m/s).
        DopplerOutput = 'Frequency'

        %OperatingFrequency     Signal carrier frequency (Hz)
        %   Specify the signal carrier frequency (in Hz) as a scalar. This
        %   property applies when you set the DopplerOutput property to
        %   'Speed'. The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 3e8;

        %PropagationSpeed   Propagation speed (m/s)
        %   Specify the propagation speed (in m/s) of the signal as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('lightspeed');
    end
    
    properties(Nontunable)
        %DopplerResolution Doppler resolution (Hz)
        %   Specify the spacing between Doppler samples as a finite numeric
        %   scalar in Hz. The default value of this property is 1 Hz.
        DopplerResolution  = 1;
        
        %DopplerOffset Doppler offset (Hz)
        %   Specify the offset to apply to the Doppler-axis in Hz. This
        %   property is a numeric scalar with default value 0 Hz.
        DopplerOffset  = 0;
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
        %Name Caption to display on the scope window
        %   Specify the caption to display on the scope window as any
        %   string. The default value of this property is 'Doppler Time
        %   Intensity Scope'. This property is tunable.
        Name = 'Doppler Time Intensity Scope';
        
        %Position Scope window position in pixels
        %   Specify the size and location of the scope window in pixels, as
        %   a four-element double vector of the form: [left bottom width
        %   height]. The default value of obj property is dependent on the
        %   screen resolution, and is such that the window is positioned in
        %   the center of the screen, with a width and height of 800 and
        %   450 pixels respectively. This property is tunable.
        Position = [560 375 800 450];
    end

    properties (Access = private)
       cScope
       cFFT
    end

    properties (Constant, Hidden)
        DopplerOutputSet = matlab.system.StringSet({'Frequency','Speed'});
        IntensityUnitsSet = matlab.system.StringSet({'db','power',...
                'magnitude'});
    end

    methods
        function obj = DTIScope(varargin)
            setProperties(obj, nargin, varargin{:});
            
            % matlabshared.scopes.MatrixViewer - Visualization
            obj.cScope = matlabshared.scopes.MatrixViewer(...
                'XDataMode','Offset and resolution',...
                'YDataMode','Span and resolution',...
                'AxisOrigin','Lower left corner',...
                'ColorBarLocation','northoutside',...            
                'YLabel',getString(message('phased:scopes:timeLabel')),...
                'Colormap','parula',...
                'Title',getString(message('phased:scopes:dtiLabel')));                   
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
        
        function set.PropagationSpeed(obj,value)
            validateattributes( value, { 'double' },...
                { 'scalar', 'positive', 'finite' }, '', 'PropagationSpeed');
            obj.PropagationSpeed = value;
        end

        function set.OperatingFrequency(obj,value)
            validateattributes( value, { 'double' },...
                { 'scalar', 'positive', 'finite' }, '', 'OperatingFrequency');
            obj.OperatingFrequency = value;
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
    end
    
    methods (Access = protected)
        %% Common functions 
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;

            if ~obj.IQDataInput
                if strcmp(prop,'DopplerFFTLength')
                    flag = true;
                end
            end

            if strcmp(obj.DopplerOutput,'Frequency') && ...
                    strcmp(prop,'OperatingFrequency') 
                flag = true;
            end

            if strcmp(obj.DopplerOutput,'Frequency') && ...
                    strcmp(prop,'PropagationSpeed')
                flag = true;
            end

        end

        function setupImpl(obj, ~)
            
            % Doppler Shift in Hz to Speed
            if strcmp(obj.DopplerOutput,'Speed')
                XResolution = dop2speed(obj.DopplerResolution ,...
                    obj.PropagationSpeed/obj.OperatingFrequency)/2;
                XOffset = dop2speed(obj.DopplerOffset ,...
                    obj.PropagationSpeed/obj.OperatingFrequency)/2;
                XLabel = getString(message('phased:scopes:msLabel'));
            else
                XResolution = obj.DopplerResolution ;
                XOffset = obj.DopplerOffset ;
                XLabel = getString(message('phased:scopes:HzLabel'));
            end
            
            if obj.IQDataInput
                % FFT
                obj.cFFT = dsp.FFT('FFTLengthSource','Property',...
                    'FFTLength',obj.DopplerFFTLength);
            end
            
            % MatrixViewer
            scope = obj.cScope;
            
            % Set properties   
            scope.XResolution = XResolution;
            scope.XOffset = XOffset;
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
            
            scope.XLabel = XLabel;
            scope.Name = obj.Name;
            scope.Position = obj.Position;

            % Set data cursor labels
            setCursorDataLabels(scope,["Doppler","Time","Intensity"]);
        end

        function validateInputsImpl(~,x)
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
            
            sz_x = size(x);
            cond = sz_x(1) < 2;
            if cond
                coder.internal.errorIf(cond,...
                    'phased:scopes:InvalidInputDimension');
            end
        end
        
        function stepImpl(obj,data) 
            
            if obj.IQDataInput         % IQDataInput - True               
                % FFT
                out = obj.cFFT(data.');
                out = fftshift(abs(out),1);
                dop = mean(out,2);                       
            else                        % IQDataInput - false
                dop = abs(data);
            end  

            switch obj.IntensityUnits
                case 'magnitude'
                    data = dop;
                case 'power'
                    data = (dop).^2;
                case 'db'
                    data = mag2db(dop);
            end
            
            % Visualization (Launch Scope)
            obj.cScope(data.');  
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            reset(obj.cScope);
            reset(obj.cFFT);
        end

        function releaseImpl(obj)
            release(obj.cScope);
            if obj.IQDataInput
                release(obj.cFFT);
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
                s.cFFT = obj.cFFT;
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
                obj.cFFT = dsp.FFT.loadobj(s.cFFT);
                s = rmfield(s,'cScope');
                s = rmfield(s,'cFFT');
            end
        end
 
        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end
        
        function flag = isInputComplexityLockedImpl(obj,~)
            flag = false;  
            if obj.IQDataInput 
                flag = true;
            end
        end 

        function num = getNumInputsImpl(~)
            num = 1;
        end

        function num = getNumOutputsImpl(~)
            % Define total number of outputs for system with optional
            % outputs
            num = 0;
        end
    end
    
    methods(Access = protected, Static,Hidden)
 
        function group = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
            
            DopplerSection = matlab.system.display.Section(...
                'Title','Doppler Settings',...
                'PropertyList',{...
                'DopplerOutput',...
                'PropagationSpeed',...
                'OperatingFrequency',...
                'DopplerFFTLength'});
            
            ScopeSection = matlab.system.display.Section(...
                'Title','Scope Settings',...
                'PropertyList',{...
                'Name',...
                'Position',...
                'IQDataInput',...
                'DopplerResolution',...
                'DopplerOffset',...
                'TimeResolution',...
                'TimeSpan',...
                'IntensityUnits'...     
                });
            
            group = [ScopeSection,DopplerSection];
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

%[EOF]