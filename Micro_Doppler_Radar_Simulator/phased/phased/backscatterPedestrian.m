classdef backscatterPedestrian < handle
%backscatterPedestrian   Backscatter pedestrian model for radar
%   H = backscatterPedestrian creates a backscatter pedestrian, H,
%   for radar simulation. This object simulates the signal reflected off
%   the pedestrian while the pedestrian is in motion.
%
%   H = backscatterPedestrian(Name,Value) returns a backscatter
%   pedestrian model, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   move method syntax:
%
%   [BPPOS,BPVEL,BPAX] = move(H,T,ANGH) returns the pedestrian's body
%   segments position, velocity, and orientation axes at the current time
%   in BPPOS, BPVEL, and BPAX, respectively. It then simulates the walking
%   motion in the next duration, specified in T (in seconds). ANGH
%   specifies the current heading angle (in degrees) as a scalar. The angle
%   is in the xy-plane.
%
%   The pedestrian model includes 16 body segments: left and right feet,
%   left and right lower legs, left and right upper legs, left and right
%   hip, left and right lower arms, left and right upper arms, left and
%   right shoulders, neck, and head. Therefore BPPOS is a 3x16 matrix with
%   each column representing the position of the corresponding body
%   segments in the [x;y;z] format (in meters). BPVEL is also a 3x16 matrix
%   whose columns are velocities of corresponding body segments in the
%   [x;y;z] form (in m/s). BPAX is a 3x3x16 array whose pages are
%   orientation axes of the corresponding body segments. The three columns
%   represents the 3 axes and each column is in [x;y;z] format.
%
%   reflect method syntax:
%
%   Y = reflect(H,X,ANG) returns the reflected signal Y off the pedestrian
%   target due to the input signal X.
%
%   The human model consists of 16 body segments: left and right feet, left
%   and right lower legs, left and right upper legs, left and right hip,
%   left and right lower arms, left and right upper arms, left and right
%   shoulders, neck, and head. Each body segment is represented by a
%   cylinder. Therefore, X is a 16-column matrix whose columns are incident
%   signals to each body segment.
% 
%   ANG is a 2x16 matrix representing the signal's incident direction to
%   each body segment. Each column of ANG specifies the incident direction
%   of the corresponding signal in the form of an [AzimuthAngle;
%   ElevationAngle] pair (in degrees).
%
%   Y is a column vector containing the combined reflected signal from all
%   body segments. The number of rows in Y is the same as the number of
%   rows in X.
%
%   plot method syntax:
%
%   HF = plot(H) returns a plot of the position of the backscatter
%   pedestrian object H at the current time.
%
%   HF is the figure handle.
%
%   backscatterPedestrian methods:
%
%   move     - Pedestrian motion (see above)
%   reflect  - Signal reflection off pedestrian (see above)
%   plot     - Plot pedestrian's position (see above) 
%   release  - Allow property value and input characteristics changes
%   clone    - Create pedestrian target object with same property values
%   reset    - Reset internal states of the pedestrian target
%
%   backscatterPedestrian properties:
%
%   Height              - Pedestrian height
%   WalkingSpeed        - Pedestrian walking speed 
%   PropagationSpeed    - Propagation speed
%   OperatingFrequency  - Operating frequency 
%   InitialPosition     - Initial position 
%   InitialHeading      - Initial heading 
%
%   % Example:
%   %   Compute the reflected signal from a pedestrian moving along
%   %   x axis. Assume the radar works at 24 GHz and the signal has
%   %   a 300 MHz bandwidth. The signal is captured at the moment 
%   %   the pedestrian starts to move and 1 second into the
%   %   movement. Assume the radar is at the origin.
%
%   c = 3e8; bw = 3e8; fs = bw; fc = 24e9;
%   wav = phased.LinearFMWaveform('SampleRate',fs,'SweepBandwidth',bw);
%   x = wav();
%
%   chan = phased.FreeSpace('OperatingFrequency',fc,'SampleRate',fs,...
%          'TwoWayPropagation',true);
%   ped = backscatterPedestrian(...
%          'OperatingFrequency',fc,'InitialPosition',[100;0;0]);
%   rpos = [0;0;0];
%
%   % time 0
%   [bppos,bpvel,bpax] = move(ped,1,0);
%   plot(ped); 
%   xp = chan(repmat(x,1,16),rpos,bppos,[0;0;0],bpvel);
%   [~,ang] = rangeangle(rpos,bppos,bpax);
%   y0 = reflect(ped,xp,ang);
%
%   % time 1
%   [bppos,bpvel,bpax] = move(ped,1,0);
%   plot(ped); 
%   xp = chan(repmat(x,1,16),rpos,bppos,[0;0;0],bpvel);
%   [~,ang] = rangeangle(rpos,bppos,bpax);
%   y1 = reflect(ped,xp,ang);
%
%   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(wav));
%   ymf = mf([y0 y1]);
%   t = (0:size(ymf,1)-1)/fs;
%   figure
%   plot(t,abs(ymf));
%   xlabel('Time (s)'); ylabel('Magnitude'); title('Pedestrian Return')
%
%   See also phased, phased.BackscatterRadarTarget, phased.Platform.
    
%   Copyright 2018-2019 The MathWorks, Inc.

%   Reference
%   [1] Victor Chen, The Micro-Doppler Effect in Radar, Artech House, 2011
%   
%   [2] Boulic, Ronan, et al. A Global Human Walking Model with Real-time
%   Kinematic Personification, The Visual Computer: International Journal
%   of Computer Graphics, Vol. 6, Issue 6, Dec 1990

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties 
        %Height     Pedestrian height (m)
        %   Specify the height (in meters) of the pedestrian as a positive
        %   scalar. The default value is 1.65.
        Height(1,1) double {mustBePositive,mustBeFinite,mustBeNonempty}= 1.65
        %WalkingSpeed   Pedestrian walking speed (m/s)
        %   Specify the walking speed (in m/s) of the pedestrian as a
        %   nonnegative scalar. The default value is 1.4.
        %
        %   The motion model limits the walking speed to about 1.4 times 
        %   pedestrian's height.
        WalkingSpeed(1,1) double {mustBePositive,mustBeFinite,mustBeNonempty} = 1.4
        %PropagationSpeed Propagation speed (m/s)
        %   Specify the wave propagation speed (in m/s) in free space as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed(1,1) double {mustBeNonnegative,mustBeFinite,mustBeNonempty} = physconst('lightspeed')
        %OperatingFrequency Signal carrier frequency (Hz)
        %   Specify the carrier frequency (in Hz) of the narrowband signal
        %   as a scalar. The default value of this property is 3e8 (300
        %   MHz).
        OperatingFrequency = 300e6
        %InitialPosition    Initial position (m)
        %   Specify the initial position of the pedestrian as a 3x1 vector
        %   in the form of [x; y; z] (in meters). The default value of this
        %   property is [0; 0; 0].
        InitialPosition = [0;0;0]
        %InitialHeading    Initial heading (deg)
        %   Specify the initial heading (in degrees) of the pedestrian as a
        %   scalar. The heading is measured from x-axis towards y-axis in
        %   the xy-plane. The default value is 0.
        InitialHeading = 0
    end
    
    properties (Access=private)
        cWalkingHuman
        cHumanTarget
        pIsInitialized = false
        pIsMotionInitialized = false
        
        pBodyPartsJointsInitialPosition
        pBodyPartsPosition
        pBodyPartsVelocity
        pBodyPartsOrientationAxes
        pBodyPartsJointsPosition
        
        pFigInitialized = false
        pLineObjInitialized = false
        pHeadObjInitialized = false
        pFigHandle
        pLineObjHandle
        pHeadObjHandle
        pFigAxes
        pFigAxesNV
    end
    
    properties (Constant)
        pMaxHeight = 3 
        pMinHeight = 0.3 
    end
    
    methods
        function set.Height(obj,value)
            sigdatatypes.validateDistance(value,'','Height',{'double'},{'scalar','positive'});
            
            % Check heights are within reasonable limit
            cond = (value > obj.pMaxHeight) ||(value < obj.pMinHeight); % Height: 0.3 m to 3 m (~1ft to ~10 ft)
            if cond
                coder.internal.errorIf(cond,...
                    'phased:target:valueOutsideExpectedRange','Height',...
                    sprintf('%0.1f',obj.pMinHeight),sprintf('%0.1f',obj.pMaxHeight),'m');
            end
            obj.Height = value;
        end
        function set.OperatingFrequency(obj,value)
            className = class(obj);
            errStr = 'backscatterPedestrian';
            if strcmp(className(1),'p')
                errStr = 'BackscatterPedestrian';
            end
            validateattributes(value,{'double'},{'scalar','finite',...
                'positive'},errStr,'OperatingFrequency');
            obj.OperatingFrequency = value;
        end
        function set.PropagationSpeed(obj,value)
            className = class(obj);
            errStr = 'backscatterPedestrian';
            if strcmp(className(1),'p')
                errStr = 'BackscatterPedestrian';
            end
            sigdatatypes.validateSpeed(value,errStr,...
                'PropagationSpeed',{'scalar','positive'});
            obj.PropagationSpeed = value;
        end
        function set.InitialPosition(obj,val)
            className = class(obj);
            errStr = 'backscatterPedestrian';
            if strcmp(className(1),'p')
                errStr = 'BackscatterPedestrian';
            end
            sigdatatypes.validate3DCartCoord(val,...
                errStr,'InitialPosition',{'size',[3 1]});
            obj.InitialPosition = val;
        end
        function set.InitialHeading(obj,val)
            className = class(obj);
            errStr = 'backscatterPedestrian';
            if strcmp(className(1),'p')
                errStr = 'BackscatterPedestrian';
            end
            sigdatatypes.validateAngle(val,...
                errStr,'InitialHeading',{'scalar'});
            obj.InitialHeading = val;
        end
    end
    
    methods 
        function obj = backscatterPedestrian(varargin)
            defaultHeight = 1.65;
            defaultWalkingSpeed = 1.4;
            defaultPropagationSpeed = physconst('lightspeed');
            defaultOperatingFrequency = 3e8;
            defaultInitialPosition = zeros(3,1);
            defaultInitialHeading = 0;
            
            if coder.target('MATLAB')
                p = inputParser;
                p.CaseSensitive = true;
                p.FunctionName = 'WalkingHumanTarget';
                p.PartialMatching = false;
                addParameter(p,'Height',defaultHeight);
                addParameter(p,'WalkingSpeed',defaultWalkingSpeed);
                addParameter(p,'PropagationSpeed',defaultPropagationSpeed);
                addParameter(p,'OperatingFrequency',defaultOperatingFrequency);
                addParameter(p,'InitialPosition',defaultInitialPosition);
                addParameter(p,'InitialHeading',defaultInitialHeading);
                parse(p,varargin{:});

                fn = fieldnames(p.Results);
                for m = 1:numel(fn)
                    obj.(fn{m}) = p.Results.(fn{m});
                end
            else
                 parms = struct('Height',uint32(0), ...
                                'WalkingSpeed',uint32(0), ...
                                'PropagationSpeed',uint32(0), ...
                                'OperatingFrequency',uint32(0), ...
                                'InitialPosition',uint32(0), ...
                                'InitialHeading',uint32(0));
                 pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
                 obj.Height = eml_get_parameter_value(pstruct.Height,defaultHeight,varargin{:});
                 obj.WalkingSpeed = eml_get_parameter_value(pstruct.WalkingSpeed,defaultWalkingSpeed,varargin{:});
                 obj.PropagationSpeed = eml_get_parameter_value(pstruct.PropagationSpeed,defaultPropagationSpeed,varargin{:});
                 obj.OperatingFrequency = eml_get_parameter_value(pstruct.OperatingFrequency,defaultOperatingFrequency,varargin{:});
                 obj.InitialPosition = eml_get_parameter_value(pstruct.InitialPosition,defaultInitialPosition,varargin{:});
                 obj.InitialHeading = eml_get_parameter_value(pstruct.InitialHeading,defaultInitialHeading,varargin{:});
            end
            [obj.pBodyPartsJointsInitialPosition,~] = getHumanBodyParts(obj.Height);
        end
        
        function release(obj)
        %release    Release backscatterPedestrian object.
            if obj.pIsInitialized
                release(obj.cWalkingHuman);
                release(obj.cHumanTarget);
                obj.pIsInitialized = false;
            end
            obj.pIsMotionInitialized = false;
        end
        
        function reset(obj)
        %release    Reset backscatterPedestrian object.
            reset(obj.cWalkingHuman);
            reset(obj.cHumanTarget);
        end
        
        function [bodypartspos,bodypartsvel,bodypartsax] = move(obj,dt,heading)
        %move   Pedestrian motion
        %   [BPPOS,BPVEL,BPAX] = move(H,T,ANGH) returns the pedestrian's
        %   body segments position, velocity, and orientation axes at the
        %   current time in BPPOS, BPVEL, and BPAX, respectively. It then
        %   simulates the walking motion in the next duration, specified in
        %   T (in seconds). ANGH specifies the current heading angle (in
        %   degrees) as a scalar. The angle is in the xy-plane.
        %
        %   The pedestrian model includes 16 body segments: left and right
        %   feet, left and right lower legs, left and right upper legs,
        %   left and right hip, left and right lower arms, left and right
        %   upper arms, left and right shoulders, neck, and head. Therefore
        %   BPPOS is a 3x16 matrix with each column represents the position
        %   of the corresponding body segments in the [x;y;z] format (in
        %   meters). BPVEL is also a 3x16 matrix whose columns are
        %   velocities of corresponding body segments in the [x;y;z] form
        %   (in m/s). BPAX is a 3x3x16 array whose pages are orientation
        %   axes of the corresponding body segments. The three columns
        %   represents the 3 axes and each column is in [x;y;z] format.
        %
        %   % Example:
        %   %   Model a quarter-circle trajectory of a pedestrian.
        %
        %   ped = backscatterPedestrian;
        %   dt = 0.003;
        %   N = 3000;
        %   angstep = 90/N;
        %   ppos = zeros(3,16,N);
        %   for m = 1:N
        %       ppos(:,:,m) = move(ped,dt,angstep*m);
        %   end
        %   
        %   plotxpos = permute(ppos(1,1:2,:),[3 2 1]);
        %   plotypos = permute(ppos(2,1:2,:),[3 2 1]);
        %   plotzpos = permute(ppos(3,1:2,:),[3 2 1]);
        %   plot3(plotxpos,plotypos,plotzpos);
        %   title('Pedestrian Trajectory');
        %   xlabel('X (m)'); ylabel('Y (m)');
        %   legend('Left foot','Right foot');
        %   view(-8,76);
        %
        %   See also phased, backscatterPedestrian/reflect.
        
            % for the sake of performance, rely on input validation in
            % WalkingHuman
            
            if ~obj.pIsInitialized
                callSetup(obj);
            end
            
            % Check for nontunable properties 
            moveMethodNontunablePropertiesCheck(obj)
            
            [bodypartspos,bodypartsvel,bodypartsax,bodypartsjpos] = step(obj.cWalkingHuman,dt,heading);
            obj.pBodyPartsPosition = bodypartspos;
            obj.pBodyPartsVelocity = bodypartsvel;
            obj.pBodyPartsOrientationAxes = bodypartsax;
            obj.pBodyPartsJointsPosition = bodypartsjpos;
            
            obj.pIsMotionInitialized = true;
        end
        
        function y = reflect(obj,x,ang)
        %reflect    Signal reflection off pedestrian
        %   Y = reflect(H,X,ANG) returns the reflected signal Y off the
        %   pedestrian target due to the input signal X.
        %
        %   The human model consists of 16 body segments: left and right
        %   feet, left and right lower legs, left and right upper legs,
        %   left and right hip, left and right lower arms, left and right
        %   upper arms, left and right shoulders, neck, and head. Each body
        %   segment is represented by a cylinder. Therefore, X is a
        %   16-column matrix whose columns are incident signals to each
        %   body segment.
        % 
        %   ANG is a 2x16 matrix representing the signal's incident
        %   direction to each body segment. Each column of ANG specifies
        %   the incident direction of the corresponding signal in the form
        %   of an [AzimuthAngle; ElevationAngle] pair (in degrees).
        %
        %   Y is a column vector containing the combined reflected signal
        %   from all body segments. The number of rows in Y is the same as
        %   the number of rows in X.
        %
        %   % Example:
        %   %   Compute the reflected signal from a pedestrian moving along
        %   %   x axis. Assume the radar works at 24 GHz and the signal has
        %   %   a 300 MHz bandwidth. The signal is captured at the moment 
        %   %   the pedestrian starts to move and 1 second into the
        %   %   movement. Assume the radar is at the origin.
        %
        %   c = 3e8; bw = 3e8; fs = bw; fc = 24e9;
        %   wav = phased.LinearFMWaveform('SampleRate',fs,'SweepBandwidth',bw);
        %   x = wav();
        %
        %   chan = phased.FreeSpace('OperatingFrequency',fc,'SampleRate',fs,...
        %          'TwoWayPropagation',true);
        %   ped = backscatterPedestrian(...
        %          'OperatingFrequency',fc,'InitialPosition',[100;0;0]);
        %   rpos = [0;0;0];
        %
        %   % time 0
        %   [bppos,bpvel,bpax] = move(ped,1,0);
        %   xp = chan(repmat(x,1,16),rpos,bppos,[0;0;0],bpvel);
        %   [~,ang] = rangeangle(rpos,bppos,bpax);
        %   y0 = reflect(ped,xp,ang);
        %
        %   % time 1
        %   [bppos,bpvel,bpax] = move(ped,1,0);
        %   xp = chan(repmat(x,1,16),rpos,bppos,[0;0;0],bpvel);
        %   [~,ang] = rangeangle(rpos,bppos,bpax);
        %   y1 = reflect(ped,xp,ang);
        %
        %   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(wav));
        %   ymf = mf([y0 y1]);
        %   t = (0:size(ymf,1)-1)/fs;
        %   plot(t,abs(ymf));
        %   xlabel('Time (s)'); ylabel('Magnitude'); title('Pedestrian Return')
        %
        %   See also phased, backscatterPedestrian/move.
        
            % for the sake of performance, rely on input validation in
            % BackscatterHumanTarget
            
            % Check for nontunable properties 
            reflectMethodNontunablePropertiesCheck(obj)
            
            cond = ~obj.pIsMotionInitialized;
            if cond
                coder.internal.errorIf(cond,'phased:target:invalidCallSequence','move','reflect');
            end
            
            % Coupling with angle computation, could be an issue, when
            % multipath/curved path are supported.
            
            y = step(obj.cHumanTarget,x,ang);
        end

        function varargout = plot(obj,varargin)
        %plot    Plot pedestrian's position
        %   HF = plot(H) returns a plot of the position of the backscatter
        %   pedestrian object H at the current time.
        %
        %   HF is the figure handle.
        %
        %   HF = plot(...,'Parent',HAX) specifies the target axes HAX
        %   for the plotting of the pedestrian object H. 
        %
        %   % Example:
        %   %   Model a quarter-circle trajectory of a pedestrian.
        %
        %   ped = backscatterPedestrian;
        %   dt = 0.003;
        %   N = 3000;
        %   angstep = 90/N;
        %   for m = 1:N
        %       move(ped,dt,angstep*m);
        %       plot(ped); 
        %   end
        %
        %   See also phased, backscatterPedestrian/move.
        
        narginchk(1,3);
        nargoutchk(0,1);
        
        % Not support for codegen
        if ~isempty(coder.target)
            coder.internal.assert(false, ...
                'phased:Waveform:CodegenNotSupported','plot');
        end
        
        % Validate inputs to method
        parseAndValidatePlotInputs(obj,varargin{:});
        
        % Check to see if a axes handle has been provided
        if ~isempty(obj.pFigAxesNV)
            obj.pFigAxes = obj.pFigAxesNV;
            obj.pFigInitialized = true;
            obj.pFigHandle = obj.pFigAxes.Parent;
        end
        
        % Check to see if figure was deleted
        if ~ishandle(obj.pFigHandle)
            obj.pFigInitialized = false;
        end
        
        % Check to see if body was deleted
        if ~ishandle(obj.pLineObjHandle)
            obj.pLineObjInitialized = false;
        else
            if obj.pFigInitialized
                if obj.pLineObjHandle.Parent ~= obj.pFigAxes % Verify that axes have remained the same
                    obj.pLineObjInitialized = false;
                end
            end
        end
        
        % Check to see if head was deleted
        if ~ishandle(obj.pHeadObjHandle)
            obj.pHeadObjInitialized = false;
            delete(obj.pHeadObjHandle(ishandle(obj.pHeadObjHandle)));
        else
            if obj.pFigInitialized
                if obj.pHeadObjHandle.Parent ~= obj.pFigAxes % Verify that axes have remained the same
                    obj.pHeadObjInitialized = false;
                end
            end
        end
        
        % Check to see if motion is initialized
        cond = ~obj.pIsMotionInitialized;
        if cond
            bppos = obj.pBodyPartsJointsInitialPosition-...
                min(obj.pBodyPartsJointsInitialPosition(3,:)); % Use initial position
        else % Motion is initialized
            bppos = obj.pBodyPartsJointsPosition; % Use current position
        end
        
        % Setup indices for line segments
        idxBody = [2 4 6 8 17 16 15 16 14 12 10 12 14 16 13 11 9 11 13 ...
            16 17 7 5 3 1];
        idxHead = 15;
%         colorVal = [0 0.4470 0.7410];
        colorVal = [0 0 0];
        
        % Plot
        x = bppos(1,:);
        y = bppos(2,:);
        z = bppos(3,:);
        if obj.pFigInitialized
            if obj.pLineObjInitialized
                % Update body
                obj.pLineObjHandle.XData = x(idxBody);
                obj.pLineObjHandle.YData = y(idxBody);
                obj.pLineObjHandle.ZData = z(idxBody);
            else
                % Create body
                obj.pLineObjHandle = line(x(idxBody),y(idxBody),z(idxBody),...
                    'Color',colorVal,'LineStyle','-','LineWidth',2,...
                    'Marker','o','MarkerSize',4,'MarkerFaceColor',k);
                axis equal;
                grid on;
                view(36,15);
                hold on;
                obj.pLineObjInitialized = true;
                
                % Setup new figure axes in case of deletion of both
                % body and head
                obj.pFigAxes = gca;
                obj.pFigAxes.XLabel.String = 'X (m)';
                obj.pFigAxes.YLabel.String = 'Y (m)';
                obj.pFigAxes.ZLabel.String = 'Z (m)';
                obj.pFigAxes.Title.String = 'Pedestrian Trajectory';
                obj.pFigAxes.ZLim = [-0.1 obj.Height+0.1];
            end
            
            if obj.pHeadObjInitialized
                % Update head location
                obj.pHeadObjHandle.XData = x(idxHead); 
                obj.pHeadObjHandle.YData = y(idxHead); 
                obj.pHeadObjHandle.ZData = z(idxHead); 
            else
                % Create head
                obj.pHeadObjHandle = line(x(idxHead),y(idxHead),z(idxHead),...
                    'Color',colorVal,'LineStyle','none',...
                    'Marker','o','MarkerSize',15,'MarkerFaceColor',k);
                obj.pHeadObjInitialized = true;
            end
            
            % Update axes
            obj.pFigAxes.XLim = [min(x(idxHead))-1 max(x(idxHead))+1];
            obj.pFigAxes.YLim = [min(y(idxHead))-1 max(y(idxHead))+1];
            drawnow limitrate;
        else
            % Create figure
            obj.pFigHandle = figure('Name','Pedestrian Trajectory');
            obj.pFigInitialized = true;
            
            % Create body
            obj.pLineObjHandle = line(x(idxBody),y(idxBody),z(idxBody),...
                'Color',colorVal,'LineStyle','-','LineWidth',2,...
                'Marker','o','MarkerSize',4,'MarkerFaceColor',k);
            axis equal;
            grid on;
            view(36,15);
            obj.pLineObjInitialized = true;
            hold on;
            
            % Create head
            obj.pHeadObjHandle = line(x(idxHead),y(idxHead),z(idxHead),...
                'Color',colorVal,'LineStyle','none',...
                'Marker','o','MarkerSize',15,'MarkerFaceColor',k);
            obj.pHeadObjInitialized = true;
            
            % Update axes
            obj.pFigAxes = gca;
            obj.pFigAxes.XLim = [min(x(idxHead))-1 max(x(idxHead))+1];
            obj.pFigAxes.YLim = [min(y(idxHead))-1 max(y(idxHead))+1];
            obj.pFigAxes.ZLim = [-0.1 obj.Height+0.1];
            obj.pFigAxes.XLabel.String = 'X (m)';
            obj.pFigAxes.YLabel.String = 'Y (m)';
            obj.pFigAxes.ZLabel.String = 'Z (m)';
            obj.pFigAxes.Title.String = 'Pedestrian Trajectory';
            
            drawnow;
        end
        hFig = obj.pFigHandle;
        
        % Define output handle if requested
        if nargout == 1
            varargout{1} = hFig;
        end
    end
        
        function newobj = clone(obj)
        % clone Create object with same property values
        %    C = clone(OBJ) creates another instance of the object, OBJ,
        %    with the same property values.
            s = saveobj(obj);
            className = class(obj); 
            if strcmp(className(1),'p')
                newobj = phased.BackscatterPedestrian;
            else
                newobj = backscatterPedestrian;
            end
            loadPrivateProps(newobj,s);
        end
    end
    
    methods (Access=protected)
        function callSetup(obj)
            obj.setup();
            obj.pIsInitialized = true;
        end
        function setup(obj)
            validateProperties(obj);
            
            obj.cWalkingHuman = phased.internal.WalkingHuman(...
                'Height',obj.Height,'WalkingSpeed',obj.WalkingSpeed,...
                'InitialPosition',obj.InitialPosition,'HeadingSource','Input port',...
                'InitialHeading',obj.InitialHeading,...
                'OutputJointsPosition',true);
            obj.cHumanTarget = phased.internal.BackscatterHumanTarget(...
                'PropagationSpeed',obj.PropagationSpeed,'Height',obj.Height,...
                'OperatingFrequency',obj.OperatingFrequency);
        end
        
        function validateProperties(obj)
            H = obj.Height;
            Ht = 0.245*H+0.246*H;
            V = obj.WalkingSpeed;
            rv = V/Ht;
            cond = (rv > 3);
            if cond
                coder.internal.errorIf(cond,'phased:target:invalidWalkingSpeed',...
                    sprintf('%5.2f',3*Ht),sprintf('%5.2f',H));
            end
        end
        
        function moveMethodNontunablePropertiesCheck(obj)
            className = class(obj);
            errStr = 'backscatterPedestrian';
            if strcmp(className(1),'p')
                errStr = 'BackscatterPedestrian';
            end
            if obj.pIsInitialized
                cond = obj.Height ~= obj.cWalkingHuman.Height;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','Height',errStr);
                end
                cond = obj.WalkingSpeed ~= obj.cWalkingHuman.WalkingSpeed;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','WalkingSpeed',errStr);
                end
                cond = obj.InitialPosition ~= obj.cWalkingHuman.InitialPosition;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','InitialPosition',errStr);
                end
                cond = obj.InitialHeading ~= obj.cWalkingHuman.InitialHeading;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','InitialHeading',errStr);
                end
            end
        end
        
        function reflectMethodNontunablePropertiesCheck(obj)
            className = class(obj);
            errStr = 'backscatterPedestrian';
            if strcmp(className(1),'p')
                errStr = 'BackscatterPedestrian';
            end
            if obj.pIsInitialized
                cond = obj.PropagationSpeed~= obj.cHumanTarget.PropagationSpeed;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','PropagationSpeed',errStr);
                end
                cond = obj.OperatingFrequency ~= obj.cHumanTarget.OperatingFrequency;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','OperatingFrequency',errStr);
                end
            end
        end
        
        function s = saveobj(obj)
            s.Height = obj.Height;
            s.WalkingSpeed = obj.WalkingSpeed;
            s.PropagationSpeed = obj.PropagationSpeed;
            s.OperatingFrequency = obj.OperatingFrequency;
            s.InitialPosition = obj.InitialPosition;
            s.InitialHeading = obj.InitialHeading;
            
            s.cWalkingHuman = saveobj(obj.cWalkingHuman);
            s.cHumanTarget = saveobj(obj.cHumanTarget);
            s.pIsInitialized = obj.pIsInitialized;
            s.pIsMotionInitialized = obj.pIsMotionInitialized;
            s.pBodyPartsJointsInitialPosition = obj.pBodyPartsJointsInitialPosition;
            s.pBodyPartsPosition = obj.pBodyPartsPosition;
            s.pBodyPartsVelocity = obj.pBodyPartsVelocity;
            s.pBodyPartsOrientationAxes = obj.pBodyPartsOrientationAxes;
            s.pBodyPartsJointsPosition = obj.pBodyPartsJointsPosition;
            
            s.pFigInitialized = obj.pFigInitialized;
            s.pLineObjInitialized = obj.pLineObjInitialized;
            s.pHeadObjInitialized = obj.pHeadObjInitialized;
            s.pFigHandle = obj.pFigHandle; 
            s.pLineObjHandle = obj.pLineObjHandle;
            s.pHeadObjHandle = obj.pHeadObjHandle;
            s.pFigAxes = obj.pFigAxes;
            s.pFigAxesNV = obj.pFigAxesNV; 
        end
        
        function loadPrivateProps(obj,s)
            obj.cWalkingHuman = phased.internal.WalkingHuman.loadobj(s.cWalkingHuman);
            s = rmfield(s,'cWalkingHuman');
            obj.cHumanTarget = phased.internal.BackscatterHumanTarget.loadobj(s.cHumanTarget);
            s = rmfield(s,'cHumanTarget');
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
    end
    
    methods (Hidden)
        function parseAndValidatePlotInputs(obj,varargin)
            % Parse
            defaultAxes = [];
            if coder.target('MATLAB')
                % MATLAB parser
                pstruct = inputParser;
                pstruct.FunctionName = '';
                pstruct.addParameter('Parent',defaultAxes);
                pstruct.parse(varargin{:});
                
                % Get values
                obj.pFigAxesNV = pstruct.Results.Parent;
            else
                % Codegen parser
                poptions = struct( ...
                    'CaseSensitivity',false, ...
                    'PartialMatching','none', ...
                    'StructExpand',false, ...
                    'IgnoreNulls',true);
                
                parms = {'Parent'};
                pstruct = coder.internal.parseParameterInputs(parms,poptions,varargin{:});
                
                % Get values
                obj.pFigAxesNV = coder.internal.getParameterValue(pstruct.Parent,defaultAxes,varargin{:});
            end
            
            % Is this an axes handle?
            if ~isempty(obj.pFigAxesNV)
                if ~ishandle(obj.pFigAxesNV)
                    coder.internal.assert(false, ...
                        'phased:target:invalidAxesHandle','Parent');
                end
                
                % Is this a valid axes handle?
                if ~strcmpi(get(obj.pFigAxesNV,'Type'),'axes')
                    coder.internal.assert(false, ...
                        'phased:target:invalidAxesHandle','Parent');
                end
            end
        end
    end
    
    methods (Hidden,Static)
        % Codegen related
        function props = matlabCodegenNontunableProperties(~)
            props = {'Height','WalkingSpeed','PropagationSpeed',...
                'OperatingFrequency','InitialPosition','InitialHeading',...
                'cWalkingHuman','cHumanTarget'};
        end
        function names = matlabCodegenOnceNames()
            names = {'pIsInitialized', 'callSetup', 'callSetup', 'setup', 'setup'};
        end
        % load
        function obj = loadobj(s)
            obj = backscatterPedestrian;
            loadPrivateProps(obj,s);
        end
   end  
end