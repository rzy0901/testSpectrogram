classdef (Hidden) AbstractMovingBicyclist < phased.internal.AbstractSampleRateEngine & matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%MovingBicyclist   Moving bicyclist
%   H = phased.MovingBicyclist creates a bicyclist motion model, H.
%   This object simulates the motion of a bicycle.
%
%   H = phased.MovingBicyclist(Name,Value) returns a moving bicyclist,
%   H, with the specified property Name set to the specified Value. You
%   can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   The bicyclist is derived from a multi-scatterer model developed for a
%   77 GHz radar system. It is composed of 5 primary parts: bicycle frame
%   and upper-body of rider, pedals, legs, front wheel, and rear wheel.
%   Each part is composed of individual point scatterers. The total number
%   of point scatterers in the bicycle target model is dependent on the
%   number of spokes defined. The bicycle model is the aggregate of all of
%   the point scatterers.
%
%   Step method syntax:
%
%   [BBPOS,BBVEL,BBAX] = step(H,T) returns the bicyclist's scatterers'
%   position, velocity, and orientation axes at the current time in BBPOS,
%   BBVEL, and BBAX, respectively. The method then simulates the motion of
%   the bicyclist in the next duration, specified in T (in seconds). This
%   syntax applies when the InputSource property is set to the default
%   'Property'. The bicyclist heading, speed, and coasting can be updated
%   using the Heading, Speed, and Coasting properties, respectively.
%
%   BBPOS is a 3xN matrix with each column representing the position of a
%   point scatterer in the [x;y;z] format (in meters), where N is the
%   number of total point scatterers. BBVEL is also a 3xN matrix whose
%   columns are velocities of the corresponding point scatterers in [x;y;z]
%   form (in m/s). BBAX is a 3x3 matrix that is the orientation axes of the
%   corresponding bicyclist in the [x;y;z] format.
%
%   The Heading and Speed properties are tunable, but the Coast property is
%   nontunable. To change the coasting of the bicyclist target, the syntax
%   described next is recommended. 
%
%   [BBPOS,BBVEL,BBAX] = step(H,T,ANGH,SPEED,COAST) also specifies:
%   ANGH  - Heading angle (in degrees)
%   SPEED - Current bicyclist speed (m/s)
%   COAST - Current coasting state of bicyclist (logical)
%   This syntax applies when you set the InputSource property to
%   'Input port'.
%
%   System objects may be called directly like a function instead of
%   using the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   getNumScatterers method syntax:
%
%   N = getNumScatterers(H) returns the number of scatterers in the
%   bicyclist model.
%
%   N is a scalar output of the number of scatterers.
%
%   MovingBicyclist methods:
%
%   step             - Return scatterers' position, velocity, and axes
%   getNumScatterers - Return number of scatterers in the bicyclist model 
%   release          - Allow property value and input changes
%   clone            - Create a system object with same property values
%   isLocked         - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>           - Reset the moving bicyclist to its initial position
%
%   MovingBicyclist properties:
%
%   NumWheelSpokes         - Number of wheel spokes
%   GearTransmissionRatio  - Gear transmission ratio
%   InitialPosition        - Initial position (m)
%   InputSource            - Source of heading, speed, and coast inputs
%   Heading                - Heading direction (deg)
%   Speed                  - Bicyclist speed (m/s)
%   Coast                  - Coast bicyclist
%   InitialHeading         - Initial heading direction (deg)
%   InitialSpeed           - Initial bicyclist speed (m/s)
%
%   % Example 1:
%   %   Model a bicyclist that accelerates from 0 to 1.5 m/s at a rate
%   %   of 0.25 m/s^2. After the bicyclist reaches 1.5 m/s, the bicyclist
%   %   coasts. Modify behavior using inputs to the step method.
%   
%   % Define simulation
%   dt = 0.003;
%   M = 3000;
%   tSim = 0:dt:M*dt-dt;
%   acceleration = 0.25;
%   speed = 0;
%   coastingSpeed = 1.5;
%   coast = false;
%   
%   % Initialize bicyclist object
%   bicycle = phased.internal.MovingBicyclist(...
%       'InputSource','Input port',...
%       'InitialSpeed',speed,...
%       'InitialHeading',0);
%   
%   % Build movement trajectory
%   for m = 1:M
%       if speed < coastingSpeed
%           % When the speed is less than the coasting speed, continue to
%           % accelerate
%           speed = acceleration*tSim(m);
%       else
%           % When the speed is equal to the coasting speed, coast the
%           % bicyclist
%           coast = true;
%       end
%       ppos = step(bicycle,dt,0,speed,coast);
%       xAxisLims = [min(ppos(1,:)) max(ppos(1,:))];
%       
%       % Plot
%       if m == 1
%           ph = scatter3(...
%               ppos(1,:),...
%               ppos(2,:),...
%               ppos(3,:),...
%               20,'filled');
%           phAx = gca;
%           axis equal;
%           xlim(xAxisLims)
%           ylim([-2 2]);
%           zlim([-0.02 1.5]);
%           xlabel('X (m)');
%           zlabel('Z (m)');
%           title('Bicyclist Trajectory');
%           view(0,0); % x-z view
%       else
%           ph.XData = ppos(1,:);
%           ph.YData = ppos(2,:);
%           ph.ZData = ppos(3,:);
%           phAx.XLim = xAxisLims;
%       end
%       drawnow limitrate;
%   end
% 
%   % Example 2:
%   %   Model the motion of a bicyclist riding in a quarter circle. Modify
%   %   heading using the Heading property.
%   
%   % Define simulation
%   dt = 0.003;
%   M = 3000;
%   angstep = 90/M;
%   
%   % Initialize bicyclist object
%   bicycle = phased.internal.MovingBicyclist(...
%       'InputSource','Property');
%   
%   % Build movement trajectory
%   for m = 1:M
%       bicycle.Heading = angstep*m;
%       ppos = step(bicycle,dt);
%       
%       % Plot
%       if m == 1
%           ph = scatter3(ppos(1,:),...
%               ppos(2,:),...
%               ppos(3,:),...
%               20,'filled');
%           axis equal;
%           xlim([-1,25]);
%           ylim([-1 25]);
%           zlim([0 2]);
%           xlabel('X (m)');
%           ylabel('Y (m)');
%           title('Bicyclist Trajectory');
%           view(-40,20);
%       else
%           ph.XData = ppos(1,:);
%           ph.YData = ppos(2,:);
%           ph.ZData = ppos(3,:);
%       end
%       drawnow limitrate;
%   end
%
%   See also phased, backscatterPedestrian, phased.Platform.

%   Copyright 2019 The MathWorks, Inc.

%   Reference
%   [1] Stolz, M. et al. "Multi-Target Reflection Point Model of Cyclists
%   for Automotive Radar." 2017 European Radar Conference (EURAD),
%   Nuremberg, 2017, pp. 94-97.
%
%   [2] Chen, V., D. Tahmoush, and W. J. Miceli. Radar Micro-Doppler
%   Signatures: Processing and Applications. The Institution of Engineering
%   and Technology: London, 2014.


%   Copyright 2019 The MathWorks, Inc.

%   Reference
%   [1] Stolz, M. et al. "Multi-Target Reflection Point Model of Cyclists
%   for Automotive Radar." 2017 European Radar Conference (EURAD),
%   Nuremberg, 2017, pp. 94-97.
%
%   [2] Chen, V., D. Tahmoush, and W. J. Miceli. Radar Micro-Doppler
%   Signatures: Processing and Applications. The Institution of Engineering
%   and Technology: London, 2014.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %NumWheelSpokes         Number of wheel spokes
        %   Specify an integer number of wheel spokes for the bicycle. The
        %   default is set to 20 spokes.
        % 
        %   The model limits the number of spokes to a maximum of 50. The
        %   minimum number of spokes is 3.
        NumWheelSpokes = 20
        %GearTransmissionRatio  Gear transmission ratio
        %   The gear transmission ratio determines the number of wheel
        %   rotations to pedal rotations. This ratio should be in the range
        %   of 0.5 - 6. The default is set to 1.5.
        GearTransmissionRatio = 1.5
        %InitialPosition        Initial position (m)
        %   Specify the initial position of the bicyclist as a 3x1 vector
        %   in the form of [x; y; z] (in meters). The default value of this
        %   property is [0; 0; 0].
        InitialPosition = [0;0;0]
    end
    
    properties (Nontunable)
        %InputSource           Source of heading, speed, and coast values
        % Specify the inputs source as one of 'Property' | 'Input port',
        % where the default is 'Property'. When you set the InputSource
        % property to 'Property', the Heading, Speed, and Coast properties
        % are used. If InputSource is set to 'Input port', initial values
        % are used based on the 'InitialHeading' and 'InitialSpeed' values.
        % Successive values for the heading, speed, and coast are specified
        % at the input.
        InputSource = 'Property'
    end
    
    properties
        %Heading                Heading direction (deg)
        %   Specify the heading (in degrees) of the bicyclist as a scalar.
        %   The heading is measured from the x-axis towards the y-axis in
        %   the xy-plane. The default value is 0.
        Heading(1,1) double = 0
        %Speed                  Bicyclist speed (m/s)
        %   Specify the speed (in m/s) of the bicyclist as a nonnegative
        %   scalar. The default value is 4.
        %
        %   The motion model limits the speed to 60 m/s.
        Speed(1,1) double = 4
    end
    
    properties (Nontunable)
        %Coast                  Coast bicyclist
        %   A logical that controls the coasting of the bicyclist. If set
        %   to true, the bicyclist does not pedal. If set to false, the
        %   bicyclist pedals. Defaults to false.
        Coast(1,1) logical = false
    end
    
    
    properties (Nontunable)
        %InitialHeading         Initial heading direction (deg)
        %   Specify the initial heading (in degrees) of the bicyclist as a
        %   scalar. The heading is measured from x-axis towards y-axis in
        %   the xy-plane. The default value is 0.
        InitialHeading(1,1) double = 0
    end
    
    properties 
        %InitialSpeed           Initial bicyclist speed (m/s)
        %   Specify the initial speed (in m/s) of the bicyclist as a
        %   nonnegative scalar. The default value is 4.
        %
        %   The motion model limits the speed to 60 m/s.
        InitialSpeed(1,1) double = 4
    end

    properties(Constant, Hidden)
        InputSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    properties (Access=protected, Nontunable)
        pNumSamplesPerCycle
        pCycleSpeed
    end
    
    properties (Constant, Hidden)
        pPointScattererSpacing = 0.1
        pWheelDiameter = 0.75
        pPedalDiameter = 0.2
        pLengthUpperLeg = 0.4
        pLengthLowerLeg = 0.4
    end
    
    properties (Access = protected)
        pFrameAndRiderUpperBodyInitialPosition
        pPedalsInitialPosition
        pLegsInitialPosition
        pFrontWheelInitialPosition
        pRearWheelInitialPosition
        pWheelInitialPosition
        pUpperLegVec
        pLowerLegVec

        pFrameAndRiderUpperBodyPosition
        pPedalsTrajectory
        pLegsTrajectory
        pFrontWheelPosition
        pRearWheelPosition
        pWheelPosition
        pMounts
        pNumScatterers
        
        pPreviousDt = 0 
        pCurrentTime
        pCurrentPedalIdx = 1
        pHeading
        pSpeed
        pInitialSpeed
        pCoast = false
        pBicyclistInitialPosition
        pBicyclistPreviousPosition
        pPosition
        pOrientationAxes
        pBetaSample
    end

    methods            
        function set.NumWheelSpokes(obj,val)
            sigdatatypes.checkFiniteNonNegIntScalar(obj,'NumWheelSpokes', val)
            obj.NumWheelSpokes = val;
        end
        function set.GearTransmissionRatio(obj,val)
            sigdatatypes.checkFinitePosDblScalar(obj,'GearTransmissionRatio', val)
            obj.GearTransmissionRatio = val;
        end
        function set.InitialPosition(obj,val)
            sigdatatypes.validate3DCartCoord(val,...
                'phased.internal.AbstractMovingBicyclist','InitialPosition',{'size',[3 1]});
            obj.InitialPosition = val;
        end
        function set.Heading(obj,val)
            sigdatatypes.validateAngle(val,...
                'phased.internal.AbstractMovingBicyclist','Heading',{'scalar'});
            obj.Heading = val;
        end
        function set.Speed(obj,val)
            sigdatatypes.validateSpeed(val,'phased.internal.AbstractMovingBicyclist','Speed',{'scalar','real'});
            obj.Speed = val;
        end
        function set.Coast(obj,val)
            obj.Coast = val;
        end
        function set.InitialHeading(obj,val)
            sigdatatypes.validateAngle(val,...
                'phased.internal.AbstractMovingBicyclist','InitialHeading',{'scalar'});
            obj.InitialHeading = val;
        end
        function set.InitialSpeed(obj,val)
            sigdatatypes.validateSpeed(val,'phased.internal.AbstractMovingBicyclist','InitialSpeed',{'scalar','real'});
            obj.InitialSpeed = val;
        end
    end
    
    methods (Access=protected)
        function obj = AbstractMovingBicyclist(varargin)
            setProperties(obj, nargin, varargin{:});
            
            % Perform initialization here for performance and convenience 
            setupBicyclistScattererPoints(obj)
        end

        function validatePropertiesImpl(obj)
            % Number of wheel spokes
            cond = (obj.NumWheelSpokes < 3) ||(obj.NumWheelSpokes > 50); % Wheel spokes range = 3 - 50
            if cond
                coder.internal.errorIf(cond,'phased:target:invalidBicyclistNumWheelSpokes','NumWheelSpokes');
            end
            
            % Gear transmission ratio
            cond = (obj.GearTransmissionRatio < 0.5) || (obj.GearTransmissionRatio > 6); % Gear transmission ratio range = 0.5 - 6
            if cond
                coder.internal.errorIf(cond,'phased:target:invalidGearTransmissionRatio','GearTransmissionRatio');
            end
        end
        
        function validateBicyclistSpeed(obj,varargin)
            % Bicyclist speed
            cond = (obj.Speed > 60); % Don't let the bicyclist go faster than 60 m/s ~ 134 mph
            if cond
                coder.internal.errorIf(cond,'phased:target:invalidBicyclistSpeed','property','Speed');
            end
        end
        
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and system block dialog
            inputPortChk = strcmpi(prop,'InitialSpeed') || ...
                strcmpi(prop,'InitialHeading');
            propChk = strcmpi(prop,'Speed') || ...
                strcmpi(prop,'Heading') || ...
                strcmpi(prop,'Coast');
            flag = false; 
            if inputPortChk
                if strcmp(obj.InputSource,'Property')
                    flag = true;
                end
            elseif propChk
                if strcmp(obj.InputSource,'Input port')
                    flag = true;
                end
            end
        end
        
        function flag = isInputSizeMutableImpl(obj,index) %#ok<INUSD>
            % Return false if input size cannot change
            % between calls to the System object
            flag = false;
        end
        
        function flag = isInputComplexityMutableImpl(obj,index) %#ok<INUSD>
            % Return false if input complexity cannot change
            % between calls to the System object
            flag = false;
        end
        
        function flag = isInputDataTypeMutableImpl(obj,index) %#ok<INUSD>
            % Return false if input data type cannot change
            % between calls to the System object
            flag = false;
        end
        
        function [out,out2,out3] = getOutputSizeImpl(obj)
            % Return size for each output port
            nPtsUpper = ceil(obj.pLengthUpperLeg/obj.pPointScattererSpacing);
            nPtsLower = ceil(obj.pLengthLowerLeg/obj.pPointScattererSpacing)-1;
            upperLegVec = (obj.pLengthUpperLeg/nPtsUpper).*(1:nPtsUpper);
            lowerLegVec = (obj.pLengthLowerLeg/(nPtsLower+1)).*(1:nPtsLower);
            [frameAndRiderUpperBody, pedals, legs, ...
                frontWheel, rearWheel, ~] = ...
                phased.internal.BackscatterBicyclist.getBicyclistParts(obj.pWheelDiameter,...
                obj.NumWheelSpokes,obj.pPointScattererSpacing,...
                upperLegVec,lowerLegVec);
            BicyclistScatterers = [frameAndRiderUpperBody pedals legs frontWheel rearWheel];
            N = size(BicyclistScatterers,2);
            out = [3 N];
            out2 = [3 N];
            out3 = [3 3];
        end
        
        function [out,out2,out3] = getOutputDataTypeImpl(obj) %#ok<MANU>
            % Return data type for each output port
            out = "double";
            out2 = "double";
            out3 = "double";
        end
        
        function [out,out2,out3] = isOutputComplexImpl(obj) %#ok<MANU>
            % Return true for each output port with complex data
            out = false;
            out2 = false;
            out3 = false;
        end
        
        function [out,out2,out3] = isOutputFixedSizeImpl(obj) %#ok<MANU>
            % Return true for each output port with fixed size
            out = true;
            out2 = true;
            out3 = true;
        end
        
        function setupImpl(obj)
            % Setup
            setupBicylistMotionCharacteristics(obj)
            
            % Simulation Parameters
            obj.pCurrentTime = 0;
            obj.pPreviousDt = 0; 
            nspercyc = obj.pNumSamplesPerCycle;
            
            % Set positions 
            obj.pFrameAndRiderUpperBodyPosition = obj.pFrameAndRiderUpperBodyInitialPosition;
            obj.pFrontWheelPosition = obj.pFrontWheelInitialPosition;
            obj.pRearWheelPosition = obj.pRearWheelInitialPosition;
            obj.pWheelPosition = obj.pWheelInitialPosition;
            
            % Initialize trajectory of pedals and legs
            obj.pPedalsTrajectory = zeros(3,size(obj.pPedalsInitialPosition,2),nspercyc);
            obj.pLegsTrajectory = zeros(3,size(obj.pLegsInitialPosition,2),nspercyc);
            obj.pPedalsTrajectory(:,:,1) = obj.pPedalsInitialPosition;
            obj.pLegsTrajectory(:,:,1) = obj.pLegsInitialPosition;
            
            % Compute general trajectory
            loadTrajectory(obj);
            
            % Assign initial values
            if strcmp(obj.InputSource,'Property')
                obj.pHeading = obj.Heading;
                obj.pSpeed = obj.Speed;
                obj.pCoast = obj.Coast;
                obj.pInitialSpeed = obj.InitialSpeed; % Placeholder
            else
                obj.pHeading = obj.InitialHeading;
                obj.pSpeed = obj.InitialSpeed;
                obj.pInitialSpeed = obj.InitialSpeed; % Used to check for changes to a "nontunable" property; this is done to support codegen
            end
        end
        
        function resetImpl(obj)
            obj.pCurrentTime = 0;
            obj.pPreviousDt = 0;
            obj.pCurrentPedalIdx = 1;
            obj.pPosition = obj.InitialPosition;
            if strcmp(obj.InputSource,'Property')
                obj.pHeading = obj.Heading;
                obj.pSpeed = obj.Speed;
                obj.pCoast = obj.Coast;
                obj.pInitialSpeed = obj.InitialSpeed; % Placeholder
            else
                obj.pHeading = obj.InitialHeading;
                obj.pSpeed = obj.InitialSpeed;
                obj.pInitialSpeed = obj.InitialSpeed; % Used to check for changes to a "nontunable" property; this is done to support codegen
            end
            obj.pOrientationAxes = rotz(obj.pHeading);
            obj.pBicyclistPreviousPosition = obj.pBicyclistInitialPosition;     
            obj.pFrameAndRiderUpperBodyPosition = obj.pFrameAndRiderUpperBodyInitialPosition;
            obj.pFrontWheelPosition = obj.pFrontWheelInitialPosition;
            obj.pRearWheelPosition = obj.pRearWheelInitialPosition;
            obj.pWheelPosition = obj.pWheelInitialPosition;
        end
        
        function loadObjectImpl(obj,s,~)
            % Set properties in object obj to values in structure s
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function [bicyclistpos,bicyclistvel,bicyclistaxes] = pedaling(obj,dt,varargin)
            moveMethodNontunablePropertiesCheck(obj); % Check to make sure speed has not changed

            % Pedal
            pos = obj.pPosition;
            ax = obj.pOrientationAxes;
            if obj.pCurrentTime == 0
                % This is time == 0 
                rtsamps = 0;
                ptidx = obj.pCurrentPedalIdx;
                ptidxBounds = [1 1];
                bicyclistpos = obj.getBicyclistPosition(rtsamps,ptidx,ptidxBounds,ax,pos);
                bicyclistvel = zeros(size(bicyclistpos));
            else
                % This is some time > 0 
                % Wheel index
                Wc = obj.pCycleSpeed/obj.pSpeed; % sec/cycle 
                rt = mod(obj.pPreviousDt,Wc)/Wc;
                rtsamps = rt*obj.pNumSamplesPerCycle;
                
                % Pedal index
                if obj.pCoast
                    ptidx = obj.pCurrentPedalIdx;
                else
                    ptidx = max(mod(rtsamps+obj.pCurrentPedalIdx,obj.pNumSamplesPerCycle),1);
                end
                ptidxBounds = [floor(ptidx) min(ceil(ptidx),obj.pNumSamplesPerCycle)];
                
                % Set indices 
                obj.pCurrentPedalIdx = ptidx;
                
                % Calculate velocity
                bicyclistpos_old = obj.pBicyclistPreviousPosition;
                bicyclistpos = obj.getBicyclistPosition(rtsamps,ptidx,ptidxBounds,ax,pos);
                bicyclistvel = (bicyclistpos-bicyclistpos_old)/obj.pPreviousDt;
            end
            
            % Set previous position for future use
            obj.pPreviousDt = dt; 
            obj.pBicyclistPreviousPosition = bicyclistpos;
            
            % Set axes for output
            bicyclistaxes = ax;
            
            % Assign inputs to private properties
            if strcmp(obj.InputSource,'Property')
                obj.pHeading = obj.Heading;
                obj.pSpeed = obj.Speed;
                obj.pCoast = obj.Coast;
            else
                obj.pHeading = varargin{1};
                obj.pSpeed = varargin{2};
                obj.pCoast = logical(varargin{3});
            end
        end
        
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj
            
            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            
            if isLocked(obj)
                s.pNumSamplesPerCycle = obj.pNumSamplesPerCycle;
                s.pCycleSpeed = obj.pCycleSpeed;

                s.pFrameAndRiderUpperBodyInitialPosition = obj.pFrameAndRiderUpperBodyInitialPosition;
                s.pPedalsInitialPosition = obj.pPedalsInitialPosition;
                s.pLegsInitialPosition = obj.pLegsInitialPosition;
                s.pFrontWheelInitialPosition = obj.pFrontWheelInitialPosition;
                s.pRearWheelInitialPosition = obj.pRearWheelInitialPosition;
                s.pWheelInitialPosition = obj.pWheelInitialPosition; 
                s.pUpperLegVec = obj.pUpperLegVec;
                s.pLowerLegVec = obj.pLowerLegVec; 
                
                s.pFrameAndRiderUpperBodyPosition = obj.pFrameAndRiderUpperBodyPosition;
                s.pPedalsTrajectory = obj.pPedalsTrajectory;
                s.pLegsTrajectory = obj.pLegsTrajectory;
                s.pFrontWheelPosition = obj.pFrontWheelPosition;
                s.pRearWheelPosition = obj.pRearWheelPosition;
                s.pWheelPosition = obj.pWheelPosition; 
                s.pMounts = obj.pMounts;
                s.pNumScatterers = obj.pNumScatterers; 
                
                s.pPreviousDt = obj.pPreviousDt; 
                s.pCurrentTime = obj.pCurrentTime;
                s.pCurrentPedalIdx = obj.pCurrentPedalIdx;
                s.pHeading = obj.pHeading;
                s.pSpeed = obj.pSpeed;
                s.pInitialSpeed = obj.pInitialSpeed; 
                s.pCoast = obj.pCoast;
                s.pBicyclistInitialPosition = obj.pBicyclistInitialPosition;
                s.pBicyclistPreviousPosition = obj.pBicyclistPreviousPosition;
                s.pPosition = obj.pPosition;
                s.pOrientationAxes = obj.pOrientationAxes;
                s.pBetaSample = obj.pBetaSample; 
            end
        end
    end
    
    methods (Access=protected)
        function setupBicyclistScattererPoints(obj)
            nPtsUpper = ceil(obj.pLengthUpperLeg/obj.pPointScattererSpacing);
            nPtsLower = ceil(obj.pLengthLowerLeg/obj.pPointScattererSpacing)-1;
            obj.pUpperLegVec = (obj.pLengthUpperLeg/nPtsUpper).*(1:nPtsUpper);
            obj.pLowerLegVec = (obj.pLengthLowerLeg/(nPtsLower+1)).*(1:nPtsLower);
            
            % Create bicycle 
            [frameAndRiderUpperBody, pedals, legs, ...
                frontWheel, rearWheel, mounts] = ...
                phased.internal.BackscatterBicyclist.getBicyclistParts(obj.pWheelDiameter,...
                obj.NumWheelSpokes, obj.pPointScattererSpacing,...
                obj.pUpperLegVec, obj.pLowerLegVec);
            
            % Set initial positions
            obj.pFrameAndRiderUpperBodyInitialPosition = frameAndRiderUpperBody;
            obj.pPedalsInitialPosition = pedals;
            obj.pLegsInitialPosition = legs;
            obj.pFrontWheelInitialPosition = frontWheel;
            obj.pRearWheelInitialPosition = rearWheel;
            obj.pWheelInitialPosition = bsxfun(@minus, obj.pFrontWheelInitialPosition, mounts.FrontWheel);
            obj.pBicyclistInitialPosition = [frameAndRiderUpperBody pedals legs frontWheel rearWheel];
            obj.pMounts = mounts; 
            obj.pNumScatterers = size(obj.pBicyclistInitialPosition,2);
        end
        
        function setupBicylistMotionCharacteristics(obj)
            obj.pNumSamplesPerCycle = 1e5; 
            obj.pCycleSpeed = 2*pi*(obj.pWheelDiameter/2)*obj.GearTransmissionRatio;
            dt = 1/obj.pNumSamplesPerCycle;
            obj.pBetaSample = 2*pi*obj.GearTransmissionRatio*dt;
        end  
    end
    
    methods(Hidden)
        function loadTrajectory(obj)
            % % Code to save angles of legs 
            % % Get rotational matrix for pedals
            % dt = 1/obj.pNumSamplesPerCycle;
            % rotmatPedal = obj.getRotationalMatrix(2*pi*dt);
            % 
            % % Rotate pedals
            % nspercyc = obj.pNumSamplesPerCycle;
            % pedalsTrajectory = obj.pPedalsTrajectory; 
            % pedalsTrajectory(:,:,1) = pedalsTrajectory(:,:,1) - obj.pMounts.Pedals; 
            % for nt = 2:nspercyc
            %     pedalsTrajectory(:,:,nt) = rotmatPedal*pedalsTrajectory(:,:,nt-1);
            % end
            % 
            % % Translate pedals
            % pedalsTrajectory = pedalsTrajectory + obj.pMounts.Pedals; 
            % 
            % % Update feet mounting positions
            % mountsFeetMat = [pedalsTrajectory(:,obj.pMounts.IdxFeet(1),2:end) ...
            %     pedalsTrajectory(:,obj.pMounts.IdxFeet(2),2:end)];
            % 
            % % Initialize angles and scatterers of legs 
            % angUpperL = zeros(1,nspercyc-1); 
            % angLowerL = zeros(1,nspercyc-1);
            % angUpperR = zeros(1,nspercyc-1); 
            % angLowerR = zeros(1,nspercyc-1);
            % 
            % % Calculate angle of legs 
            % mounts = obj.pMounts; 
            % spacing = obj.pPointScattererSpacing;
            % for nt = 1:nspercyc-1
            %     mounts.Feet = mountsFeetMat(:,:,nt);
            %     
            %     [angUpperL(:,nt),angLowerL(:,nt),...
            %         angUpperR(:,nt),angLowerR(:,nt)] = ...
            %         phased.internal.BackscatterBicyclist.getBicyclistLegAngles(spacing,mounts,obj.pUpperLegVec,obj.pLowerLegVec);
            % end
            % 
            % % Save angles
            % saveDir = fullfile(matlabroot,'toolbox','phased','phased');
            % save(fullfile(saveDir,'BicyclistLegsAngles.mat'),'angUpperL','angLowerL','angUpperR','angLowerR');
            
            % Get rotational matrix for pedals
            dt = 1/obj.pNumSamplesPerCycle;
            rotmatPedal = obj.getRotationalMatrix(2*pi*dt);
            
            % Rotate pedals
            nspercyc = obj.pNumSamplesPerCycle;
            obj.pPedalsTrajectory(:,:,1) = bsxfun(@minus, obj.pPedalsTrajectory(:,:,1), obj.pMounts.Pedals);
            for nt = 2:nspercyc
                obj.pPedalsTrajectory(:,:,nt) = rotmatPedal*obj.pPedalsTrajectory(:,:,nt-1);
            end
            
            % Translate pedals
            obj.pPedalsTrajectory = bsxfun(@plus, obj.pPedalsTrajectory, obj.pMounts.Pedals);
                
            % Update feet mounting positions
            mounts = obj.pMounts;
            spacing = obj.pPointScattererSpacing;
                
            % Move legs
            nPtsUpper = length(obj.pUpperLegVec);
            if coder.target('MATLAB')
                savedAngles = load('BicyclistLegsAngles.mat','angUpperL','angLowerL','angUpperR','angLowerR');     
                
                % Calculate leg position
                [leftLeg, rightLeg] = phased.internal.BackscatterBicyclist.getBicyclistLegs(...
                        obj.pUpperLegVec,obj.pLowerLegVec,...
                        savedAngles.angUpperL,...
                        savedAngles.angLowerL,...
                        savedAngles.angUpperR,...
                        savedAngles.angLowerR);
            else
                % Initialize 
                mountsFeetMat = [obj.pPedalsTrajectory(:,obj.pMounts.IdxFeet(1),2:end) ...
                    obj.pPedalsTrajectory(:,obj.pMounts.IdxFeet(2),2:end)];
                angUpperL = zeros(1,nspercyc-1);
                angLowerL = zeros(1,nspercyc-1);
                angUpperR = zeros(1,nspercyc-1);
                angLowerR = zeros(1,nspercyc-1);
                
                for nt = 1:nspercyc-1
                    mounts.Feet = mountsFeetMat(:,:,nt);
                    
                    % Calculate angle of legs 
                    [angUpperL(:,nt),angLowerL(:,nt),...
                        angUpperR(:,nt),angLowerR(:,nt)] = ...
                        phased.internal.BackscatterBicyclist.getBicyclistLegAngles(...
                        spacing,mounts,obj.pUpperLegVec,obj.pLowerLegVec);
                end
                
                % Calculate leg position
                [leftLeg, rightLeg] = ...
                    phased.internal.BackscatterBicyclist.getBicyclistLegs(...
                    obj.pUpperLegVec,obj.pLowerLegVec,...
                    angUpperL,angLowerL,angUpperR,angLowerR);
            end
                            
            % Shift legs into position about hips
            leftLeg(:,1:nPtsUpper,:) = bsxfun(@plus, leftLeg(:,1:nPtsUpper,:), mounts.Hip(:,1));
            leftLeg(:,(nPtsUpper+1):end,:) = bsxfun(@plus, leftLeg(:,(nPtsUpper+1):end,:), leftLeg(:,nPtsUpper,:));
            rightLeg(:,1:nPtsUpper,:) = bsxfun(@plus, rightLeg(:,1:nPtsUpper,:), mounts.Hip(:,2));
            rightLeg(:,(nPtsUpper+1):end,:) = bsxfun(@plus, rightLeg(:,(nPtsUpper+1):end,:), rightLeg(:,nPtsUpper,:));
            
            % Assign legs to object 
            obj.pLegsTrajectory(:,:,2:nspercyc) = [leftLeg rightLeg];
        end
        
        function bicyclistpos = getBicyclistPosition(obj,rtsamps,ptidx,ptidxBounds,ax,pos)
            % Calculate rotation matrices for wheels 
            beta = rtsamps*obj.pBetaSample; 
            rotmat = obj.getRotationalMatrix(beta);
            
            % Rotate wheels
            obj.pWheelPosition = rotmat*obj.pWheelPosition;
            obj.pFrontWheelPosition = bsxfun(@plus, obj.pWheelPosition, obj.pMounts.FrontWheel);
            obj.pRearWheelPosition = bsxfun(@plus, obj.pWheelPosition, obj.pMounts.RearWheel);
            
            % manual interp1 to avoid overhead
            alpha = ptidx-ptidxBounds(1);
            traj = [obj.pPedalsTrajectory(:,:,ptidxBounds) obj.pLegsTrajectory(:,:,ptidxBounds)];
            legPedealInterpTraj = (1-alpha)*traj(:,:,1)+...
                alpha*traj(:,:,2);
                
            % Get current bicycle position and orient as appropriate 
            bicyclistpos = [...
                obj.pFrameAndRiderUpperBodyPosition, ...
                legPedealInterpTraj,...
                obj.pFrontWheelPosition, ...
                obj.pRearWheelPosition];
            bicyclistpos = bsxfun(@plus, ax*bicyclistpos, pos);
        end
        
        function moveMethodNontunablePropertiesCheck(obj)
            cond = obj.InitialSpeed ~= obj.pInitialSpeed;
            if cond
                coder.internal.errorIf(cond,'phased:target:nontunableProperty','InitialSpeed','MovingBicyclist');
            end
        end
    end
    
    methods (Hidden,Static)
        function rotmat = getRotationalMatrix(beta)   
            % Calculate rotation matrices for all instances
            % For a single instance, the rotation matrix is
            % rotmatWheel = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
            cosBeta = cos(beta);
            sinBeta = sin(beta);
            rotmat = [cosBeta 0 sinBeta; 0 1 0; -sinBeta 0 cosBeta];   
        end
    end
end
