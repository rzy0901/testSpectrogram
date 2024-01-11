classdef (Sealed) MovingBicyclist < phased.internal.AbstractMovingBicyclist
%This class is for internal use only. It may be removed in the future.

%MovingBicyclist   Moving bicyclist
%   H = phased.internal.MovingBicyclist creates a bicyclist motion model,
%   H. This object simulates the motion of a bicycle.
%
%   H = phased.internal.MovingBicyclist(Name,Value) returns a moving
%   bicyclist, H, with the specified property Name set to the specified
%   Value. You can specify additional name-value pair arguments in any
%   order as (Name1,Value1,...,NameN,ValueN).
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

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods 
        function obj = MovingBicyclist(varargin)
            obj@phased.internal.AbstractMovingBicyclist(varargin{:});
        end
        
        function N = getNumScatterers(obj)
            %getNumScatterers    Return number of scatterers in the model 
            %   N = getNumScatterers(H) returns the number of scatterers in
            %   the bicyclist model.
            %
            %   N is a scalar output of the number of scatterers.
            %
            %   % Example 1:
            %   %   Create a default bicyclist object. Determine the number
            %   %   of scatterers in the model.
            %
            %   % Initialize bicyclist object
            %   bicycle = backscatterBicyclist;
            %   N = bicycle.getNumScatterers
            %
            %   See also phased, backscatterBicyclist/move.
            
            % Reference
            %   [1] Stolz, M. et al. "Multi-Target Reflection Point Model
            %   of Cyclists for Automotive Radar." 2017 European Radar
            %   Conference (EURAD), Nuremberg, 2017, pp. 94-97.
            
            N = obj.pNumScatterers;
        end
    end
    
    methods (Access=protected)
        function num = getNumInputsImpl(obj)
            if strcmpi(obj.InputSource,'Input port')
                num = 4;
            else
                num = 1;
            end
        end
        
        function validateStepInputs(obj,dt,heading,speed,coast)
            % Validate inputs to the step method at initialization
            cond =  ~isscalar(dt);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeScalar','T');
            end
            cond =  ~isa(dt,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','T','double');
            end
            cond =  ~isreal(dt);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:Platform:MustBeReal', 'T');
            end
            if strcmp(obj.InputSource,'Input port')
                cond =  ~isscalar(heading);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeScalar','ANGH');
                end
                cond =  ~isa(heading,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','ANGH','double');
                end
                cond =  ~isreal(heading) || isnan(heading) || isinf(heading) || isempty(heading);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:Platform:MustBeReal', 'ANGH');
                end
                sigdatatypes.validateAngle(heading,...
                    '','ANGH',{'scalar'});
                
                % Check speed
                cond =  ~isscalar(speed);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeScalar','SPEED');
                end
                cond =  ~isa(speed,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','SPEED','double');
                end
                cond =  ~isreal(speed) || isnan(speed) || isinf(speed) ...
                    || isempty(speed) || (speed < 0) || (speed > 60); % Speed range = 0 - 60 m/s
                if cond
                    coder.internal.errorIf(cond,'phased:target:invalidBicyclistSpeed','input','SPEED');
                end
                
                % Check coast
                cond =  ~isscalar(coast);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeScalar','COAST');
                end
                cond =  ~(isa(coast,'numeric') || isa(coast,'logical'));
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','COAST',{'numeric','logical'});
                end
                cond =  ~isreal(coast) || isnan(coast) || ...
                    isinf(coast) || isempty(coast);
                if cond
                    coder.internal.errorIf(cond,'phased:target:invalidCoast','input','COAST');
                end
            end
        end

        function [scattererspos,scatterersvel,scatterersaxes] = stepImpl(obj,dt,varargin) 
            cond = (dt<0);
            if cond
                coder.internal.errorIf(cond,'phased:step:expectedNonnegative','T');
            end
            
            if strcmpi(obj.InputSource,'Input port')
                validateStepInputs(obj,dt,varargin{:}); 
            else
                validateBicyclistSpeed(obj); 
            end
            
            [scattererspos,scatterersvel,scatterersaxes] = pedaling(obj,dt,varargin{:});
            
            obj.pCurrentTime = obj.pCurrentTime+dt;
            phi = obj.pHeading;
            v = obj.pSpeed;
            obj.pPosition = obj.pPosition+v*dt*[cosd(phi);sind(phi);0];
            obj.pOrientationAxes = rotz(phi);
        end
    end   
end
