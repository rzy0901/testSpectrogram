classdef (Sealed,StrictDefaults) Platform < phased.internal.AbstractPlatform
%Platform Motion platform
%   H = phased.Platform creates a platform System object, H. This object
%   models the motion of one or more platforms in space.
%
%   H = phased.Platform(Name,Value) creates a platform object, H, with the
%   specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   H = phased.Platform(P,V,Name,Value) creates a platform object, H,
%   with the InitialPosition property set to P, the Velocity property set
%   to V and other specified property Names set to the specified Values.
%   P and V are value-only arguments. To specify a value-only argument,
%   you must also specify all preceding value-only arguments. You can
%   specify Name-Value pair arguments in any order.
%
%   Step method syntax:
%
%   [POS,VEL] = step(H,T) returns the current positions, POS (in meters),
%   and the current velocities, VEL (in m/s), of the platforms.
%
%   If you set the MotionModel property to 'Velocity' and the
%   VelocitySource property to 'Property', the velocity is assumed to be a
%   constant within each time step T. The object updates each position
%   using the equation
%
%       POS = POS + VEL*T
%
%   where T specifies the elapsed time (in seconds) for the current step.
%
%   If you set the MotionModel property to 'Acceleration' and the
%   AccelerationSource property to 'Property', the acceleration is assumed
%   to be a constant within the time step T. The object updates each
%   position using the equation
%
%       POS = POS + VEL*T + 1/2*ACL*T^2
%       VEL = VEL + ACL*T
%
%   where ACL is the acceleration (in m/s^2).
%
%   If you set the MotionModel property to 'Custom'. The object performs a
%   piecewise cubic interpolation on the waypoints to derive the position
%   and velocity at each time step.
%   
%   [POS,VEL] = step(H,T,V) specifies the current velocity in V. V is a
%   3xN matrix where N is the number of platforms to model. This syntax
%   applies when you set the MotionModel property to 'Velocity' and the
%   VelocitySource property to 'Input port'.
%
%   [POS,VEL] = step(H,T,A) specifies the current acceleration in A. A is a
%   3xN matrix where N is the number of platforms to model. This syntax
%   applies when you set the MotionModel property to 'Acceleration' and the
%   AccelerationSource property to 'Input port'.
%
%   [POS,VEL,LAXES] = step(...) returns the additional output LAXES as the
%   platform's orientation axes when you set the OrientationAxesOutputPort
%   property to true.
%    
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   Platform methods:
%
%   step     - Output current position, velocity and orientation axes of 
%              the platform (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a platform object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset the platform to its initial position
%
%   Platform properties:
%
%   MotionModel               - Model of object motion
%   InitialPosition           - Initial position
%   InitialVelocity           - Initial velocity
%   VelocitySource            - Source of velocity
%   Velocity                  - Velocity 
%   AccelerationSource        - Source of acceleration
%   Acceleration              - Acceleration
%   CustomTrajectory          - Custom trajectory
%   ScanMode                  - Mechanical scanning mode
%   InitialScanAngle          - Initial scan angle
%   AzimuthSpan               - Azimuth scan angle span
%   AzimuthScanRate           - Azimuth scan rate
%   InitialOrientationAxes    - Initial orientation axes 
%   OrientationAxesOutputPort - Enable orientation axes output
%
%   % Examples:
%
%   % Example 1:
%   %   Define a platform at the origin with a velocity of [100; 100; 0] in
%   %   meters per second. Simulate the motion of the platform for 2 steps,
%   %   assuming the time elapsed for each step is 1 second.
%
%   platform = phased.Platform([0; 0; 0],[100; 100; 0]);
%   T = 1;
%   [pos,v] = platform(T)
%   [pos,v] = platform(T)
%
%   % Example 2:
%   %   Simulates the trajectory of an object for 15 seconds. The object
%   %   starts from the origin and has an acceleration of [1;0;0] m/s^2.
%   
%   platform = phased.Platform('MotionModel','Acceleration',...
%           'Acceleration',[1;0;0]);
%   dt = 1;  N = 15;  t = (0:N-1)*dt;
%   pos = zeros(3,N); vel = zeros(3,N);
% 
%   for m = 1:N
%       [pos(:,m), vel(:,m)] = platform(dt);
%   end
% 
%   ax = plotyy(t,pos(1,:),t,vel(1,:)); xlabel(ax(1),'Time (s)');
%   ylabel(ax(1),'Position (m)'); ylabel(ax(2),'Velocity (m/s)');
%
%   % Example 3:
%   %   Simulate a circular trajectory. Assume the initial orientation is
%   %   along x-axes, plot both the trajectory and the orientation along
%   %   the way.
%
%   platform = phased.Platform('VelocitySource','Input port',...
%       'OrientationAxesOutputPort',true);
%   dt = 1; N = 360;
%   pos = zeros(3,N); vel = zeros(3,N); naxes = zeros(3,N);
%  
%   for m = 1:N
%       [pos(:,m),vel(:,m), ax] = platform(dt,...
%               [cosd((m-1)*dt);sind((m-1)*dt);0]);
%       naxes(:,m) = ax(:,1);
%   end
% 
%   plot(pos(1,:),pos(2,:)); hold on;
%   quiver(pos(1,1:20:end),pos(2,1:20:end),...
%           naxes(1,1:20:end),naxes(2,1:20:end));
%   xlabel('X'); ylabel('Y'); axis square;
%   legend({'Trajectory','Orientation'});
%
%   See also phased, phased.Radiator, phased.Collector, global2localcoord,
%   local2globalcoord, rangeangle.

%   Copyright 2010-2018 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods

        function obj = Platform(varargin)
            %Platform   Construct the Platform class.
            obj@phased.internal.AbstractPlatform(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        
        function t = getElapsedTime(obj,varargin)
            t = varargin{1};
            sigdatatypes.validateDuration(t,'step','T');
            if obj.pIsVelocityViaInput
                obj.pVelocity = varargin{2};
            elseif obj.pIsAccelerationViaInput
                obj.pAcceleration = varargin{2};
            end
        end
        
        function validateInputsImpl(obj, varargin)
            numOfPlatforms = size(obj.InitialPosition,2);                        

            t = varargin{1};
            cond =  ~isscalar(t);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeScalar','T');
            end
            cond =  ~isa(t,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','T','double');
            end
            cond =  ~isreal(t);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:Platform:MustBeReal', 'T');
            end
            if strcmp(obj.MotionModel,'Velocity') 
                if strcmp(obj.VelocitySource,'Input port')
                    v = varargin{2};
                    sigdatatypes.validate3DCartCoord(v,...
                        'step','V',{'ncols',numOfPlatforms});
                end
            else % MotionModel is Acceleration
                if strcmp(obj.AccelerationSource,'Input port')
                    a = varargin{2};
                    sigdatatypes.validate3DCartCoord(a,...
                        'step','A',{'ncols',numOfPlatforms});
                end
            end
                        
            
        end
        
        function num = getNumInputsImpl(obj)
            num = 1;
            if (strcmp(obj.MotionModel,'Velocity') && ...
                    strcmp(obj.VelocitySource,'Input port')) || ...
                    (strcmp(obj.MotionModel,'Acceleration') && ...
                    strcmp(obj.AccelerationSource,'Input port'))
                num = num+1;
            end
        end
        
        function loadObjectImpl(obj,s,~)
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
    end
    
    methods (Static, Hidden, Access = protected)        
        function groups = getPropertyGroupsImpl
            dInitialOrientationAxes = matlab.system.display.internal.Property(...
                'InitialOrientationAxes','Alias', 'OrientationAxes');
            groups = matlab.system.display.Section('phased.Platform', ...
                            'DependOnPrivatePropertyList',...
                            {'Velocity','Acceleration'});            
            for m = 1:numel(groups.PropertyList)
                if strcmp(groups.PropertyList{m},'InitialOrientationAxes')
                    groups.PropertyList{m} = dInitialOrientationAxes;
                end
            end
        end
    end
    
    methods (Static,Hidden)
        function a = getAlternateBlock
            a = 'phasedenvlib/Motion Platform';
        end
    end
end

