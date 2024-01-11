classdef (Sealed, StrictDefaults) SimulinkPlatform < phased.internal.AbstractPlatform & ...
        matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.

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
%   Platform methods:
%
%   step     - Output current position, velocity and orientation axes of 
%              the platform (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a platform object with same property values
%   isLocked - Locked status (logical)
%   reset    - Reset the platform to its initial position
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
%   plat = phased.Platform([0; 0; 0],[100; 100; 0]);
%   T = 1;
%   [pos,v] = step(plat,T)
%   [pos,v] = step(plat,T)
%
%   % Example 2:
%   %   Simulates the trajectory of an object for 15 seconds. The object
%   %   starts from the origin and has an acceleration of [1;0;0] m/s^2.
%   
%   plat = phased.Platform('MotionModel','Acceleration',...
%           'Acceleration',[1;0;0]);
%   dt = 1;  N = 15;  t = (0:N-1)*dt;
%   pos = zeros(3,N); vel = zeros(3,N);
% 
%   for m = 1:N
%       [pos(:,m), vel(:,m)] = step(plat,dt);
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
%   plat = phased.Platform('VelocitySource','Input port',...
%       'OrientationAxesOutputPort',true);
%   dt = 1; N = 360;
%   pos = zeros(3,N); vel = zeros(3,N); naxes = zeros(3,N);
%  
%   for m = 1:N
%       [pos(:,m),vel(:,m), ax] = step(plat,dt,...
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

%   Copyright 2016-2017 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %ElapsedTimeSource Source of elapsed simulation time
        %   Specify the source of elapsed simulation time between
        %   simulation steps using one of 'Auto' | 'Derive from reference
        %   input port', where the default is 'Auto'. If you set the
        %   ElapsedTimeSource property to 'Auto', the elapsed time is
        %   derived from Simulink time engine. If you set the
        %   ElapsedTimeSource property to 'Derive from reference input
        %   port', the elapsed time is computed using the number of samples
        %   in the reference input signal and the signal's sample rate.
        ElapsedTimeSource = 'Auto'
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
        %SampleRate     Sample rate (Hz)
        %   Specify the matched filter coefficients sample rate (in Hz) as
        %   a positive scalar. The default value of this property is 1e6.
        %   This property applies when you set SpectrumWindow property to a
        %   value other than 'None'.
        SampleRate = 1e6;
    end
    
    properties (Access = protected)
        pSampleRate
    end
    
    
    properties(Constant, Hidden)
        ElapsedTimeSourceSet = matlab.system.StringSet(...
            {'Auto','Derive from reference input port'});
    end
    
    methods

        function obj = SimulinkPlatform(varargin)
            %Platform   Construct the Platform class.
            obj@phased.internal.AbstractPlatform(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        function setupImpl(obj,varargin)
            setupImpl@phased.internal.AbstractPlatform(obj,varargin{:});
            if ~strcmp(obj.ElapsedTimeSource,'Auto')
                if obj.SampleRateFromInputCheckbox
                    fs = getSampleRateInSimulation(obj);
                    cond = ~isscalar(fs) || (fs<=0);
                    if cond
                        coder.internal.errorIf(cond,...
                             'phased:phased:invalidSampleTime');
                    end
                else
                    fs = obj.SampleRate;
                end
                obj.pSampleRate = fs;

                X = varargin{1};
                obj.pNumInputChannels = getNumChannels(obj,X);
                obj.pValidatedNumInputChannels = getNumChannels(obj,X);
            else
                sts = getSampleTime(obj);
                cond = strcmp(sts.Type,'Controllable');
                if cond
                    coder.internal.errorIf(cond,...
                        'phased:Platform:invalidPlatformSampleTime',...
                        'ElapsedTimeSource','Auto','Controllable',...
                        'Derive from reference input port');
                end
            end
        end
        
        function validateInputsImpl(obj, varargin)
            numOfPlatforms = size(obj.InitialPosition,2);                        
            if ~strcmp(obj.ElapsedTimeSource,'Auto')
                x = varargin{1};
                if ~isempty(x)  % g1398193
                    validateNumChannels(obj,x);
                end
                
                cond =  ~isa(x,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','Ref','double');
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
                        
            else  % Auto
                if strcmp(obj.MotionModel,'Velocity') 
                    if strcmp(obj.VelocitySource,'Input port')
                        v = varargin{1};
                        sigdatatypes.validate3DCartCoord(v,...
                            'step','V',{'ncols',numOfPlatforms});
                    end
                else % MotionModel is Acceleration
                    if strcmp(obj.AccelerationSource,'Input port')
                        a = varargin{1};
                        sigdatatypes.validate3DCartCoord(a,...
                            'step','A',{'ncols',numOfPlatforms});
                    end
                end
                        
            end
            
        end
        
        function t = getElapsedTime(obj,varargin)
            if strcmp(obj.ElapsedTimeSource,'Auto')
                sts = getSampleTime(obj);
                if strcmp(sts.Type,'Discrete') || strcmp(sts.Type,'Inherited')
                    t = sts.SampleTime; 
                else
                    coder.internal.errorIf(true, ...
                        'phased:phased:invalidSampleTimeType',sts.Type);
                end
                cond = (t<=0); % scalar guaranteed in input validation
                if cond
                    coder.internal.errorIf(cond,...
                         'phased:phased:invalidSampleTime');
                end
                if obj.pIsVelocityViaInput
                    obj.pVelocity = varargin{1};
                elseif obj.pIsAccelerationViaInput
                    obj.pAcceleration = varargin{1};
                end
            else
                x = varargin{1};
                cond =  ~isa(x,'double') ;
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','Ref','double');
                end
                t = size(x,1)/obj.pSampleRate;
                if obj.pIsVelocityViaInput
                    obj.pVelocity = varargin{2};
                elseif obj.pIsAccelerationViaInput
                    obj.pAcceleration = varargin{2};
                end
            end
        end
        
        function flag = isInputSizeLockedImpl(obj,ind)
            flag = isInputSizeLockedImpl@phased.internal.AbstractPlatform(obj,ind);
            if ~strcmp(obj.ElapsedTimeSource,'Auto') && ind==1
                flag = false;
            end
        end

        function num = getNumInputsImpl(obj)
            if strcmp(obj.ElapsedTimeSource,'Auto')
                num = 0;
            else
                num = 1;
            end
            if (strcmp(obj.MotionModel,'Velocity') && ...
                    strcmp(obj.VelocitySource,'Input port')) || ...
                    (strcmp(obj.MotionModel,'Acceleration') && ...
                    strcmp(obj.AccelerationSource,'Input port'))
                num = num+1;
            end
                
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractPlatform(obj, prop);
            if strcmp(obj.ElapsedTimeSource,'Auto') && ...
                    (strcmp(prop,'SampleRateFromInputCheckbox') || ...
                    strcmp(prop,'SampleRate'))
                flag = true;
            elseif obj.SampleRateFromInputCheckbox && strcmp(prop,'SampleRate')
                flag = true;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractPlatform(obj);
            if isLocked(obj)
                s.pSampleRate = obj.pSampleRate;
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
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:PlatformTitle')),...
              'Text',getString(message('phased:library:block:PlatformDesc')));
        end
        
        function groups = getPropertyGroupsImpl
            dInitialOrientationAxes = matlab.system.display.internal.Property(...
                'InitialOrientationAxes','Alias', 'OrientationAxes');
            groups = matlab.system.display.Section('phased.internal.SimulinkPlatform', ...
                            'DependOnPrivatePropertyList',...
                            {'Velocity','Acceleration'});     
            groups.PropertyList = groups.PropertyList([4:end 1:3]);
            for m = 1:numel(groups.PropertyList)
                if strcmp(groups.PropertyList{m},'InitialOrientationAxes')
                    groups.PropertyList{m} = dInitialOrientationAxes;
                end
            end
        end
    end
    
    methods (Access = protected)
        function varargout = getOutputNamesImpl(obj) 
            if obj.OrientationAxesOutputPort
                varargout = {'Pos','Vel','LAxes'};
            else
                varargout = {'Pos','Vel'};
            end
        end
        
        function varargout = getInputNamesImpl(obj)
            idx = 1;
            if ~strcmp(obj.ElapsedTimeSource,'Auto')
                varargout{idx} = 'Ref';
                idx = idx+1;
            end
            if strcmp(obj.MotionModel,'Velocity') && ...
                    strcmp(obj.VelocitySource,'Input port')
                varargout{idx} = 'Vel';
            elseif strcmp(obj.MotionModel,'Acceleration') && ...
                    strcmp(obj.AccelerationSource,'Input port')
                varargout{idx} = 'Acl';
            end
        end
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = 'Platform';
        end
        function varargout = getOutputSizeImpl(obj) 
            sz = size(obj.InitialPosition);
            axSz = [3 3 sz(2)];
            varargout = {sz,sz,axSz};
        end
        function varargout = isOutputFixedSizeImpl(obj) %#ok<MANU>
            varargout = {true, true, true};
            
        end
        function varargout = getOutputDataTypeImpl(obj) %#ok<MANU>
            varargout = {'double','double','double'};
        end
        function varargout = isOutputComplexImpl(obj) %#ok<MANU>
            varargout = {false, false,false};
        end    
    end
    
end

