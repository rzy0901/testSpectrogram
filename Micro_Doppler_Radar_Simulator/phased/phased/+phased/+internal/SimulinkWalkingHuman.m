classdef (Sealed) SimulinkWalkingHuman < phased.internal.AbstractWalkingHuman & matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%WalkingHuman   Walking human
%   H = phased.WalingHuman creates a walking human motion model, H. This
%   object simulates the motion of a pedestrian
%
%   H = phased.WalkingHuman(Name,Value) returns a walking human, H, with
%   the specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   [BPPOS,BPVEL,BPAX] = step(H,T) returns the pedestrian's body parts
%   position, velocity, and orientation axes at the current time in BPPOS,
%   BPVEL, and BPAX, respectively. It then simulates the walking motion in
%   the next duration, specified in T (in seconds).
%
%   The pedestrian model includes 16 body parts: left and right feet, left
%   and right lower legs, left and right upper legs, left and right hip,
%   left and right lower arms, left and right upper arms, left and right
%   shoulders, neck, and head. Therefore BPPos is a 3x16 matrix with each
%   column represents the position of the corresponding body parts in the
%   [x;y;z] format (in meters). BPVel is also a 3x16 matrix whose columns
%   are velocities of corresponding body parts in the [x;y;z] form (in
%   m/s). BPAX is a 3x3x16 array whose pages are orientation axes of the
%   corresponding body parts. The three columns represents the 3 axes and
%   each column is in [x;y;z] format.
%
%   [BPPOS,BPVEL,BPAX] = step(H,T,ANGH) also specifies the current heading
%   angle, ANGH (in degree), as a scalar. The angle is in the xy-plane.
%   This syntax applies when you set the HeadingSource property to 'Input
%   port'.
%
%   [...,BJPOS] = step(...) also returns body joints positions in BJPOS.
%   There are 17 joints positions returned: left and right toes, left and
%   right ankles, left and right knees, left and right hips, left and right
%   hands, left and right elbows, left and right shoulders, neck, head, and
%   center between the hips. Thus BJPOS is a 3x17 matrix whose columns are
%   positions of corresponding body joints in the [x;y;z] format.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   WalkingHuman methods:
%
%   step     - Output current body parts position, velocity and orientation 
%              axes of the pedestrian (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a walking human object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset the walking human to its initial position
%
%   WalkingHuman properties:
%
%   Height                  - Pedestrian height 
%   WalkingSpeed            - Pedestrian walking speed 
%   InitialPosition         - Initial position 
%   OutputJointsPosition    - Whether to output body joints position
%   HeadingSource           - Source of heading direction
%   Heading                 - Heading direction 
%   InitialHeading          - Initial heading direction 
%
%   % Example:
%   %   Model the motion of a pedestrian to walk along a circle. 
%
%   ped = phased.WalkingHuman('HeadingSource','Input port',...
%                             'OutputJointsPosition',false);
%   dt = 0.003;
%   N = 3000;
%   ppos = zeros(3,16,N);
%   pax = zeros(3,3,16,N);
%   for m = 1:N
%       [ppos(:,:,m),~,pax(:,:,:,m)] = step(ped,dt,0.03*m);
%   end
%     
%   % plot
%   for m = 1:N
%       if m == 1
%           ph = scatter3(ppos(1,:,m),ppos(2,:,m),ppos(3,:,m),20);
%           phax = line(...
%             [ppos(1,:,m);ppos(1,:,m)+0.1*shiftdim(pax(1,1,:,m),1)],...
%             [ppos(2,:,m);ppos(2,:,m)+0.1*shiftdim(pax(2,1,:,m),1)],...
%             [ppos(3,:,m);ppos(3,:,m)+0.1*shiftdim(pax(3,1,:,m),1)],...
%             'Color','r');
%           axis equal
%           xlim([-1 20]);
%           ylim([-1 10]);
%           zlim([0 2]);
%           view(36,15);
%       else
%           ph.XData = ppos(1,:,m);
%           ph.YData = ppos(2,:,m);
%           ph.ZData = ppos(3,:,m);
%           lx = [ppos(1,:,m);ppos(1,:,m)+0.1*shiftdim(pax(1,1,:,m),1)];
%           ly = [ppos(2,:,m);ppos(2,:,m)+0.1*shiftdim(pax(2,1,:,m),1)];
%           lz = [ppos(3,:,m);ppos(3,:,m)+0.1*shiftdim(pax(3,1,:,m),1)];
%           for lidx = 1:16
%               phax(lidx).XData = lx(:,lidx).';
%               phax(lidx).YData = ly(:,lidx).';
%               phax(lidx).ZData = lz(:,lidx).';
%           end
%       end
%       drawnow;
%   end
%
%   See also phased, phased.Platform.
    
%   Copyright 2018 The MathWorks, Inc.

%   Reference
%   [1] Victor Chen, The Micro-Doppler Effect in Radar, Artech House, 2011
%   
%   [2] Boulic, Ronan, et al. A Global Human Walking Model with Real-time
%   Kinematic Personification, The Visual Computer: International Journal
%   of Computer Graphics, Vol. 6, Issue 6, Dec 1990

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods 
        function obj = SimulinkWalkingHuman(varargin)
            obj@phased.internal.AbstractWalkingHuman(varargin{:});
        end
        
    end
    
    methods (Access=protected)
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = isInactivePropertyImpl@phased.internal.AbstractWalkingHuman(obj,prop);
        end

        function validateInputsImpl(obj,heading)
            % Validate inputs to the step method at initialization
            if strcmp(obj.HeadingSource,'Input port')
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
                cond =  ~isreal(heading);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:Platform:MustBeReal', 'ANGH');
                end
            end
            
        end

        function num = getNumInputsImpl(obj)
            % Define total number of inputs for system with optional inputs
            if strcmp(obj.HeadingSource,'Property')
                num = 0;
            else
                num = 1;
            end
        end

        function icon = getIconImpl(obj) %#ok<MANU>
            % Define icon for System block
            icon = "Walking Human"; 
            
        end

        function name = getInputNamesImpl(obj) %#ok<MANU>
            % Return input port names for System block
            name = 'angh';
        end

        function [name,name2,name3,name4] = getOutputNamesImpl(obj) %#ok<MANU>
            % Return output port names for System block
            name = 'pos';
            name2 = 'vel';
            name3 = 'ax';
            name4 = 'jtpos';
        end

        function setupImpl(obj)
            setupImpl@phased.internal.AbstractWalkingHuman(obj);
        end
        
        function resetImpl(obj)
            obj.pCurrentTime = 0;
            obj.pPosition = obj.InitialPosition;
            if strcmp(obj.HeadingSource,'Property')
                obj.pHeading = obj.Heading;
            else
                obj.pHeading = obj.InitialHeading;
            end
            obj.pOrientationAxes = rotz(obj.pHeading);
            obj.pPreviousBodyPartsPosition = obj.getBodyPartsPosition(...
                obj.pBodyPartsInitialJointPosition);
        end
        
        function [bodypartspos,bodypartsvel,bodypartsax,bodyjointspos] = stepImpl(obj,heading)
            T = getCurrentTime(obj);
            if T==0
                dt = 0;
            else
                dt = T-obj.pCurrentTime;
                phi = obj.pHeading;
                v = obj.WalkingSpeed;
                obj.pPosition = obj.pPosition+v*dt*[cosd(phi);sind(phi);0];
                obj.pOrientationAxes = rotz(phi);
            end
            obj.pCurrentTime = T;
            [bodypartspos,bodypartsvel,bodypartsax,bodyjointspos] = walking(obj,dt,heading);
        end
        
        function loadObjectImpl(obj,s,~)
            % Set properties in object obj to values in structure s
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractWalkingHuman(obj);
            
        end

    end
    
end
