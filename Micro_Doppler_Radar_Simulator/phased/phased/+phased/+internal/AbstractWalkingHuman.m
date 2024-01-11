classdef (Hidden) AbstractWalkingHuman < phased.internal.AbstractSampleRateEngine & matlab.system.mixin.Propagates 
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

    properties (Nontunable)
        %Height     Pedestrian height (m)
        %   Specify the height (in meters) of the pedestrian as a positive
        %   scalar. The default value is 1.65.
        Height = 1.65
        %WalkingSpeed   Pedestrian walking speed (m/s)
        %   Specify the walking speed (in m/s) of the pedestrian as a
        %   nonnegative scalar. The default value is 1.4.
        %
        %   The motion model limits the walking speed to about 1.4 times of
        %   pedestrian's height.
        WalkingSpeed = 1.4
        %InitialPosition    Initial position (m)
        %   Specify the initial position of the pedestrian as a 3x1 vector
        %   in the form of [x; y; z] (in meters). The default value of this
        %   property is [0; 0; 0].
        InitialPosition = [0;0;0]
    end
    
    properties (Nontunable,Logical)
        %OutputJointsPosition   Whether to output body joints position
        %   Set this property to true to output body joints position at
        %   each time step when walking. Set this property to false to
        %   not output body joints position at each time step. The default
        %   value is false.
        OutputJointsPosition = false
    end
    
    properties (Nontunable)
        %HeadingSource  Source of heading direction
        % Specify the source of heading direction as one of 'Property' |
        % 'Input port', where the default is 'Property'. When you set the
        % HeadingSource property to 'Property', the walking is along a
        % straight line and the direction is specified in the Heading
        % property. If you set the HeadingSource property to 'Input port',
        % the initial heading direction is specified in the InitialHeading
        % property and then the successive heading at each time step is
        % specified via the input.
        HeadingSource = 'Property'
    end
    
    properties (Nontunable)
        %Heading    Heading direction (deg)
        %   Specify the initial heading (in degrees) of the pedestrian as a
        %   scalar. The heading is measured from x-axis towards y-axis in
        %   the xy-plane. The default value is 0.
        Heading = 0;
        %InitialHeading    Initial heading direction (deg)
        %   Specify the initial heading (in degrees) of the pedestrian as a
        %   scalar. The heading is measured from x-axis towards y-axis in
        %   the xy-plane. The default value is 0.
        InitialHeading = 0;
    end
    
    properties(Constant, Hidden)
        HeadingSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    properties (Access=protected, Nontunable)
        
        pRelativeSpeed
        pRelativeCycleLength
        pFootOpenAngle
        pDurationCycle
        pDurationSupport
        pDurationBalance
        pDurationDoubleSupport
        pRelativeDurationSupport
        
        pNumSamplesPerCycle = 1000
        pCycleTime
        
        pBodyPartsInitialJointPosition
        pBodyPartsLength
        
        pApplyVirtualFloor = true
    end
    
    properties (Access = protected)
        pLeftToeTrajectory
        pRightToeTrajectory
        pLeftAnkleTrajectory
        pRightAnkleTrajectory
        pLeftKneeTrajectory
        pRightKneeTrajectory
        pLeftHipTrajectory
        pRightHipTrajectory
        pLeftHandTrajectory
        pRightHandTrajectory
        pLeftElbowTrajectory
        pRightElbowTrajectory
        pLeftShoulderTrajectory
        pRightShoulderTrajectory
        pHeadTrajectory
        pNeckTrajectory
        pOriginTrajectory
        
        pBodyJointsTrajectory
        pBodyPartsLocalAxes
        pPreviousBodyPartsPosition
        pCurrentTime
        pHeading
        pPosition
        pOrientationAxes
    end
    
    methods
        function set.Height(obj,val)
            sigdatatypes.validateDistance(val,'','Height',{'double','single'},{'scalar','positive'});
            obj.Height = val;
        end
        function set.WalkingSpeed(obj,val)
            sigdatatypes.validateSpeed(val,'','WalkingSpeed',{'double','single'},{'scalar'});
            obj.WalkingSpeed = val;
        end
        function set.InitialHeading(obj,val)
            sigdatatypes.validateAngle(val,...
                'BackscatterPedestrianTarget','InitialHeading',{'scalar'});
            obj.InitialHeading = val;
        end
    end
    
    methods (Access=protected)
        function obj = AbstractWalkingHuman(varargin)
            setProperties(obj, nargin, varargin{:});
        end
        
        function validatePropertiesImpl(obj)
            H = obj.Height;
            Ht = 0.245*H+0.246*H;
            V = obj.WalkingSpeed;
            rv = V/Ht;
            cond = (rv > 3);
%             if cond
%                 coder.internal.errorIf(cond,'phased:target:invalidWalkingSpeed',...
%                     num2str(V),num2str(H));
%             end
                
        end

        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            if strcmp(prop,'InitialHeading')
                if strcmp(obj.HeadingSource,'Property')
                    flag = true;
                else
                    flag = false;
                end
            elseif strcmp(prop,'Heading')
                if strcmp(obj.HeadingSource,'Property')
                    flag = false;
                else
                    flag = true;
                end
            else
                flag = false;
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

        function num = getNumOutputsImpl(obj) 
            % Define total number of outputs for system with optional
            % outputs
            if obj.OutputJointsPosition
                num = 4;
            else
                num = 3;
            end
        end

        function [out,out2,out3,out4] = getOutputSizeImpl(obj) %#ok<MANU>
            % Return size for each output port
            out = [3 16];
            out2 = [3 16];
            out3 = [3 3 16];
            out4 = [3 17];

            % Example: inherit size from first input port
            % out = propagatedInputSize(obj,1);
        end

        function [out,out2,out3,out4] = getOutputDataTypeImpl(obj) %#ok<MANU>
            % Return data type for each output port
            out = "double";
            out2 = "double";
            out3 = "double";
            out4 = "double";

        end

        function [out,out2,out3,out4] = isOutputComplexImpl(obj) %#ok<MANU>
            % Return true for each output port with complex data
            out = false;
            out2 = false;
            out3 = false;
            out4 = false;

        end

        function [out,out2,out3,out4] = isOutputFixedSizeImpl(obj) %#ok<MANU>
            % Return true for each output port with fixed size
            out = true;
            out2 = true;
            out3 = true;
            out4 = true;

        end
        
        function setupImpl(obj)
            setupHumanScattererPoints(obj);
            setupWalkingCharacteristics(obj);
            
            nspercyc = obj.pNumSamplesPerCycle;
            tcyc = (0:nspercyc-1)/nspercyc; % 1 cycle
            obj.pCycleTime = tcyc(:);

            RV = obj.pRelativeSpeed;
            RDs = obj.pRelativeDurationSupport;
            
            % compute joint trajectory
            [transvert,translat,transfb] = obj.getTranslationTrajectory(RV,RDs,tcyc);
            [rotfb,rotlr,rottor] = obj.getPelvisRotationTrajectory(RV,tcyc);
            [flexhip,flexknee,flexankle] = obj.getFlexTrajectory(RV,RDs,tcyc);
            [thoraxtor,flexshoulder,flexelbow] = obj.getUpperBodyTrajectory(RV,tcyc);
            
            % compute body part trajectory
            [ltoe,rtoe,ltoeax,rtoeax] = flexAnkle(obj,flexankle);
            [ltoe,rtoe,lankle,rankle,ltoeax,rtoeax,lankleax,rankleax] = ...
                flexKnee(obj,flexknee,ltoe,rtoe,ltoeax,rtoeax);
            [ltoe,rtoe,lankle,rankle,lknee,rknee,ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax] = ...
                flexHip(obj,flexhip,ltoe,rtoe,lankle,rankle,ltoeax,rtoeax,lankleax,rankleax);
            [ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax,lhipax,rhipax] = ...
                rotateLowerBody(obj,rotlr,rottor,ltoe,rtoe,lankle,rankle,lknee,rknee,...
                ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax);
            [lhand,rhand,lhandax,rhandax] = flexElbow(obj,flexelbow);
            [lhand,rhand,lelbow,relbow,lhandax,rhandax,lelbowax,relbowax] = ...
                flexShoulder(obj,flexshoulder,lhand,rhand,lhandax,rhandax);
            [lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head,...
                lhandax,rhandax,lelbowax,relbowax,lshoulderax,rshoulderax,neckax,headax] = ...
                rotateUpperBody(obj,rotfb,thoraxtor,lhand,rhand,lelbow,relbow,...
                lhandax,rhandax,lelbowax,relbowax);
            [ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head,lclorigin,lclax] = translateBody(obj,...
                transfb,translat,transvert,ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head);
            
            obj.pLeftToeTrajectory = ltoe;
            obj.pRightToeTrajectory = rtoe;
            obj.pLeftAnkleTrajectory = lankle;
            obj.pRightAnkleTrajectory = rankle;
            obj.pLeftKneeTrajectory = lknee;
            obj.pRightKneeTrajectory = rknee;
            obj.pLeftHipTrajectory = lhip;
            obj.pRightHipTrajectory = rhip;
            obj.pLeftHandTrajectory = lhand;
            obj.pRightHandTrajectory = rhand;
            obj.pLeftElbowTrajectory = lelbow;
            obj.pRightElbowTrajectory = relbow;
            obj.pLeftShoulderTrajectory = lshoulder;
            obj.pRightShoulderTrajectory = rshoulder;
            obj.pHeadTrajectory = head;
            obj.pNeckTrajectory = neck;
            obj.pOriginTrajectory = lclorigin;
            
            bodyjointstraj = cat(3,ltoe,rtoe,lankle,rankle,lknee,rknee,...
                lhip,rhip,lhand,rhand,lelbow,relbow,lshoulder,rshoulder,head,neck,lclorigin);
            if obj.pApplyVirtualFloor
                bodyjointstraj(3,:,:) = bsxfun(@minus,bodyjointstraj(3,:,:),min(bodyjointstraj(3,:,:),[],3));
            end
            
            torsoax = neckax;
            bodypartsax = cat(4,ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax,lhipax,rhipax,...
                lhandax,rhandax,lelbowax,relbowax,lshoulderax,rshoulderax,headax,torsoax,lclax);
            
            obj.pBodyJointsTrajectory = bodyjointstraj;
            obj.pBodyPartsLocalAxes = bodypartsax(:,:,:,1:end-1); % no need for fixed lclax
            
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
        
        function [bodypartspos,bodypartsvel,bodypartsax,bodyjointspos] = walking(obj,dt,heading)
            pos = obj.pPosition;
            ax = obj.pOrientationAxes;
            if obj.pCurrentTime == 0
                bodyjointspos = bsxfun(@plus,ax*obj.pBodyPartsInitialJointPosition,pos);
                bodypartspos = obj.getBodyPartsPosition(bodyjointspos);
                bodypartsvel = zeros(size(bodypartspos));
                numbodyparts = size(bodypartspos,2);
                bodypartsax = repmat(ax,1,1,numbodyparts);
            else
                traj = obj.pBodyJointsTrajectory;
                bodypartspos_old = obj.pPreviousBodyPartsPosition;
                % use linear interpolation for position and velocity
                % instead of nearest neighbor to get better velocity
                % estimate to avoid upsample aliasing side effect
                Dc = obj.pDurationCycle;
                Nspc = obj.pNumSamplesPerCycle;
                rt = mod(obj.pCurrentTime,Dc)/Dc; % always less than 1
                rtidx = floor(rt*Nspc)+1;  % rtidx start with 1
                rtidxbound = [rtidx rtidx+1]-1;  % Bound start with 0 to align with period
                interp_idx = mod(rtidxbound,Nspc)+1; % used to assess cached trajctory
                
                % manual interp1 to avoid overhead
                alpha = Nspc*rt-(rtidx-1);
                interp_traj = (1-alpha)*traj(:,interp_idx(1),:)+...
                    alpha*traj(:,interp_idx(2),:);
                
                bodyjointspos = bsxfun(@plus,ax*squeeze(interp_traj),pos);
                bodypartspos = obj.getBodyPartsPosition(bodyjointspos);
                bodypartsvel = (bodypartspos-bodypartspos_old)/dt;
                bodypartsax = squeeze(obj.pBodyPartsLocalAxes(:,:,rtidx,:));
                for m = 1:size(bodypartsax,3)
                    bodypartsax(:,:,m) = ax*bodypartsax(:,:,m);
                end
            end
            obj.pPreviousBodyPartsPosition = bodypartspos;
            
            if ~strcmp(obj.HeadingSource,'Property')
                obj.pHeading = heading;
            end
        end
        
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            
            if isLocked(obj)
                s.pRelativeSpeed = obj.pRelativeSpeed;
                s.pRelativeCycleLength = obj.pRelativeCycleLength;
                s.pFootOpenAngle = obj.pFootOpenAngle;
                s.pDurationCycle = obj.pDurationCycle;
                s.pDurationSupport = obj.pDurationSupport;
                s.pDurationBalance = obj.pDurationBalance;
                s.pDurationDoubleSupport = obj.pDurationDoubleSupport;
                s.pRelativeDurationSupport = obj.pRelativeDurationSupport;

                s.pNumSamplesPerCycle = obj.pNumSamplesPerCycle;
                s.pCycleTime = obj.pCycleTime;

                s.pBodyPartsInitialJointPosition = obj.pBodyPartsInitialJointPosition;
                s.pBodyPartsLength = obj.pBodyPartsLength;
                
                s.pApplyVirtualFloor = obj.pApplyVirtualFloor;
                
                s.pLeftToeTrajectory = obj.pLeftToeTrajectory;
                s.pRightToeTrajectory = obj.pRightToeTrajectory;
                s.pLeftAnkleTrajectory = obj.pLeftAnkleTrajectory;
                s.pRightAnkleTrajectory = obj.pRightAnkleTrajectory;
                s.pLeftKneeTrajectory = obj.pLeftKneeTrajectory;
                s.pRightKneeTrajectory = obj.pRightKneeTrajectory;
                s.pLeftHipTrajectory = obj.pLeftHipTrajectory;
                s.pRightHipTrajectory = obj.pRightHipTrajectory;
                s.pLeftHandTrajectory = obj.pLeftHandTrajectory;
                s.pRightHandTrajectory = obj.pRightHandTrajectory;
                s.pLeftElbowTrajectory = obj.pLeftElbowTrajectory;
                s.pRightElbowTrajectory = obj.pRightElbowTrajectory;
                s.pLeftShoulderTrajectory = obj.pLeftShoulderTrajectory;
                s.pRightShoulderTrajectory = obj.pRightShoulderTrajectory;
                s.pHeadTrajectory = obj.pHeadTrajectory;
                s.pNeckTrajectory = obj.pNeckTrajectory;
                s.pOriginTrajectory = obj.pOriginTrajectory;

                s.pBodyJointsTrajectory = obj.pBodyJointsTrajectory;
                s.pBodyPartsLocalAxes = obj.pBodyPartsLocalAxes;
                s.pPreviousBodyPartsPosition = obj.pPreviousBodyPartsPosition;
                s.pCurrentTime = obj.pCurrentTime;
                s.pHeading = obj.pHeading;
                s.pPosition = obj.pPosition;
                s.pOrientationAxes = obj.pOrientationAxes;
            end

        end

    end
    
    methods (Access=protected)
        function setupHumanScattererPoints(obj)
            [bodypartspos,bodypartslen] = phased.internal.getHumanBodyParts(obj.Height);
            
            if obj.pApplyVirtualFloor
                bodypartspos = phased.internal.AbstractWalkingHuman.applyVirtualFloor(bodypartspos);
            end
            
            obj.pBodyPartsInitialJointPosition = bodypartspos;
            obj.pBodyPartsLength = bodypartslen;
            
        end
        
        function setupWalkingCharacteristics(obj)
            hipheight = obj.pBodyPartsLength(15);
            RV = obj.WalkingSpeed/hipheight;
            obj.pRelativeSpeed = RV;
            obj.pRelativeCycleLength = 1.346*sqrt(RV);
            obj.pFootOpenAngle = -1.4*RV+8.5;%8.5
            
            Dc = obj.pRelativeCycleLength/RV;
            obj.pDurationCycle = Dc;
            obj.pDurationSupport = 0.752*Dc-0.143;
            obj.pDurationBalance = 0.248*Dc+0.143;
            obj.pDurationDoubleSupport = 0.252*Dc-0.143;
            obj.pRelativeDurationSupport = obj.pDurationSupport/Dc;
        end
        
        function [ltoe,rtoe,ltoeax,rtoeax] = flexAnkle(obj,angAnkle)
            nt = obj.pNumSamplesPerCycle;
            bodypartslen = obj.pBodyPartsLength;
            ltoe = zeros(3,nt);
            rtoe = zeros(3,nt);
            ltoeax = zeros(3,3,nt);
            rtoeax = zeros(3,3,nt);
            footlen = bodypartslen(1);
            ankleheight = bodypartslen(11);
            footvec = [footlen;0;-ankleheight];  % relative to ankle
            
            angAnkle_l = fftshift(angAnkle);
            angAnkle_r = angAnkle;
            
            for m = 1:nt
                rotmat_l = roty(-angAnkle_l(m));
                ltoe(:,m) = rotmat_l*footvec;
                ltoeax(:,:,m) = rotmat_l;
                rotmat_r = roty(-angAnkle_r(m));
                rtoe(:,m) = rotmat_r*footvec;
                rtoeax(:,:,m) = rotmat_r;
            end
            
            hiplen = bodypartslen(4);
            upperleglen = bodypartslen(3);
            lowerleglen = bodypartslen(2);
            
            ltoe = bsxfun(@plus,ltoe,[0;hiplen;-(upperleglen+lowerleglen)]);   % relative to hip
            rtoe = bsxfun(@plus,rtoe,[0;-hiplen;-(upperleglen+lowerleglen)]); % relative to hip
        end
        
        function [ltoe,rtoe,lankle,rankle,ltoeax,rtoeax,lankleax,rankleax] = ...
                flexKnee(obj,angKnee,ltoe,rtoe,ltoeax,rtoeax)
            nt = obj.pNumSamplesPerCycle;
            bodypartslen = obj.pBodyPartsLength;
            lankle = zeros(3,nt);
            rankle = zeros(3,nt);
            lankleax = zeros(3,3,nt);
            rankleax = zeros(3,3,nt);
            lowerleglen = bodypartslen(2);
            anklevec = [0;0;-lowerleglen];      % relative to knee
            
            upperleglen = bodypartslen(3);
            ltoe(2,:) = 0;
            ltoe(3,:) = ltoe(3,:)+upperleglen;  % relative to knee
            rtoe(2,:) = 0;
            rtoe(3,:) = rtoe(3,:)+upperleglen;  % relative to knee
            
            angKnee_l = fftshift(angKnee);
            angKnee_r = angKnee;
            
            for m = 1:nt
                rotmatl = roty(angKnee_l(m));
                lankle(:,m) = rotmatl*anklevec;
                ltoe(:,m) = rotmatl*ltoe(:,m);
                lankleax(:,:,m) = rotmatl;
                ltoeax(:,:,m) = rotmatl*ltoeax(:,:,m);
                rotmatr = roty(angKnee_r(m));
                rankle(:,m) = rotmatr*anklevec;
                rtoe(:,m) = rotmatr*rtoe(:,m);
                rankleax(:,:,m) = rotmatr;
                rtoeax(:,:,m) = rotmatr*rtoeax(:,:,m);
            end
            
            hiplen = bodypartslen(4);
            
            ltoe = bsxfun(@plus,ltoe,[0;hiplen;-upperleglen]);       % relative to hip
            rtoe = bsxfun(@plus,rtoe,[0;-hiplen;-upperleglen]);     % relative to hip
            lankle = bsxfun(@plus,lankle,[0;hiplen;-upperleglen]);   % relative to hip
            rankle = bsxfun(@plus,rankle,[0;-hiplen;-upperleglen]); % relative to hip
        end
        
        function [ltoe,rtoe,lankle,rankle,lknee,rknee,...
                ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax] = ...
                flexHip(obj,angHip,ltoe,rtoe,lankle,rankle,ltoeax,rtoeax,lankleax,rankleax)
            nt = obj.pNumSamplesPerCycle;
            bodypartslen = obj.pBodyPartsLength;
            lknee = zeros(3,nt);
            rknee = zeros(3,nt);
            lkneeax = zeros(3,3,nt);
            rkneeax = zeros(3,3,nt);
            upperleglen = bodypartslen(3);
            kneevec = [0;0;-upperleglen];       % relative to hip
            
            ltoe(2,:) = 0;                      % relative to hip
            rtoe(2,:) = 0;                      % relative to hip
            
            lankle(2,:) = 0;                      % relative to hip
            rankle(2,:) = 0;                      % relative to hip

            angHip_l = fftshift(-angHip);
            angHip_r = -angHip;
            
            for m = 1:nt
                rotmatl = roty(angHip_l(m));
                lknee(:,m) = rotmatl*kneevec;
                lankle(:,m) = rotmatl*lankle(:,m);
                ltoe(:,m) = rotmatl*ltoe(:,m);
                lkneeax(:,:,m) = rotmatl;
                lankleax(:,:,m) = rotmatl*lankleax(:,:,m);
                ltoeax(:,:,m) = rotmatl*ltoeax(:,:,m);
                rotmatr = roty(angHip_r(m));
                rknee(:,m) = rotmatr*kneevec;
                rankle(:,m) = rotmatr*rankle(:,m);
                rtoe(:,m) = rotmatr*rtoe(:,m);
                rkneeax(:,:,m) = rotmatr;
                rankleax(:,:,m) = rotmatl*rankleax(:,:,m);
                rtoeax(:,:,m) = rotmatl*rtoeax(:,:,m);
            end
            
            hiplen = bodypartslen(4);
            
            ltoe = bsxfun(@plus,ltoe,[0;hiplen;0]);       % relative to hip
            rtoe = bsxfun(@plus,rtoe,[0;-hiplen;0]);     % relative to hip
            lankle = bsxfun(@plus,lankle,[0;hiplen;0]);   % relative to hip
            rankle = bsxfun(@plus,rankle,[0;-hiplen;0]); % relative to hip
            lknee = bsxfun(@plus,lknee,[0;hiplen;0]);     % relative to hip
            rknee = bsxfun(@plus,rknee,[0;-hiplen;0]);   % relative to hip
        end
        
        function [ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax,lhipax,rhipax] = ...
                rotateLowerBody(obj,rotlr,rottor,ltoe,rtoe,lankle,rankle,lknee,rknee,...
                ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax)
            nt = obj.pNumSamplesPerCycle;
            bodypartslen = obj.pBodyPartsLength;
            lhip = zeros(3,nt);
            rhip = zeros(3,nt);
            lhipax = zeros(3,3,nt);
            rhipax = zeros(3,3,nt);
            
            hiplen = bodypartslen(4);
            lhipvec = [0;hiplen;0];
            rhipvec = [0;-hiplen;0];
            
            for m = 1:nt
                rotmat = rotz(-rottor(m))*rotx(rotlr(m));
                lhip(:,m) = rotmat*lhipvec;
                rhip(:,m) = rotmat*rhipvec;
                
                lknee(:,m) = rotmat*lknee(:,m);
                rknee(:,m) = rotmat*rknee(:,m);
                
                lankle(:,m) = rotmat*lankle(:,m);
                rankle(:,m) = rotmat*rankle(:,m);
                
                ltoe(:,m) = rotmat*ltoe(:,m);
                rtoe(:,m) = rotmat*rtoe(:,m);
                
                lhipax(:,:,m) = rotmat;
                lkneeax(:,:,m) = rotmat*lkneeax(:,:,m);
                lankleax(:,:,m) = rotmat*lankleax(:,:,m);
                ltoeax(:,:,m) = rotmat*ltoeax(:,:,m);
                
                rhipax(:,:,m) = rotmat;
                rkneeax(:,:,m) = rotmat*rkneeax(:,:,m);
                rankleax(:,:,m) = rotmat*rankleax(:,:,m);
                rtoeax(:,:,m) = rotmat*rtoeax(:,:,m);
                
            end
                
            % ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip all relative to
            % hip
        end
        
        function [lhand,rhand,lhandax,rhandax] = flexElbow(obj,angElbow)
            nt = obj.pNumSamplesPerCycle;
            bodypartslen = obj.pBodyPartsLength;
            lhand = zeros(3,nt);
            rhand = zeros(3,nt);
            lhandax = zeros(3,3,nt);
            rhandax = zeros(3,3,nt);
            
            lowerarmlen = bodypartslen(5);
            handvec = [0;0;-lowerarmlen];  % relative to elbow
            
            angElbow_l = fftshift(angElbow); 
            angElbow_r = angElbow;
            
            for m = 1:nt
                rotmatl = roty(-angElbow_l(m));
                lhand(:,m) = rotmatl*handvec;
                lhandax(:,:,m) = rotmatl;
                rotmatr = roty(-angElbow_r(m));
                rhand(:,m) = rotmatr*handvec;
                rhandax(:,:,m) = rotmatr;
            end
            
            shoulderlen = bodypartslen(7);
            upperarmlen = bodypartslen(6);
            torsolen = bodypartslen(10);
            
            lhand = bsxfun(@plus,lhand,[0;shoulderlen;torsolen-upperarmlen]);   % relative to hip
            rhand = bsxfun(@plus,rhand,[0;-shoulderlen;torsolen-upperarmlen]); % relative to hip
            
        end
        
        function [lhand,rhand,lelbow,relbow,lhandax,rhandax,lelbowax,relbowax] = ...
                flexShoulder(obj,angShoulder,lhand,rhand,lhandax,rhandax)
            nt = obj.pNumSamplesPerCycle;
            bodypartslen = obj.pBodyPartsLength;
            lelbow = zeros(3,nt);
            relbow = zeros(3,nt);
            lelbowax = zeros(3,3,nt);
            relbowax = zeros(3,3,nt);
            
            upperarmlen = bodypartslen(6);
            elbowvec = [0;0;-upperarmlen];    % relative to shoulder
            
            torsolen = bodypartslen(10);
            lhand(2,:) = 0;
            lhand(3,:) = lhand(3,:)-torsolen; % relative to shoulder
            rhand(2,:) = 0;
            rhand(3,:) = rhand(3,:)-torsolen; % relative to shoulder
            
            angShoulder_l = fftshift(angShoulder); 
            angShoulder_r = angShoulder;
            
            for m = 1:nt
                rotmatl = roty(-angShoulder_l(m));
                lelbow(:,m) = rotmatl*elbowvec;
                lhand(:,m) = rotmatl*lhand(:,m);
                lelbowax(:,:,m) = rotmatl;
                lhandax(:,:,m) = rotmatl*lhandax(:,:,m);
                
                rotmatr = roty(-angShoulder_r(m));
                relbow(:,m) = rotmatr*elbowvec;
                rhand(:,m) = rotmatr*rhand(:,m);
                relbowax(:,:,m) = rotmatr;
                rhandax(:,:,m) = rotmatr*rhandax(:,:,m);
            end
            
            shoulderlen = bodypartslen(7);
            
            lhand = bsxfun(@plus,lhand,[0;shoulderlen;torsolen]);     % relative to hip
            rhand = bsxfun(@plus,rhand,[0;-shoulderlen;torsolen]);   % relative to hip
            lelbow = bsxfun(@plus,lelbow,[0;shoulderlen;torsolen]);   % relative to hip
            relbow = bsxfun(@plus,relbow,[0;-shoulderlen;torsolen]); % relative to hip
            
        end
        
        function [lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head,...
                lhandax,rhandax,lelbowax,relbowax,lshoulderax,rshoulderax,neckax,headax] = ...
                rotateUpperBody(obj,rotfb,thoraxtor,lhand,rhand,lelbow,relbow,...
                lhandax,rhandax,lelbowax,relbowax)
            nt = obj.pNumSamplesPerCycle;
            bodypartslen = obj.pBodyPartsLength;
            lshoulder = zeros(3,nt);
            rshoulder = zeros(3,nt);
            head = zeros(3,nt);
            neck = zeros(3,nt);
            lshoulderax = zeros(3,3,nt);
            rshoulderax = zeros(3,3,nt);
            neckax = zeros(3,3,nt);
            headax = zeros(3,3,nt);
            
            headlen = bodypartslen(8);
            torsolen = bodypartslen(10);
            headvec = [0;0;torsolen+headlen];         % relative to hip
            neckvec = [0;0;torsolen];                 % relative to hip
            
            shoulderlen = bodypartslen(7);
            lshouldervec = [0;shoulderlen;torsolen];  % relative to hip
            rshouldervec = [0;-shoulderlen;torsolen]; % relative to hip
            
            for m = 1:nt
                rotmat = rotz(-thoraxtor(m))*roty(-rotfb(m));
                head(:,m) = rotmat*headvec;
                neck(:,m) = rotmat*neckvec;
                lshoulder(:,m) = rotmat*lshouldervec;
                rshoulder(:,m) = rotmat*rshouldervec;
                lelbow(:,m) = rotmat*lelbow(:,m);
                relbow(:,m) = rotmat*relbow(:,m);
                lhand(:,m) = rotmat*lhand(:,m);
                rhand(:,m) = rotmat*rhand(:,m);
                
                headax(:,:,m) = rotmat;
                neckax(:,:,m) = rotmat;
                lshoulderax(:,:,m) = rotmat;
                rshoulderax(:,:,m) = rotmat;
                lelbowax(:,:,m) = rotmat*lelbowax(:,:,m);
                relbowax(:,:,m) = rotmat*relbowax(:,:,m);
                lhandax(:,:,m) = rotmat*lhandax(:,:,m);
                rhandax(:,:,m) = rotmat*rhandax(:,:,m);
            end
            
            % lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head all
            % relative to hip
            
        end
        
        function [ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head,lclorigin,lclax] = ...
                translateBody(obj,transfb,translat,transvert,...
                ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head) %#ok<INUSL>
            
            transmat = [transfb;translat;transvert];

            lclorigin = transmat;
            head = head+transmat;
            neck = neck+transmat;
            lshoulder = lshoulder+transmat;
            rshoulder = rshoulder+transmat;
            lelbow = lelbow+transmat;
            relbow = relbow+transmat;
            lhand = lhand+transmat;
            rhand = rhand+transmat;
            lhip = lhip+transmat;
            rhip = rhip+transmat;
            lknee = lknee+transmat;
            rknee = rknee+transmat;
            lankle = lankle+transmat;
            rankle = rankle+transmat;
            ltoe = ltoe+transmat;
            rtoe = rtoe+transmat;
            
            lclax = repmat(eye(3),1,1,size(transmat,2));
        end
        
    end
    
    methods (Hidden, Static)
        
        function bodypartspos_floor = applyVirtualFloor(bodypartspos)
            bodypartspos_floor = bodypartspos;
            bodypartspos_floor(3,:) = bsxfun(@minus,bodypartspos(3,:),min(bodypartspos(3,:)));
        end

        function [transvert,translat,transfb] = getTranslationTrajectory(RV,RDs,t)
            % vertical
            Av = 0.015*RV;
            transvert = -Av+Av*sin(2*pi*(2*t-0.35));
            
            % lateral
            if RV > 0.5
                Al = -0.032;
            else
                Al = -0.128*RV.^2+0.128*RV;
            end
            translat = Al*sin(2*pi*(t-0.1));
            
            % forward/backward
            if RV > 0.5
                Afb = -0.021;
            else
                Afb = -0.084*RV.^2+0.084*RV;
            end
            phifb = 0.625-RDs;
            transfb = Afb*sin(2*pi*(2*t+2*phifb));
        end
        
        function [rotfb,rotlr,rottor] = getPelvisRotationTrajectory(RV,t)
            % forward/backward
            if RV > 0.5
                Afb = 2;
            else
                Afb = -8*RV.^2+8*RV;
            end
            rotfb = -Afb+Afb*sin(2*pi*(2*t-0.1));
            
            % left/right
            Alr = 1.66*RV;
            rotlr = zeros(size(t));
            idx = t<0.15;
            rotlr(idx) = -Alr+Alr*cos(2*pi*(10/3*t(idx)));
            idx = (t>=0.15)&(t<0.5);
            rotlr(idx) = -Alr-Alr*cos(2*pi*(10/7*(t(idx)-0.15)));
            idx = (t>=0.5)&(t<0.65);
            rotlr(idx) = Alr-Alr*cos(2*pi*(10/3*(t(idx)-0.5)));
            idx = t>=0.65;
            rotlr(idx) = Alr+Alr*cos(2*pi*(10/7*(t(idx)-0.65)));
            
            % torsion
            Ator = 4*RV;
            rottor = -Ator*cos(2*pi*t);
        end
        
        function [flexhip,flexknee,flexankle] = getFlexTrajectory(RV,RDs,t)
            
            % gait index based on relative speed, see [1,2] for details
            if RV<0.5         % A=0.5
                gaitidx = 1;
            elseif RV < 1.3   % B=1.3
                gaitidx = 2;
            else              % C=3
                gaitidx = 3;
            end
            
            nt = numel(t);
            t3 = [t-1 t t+1];  % span 3 cycles
            
            % hip
            switch gaitidx
                case 1
                    xhip = [-0.1 0.5 0.9];
                    yhip = [50*RV -30*RV 50*RV];
                    x3hip = [xhip(1:2)-1 xhip xhip(2:3)+1];
                    y3hip = [yhip(1:2) yhip yhip(2:3)];
                case 2
                    xhip = [-0.1 0.5 0.9];
                    yhip = [25 -15 25];
                    x3hip = [xhip(1:2)-1 xhip xhip(2:3)+1];
                    y3hip = [yhip(1:2) yhip yhip(2:3)];
                case 3
                    xhip = [0.2/1.7*(RV-1.3)-0.1 0.5 0.9];
                    yhip = [5/1.7*(RV-1.3)+25 -15 6/1.7*(RV-1.3)+25];
                    x3hip = [xhip-1 xhip xhip+1];
                    y3hip = [yhip yhip yhip];
            end
            temphip = pchip(x3hip,y3hip,t3);
            flexhip = temphip(nt+1:2*nt);
            
            % knee
            switch gaitidx
                case 1
                    xknee = [0.17 0.4 0.75 1];
                    yknee = [3 3 140*RV 3];
                case 2
                    xknee = [0.17 0.4 0.75 1];
                    yknee = [3 3 70 3];
                case 3
                    xknee = [-0.05/1.7*(RV-1.3)+0.17 0.4 -0.05/1.7*(RV-1.3)+0.75 -0.03/1.7*(RV-1.3)+1];
                    yknee = [22/1.7*(RV-1.3)+3 3 -5/1.7*(RV-1.3)+70 3/1.7*(RV-1.3)+3];
            end
            x3knee = [xknee-1 xknee xknee+1];
            y3knee = [yknee yknee yknee];
            tempknee = pchip(x3knee,y3knee,t3);
            flexknee = tempknee(nt+1:2*nt);
                
            % ankle
            switch gaitidx
                case 1
                    xankle = [0 0.08 0.5 RDs 0.85];
                    yankle = [-3 -30*RV-3 22*RV-3 -34*RV-3 -3];
                case 2
                    xankle = [0 0.08 0.5 RDs 0.85];
                    yankle = [-3 -18 8 -20 -3];
                case 3
                    xankle = [0 0.08 -0.1/1.7*(RV-1.3)+0.5 RDs 0.85];
                    yankle = [5/1.7*(RV-1.3)-3 4/1.7*(RV-1.3)-18 -3/1.7*(RV-1.3)+8 -8/1.7*(RV-1.3)-20 5/1.7*(RV-1.3)-3];
            end
            x3ankle = [xankle-1 xankle xankle+1];
            y3ankle = [yankle yankle yankle];
            tempankle = pchip(x3ankle,y3ankle,t3);
            flexankle = tempankle(nt+1:2*nt);
            
        end
        
        function [thoraxtor,flexshoulder,flexelbow] = getUpperBodyTrajectory(RV,t)
            
            % gait index
            if RV<0.5         % A=0.5
                gaitidx = 1;
            elseif RV < 1.3   % B=1.3
                gaitidx = 2;
            else              % C=3
                gaitidx = 3;
            end
            
            nt = numel(t);
            t3 = [t-1 t t+1];  % span 3 cycles
            
            % Torsion of the thorax
            x = [0.1 0.4 0.6 0.9];
            y = [4/3*RV -4.5/3*RV -4/3*RV 4.5/3*RV];
            x3 = [x-1 x x+1];
            y3 = [y y y];
            temptor = pchip(x3,y3,t3);
            thoraxtor = temptor(nt+1:2*nt);
            
            % shoulder
            As = 9.88*RV;
            restval = 3;
            flexshoulder = restval-As./2-As.*cos(2*pi*t);
            
            % elbow
            switch gaitidx
                case 1
                    xelbow = [0.05 0.5 0.9];
                    yelbow = [6*RV+3 34*RV+3 10*RV+3];
                case 2
                    xelbow = [0.05 0.01/0.8*(RV-0.5)+0.5 0.9];
                    yelbow = [8/0.8*(RV-0.5)+6 24/0.8*(RV-0.5)+20 9/0.8*(RV-0.5)+8];
                case 3
                    xelbow = [0.05 0.04/1.7*(RV-1.3)+0.51 -0.1/1.7*(RV-1.3)+0.9];
                    yelbow = [-6/1.7*(RV-1.3)+14  26/1.7*(RV-1.3)+44 -6/1.7*(RV-1.3)+17];
            end
            x3elbow = [xelbow-1 xelbow xelbow+1];
            y3elbow = [yelbow yelbow yelbow];
            tempelbow = pchip(x3elbow,y3elbow,t3);
            flexelbow = tempelbow(nt+1:2*nt);
        end
        
        function bodypartspos = getBodyPartsPosition(bodyjointspos)
            ltoe = bodyjointspos(:,1);
            rtoe = bodyjointspos(:,2);
            lankle = bodyjointspos(:,3);
            rankle = bodyjointspos(:,4);
            lknee = bodyjointspos(:,5);
            rknee = bodyjointspos(:,6);
            lhip = bodyjointspos(:,7);
            rhip = bodyjointspos(:,8);
            lhand = bodyjointspos(:,9);
            rhand = bodyjointspos(:,10);
            lelbow = bodyjointspos(:,11);
            relbow = bodyjointspos(:,12);
            lshoulder = bodyjointspos(:,13);
            rshoulder = bodyjointspos(:,14);
            head = bodyjointspos(:,15);
            neck = bodyjointspos(:,16);
            lclorigin = bodyjointspos(:,17);
            
            lfootc = (ltoe+lankle)/2;
            rfootc = (rtoe+rankle)/2;
            llowerlegc = (lankle+lknee)/2;
            rlowerlegc = (rankle+rknee)/2;
            lupperlegc = (lknee+lhip)/2;
            rupperlegc = (rknee+rhip)/2;
            lhipc = (lhip+lclorigin)/2;
            rhipc = (rhip+lclorigin)/2;
            llowerarmc = (lhand+lelbow)/2;
            rlowerarmc = (rhand+relbow)/2;
            lupperarmc = (lelbow+lshoulder)/2;
            rupperarmc = (relbow+rshoulder)/2;
            lshoulderc = (lshoulder+neck)/2;
            rshoulderc = (rshoulder+neck)/2;
            headc = (head+neck)/2;
            torsoc = (neck+lclorigin)/2;
            
            bodypartspos = [lfootc rfootc llowerlegc rlowerlegc lupperlegc rupperlegc ...
                lhipc rhipc llowerarmc rlowerarmc lupperarmc rupperarmc ...
                lshoulderc rshoulderc headc torsoc];
        end
        
    end
    
    
end
