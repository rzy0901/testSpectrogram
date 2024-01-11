classdef (Sealed) WalkingHuman < phased.internal.AbstractWalkingHuman
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

    properties (Access=private)
        pDuration
    end
    
    methods 
        function obj = WalkingHuman(varargin)
            obj@phased.internal.AbstractWalkingHuman(varargin{:});
        end
        
        function [phs,phl] = plot(obj)
            if isempty(obj.pBodyPartsInitialJointPosition)
                setupHumanScattererPoints(obj);
            end
            pts = obj.pBodyPartsInitialJointPosition;
            ptsc = [1 2 1 2 1 2 1 2 3 4 3 4 3 4 5 5 5];
            phs = scatter3(pts(1,:),pts(2,:),pts(3,:),20,ptsc);
            phl = line(...
                [pts(1,[1 3 5 7 2 4 6 8 9 11 13 10 12 14 17 15]);...
                pts(1,[3 5 7 17 4 6 8 17 11 13 16 12 14 16 16 16])],...
                [pts(2,[1 3 5 7 2 4 6 8 9 11 13 10 12 14 17 15]);...
                pts(2,[3 5 7 17 4 6 8 17 11 13 16 12 14 16 16 16])],...
                [pts(3,[1 3 5 7 2 4 6 8 9 11 13 10 12 14 17 15]);...
                pts(3,[3 5 7 17 4 6 8 17 11 13 16 12 14 16 16 16])],'Color','b');
            axis equal
            view(90,0)
        end
        
        function [phs,phl] = animate(obj,nsteps)
            if isempty(obj.pBodyPartsInitialJointPosition)
                setupHumanScattererPoints(obj);
                setupWalkingCharacteristics(obj);
            end
            
            nspercyc = obj.pNumSamplesPerCycle;
            tcyc = (0:nspercyc-1)/nspercyc; % 1 cycle
            Dc = obj.pDurationCycle;
            dt = Dc/nspercyc;
            ns = nspercyc*nsteps;
            t = (0:ns-1)*dt;
            
            pos = [obj.WalkingSpeed*t;zeros(2,ns)];
            
            RV = obj.pRelativeSpeed;
            RDs = obj.pRelativeDurationSupport;
            
            % compute joint trajectory
            [transvert,translat,transfb] = phased.internal.WalkingHuman.getTranslationTrajectory(RV,RDs,tcyc);
            [rotfb,rotlr,rottor] = phased.internal.WalkingHuman.getPelvisRotationTrajectory(RV,tcyc);
            [flexhip,flexknee,flexankle] = phased.internal.WalkingHuman.getFlexTrajectory(RV,RDs,tcyc);
            [thoraxtor,flexshoulder,flexelbow] = phased.internal.WalkingHuman.getUpperBodyTrajectory(RV,tcyc);
            
            % compute body part trajectory
            [ltoe,rtoe,ltoeax,rtoeax] = flexAnkle(obj,flexankle);
            [ltoe,rtoe,lankle,rankle,ltoeax,rtoeax,lankleax,rankleax] = ...
                flexKnee(obj,flexknee,ltoe,rtoe,ltoeax,rtoeax);
            [ltoe,rtoe,lankle,rankle,lknee,rknee,ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax] = ...
                flexHip(obj,flexhip,ltoe,rtoe,lankle,rankle,ltoeax,rtoeax,lankleax,rankleax);
            [ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip] = ...
                rotateLowerBody(obj,rotlr,rottor,ltoe,rtoe,lankle,rankle,lknee,rknee,...
                ltoeax,rtoeax,lankleax,rankleax,lkneeax,rkneeax);
            [lhand,rhand,lhandax,rhandax] = flexElbow(obj,flexelbow);
            [lhand,rhand,lelbow,relbow,lhandax,rhandax,lelbowax,relbowax] = ...
                flexShoulder(obj,flexshoulder,lhand,rhand,lhandax,rhandax);
            [lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head] = ...
                rotateUpperBody(obj,rotfb,thoraxtor,lhand,rhand,lelbow,relbow,...
                lhandax,rhandax,lelbowax,relbowax);
            [ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head,lclorigin] = translateBody(obj,...
                transfb,translat,transvert,ltoe,rtoe,lankle,rankle,lknee,rknee,lhip,rhip,...
                lhand,rhand,lelbow,relbow,lshoulder,rshoulder,neck,head);
            
            % combine with macro body movement
            ltoepos = repmat(ltoe,1,nsteps)+pos;
            rtoepos = repmat(rtoe,1,nsteps)+pos;
            lanklepos = repmat(lankle,1,nsteps)+pos;
            ranklepos = repmat(rankle,1,nsteps)+pos;
            lkneepos = repmat(lknee,1,nsteps)+pos;
            rkneepos = repmat(rknee,1,nsteps)+pos;
            lhippos = repmat(lhip,1,nsteps)+pos;
            rhippos = repmat(rhip,1,nsteps)+pos;
            lhandpos = repmat(lhand,1,nsteps)+pos;
            rhandpos = repmat(rhand,1,nsteps)+pos;
            lelbowpos = repmat(lelbow,1,nsteps)+pos;
            relbowpos = repmat(relbow,1,nsteps)+pos;
            lshoulderpos = repmat(lshoulder,1,nsteps)+pos;
            rshoulderpos = repmat(rshoulder,1,nsteps)+pos;
            neckpos = repmat(neck,1,nsteps)+pos;
            headpos = repmat(head,1,nsteps)+pos;
            lcloriginpos = repmat(lclorigin,1,nsteps)+pos;
            
            % form data matrix
            bodypartspos = cat(3,ltoepos,rtoepos,lanklepos,ranklepos,lkneepos,rkneepos,...
                lhippos,rhippos,lhandpos,rhandpos,lelbowpos,relbowpos,lshoulderpos,rshoulderpos,...
                headpos,neckpos,lcloriginpos);
            
            % virtual floor
            bodypartspos(3,:,:) = bodypartspos(3,:,:)-min(bodypartspos(3,:,:),[],3);
            
            [phs,phl] = plot(obj);
            xlim([-1,2*nsteps]);
            ylim([-1 1]);
            zlim([0 2]);
            view(36,15);
            
            bodysegidx = [1 3 5 7 2 4 6 8 9 11 13 10 12 14 17 15; ...
                3 5 7 17 4 6 8 17 11 13 16 12 14 16 16 16];
            
            for m = 1:ns
                currentframe = squeeze(bodypartspos(:,m,:));
                phs.XData = currentframe(1,:);
                phs.YData = currentframe(2,:);
                phs.ZData = currentframe(3,:);
                
                
                for n = 1:numel(phl)
                    phl(n).XData = currentframe(1,bodysegidx(:,n));
                    phl(n).YData = currentframe(2,bodysegidx(:,n));
                    phl(n).ZData = currentframe(3,bodysegidx(:,n));
                end
                drawnow;
            end
            
        end
    end
    
    methods (Access=protected)
        function validateInputsImpl(obj,dt,heading)
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
                num = 1;
            else
                num = 2;
            end
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractWalkingHuman(obj);
            obj.pDuration = 0; % intialize for codegen
        end


        function [bodypartspos,bodypartsvel,bodypartsax,bodyjointspos] = stepImpl(obj,dt,heading)
            cond = (dt<0);
            if cond
                coder.internal.errorIf(cond,'phased:step:expectedNonnegative','T');
            end
            [bodypartspos,bodypartsvel,bodypartsax,bodyjointspos] = walking(obj,obj.pDuration,heading);
            obj.pCurrentTime = obj.pCurrentTime+dt;
            obj.pDuration = dt;
            phi = obj.pHeading;
            v = obj.WalkingSpeed;
            obj.pPosition = obj.pPosition+v*dt*[cosd(phi);sind(phi);0];
            obj.pOrientationAxes = rotz(phi);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractWalkingHuman(obj);
            s.pDuration = obj.pDuration;
        end
        
        function loadObjectImpl(obj,s,~)
            % Set properties in object obj to values in structure s
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

    end
    

    
end
