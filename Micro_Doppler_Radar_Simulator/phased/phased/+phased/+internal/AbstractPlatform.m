classdef (Hidden) AbstractPlatform < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.Nondirect & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%   Copyright 2016-2018 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable) 

        %MotionModel    Model of object motion
        %   Specify the object motion model as one of 'Velocity' |
        %   'Acceleration' | 'Custom', where the default is 'Velocity'.
        %   When you set this property to 'Velocity', the object motion is
        %   assumed to be constant velocity during each simulation step.
        %   When you set this property to 'Acceleration', the object motion
        %   is assumed to be constant acceleration during each simulation
        %   step. When you set the property to 'Custom', the object motion
        %   is specified by a series of waypoints in the CustomTrajectory
        %   property.
        MotionModel = 'Velocity'
        %InitialPosition    Initial position (m)
        %   Specify the initial position of the platforms as a 3xN vector
        %   in the form of [x; y; z] (in meters). N is the number of
        %   platforms to model. The default value of this property is [0;
        %   0; 0].
        InitialPosition = [0; 0; 0]
        %InitialVelocity    Initial velocity (m/s)
        %   Specify the initial velocity of the platforms as a 3xN vector
        %   in the form of [Vx; Vy; Vz] (in meters/second). N is the number
        %   of platforms to model. The default value of this property is
        %   [0; 0; 0]. This property only applies when you set the
        %   MotionModel property to 'Velocity' and the VelocitySource to
        %   'Input port', or when you set the MotionModel property to
        %   'Acceleration'.
        InitialVelocity = [0; 0; 0]
        %VelocitySource     Source of velocity
        %   Specify how the velocity is specified as one of 'Property' |
        %   'Input port', where the default is 'Property'. When you set
        %   this property to 'Property', the velocity is specified by the
        %   Velocity property. When you set this property to 'Input port',
        %   the velocity is specified as an input argument. This property
        %   applies when you set the MotionModel property to 'Velocity'.
        VelocitySource = 'Property'
    end
        
    properties (Dependent)
        %Velocity   Velocity (m/s)
        %   Specify the current velocity of the platforms as a 3xN vector
        %   in the form of [Vx; Vy; Vz] (in meters/second). N is the number
        %   of platforms to model. This property is tunable. The default
        %   value of this property is [0; 0; 0]. This property applies when
        %   you set the MotionModel property to 'Velocity' and the
        %   VelocitySource to 'Property'.
        Velocity
    end
    
    properties (Nontunable)
        %AccelerationSource     Source of acceleration
        %   Specify how the acceleration is specified as one of 'Property'
        %   | 'Input port', where the default is 'Property'. When you set
        %   this property to 'Property', the acceleration is specified by
        %   the Acceleration property. When you set this property to 'Input
        %   port', the acceleration is specified as an input argument. This
        %   property applies when you set the MotionModel property to
        %   'Acceleration'.
        AccelerationSource = 'Property'
    end
    
    properties (Dependent)
        %Acceleration   Acceleration (m/s^2)
        %   Specify the current acceleration of the platforms as a 3xN
        %   vector in the form of [Ax; Ay; Az] (in meters/second^2). N is
        %   the number of platforms to model. This property is tunable. The
        %   default value of this property is [0; 0; 0]. This property
        %   applies when you set the MotionModel property to 'Acceleration'
        %   and the AccelerationSource to 'Property'.
        Acceleration
    end
    
    properties (Nontunable)
        %CustomTrajectory   Custom trajectory waypoints
        %   Specify the custom trajectory waypoints as a 4-column or
        %   7-column matrix or 3-dimensional array. The first column
        %   indicates the time incidences, (in seconds), where the platform
        %   position gets measured. Note that the first value in the time
        %   vector must be 0, indicating the start of simulation time.
        %
        %   The 2nd through 4th columns are position measurements (in
        %   meters) in x, y, and z coordinates. The 5th through 7th columns
        %   in the matrix are velocity measurements (in meters/second) in x,
        %   y, and z coordinates. If only 4 columns are specified, the
        %   velocity is derived from the position measurements.
        %
        %   When you set the CustomTrajectory property to a 3-dimensional
        %   array, the number of pages, N, indicates the number of
        %   platforms.
        CustomTrajectory = [0 0 0 0 0 0 0;1 0 0 0 0 0 0]
    end

    properties (Nontunable)
        %ScanMode   Mechanical scanning mode
        %   Specify the mechanical scanning mode for the platform as one of
        %   'None' | 'Circular' | 'Sector', where 'None' is the default.
        %   When you set the ScanMode property to 'Circular', the platform
        %   scan clockwise 360 degrees continuously along azimuth direction
        %   of the platform's orientation axes. When you set the ScanMode
        %   property to 'Sector', the platform starts clockwise scanning
        %   along the azimuth angles in the platform's orientation axes
        %   within a range specified by the AzimuthSpan property. When the
        %   platform reaches the boundary of the span, it reverses the
        %   direction and scan back.
        %
        %   Note that all scanning happens within the orientation axes of
        %   the platform. 
        ScanMode = 'None'
        %InitialScanAngle   Initial scan angle (deg)
        %   Specify the initial scan angle (in degrees) as a 1xN vector
        %   where N is the number of platforms. The default value is 0.
        %   This property applies when you set the ScanMode to either
        %   'Circular' or 'Sector'.
        %
        %   The scanning occurs in the local coordinate system of the
        %   platform. The InitialOrientationAxes specifies the original
        %   local coordinate system. At the beginning of the simulation,
        %   the orientation axes specified in the orientation axes is
        %   rotated to the corresponding angle specified in the
        %   InitialScanAngle property.
        InitialScanAngle = 0
        %AzimuthSpan    Azimuth scan angle span (deg)
        %   Specify the azimuth angle span as an Nx2 matrix where N is the
        %   number of platforms. Each row of the matrix specifies the scan
        %   range of the corresponding platform in the form of
        %   [ScanAngleLowerBound ScanAngleHigherBound] (in degrees). The
        %   default value is [-60 60]. This property applies when you set
        %   the ScanMode to 'Sector'.
        AzimuthSpan = [-60 60]
        %AzimuthScanRate    Azimuth scan rate (deg/s)
        %   Specify the azimuth scan rate (in degrees/second) as a 1xN
        %   vector where N is the number of platforms. Each entry in the
        %   vector is the azimuth scan rate for the corresponding platform.
        %   The default value is 10 degrees/second. This property applies
        %   when you set the ScanMode property to 'Circular' or 'Sector'.
        AzimuthScanRate = 10
        %ElevationRate = 10  % raster
        %ElevationStep       % raster
    end
    
    properties (Hidden, Dependent)
        %OrientationAxes   Orientation axes
        %   Specify the 3 axes that defines the local (x, y, z) coordinate
        %   system for the N platforms as a 3x3xN matrix (each column
        %   corresponds to an axis). The 3 axes of each platform must be
        %   orthonormal. It is sufficient to specify only one set of axes
        %   as a 3x3 matrix when the local coordinate systems of the N
        %   platforms coincide. The default value of this property is [1 0
        %   0; 0 1 0; 0 0 1]. This property applies when you set the
        %   OrientationAxesOutputPort property to true.
        OrientationAxes 
    end
    
    properties (Nontunable)
        %InitialOrientationAxes   Initial orientation axes
        %   Specify the 3 axes that defines the local (x, y, z) coordinate
        %   system for the N platforms as a 3x3xN matrix (each column
        %   corresponds to an axis). The 3 axes of each platform must be
        %   orthonormal. It is sufficient to specify only one set of axes
        %   as a 3x3 matrix when the local coordinate systems of the N
        %   platforms coincide. The default value of this property is [1 0
        %   0; 0 1 0; 0 0 1]. 
        InitialOrientationAxes = eye(3)
    end
    
    properties (Nontunable, Logical) 
        %OrientationAxesOutputPort  Enable orientation axes output
        %   Set this property to true to output the orientation axes of the
        %   platform. Set this property to false to not output the
        %   orientation axes of the platform. The default value of this
        %   property is false.
        OrientationAxesOutputPort = false;
    end
    
    properties(Constant, Hidden)
        MotionModelSet = matlab.system.StringSet(...
            {'Velocity','Acceleration','Custom'});
        VelocitySourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        AccelerationSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        ScanModeSet = matlab.system.StringSet({'None','Circular',...
            'Sector'});  % 'Raster','Helical'
    end
    
    properties (Access = protected)
        pPosition
        pVelocity  = [0; 0; 0]
        pAcceleration = [0; 0; 0]
        pOrientationAxes  
        pPreviousVelocity
        pScanAngle
        pCurrentScanDirection
        pPPx
        pPPy
        pPPz
        pPPvx
        pPPvy
        pPPvz
        pCurrentTime
        pTrajectoryMaxTime
    end
    
    properties (Access = protected, Nontunable)
        pMotionModelIndex 
        pIsVelocityViaInput = false
        pIsAccelerationViaInput = false
        pScanModeIndex
    end
    
    methods
        function set.InitialPosition(obj,val)
            sigdatatypes.validate3DCartCoord(val,...
                'phased.Platform','InitialPosition');
            obj.InitialPosition = val;
        end
        
        function set.InitialVelocity(obj,val)
            sigdatatypes.validate3DCartCoord(val,...
                'phased.Platform','InitialVelocity');
            obj.InitialVelocity = val;
        end
        
        function set.Velocity(obj,val)
            sigdatatypes.validate3DCartCoord(val,...
                'phased.Platform','Velocity');
            obj.pVelocity = val;
        end
        
        function val = get.Velocity(obj)
            val = obj.pVelocity;
        end
              
        function set.Acceleration(obj,val)
            sigdatatypes.validate3DCartCoord(val,...
                'phased.Platform','Acceleration');
            obj.pAcceleration = val;
        end
        
        function val = get.Acceleration(obj)
            val = obj.pAcceleration;
        end
        
        function set.InitialOrientationAxes(obj,val)
            validateattributes(val,{'double'},{'finite','nonnan','nonempty','real',...
                                '3d'},'phased.Platform','InitialOrientationAxes');

            cond =  size(val,1) ~= 3 || size(val,2) ~= 3;
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:Platform:OrientAxisCheckDims','InitialOrientationAxes');
            end
            %Check that all Orientation axis are orthonormal
            numOrient = size(val,3);
            for orientId = 1:numOrient
                curVal = val(:,:,numOrient);
                cond =  norm(curVal'*curVal-eye(3)) > sqrt(eps);
                if cond
                    coder.internal.errorIf(cond, ...
                      'phased:Platform:NonOrthonormal');
                end
            end
            obj.InitialOrientationAxes = val;
        end
        
        function set.OrientationAxes(obj,val)
            obj.InitialOrientationAxes = val;
        end
        
        function val = get.OrientationAxes(obj)
            val = obj.InitialOrientationAxes;
        end
        
        function set.InitialScanAngle(obj,val)
            sigdatatypes.validateAngle(val,'phased.Platform','ScanAngle',...
                {'row','>=',-180,'<=',180});
            obj.InitialScanAngle = val;
        end
        
        function set.AzimuthSpan(obj,val)
            sigdatatypes.validateAngle(val,'phased.Platform','AzimuthSpan',...
                {'2d','ncols',2,'>=',-180,'<=',180});
            cond = any((val(:,2)-val(:,1))<=0);
            if cond
                coder.internal.errorIf(cond, ...
                  'phased:phased:expectedIncreasingRows','AzimuthSpan');
            end
            obj.AzimuthSpan = val;
        end
        
        function set.AzimuthScanRate(obj,val)
            sigdatatypes.validateSpeed(val,'phased.Platform','AzimuthScanRate',...
                {'row','positive'});
            obj.AzimuthScanRate = val;
        end
        
        function set.CustomTrajectory(obj,val)
            validateattributes(val,{'double'},{'3d','real','finite','nonempty','nonnan'},...
                'phased.Platform','Trajectory');
            cond = size(val,1)<2;
            if cond
                coder.internal.errorIf(cond,'phased:phased:expectedMoreRows','CustomTrajectory','2');
            end
            cond = (size(val,2)~=4 && size(val,2)~=7);
            if cond
                coder.internal.errorIf(cond,'phased:phased:invalidColumnDimension','CustomTrajectory','4','7');
            end
            cond = any(val(1,1,:)~=0);
            if cond
                coder.internal.errorIf(cond,'phased:Platform:StartingTimeMustBeZero','CustomTrajectory');
            end
            cond = any(diff(squeeze(val(:,1,:)),1,1)<=0);
            if cond
                coder.internal.errorIf(cond,'phased:Platform:ExpectIncreasingTime','CustomTrajectory');
            end
            obj.CustomTrajectory = val;
        end
    end
    
    methods (Access = protected, Abstract)
        T = getElapsedTime(obj,varargin)
    end
            

    methods (Access = protected)

        function obj = AbstractPlatform(varargin)
            %Platform   Construct the Platform class.
            setProperties(obj, nargin, varargin{:}, 'InitialPosition', 'Velocity');
        end
        
    end
    
    methods (Access = protected)
        
        function setupImpl(obj,varargin)
            if strcmp(obj.MotionModel,'Velocity')
                obj.pMotionModelIndex = 1;
                obj.pIsVelocityViaInput = ...
                    strcmp(obj.VelocitySource,'Input port');
                nplat = size(obj.InitialPosition,2);
            elseif strcmp(obj.MotionModel,'Acceleration')
                obj.pMotionModelIndex = 2;
                obj.pIsAccelerationViaInput = ...
                    strcmp(obj.AccelerationSource,'Input port');
                nplat = size(obj.InitialPosition,2);
            else % Custom
                obj.pMotionModelIndex = 3;
                traj = obj.CustomTrajectory;
                nplat = size(traj,3);
                tmax = squeeze(traj(end,1,:));
                obj.pTrajectoryMaxTime = tmax;
                
                if size(traj,2) == 4
                    [PPx,PPvx] = phased.internal.motionhermite(traj(:,1,1)/tmax(1),traj(:,2,1));
                    [PPy,PPvy] = phased.internal.motionhermite(traj(:,1,1)/tmax(1),traj(:,3,1));
                    [PPz,PPvz] = phased.internal.motionhermite(traj(:,1,1)/tmax(1),traj(:,4,1));
                    obj.pPPx = repmat(PPx,1,nplat);
                    obj.pPPy = repmat(PPy,1,nplat);
                    obj.pPPz = repmat(PPz,1,nplat);
                    obj.pPPvx = repmat(PPvx,1,nplat);
                    obj.pPPvy = repmat(PPvy,1,nplat);
                    obj.pPPvz = repmat(PPvz,1,nplat);
                    for m =  2:nplat
                        [obj.pPPx(m),obj.pPPvx(m)] = phased.internal.motionhermite(traj(:,1,m)/tmax(m),traj(:,2,m));
                        [obj.pPPy(m),obj.pPPvy(m)] = phased.internal.motionhermite(traj(:,1,m)/tmax(m),traj(:,3,m));
                        [obj.pPPz(m),obj.pPPvz(m)] = phased.internal.motionhermite(traj(:,1,m)/tmax(m),traj(:,4,m));
                    end
                else
                    [PPx,PPvx] = phased.internal.motionhermite(traj(:,1,1)/tmax(1),traj(:,2,1),traj(:,5,1));
                    [PPy,PPvy] = phased.internal.motionhermite(traj(:,1,1)/tmax(1),traj(:,3,1),traj(:,6,1));
                    [PPz,PPvz] = phased.internal.motionhermite(traj(:,1,1)/tmax(1),traj(:,4,1),traj(:,7,1));
                    obj.pPPx = repmat(PPx,1,nplat);
                    obj.pPPy = repmat(PPy,1,nplat);
                    obj.pPPz = repmat(PPz,1,nplat);
                    obj.pPPvx = repmat(PPvx,1,nplat);
                    obj.pPPvy = repmat(PPvy,1,nplat);
                    obj.pPPvz = repmat(PPvz,1,nplat);
                    for m = 2:nplat
                        [obj.pPPx(m),obj.pPPvx(m)] = phased.internal.motionhermite(traj(:,1,m)/tmax(m),traj(:,2,m),traj(:,5,m));
                        [obj.pPPy(m),obj.pPPvy(m)] = phased.internal.motionhermite(traj(:,1,m)/tmax(m),traj(:,3,m),traj(:,6,m));
                        [obj.pPPz(m),obj.pPPvz(m)] = phased.internal.motionhermite(traj(:,1,m)/tmax(m),traj(:,4),m,traj(:,7,m));
                    end
                end
            end
            
            if strcmp(obj.ScanMode,'None')
                obj.pScanModeIndex = 0;
            elseif strcmp(obj.ScanMode,'Circular')
                obj.pScanModeIndex = 1;
                obj.pScanAngle = [obj.InitialScanAngle;zeros(1,nplat)];
            else % Sector
                obj.pScanModeIndex = 2;
                obj.pCurrentScanDirection = -ones(1,nplat);
                obj.pScanAngle = [obj.InitialScanAngle;zeros(1,nplat)];
            end
        end
        
        function initializeOrientationAxes(obj)
            if size(obj.InitialPosition,2) ~= ...
                    size(obj.InitialOrientationAxes,3)
                obj.pOrientationAxes = repmat(obj.InitialOrientationAxes,...
                    1,1,size(obj.InitialPosition,2));
            else
                obj.pOrientationAxes = obj.InitialOrientationAxes;
            end
            if ~strcmp(obj.ScanMode,'None')
                initang = obj.InitialScanAngle;
                ax = obj.pOrientationAxes;
                for m = 1:size(ax,3)
                    ax(:,:,m) = rotz(initang(1,m))*ax(:,:,m);
                end
                obj.pOrientationAxes = ax;
            end
        end
        
        function updateOrientationAxes(obj,vel)
            current_vel = obj.pPreviousVelocity;
            ax = obj.pOrientationAxes;
            for m = 1:size(ax,3)
                ax(:,:,m) = rotvv(current_vel(:,m),vel(:,m))*ax(:,:,m);
                for n = 1:3
                    % ensure norm for each vector to be 1 to preserve
                    % orthonormal property, avoid cumulative error
                    ax(:,n,m) = ax(:,n,m)/norm(ax(:,n,m));
                end
            end
            obj.pOrientationAxes = ax;
        end
        
        function updateScanAngle(obj,t)
            scanidx = obj.pScanModeIndex;
            switch scanidx
                case 1
                    w = obj.AzimuthScanRate;
                    ang_current = obj.pScanAngle(1,:);
                    deltaang = -rem(w*t,360);
                    ang = ang_current+deltaang;
                    obj.pScanAngle(1,:) = ang;
                    pax = obj.pOrientationAxes;
                    for m = 1:size(deltaang,2)
                        pax(:,:,m) = rotz(deltaang(m))*pax(:,:,m);
                    end
                    obj.pOrientationAxes = pax;
                case 2
                    w = obj.AzimuthScanRate;
                    ang_current = obj.pScanAngle(1,:);
                    angspan = diff(obj.AzimuthSpan,1,2);
                    scandir = obj.pCurrentScanDirection;
                    if scandir > 0  % Counterclockwise
                        azlb = obj.AzimuthSpan(:,1).';
                        azhb = obj.AzimuthSpan(:,2).';
                        ang = ang_current+w*t;
                        
                        resanghb = rem(ang-azlb,angspan);
                        resanghb2 = rem(ang-azlb,angspan*2);
                        
                        noswitchidx = resanghb==resanghb2;
                        ang(noswitchidx) = resanghb(noswitchidx)+azlb(noswitchidx);
                        ang(~noswitchidx) = azhb(~noswitchidx)-resanghb(~noswitchidx);
                        scandir(~noswitchidx) = -scandir(~noswitchidx);
                        
                        deltaang = ang-ang_current;
                    else  % Clockwise
                        azlb = obj.AzimuthSpan(:,2).';
                        azhb = obj.AzimuthSpan(:,1).';
                        ang = ang_current-w*t;
                        
                        resanghb = rem(abs(ang-azlb),angspan);
                        resanghb2 = rem(abs(ang-azlb),angspan*2);
                        
                        noswitchidx = resanghb==resanghb2;
                        ang(noswitchidx) = azlb(noswitchidx)-resanghb(noswitchidx);
                        ang(~noswitchidx) = azhb(~noswitchidx)+resanghb(~noswitchidx);
                        scandir(~noswitchidx) = -scandir(~noswitchidx);
                        
                        deltaang = ang-ang_current;
                    end
                    obj.pScanAngle(1) = ang;
                    obj.pCurrentScanDirection = scandir;
                    obj.pOrientationAxes = rotz(deltaang)*obj.pOrientationAxes;
            end
            
        end
            

        function [pos,vel,ax] = outputImpl(obj,varargin) 
            
            pos = obj.pPosition;
            vel = obj.pVelocity;
            ax = obj.pOrientationAxes;
            
        end
        
        function updateImpl(obj,varargin)
            t = getElapsedTime(obj,varargin{:});
            switch obj.pMotionModelIndex
                case 1 %Velocity
                    obj.pPosition = obj.pPosition + ...
                        obj.pVelocity * t;
                case 2 % Acceleration
                    obj.pPosition = obj.pPosition + ...
                        obj.pVelocity*t + obj.pAcceleration*t^2/2;
                    obj.pVelocity = obj.pVelocity + obj.pAcceleration*t;
                case 3
                    tc = obj.pCurrentTime+t;
                    tmax = obj.pTrajectoryMaxTime;
                    traj = obj.CustomTrajectory;
                    obj.pCurrentTime = tc;
                    PPx = obj.pPPx;
                    PPy = obj.pPPy;
                    PPz = obj.pPPz;
                    PPvx = obj.pPPvx;
                    PPvy = obj.pPPvy;
                    PPvz = obj.pPPvz;
                    for m = 1:size(obj.CustomTrajectory,3)
                        posx = ppval(PPx(m),tc/tmax(m));
                        posy = ppval(PPy(m),tc/tmax(m));
                        posz = ppval(PPz(m),tc/tmax(m));
                        if size(traj,2) == 4
                            posvx = ppval(PPvx(m),tc/tmax(m))/tmax(m);
                            posvy = ppval(PPvy(m),tc/tmax(m))/tmax(m);
                            posvz = ppval(PPvz(m),tc/tmax(m))/tmax(m);
                        else
                            posvx = ppval(PPvx(m),tc/tmax(m));
                            posvy = ppval(PPvy(m),tc/tmax(m));
                            posvz = ppval(PPvz(m),tc/tmax(m));
                        end
                        obj.pPosition(:,m) = [posx;posy;posz];
                        obj.pVelocity(:,m) = [posvx;posvy;posvz];
                    end
            end
            updateScanAngle(obj,t);
            updateOrientationAxes(obj,obj.pVelocity);
            obj.pPreviousVelocity = obj.pVelocity;
        end
        
        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end

        function resetImpl(obj)
            obj.pCurrentTime = 0;
            if strcmp(obj.MotionModel,'Velocity')
                obj.pPosition = obj.InitialPosition;
                obj.pAcceleration = [0;0;0];
                if strcmp(obj.VelocitySource,'Input port')
                    obj.pVelocity = obj.InitialVelocity;
                    obj.pPreviousVelocity = obj.InitialVelocity;
                else
                    obj.pPreviousVelocity = obj.pVelocity;
                end
            elseif strcmp(obj.MotionModel,'Acceleration') 
                obj.pPosition = obj.InitialPosition;
                obj.pVelocity = obj.InitialVelocity;
                obj.pPreviousVelocity = obj.InitialVelocity;
            else %Custom
                traj = obj.CustomTrajectory;
                if ismatrix(traj)
                    obj.pPosition = traj(1,2:4).';
                else
                    obj.pPosition = squeeze(traj(1,2:4,:));
                end
                if size(traj,2) == 4
                    tmax = obj.pTrajectoryMaxTime;
                    Nptraj = numel(obj.pPPvx);
                    vx = zeros(1,Nptraj);
                    vy = zeros(1,Nptraj);
                    vz = zeros(1,Nptraj);
                    for m = 1:numel(obj.pPPvx)
                        vx(m) = ppval(obj.pPPvx(m),0)/tmax(m);
                        vy(m) = ppval(obj.pPPvy(m),0)/tmax(m);
                        vz(m) = ppval(obj.pPPvz(m),0)/tmax(m);
                    end
                    obj.pVelocity = [vx;vy;vz];
                    obj.pPreviousVelocity = [vx;vy;vz];
                else
                    if ismatrix(traj)
                        obj.pVelocity = traj(1,5:7).';
                        obj.pPreviousVelocity = traj(1,5:7).';
                    else
                        obj.pVelocity = squeeze(traj(1,5:7,:));
                        obj.pPreviousVelocity = squeeze(traj(1,5:7,:));
                    end
                end
            end
            if ~strcmp(obj.ScanMode,'None')
                obj.pScanAngle(1,:) = obj.InitialScanAngle;
            end
            initializeOrientationAxes(obj);
        end       
        
        function validatePropertiesImpl(obj)
            coder.extrinsic('num2str');   
            if strcmp(obj.MotionModel,'Velocity') 
                if strcmp(obj.VelocitySource,'Property')
                    cond = ~isequal(size(obj.InitialPosition),size(obj.Velocity));
                    if cond
                        coder.internal.errorIf(cond,'phased:Platform:ExpectedSameDim',...
                            'InitialPosition','Velocity');
                    end
                else % velocity input port
                    cond = ~isequal(size(obj.InitialPosition),size(obj.InitialVelocity));
                    if cond
                        coder.internal.errorIf(cond,'phased:Platform:ExpectedSameDim',...
                            'InitialPosition','InitialVelocity');
                    end
                end
                numOfPlatforms = size(obj.InitialPosition,2);                        
            elseif strcmp(obj.MotionModel,'Acceleration') 
                cond = ~isequal(size(obj.InitialPosition),size(obj.InitialVelocity));
                if cond
                    coder.internal.errorIf(cond,'phased:Platform:ExpectedSameDim',...
                        'InitialPosition','InitialVelocity');
                end
                if strcmp(obj.AccelerationSource,'Property')
                    cond = ~isequal(size(obj.InitialPosition),size(obj.Acceleration));
                    if cond
                        coder.internal.errorIf(cond,'phased:Platform:ExpectedSameDim',...
                            'InitialPosition','Acceleration');
                    end
                end
                numOfPlatforms = size(obj.InitialPosition,2); 
            else  % Custom
                numOfPlatforms = size(obj.CustomTrajectory,3);
            end
            
            cond =  numOfPlatforms ~= 1 && ...
                    size(obj.InitialOrientationAxes,3) ~= size(obj.InitialPosition,2) && ...
                    size(obj.InitialOrientationAxes,3) ~= 1;
            if cond
                coder.internal.errorIf(cond,'phased:Platform:OrientAxisCheckAllDims',...
                    'InitialOrientationAxes',['[3 3] or [3 3 ' coder.const(num2str(numOfPlatforms)) ']']);
            end
            
            if strcmp(obj.ScanMode,'Circular')
                cond = ~isequal(size(obj.InitialScanAngle,2),numOfPlatforms);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidColumnNumbers',...
                        'InitialScanAngle',coder.const(num2str(numOfPlatforms)));
                end
                cond = ~isequal(size(obj.AzimuthScanRate,2),numOfPlatforms);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidColumnNumbers',...
                        'AzimuthScanRate',coder.const(num2str(numOfPlatforms)));
                end
            elseif strcmp(obj.ScanMode,'Sector')
                cond = ~isequal(size(obj.InitialScanAngle,2),numOfPlatforms);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidColumnNumbers',...
                        'InitialScanAngle',coder.const(num2str(numOfPlatforms)));
                end
                cond = ~isequal(size(obj.AzimuthScanRate,2),numOfPlatforms);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidColumnNumbers',...
                        'AzimuthScanRate',coder.const(num2str(numOfPlatforms)));
                end
                cond = ~isequal(size(obj.AzimuthSpan,1),numOfPlatforms);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidRowNumbers',...
                        'AzimuthSpan',coder.const(num2str(numOfPlatforms)));
                end
                cond = any(obj.InitialScanAngle<obj.AzimuthSpan(:,1).') || ...
                    any(obj.InitialScanAngle>obj.AzimuthSpan(:,2).');
                if cond
                    coder.internal.errorIf(cond,'phased:Platform:OutOfBound',...
                        'InitialScanAngle','AzimuthSpan');
                end
            end
        end
        
        function flag = isInputComplexityLockedImpl(~,~) 
            flag = true;  % (index == 1) || (index == 2) || (index == 3)
        end
        
        function flag = isOutputComplexityLockedImpl(~,~) 
            flag = true;  % (index == 1) || (index == 2) || (index == 3)
        end
        
        function num = getNumOutputsImpl(obj) 
            if obj.OrientationAxesOutputPort
                num = 3;
            else
                num = 2;
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if strcmp(obj.MotionModel,'Acceleration') && ...
                    (strcmp(prop, 'VelocitySource') || strcmp(prop, 'Velocity'))
                flag = true;
            elseif strcmp(obj.MotionModel,'Acceleration') && ...
                    strcmp(obj.AccelerationSource,'Input port') && strcmp(prop, 'Acceleration')
                flag = true;
            elseif strcmp(obj.MotionModel,'Velocity') && ...
                    (strcmp(prop, 'AccelerationSource') || strcmp(prop, 'Acceleration'))
                flag = true;
            elseif strcmp(obj.MotionModel,'Velocity') && ...
                    strcmp(obj.VelocitySource,'Input port') && strcmp(prop, 'Velocity')
                flag = true;
            elseif strcmp(obj.MotionModel,'Velocity') && ...
                    strcmp(obj.VelocitySource,'Property') && strcmp(prop, 'InitialVelocity')
                flag = true;
            elseif strcmp(obj.MotionModel,'Custom') && ...
                    (strcmp(prop,'InitialPosition') || strcmp(prop,'VelocitySource') || ...
                    strcmp(prop,'Velocity') || strcmp(prop,'AccelerationSource') || ...
                    strcmp(prop,'InitialVelocity') || strcmp(prop,'Acceleration'))
                flag = true;
            elseif ~strcmp(obj.MotionModel,'Custom') && strcmp(prop,'CustomTrajectory')
                flag = true;
            elseif strcmp(obj.ScanMode,'None') && ...
                    (strcmp(prop,'InitialScanAngle') || strcmp(prop,'AzimuthSpan') || ...
                    strcmp(prop,'AzimuthScanRate'))
                flag = true;
            elseif strcmp(obj.ScanMode,'Circular') && strcmp(prop,'AzimuthSpan')
                flag = true;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            s.pPreviousVelocity = obj.pPreviousVelocity;
            s.pVelocity = obj.pVelocity;
            s.pAcceleration = obj.pAcceleration;
            s.pOrientationAxes = obj.pOrientationAxes;
            s.pMotionModelIndex = obj.pMotionModelIndex;
            s.pIsVelocityViaInput = obj.pIsVelocityViaInput;
            s.pIsAccelerationViaInput = obj.pIsAccelerationViaInput;
            s.pScanModeIndex = obj.pScanModeIndex;
            s.pScanAngle = obj.pScanAngle;
            s.pCurrentScanDirection = obj.pCurrentScanDirection;
            s.pPPx = obj.pPPx;
            s.pPPy = obj.pPPy;
            s.pPPz = obj.pPPz;
            s.pPPvx = obj.pPPvx;
            s.pPPvy = obj.pPPvy;
            s.pPPvz = obj.pPPvz;
            s.pCurrentTime = obj.pCurrentTime;
            s.pTrajectoryMaxTime = obj.pTrajectoryMaxTime;
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.pPosition = obj.pPosition;
            end
        end
        
        function s = loadSubObjects(obj,s) 
            if isfield(s,'isLocked')
                s = rmfield(s,'isLocked');
                if ~isfield(s,'pPreviousVelocity')
                    obj.pPreviousVelocity = s.pVelocity;
                end
                if ~isfield(s,'InitialOrientationAxes')
                    obj.InitialOrientationAxes = s.pOrientationAxes;
                end
            end
        end
        
    end
    
    methods (Access = protected)
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

function T = rotvv(v1,v2)
%rotvv    Rotation matrix for rotating one vector to align with another
%   ROTMAT = rotvv(V1,V2) returns the rotation matrix, ROTMAT, that
%   rotates the vector, V1, to be aligned with the vector, V2. Both V1 and
%   V2 are originated from origin and are specified in the form of [x;y;z].
%
%   ROTMAT is a 3x3 matrix. The rotation of the point can be achieved by
%   left-multiplying ROTMAT with the point's coordinate vector [x;y;z].
%
%   % Example:
%   %   Rotate a vector, [1;2;3] to be aligned with the y axis. 
%
%   v1 = [1;2;3];
%   v2 = [0;1;0];
%   v1r = rotvv(v1,v2)*v1

%   Copyright 2015 The MathWorks, Inc.

%   References:
%   [1] James Foley, et. al. Computer Graphics Principles and Practices in
%       C, 2nd Edition, Addison-Wesley, 1995

%#codegen
%#ok<*EMCA>

if isequal(v1,v2)
    T = eye(3);
else
    nzidx = any(v1~=0);
    v1(:,nzidx) = bsxfun(@rdivide,v1(:,nzidx),sqrt(sum(abs(v1(:,nzidx).^2))));
    nzidx = any(v2~=0);
    v2(:,nzidx) = bsxfun(@rdivide,v2(:,nzidx),sqrt(sum(abs(v2(:,nzidx).^2))));

    v1azel = phased.internal.dirvec2azel(v1);
    v2azel = phased.internal.dirvec2azel(v2);
    % first rotate v1 back to x axis and then rotate to v2
    T = rotazel(v2azel(1),v2azel(2))*rotazel(v1azel(1),v1azel(2)).';
end

end

function T = rotazel(az,el)
%rotazel    Rotation matrix for specified azimuth and elevation angles
%   ROTMAT = rotazel(AZ,EL) returns the rotation matrix, ROTMAT, that
%   rotates a point by the specified azimuth angle, AZ (in degrees), and
%   elevation angle, EL (in degrees). The point is specified in the form of
%   [x;y;z]. The azimuth angle is defined as the angle from the X-axis
%   toward the Y-axis and the elevation angle is the angle from the X-Y
%   plane toward the Z-axis.
%
%   ROTMAT is a 3x3 matrix. The rotation of the point can be achieved by
%   left-multiplying ROTMAT with the point's coordinate vector [x;y;z].
%
%   % Example:
%   %   Rotate a point, (1,0,0) by 30 degrees azimuth and 45 degrees
%   %   elevation.
%
%   p = [1;0;0];
%   p = rotazel(30,45)*p

%   Copyright 2015 The MathWorks, Inc.

%   References:
%   [1] James Foley, et. al. Computer Graphics Principles and Practices in
%       C, 2nd Edition, Addison-Wesley, 1995

%#codegen
%#ok<*EMCA>

T = rotz(az)*roty(-el);  % -el because it's clockwise around y axis

end
