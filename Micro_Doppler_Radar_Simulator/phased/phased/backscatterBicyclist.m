classdef backscatterBicyclist < handle
%backscatterBicyclist   Backscatter bicyclist model for radar
%   H = backscatterBicyclist returns a backscatter bicyclist, H, for
%   radar simulation. This object simulates the signal reflected off the
%   bicyclist while the bicyclist is in motion.
%
%   H = backscatterBicyclist(Name,Value) returns a backscatter
%   bicyclist model, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The bicyclist is derived from a multi-scatterer model developed for a
%   77 GHz radar system. It is composed of 5 primary parts: bicycle frame
%   and upper-body of rider, pedals, legs, front wheel, and rear wheel.
%   Each part is composed of individual point scatterers. The total number
%   of point scatterers in the bicycle target model is dependent on the
%   number of spokes defined. The bicycle model is the aggregate of all of
%   the point scatterers.
%
%   move method syntax: [BBPOS,BBVEL,BBAX] = move(H,T,ANGH) returns the
%   bicyclist's scatterers' position, velocity, and orientation axes at the
%   current time in BBPOS, BBVEL, and BBAX, respectively. The method then
%   simulates the motion of the bicyclist in the next duration, specified
%   in T (in seconds). ANGH specifies the current heading angle (in
%   degrees) as a scalar. The angle is in the xy-plane. The bicyclist speed
%   and coasting may be updated using the backscatterBicyclist properties
%   Speed and Coast.
%
%   [BBPOS,BBVEL,BBAX] = move(...,SPEED) specifies the translational speed
%   (in m/s) of the bicyclist as a non-negative scalar. The bicyclist
%   coasting may be updated using the backscatterBicyclist property Coast.
%
%   [BBPOS,BBVEL,BBAX] = move(...,COAST) specifies the coasting of the
%   bicyclist as a logical. If set to true, the bicyclist does not pedal.
%   If set to false, the bicyclist pedals.
%
%   BBPOS is a 3xN matrix with each column representing the position of a
%   point scatterer in the [x;y;z] format (in meters), where N is the
%   number of total point scatterers. BBVEL is also a 3xN matrix whose
%   columns are velocities of the corresponding point scatterers in [x;y;z]
%   form (in m/s). BBAX is a 3x3 matrix that is the orientation axes of the
%   corresponding bicyclist in the [x;y;z] format.
%
%   reflect method syntax:
%
%   Y = reflect(H,X,ANG) returns the reflected signal Y off the bicyclist
%   target due to the input signal X.
%
%   X is an MxN matrix where M is the number of samples in the signal
%   incident to each point scatterer, and N is the number of point
%   scatterers.
%
%   ANG is a 2xN matrix representing the signal's incident direction to
%   each point scatterer. Each column of ANG specifies the incident
%   direction of the corresponding signal in the form of an [AzimuthAngle;
%   ElevationAngle] pair (in degrees).
%
%   Y is an M-length vector containing the combined reflected signals from
%   all point scatterers. The number of rows in Y is the same as the number
%   of rows in X.
%
%   RCS values depend on the viewing angle. The RCS values for the
%   individual point scatterers are computed from the angle of view of the
%   whole bicyclist. The default value of this property is a 1x361 matrix
%   with values that were derived from measurements taken at 77 GHz.
%
%   Internal occlusions within the bicyclist are not modeled.
%
%   plot method syntax:
%
%   HF = plot(H) returns a plot of the scatterers' positions for the
%   backscatter bicyclist object H at the current time.
%
%   HF is the figure handle.
%
%   getNumScatterers method syntax:
%
%   N = getNumScatterers(H) returns the number of scatterers in the
%   bicyclist model.
%
%   N is a scalar output of the number of scatterers.
%
%   backscatterBicyclist methods:
%
%   move             - Bicyclist motion (see above)
%   reflect          - Signal reflection off bicyclist (see above)
%   plot             - Plot bicyclist scatterers' positions (see above) 
%   getNumScatterers - Return number of scatterers in the bicyclist model
%   release          - Allow property value and input changes
%   clone            - Create bicyclist object with same property values
%   reset            - Reset internal states of the bicyclist target
%
%   backscatterBicyclist properties:
%
%   NumWheelSpokes         - Number of wheel spokes
%   GearTransmissionRatio  - Gear transmission ratio
%   OperatingFrequency     - Signal carrier frequency (Hz)
%   InitialPosition        - Initial position (m)
%   InitialHeading         - Initial heading direction (deg)
%   Speed                  - Bicyclist speed (m/s)
%   Coast                  - Coast bicyclist
%   PropagationSpeed       - Propagation speed (m/s)
%   AzimuthAngles          - Azimuth angles (deg)
%   ElevationAngles        - Elevation angles (deg)
%   RCSPattern             - Radar cross section pattern (square meters)
%
%   % Example:
%   %   Compute the reflected signal from a bicyclist moving along the
%   %   x-axis. Assume the radar works at 24 GHz, and the signal has a 300
%   %   MHz bandwidth. The signal is captured at the moment the bicyclist
%   %   starts to move and 1 second into the movement. Assume the radar is
%   %   at the origin.
%   
%   % Define parameters
%   bw = 3e8;
%   fs = bw;
%   fc = 24e9;
%   
%   % Initialize bicyclist object
%   bicycle = backscatterBicyclist(...
%       'OperatingFrequency',fc,...
%       'InitialPosition',[5;0;0]);
%   
%   % Initialize radar and propagation channel
%   rpos = [0;0;0];
%   wav = phased.LinearFMWaveform(...
%       'SampleRate',fs,...
%       'SweepBandwidth',bw);
%   x = wav();
%   chan = phased.FreeSpace(...
%       'OperatingFrequency',fc,...
%       'SampleRate',fs,...
%       'TwoWayPropagation',true);
%   
%   % Time 0
%   [bbpos,bbvel,bbax] = move(bicycle,1,0);
%   N = bicycle.getNumScatterers;
%   plot(bicycle)
%   xp = chan(repmat(x,1,N),rpos,bbpos,[0;0;0],bbvel);
%   [~,ang] = rangeangle(rpos,bbpos,bbax);
%   y0 = reflect(bicycle,xp,ang);
%   
%   % Time 1
%   [bbpos,bbvel,bbax] = move(bicycle,1,0);
%   plot(bicycle)
%   xp = chan(repmat(x,1,N),rpos,bbpos,[0;0;0],bbvel);
%   [~,ang] = rangeangle(rpos,bbpos,bbax);
%   y1 = reflect(bicycle,xp,ang);
%   
%   % Plot reflected returns
%   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(wav));
%   ymf = mf([y0 y1]);
%   t = (0:size(ymf,1)-1)/fs;
%   figure()
%   plot(t,mag2db(abs(ymf)));
%   ylim([-200 0])
%   xlabel('Time (s)');
%   ylabel('Magnitude (dB)');
%
%   See also phased, backscatterPedestrian, phased.RadarTarget.

%   Copyright 2019 The MathWorks, Inc.

%   Reference
%   [1] Stolz, M. et al. "Multi-Target Reflection Point Model of Cyclists
%   for Automotive Radar." 2017 European Radar Conference (EURAD),
%   Nuremberg, 2017, pp. 94-97.
%
%   [2] Chen, V., D. Tahmoush, and W. J. Miceli. Radar Micro-Doppler
%   Signatures: Processing and Applications. The Institution of Engineering
%   and Technology: London, 2014.
%
%   [3] Belgiovane, D., and C. C. Chen. "Bicycles and  Human Riders
%   Backscattering at 77 GHz for Automotive Radar." 2016 10th European
%   Conference on Antennas and Propagation (EuCAP), Davos, 2016, pp. 1-5.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties
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
        %OperatingFrequency     Signal carrier frequency (Hz)
        %   Specify the carrier frequency (in Hz) of the narrowband signal
        %   as a scalar. The default value of this property is 77e9 (77
        %   GHz).
        OperatingFrequency = 77e9
        %InitialPosition        Initial position (m)
        %   Specify the initial position of the bicyclist as a 3x1 vector
        %   in the form of [x; y; z] (in meters). The default value of this
        %   property is [0; 0; 0].
        InitialPosition = [0;0;0]
        %InitialHeading         Initial heading direction (deg)
        %   Specify the initial heading (in degrees) of the bicyclist as a
        %   scalar. The heading is measured from the x-axis towards the
        %   y-axis in the xy-plane. The default value is 0.
        InitialHeading = 0
        %Speed                  Bicyclist speed (m/s)
        %   Specify the speed (in m/s) of the bicyclist as a nonnegative
        %   scalar. The default value is 4.
        %
        %   The motion model limits the speed to 60 m/s.
        Speed = 4
        %Coast                  Coast bicyclist
        %   A logical that controls the coasting of the bicyclist. If set
        %   to true, the bicyclist does not pedal. If set to false, the
        %   bicyclist pedals. Defaults to false.
        Coast(1,1) logical = false 
        %PropagationSpeed       Propagation speed (m/s)
        %   Specify the wave propagation speed (in m/s) in free space as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('lightspeed')
        %AzimuthAngles          Azimuth angles (deg)
        %   Specify the azimuth angles (in degrees) as a length P vector.
        %   These are the azimuth angles where the custom pattern is
        %   evaluated. P must be greater than 2. The default value of this
        %   property is -180:180.
        AzimuthAngles = -180:180
        %ElevationAngles        Elevation angles (deg)
        %   Specify the elevation angles (in degrees) as a length Q vector.
        %   These are the elevation angles where the custom pattern is
        %   evaluated. The default value of this property is 0.
        ElevationAngles = 0
        %RCSPattern             Radar cross section pattern (square meters)
        %   Specify the radar cross section pattern of the whole bicyclist
        %   (in square meters) as a QxP matrix, where Q is the number of
        %   elements presented in the ElevationAngles property, and P is
        %   the number of elements presented in the AzimuthAngles property.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can be given as a 1xP vector.
        %
        %   The individual point scatterers get an RCS value approximated
        %   by the angle of view of the whole cyclist. The default value of
        %   this property is a 1x361 matrix with values that were derived
        %   from measurements taken at 77 GHz.
        RCSPattern = backscatterBicyclist.defaultRCSPattern
    end
    
    properties (Access=private)
        cMovingBicyclist
        cBackscatterBicyclistTarget
        pIsInitialized = false
        pIsMotionInitialized = false     
        
        pBicyclistPosition
        pBicyclistVelocity
        pBicyclistAxes
        pBicyclistInitialPosition
        pNumScatterers
        
        pFigInitialized = false
        pScatterObjInitialized = false
        pFigHandle
        pScatterObjHandle
        pFigAxes
        pFigAxesNV
    end
    
    properties (Constant, Hidden)
        pPointScattererSpacing = 0.1
        pWheelDiameter = 0.75
        pLengthUpperLeg = 0.4
        pLengthLowerLeg = 0.4
    end

    methods
        function set.Speed(obj,val)
            sigdatatypes.validateSpeed(val,'backscatterBicyclist','Speed',{'scalar','real'});
            obj.Speed = val;
        end
        function set.Coast(obj,val)
            sigdatatypes.checkNumericOrLogicalScalar(obj,'Coast', val)
            obj.Coast = logical(val);
        end               
        function set.NumWheelSpokes(obj,val)
            sigdatatypes.checkFiniteNonNegIntScalar(obj,'NumWheelSpokes', val)
            obj.NumWheelSpokes = val;
        end
        function set.GearTransmissionRatio(obj,val)
            sigdatatypes.checkFinitePosDblScalar(obj,'GearTransmissionRatio', val)
            obj.GearTransmissionRatio = val;
        end
        function set.OperatingFrequency(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'nonnegative','real'},'backscatterBicyclist','OperatingFrequency');
            obj.OperatingFrequency = value;
        end
        function set.PropagationSpeed(obj,value)
            sigdatatypes.validateSpeed(value,'backscatterBicyclist',...
                'PropagationSpeed',{'scalar','positive','real'});
            obj.PropagationSpeed = value;
        end
        function set.InitialPosition(obj,val)
            sigdatatypes.validate3DCartCoord(val,...
                'backscatterBicyclist','InitialPosition',{'size',[3 1]});
            obj.InitialPosition = val;
        end
        function set.InitialHeading(obj,val)
            sigdatatypes.validateAngle(val,...
                'backscatterBicyclist','InitialHeading',{'scalar'});
            obj.InitialHeading = val;
        end
        function set.AzimuthAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'backscatterBicyclist','AzimuthAngles',...
                {'vector','>=',-180,'<=',180});
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'AzimuthAngles');
            end
            obj.AzimuthAngles = value;
        end
        function set.ElevationAngles(obj,value)
            % Allow single elevation cut, i.e., azimuth only pattern
            sigdatatypes.validateAngle(value,...
                'backscatterBicyclist','ElevationAngles',...
                {'vector','>=',-90,'<=',90});
            obj.ElevationAngles = value;
        end
        function set.RCSPattern(obj,value)
            validateattributes(value,{'double'},...
                {'real','nonnegative','finite','nonempty','2d'},...
                'backscatterBicyclist','RCSPattern');
            obj.RCSPattern  = value;
        end
    end
    
    methods
        function obj = backscatterBicyclist(varargin)
            % Assign defaults 
            defaultSpeed = 4;
            defaultNumWheelSpokes = 20;
            defaultGearTransmissionRatio = 1.5;
            defaultPropagationSpeed = physconst('lightspeed');
            defaultOperatingFrequency = 77e9;
            defaultInitialPosition = zeros(3,1);
            defaultInitialHeading = 0;
            defaultCoast = false;
            defaultAzimuthAngles = -180:180;
            defaultElevationAngles = 0;
            defaultRCSPattern = backscatterBicyclist.defaultRCSPattern;
            
            if coder.target('MATLAB')
                % Parse MATLAB properties
                p = inputParser;
                p.CaseSensitive = true;
                p.FunctionName = 'BicycleTarget';
                p.PartialMatching = false;
                addParameter(p,'Speed',defaultSpeed);
                addParameter(p,'NumWheelSpokes',defaultNumWheelSpokes);
                addParameter(p,'GearTransmissionRatio',defaultGearTransmissionRatio);
                addParameter(p,'PropagationSpeed',defaultPropagationSpeed);
                addParameter(p,'OperatingFrequency',defaultOperatingFrequency);
                addParameter(p,'InitialPosition',defaultInitialPosition);
                addParameter(p,'InitialHeading',defaultInitialHeading);
                addParameter(p,'Coast',defaultCoast);
                addParameter(p,'AzimuthAngles',defaultAzimuthAngles);
                addParameter(p,'ElevationAngles',defaultElevationAngles);
                addParameter(p,'RCSPattern',defaultRCSPattern);

                parse(p,varargin{:});
                
                fn = fieldnames(p.Results);
                for m = 1:numel(fn)
                    obj.(fn{m}) = p.Results.(fn{m});
                end
            else
                % Parse CodeGen properties 
                parms = struct(...
                    'Speed',uint32(0), ...
                    'NumWheelSpokes',uint32(0), ...
                    'GearTransmissionRatio',uint32(0), ...
                    'PropagationSpeed',uint32(0), ...
                    'OperatingFrequency',uint32(0), ...
                    'InitialPosition',uint32(0), ...
                    'InitialHeading',uint32(0),...
                    'Coast',uint32(0),...
                    'AzimuthAngles',uint32(0), ...
                    'ElevationAngles',uint32(0), ...
                    'RCSPattern',uint32(0));
                pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
                obj.Speed = eml_get_parameter_value(pstruct.Speed,defaultSpeed,varargin{:});
                obj.NumWheelSpokes = eml_get_parameter_value(pstruct.NumWheelSpokes,defaultNumWheelSpokes,varargin{:});
                obj.GearTransmissionRatio = eml_get_parameter_value(pstruct.GearTransmissionRatio,defaultGearTransmissionRatio,varargin{:});
                obj.PropagationSpeed = eml_get_parameter_value(pstruct.PropagationSpeed,defaultPropagationSpeed,varargin{:});
                obj.OperatingFrequency = eml_get_parameter_value(pstruct.OperatingFrequency,defaultOperatingFrequency,varargin{:});
                obj.InitialPosition = eml_get_parameter_value(pstruct.InitialPosition,defaultInitialPosition,varargin{:});
                obj.InitialHeading = eml_get_parameter_value(pstruct.InitialHeading,defaultInitialHeading,varargin{:});
                obj.Coast = eml_get_parameter_value(pstruct.Coast,defaultCoast,varargin{:});
                obj.AzimuthAngles = eml_get_parameter_value(pstruct.AzimuthAngles,defaultAzimuthAngles,varargin{:});
                obj.ElevationAngles = eml_get_parameter_value(pstruct.ElevationAngles,defaultElevationAngles,varargin{:});
                obj.RCSPattern = eml_get_parameter_value(pstruct.RCSPattern,defaultRCSPattern,varargin{:});
            end
            
            % Determine number of scatterers based on inputs
            nPtsUpper = ceil(obj.pLengthUpperLeg/obj.pPointScattererSpacing);
            nPtsLower = ceil(obj.pLengthLowerLeg/obj.pPointScattererSpacing)-1;
            upperLegVec = (obj.pLengthUpperLeg/nPtsUpper).*(1:nPtsUpper);
            lowerLegVec = (obj.pLengthLowerLeg/(nPtsLower+1)).*(1:nPtsLower);
            [frameAndRiderUpperBody, pedals, legs, ...
                frontWheel, rearWheel, ~] = ...
                phased.internal.BackscatterBicyclist.getBicyclistParts(obj.pWheelDiameter,...
                obj.NumWheelSpokes,obj.pPointScattererSpacing,...
                upperLegVec,lowerLegVec);
            obj.pBicyclistInitialPosition = [frameAndRiderUpperBody pedals legs frontWheel rearWheel];
            obj.pNumScatterers = size(obj.pBicyclistInitialPosition,2);
        end
        
        function release(obj)
            %release    Release backscatterBicyclist object.
            if obj.pIsInitialized
                release(obj.cMovingBicyclist);
                release(obj.cBackscatterBicyclistTarget);
                obj.pIsInitialized = false;
            end
            obj.pIsMotionInitialized = false;
            obj.pFigInitialized = false;
        end
        
        function reset(obj)
            %release    Reset backscatterBicyclist object.
            reset(obj.cMovingBicyclist);
            reset(obj.cBackscatterBicyclistTarget);
        end
        
        function [bbpos,bbvel,bbax] = move(obj,dt,heading,varargin)
            %move   Bicyclist motion
            %   [BBPOS,BBVEL,BBAX] = move(H,T,ANGH) returns the bicyclist's
            %   scatterers' position, velocity, and orientation axes at the
            %   current time in BBPOS, BBVEL, and BBAX, respectively. The
            %   method then simulates the motion of the bicyclist in the
            %   next duration, specified in T (in seconds). ANGH specifies
            %   the current heading angle (in degrees) as a scalar. The
            %   angle is in the xy-plane. The bicyclist speed and coasting
            %   may be updated using the backscatterBicyclist properties
            %   Speed and Coast.
            %
            %   [BBPOS,BBVEL,BBAX] = move(...,SPEED) specifies the
            %   translational speed (in m/s) of the bicyclist as a
            %   non-negative scalar. The bicyclist coasting may be updated
            %   using the backscatterBicyclist property Coast.
            %
            %   [BBPOS,BBVEL,BBAX] = move(...,COAST) specifies the coasting
            %   of the bicyclist as a logical. If set to true, the
            %   bicyclist does not pedal. If set to false, the bicyclist
            %   pedals.
            %
            %   BBPOS is a 3xN matrix with each column representing the
            %   position of a point scatterer in the [x;y;z] format (in
            %   meters), where N is the number of total point scatterers.
            %   BBVEL is also a 3xN matrix whose columns are velocities of
            %   the corresponding point scatterers in [x;y;z] form (in
            %   m/s). BBAX is a 3x3 matrix that is the orientation axes of
            %   the corresponding bicyclist in the [x;y;z] format.
            %
            %   % Example 1:
            %   %   Model a bicyclist that accelerates from 0 to 1.5 m/s at
            %   %   a rate of 0.25 m/s^2. After the bicyclist reaches 1.5
            %   %   m/s, the bicyclist coasts.
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
            %   bicycle = backscatterBicyclist(...
            %       'Speed',speed,...
            %       'Coast',coast);
            %
            %   % Build movement trajectory
            %   for m = 1:M
            %       if speed < coastingSpeed
            %           % When the speed is less than the coasting speed,
            %           % continue to accelerate
            %           speed = acceleration*tSim(m);
            %       else
            %           % When the speed is equal to the coasting speed,
            %           % coast the bicyclist
            %           coast = true;
            %       end
            %       ppos = move(bicycle,dt,0,speed,coast);
            %
            %       % Plot
            %       plot(bicycle);
            %   end
            %
            %   % Example 2:
            %   %   Model the motion of a bicyclist riding in a quarter
            %   %   circle.
            %
            %   % Define simulation
            %   dt = 0.003;
            %   M = 3000;
            %   angstep = 90/M;
            %
            %   % Initialize bicyclist object
            %   bicycle = backscatterBicyclist;
            %
            %   % Build movement trajectory
            %   for m = 1:M
            %       ppos = move(bicycle,dt,angstep*m);
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
            %      drawnow limitrate;
            %   end
            %
            %   See also phased, backscatterBicyclist/reflect.

            %   Reference
            %   [1] Stolz, M. et al. "Multi-Target Reflection Point Model
            %   of Cyclists for Automotive Radar." 2017 European Radar
            %   Conference (EURAD), Nuremberg, 2017, pp. 94-97.
            %
            %   [2] Chen, V., D. Tahmoush, and W. J. Miceli. Radar
            %   Micro-Doppler Signatures: Processing and Applications. The
            %   Institution of Engineering and Technology: London, 2014.
            
            % For the sake of performance, rely on input validation in
            % MovingBicyclist
            
            narginchk(3,5); 
            nargoutchk(0,3); 
            
            % Setup if not initialized
            if ~obj.pIsInitialized
                callSetup(obj);
            end
            
            % Ensure nontunable properties have not changed 
            moveMethodNontunablePropertiesCheck(obj);
            
            % Parse variable arguments
            [speed, coast] = parseVarargin(obj,varargin{:});
            
            % Step 
            [bbpos,bbvel,bbax] = step(obj.cMovingBicyclist,dt,heading,speed,coast);
            
            % Assign outputs of move to object properties 
            obj.pBicyclistPosition = bbpos;
            obj.pBicyclistVelocity = bbvel;
            obj.pBicyclistAxes = bbax;
            obj.pIsMotionInitialized = true;
        end
        
        function y = reflect(obj,x,ang)
            %reflect    Signal reflection off bicyclist
            %   Y = reflect(H,X,ANG) returns the reflected signal Y off the
            %   bicyclist target due to the input signal X.
            %
            %   X is an MxN matrix where M is the number of samples in the
            %   signal incident to each point scatterer, and N is the
            %   number of point scatterers.
            %
            %   ANG is a 2xN matrix representing the signal's incident
            %   direction to each point scatterer. Each column of ANG
            %   specifies the incident direction of the corresponding
            %   signal in the form of an [AzimuthAngle; ElevationAngle]
            %   pair (in degrees).
            %
            %   Y is an M-length vector containing the combined reflected
            %   signals from all point scatterers. The number of rows in Y
            %   is the same as the number of rows in X.
            %
            %   RCS values depend on the viewing angle. The RCS values for
            %   the individual point scatterers are computed from the angle
            %   of view of the whole bicyclist. The default value of this
            %   property is a 1x361 matrix with values that were derived
            %   from measurements taken at 77 GHz.
            %
            %   Internal occlusions within the bicyclist are not modeled.
            %
            %   % Example:
            %   %   Compute the reflected signal from a bicyclist moving
            %   %   along the x axis. Assume the radar works at 24 GHz and
            %   %   the signal has a 300 MHz bandwidth. The signal is
            %   %   captured at the moment the bicyclist starts to move and
            %   %   1 second into the movement. Assume the radar is at the
            %   %   origin.
            %
            %   % Define parameters
            %   bw = 3e8;
            %   fs = bw;
            %   fc = 24e9;
            %
            %   % Initialize bicyclist object
            %   bicycle = backscatterBicyclist(...
            %       'OperatingFrequency',fc,...
            %       'InitialPosition',[5;0;0]);
            %   N = bicycle.getNumScatterers;
            %
            %   % Initialize radar and propagation channel
            %   rpos = [0;0;0];
            %   wav = phased.LinearFMWaveform(...
            %       'SampleRate',fs,...
            %       'SweepBandwidth',bw);
            %   x = wav();
            %   chan = phased.FreeSpace(...
            %       'OperatingFrequency',fc,...
            %       'SampleRate',fs,...
            %       'TwoWayPropagation',true);
            %
            %   % Time 0
            %   [bbpos,bbvel,bbax] = move(bicycle,1,0);
            %   xp = chan(repmat(x,1,N),rpos,bbpos,[0;0;0],bbvel);
            %   [~,ang] = rangeangle(rpos,bbpos,bbax);
            %   y0 = reflect(bicycle,xp,ang);
            %
            %   % Time 1
            %   [bbpos,bbvel,bbax] = move(bicycle,1,0);
            %   xp = chan(repmat(x,1,N),rpos,bbpos,[0;0;0],bbvel);
            %   [~,ang] = rangeangle(rpos,bbpos,bbax);
            %   y1 = reflect(bicycle,xp,ang);
            %
            %   % Plot
            %   mf = phased.MatchedFilter(...
            %       'Coefficients',getMatchedFilter(wav));
            %   ymf = mf([y0 y1]);
            %   t = (0:size(ymf,1)-1)/fs;
            %   plot(t,mag2db(abs(ymf)));
            %   ylim([-200 0])
            %   xlabel('Time (s)');
            %   ylabel('Magnitude (dB)');
            %
            %   See also phased, backscatterBicyclist/move.
            
            %   Reference
            %   [1] Stolz, M. et al. "Multi-Target Reflection Point Model
            %   of Cyclists for Automotive Radar." 2017 European Radar
            %   Conference (EURAD), Nuremberg, 2017, pp. 94-97.
            %
            %   [2] Chen, V., D. Tahmoush, and W. J. Miceli. Radar
            %   Micro-Doppler Signatures: Processing and Applications. The
            %   Institution of Engineering and Technology: London, 2014.
            %
            %   [3] Belgiovane, D., and C. C. Chen. "Bicycles and  Human
            %   Riders Backscattering at 77 GHz for Automotive Radar." 2016
            %   10th European Conference on Antennas and Propagation
            %   (EuCAP), Davos, 2016, pp. 1-5.
            
            % For the sake of performance, rely on input validation in
            % BackscatterBicyclistTarget
            
            narginchk(3,3); 
            nargoutchk(0,1);
            
            % Check motion has been initialized 
            cond = ~obj.pIsMotionInitialized;
            if cond
                coder.internal.errorIf(cond,'phased:target:invalidCallSequence','move','reflect');
            end
            
            % Ensure nontunable properties have not changed 
            reflectMethodNontunablePropertiesCheck(obj)
            
            % Set the number of scatterers 
            if obj.cBackscatterBicyclistTarget.NumScatterers == -1
                obj.cBackscatterBicyclistTarget.NumScatterers = size(x,2);
            end
            
            % Reflect 
            y = step(obj.cBackscatterBicyclistTarget,x,ang);
        end
        
        function varargout = plot(obj,varargin)
            %plot    Plot bicyclist scatterers' positions
            %   HF = plot(H) returns a plot of the scatterers' positions
            %   for the backscatter bicyclist object H at the current time.
            %
            %   HF is the figure handle.
            %
            %   HF = plot(...,'Parent',HAX) specifies the target axes HAX
            %   for the plotting of the bicyclist object H. 
            %
            %   % Example:
            %   %   Model a bicyclist that accelerates from 0 to 1.5 m/s at
            %   %   a rate of 0.25 m/s^2. After the bicyclist reaches 1.5
            %   %   m/s, the bicyclist coasts.
            %
            %   % Define simulation
            %   dt = 0.003;
            %   M = 3000;
            %   tSim = 0:dt:M*dt-dt;
            %   coastingSpeed = 1.5;
            %   acceleration = 0.25;
            %
            %   % Initialize bicyclist object
            %   bicycle = backscatterBicyclist('Speed',0,'Coast',0);
            %
            %   % Build movement trajectory
            %   for m = 1:M
            %       if bicycle.Speed < coastingSpeed
            %           % When the speed is less than the coasting speed, 
            %           % continue to accelerate
            %           bicycle.Speed = acceleration*tSim(m);
            %       else
            %           % When the speed is equal to the coasting speed, 
            %           % coast the bicyclist
            %           bicycle.Coast = 1;
            %       end
            %       move(bicycle,dt,0);
            %       plot(bicycle);
            %   end
            %
            %   See also phased, backscatterBicyclist/move.
            
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
                obj.pFigInitialized = 1;
                obj.pFigHandle = obj.pFigAxes.Parent;
            end

            % Check to see if figure was deleted
            if ~ishandle(obj.pFigHandle)
                obj.pFigInitialized = false;
            end
            
            % Check to see if scatter object was deleted
            if ~ishandle(obj.pScatterObjHandle)
                obj.pScatterObjInitialized = false;
            else
                if obj.pFigInitialized
                    if obj.pScatterObjHandle.Parent ~= obj.pFigAxes % Verify that axes have remained the same
                        obj.pScatterObjInitialized = false;
                    end
                end
            end
            
            % Check to see if motion is initialized 
            cond = ~obj.pIsMotionInitialized;
            if cond
                bbpos = obj.pBicyclistInitialPosition; % Use initial position
            else % Motion is initialized
                bbpos = obj.pBicyclistPosition; % Use current position
            end

            % Plot
            colorVal = [0 0.4470 0.7410];
            x = bbpos(1,:);
            y = bbpos(2,:);
            z = bbpos(3,:);
            if obj.pScatterObjInitialized
                % Plot is initialized
                obj.pScatterObjHandle.XData = x;
                obj.pScatterObjHandle.YData = y;
                obj.pScatterObjHandle.ZData = z;
                obj.pFigAxes.XLim = [min(x)-0.1 max(x)+0.1];
                obj.pFigAxes.YLim = [min(y)-0.1 max(y)+0.1];
                drawnow limitrate;
            else
                % Plot is not initialized 
                if ~obj.pFigInitialized
                    obj.pFigHandle = figure('Name','Bicyclist Trajectory');
                    obj.pFigInitialized = true;
                end
                obj.pScatterObjHandle = line(x,y,z,...
                    'Color',colorVal,'LineStyle','none',...
                    'Marker','o','MarkerSize',4,'MarkerFaceColor',colorVal);
                axis equal;
                grid on; 
                obj.pFigAxes = gca;
                obj.pFigAxes.XLim = [min(x)-0.1 max(x)+0.1];
                obj.pFigAxes.YLim = [min(y)-0.1 max(y)+0.1];
                obj.pFigAxes.ZLim = [min(z)-0.1 max(z)+0.1];
                obj.pFigAxes.XLabel.String = 'X (m)';
                obj.pFigAxes.YLabel.String = 'Y (m)';
                obj.pFigAxes.ZLabel.String = 'Z (m)';
                obj.pFigAxes.Title.String = 'Bicyclist Trajectory';
                view(36,15);
                obj.pScatterObjInitialized = true;
                drawnow;
            end
            hFig = obj.pFigHandle; 
            
            % Define output handle if requested
            if nargout == 1
                varargout{1} = hFig; 
            end
        end
        
        function N = getNumScatterers(obj)
            %getNumScatterers    Return number of scatterers in bicyclist
            %   N = getNumScatterers(H) returns the number of scatterers in
            %   the bicyclist model.
            %
            %   N is a scalar output of the number of scatterers.
            %
            %   % Example:
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
        
        function newobj = clone(obj)
            % clone Create object with same property values
            %   C = clone(OBJ) creates another instance of the object, OBJ,
            %   with the same property values.
            s = saveobj(obj);
            newobj = backscatterBicyclist;
            loadPrivateProps(newobj,s);
        end
    end
    
    methods (Access=protected)
        function callSetup(obj)
            obj.setup();
            obj.pIsInitialized = true;
        end
        
        function setup(obj)                 
            % Move system object
            obj.cMovingBicyclist = phased.internal.MovingBicyclist(...
                'InputSource','Input port',...
                'NumWheelSpokes',obj.NumWheelSpokes,...
                'GearTransmissionRatio',obj.GearTransmissionRatio,...
                'InitialPosition',obj.InitialPosition,...
                'InitialHeading',obj.InitialHeading,...
                'InitialSpeed',obj.Speed);
            
            % Reflect system object
            obj.cBackscatterBicyclistTarget = phased.internal.BackscatterBicyclistTarget(...
                'PropagationSpeed',obj.PropagationSpeed,...
                'OperatingFrequency',obj.OperatingFrequency,...
                'RCSPatternSource','Property',...
                'AzimuthAngles',obj.AzimuthAngles,...
                'ElevationAngles',obj.ElevationAngles,...
                'RCSPattern',obj.RCSPattern);
        end
        
        function moveMethodNontunablePropertiesCheck(obj)
            if obj.pIsInitialized
                cond = obj.NumWheelSpokes ~= obj.cMovingBicyclist.NumWheelSpokes;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','NumWheelSpokes','backscatterBicyclist');
                end
                cond = obj.GearTransmissionRatio ~= obj.cMovingBicyclist.GearTransmissionRatio;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','GearTransmissionRatio','backscatterBicyclist');
                end
                cond = obj.InitialPosition ~= obj.cMovingBicyclist.InitialPosition;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','InitialPosition','backscatterBicyclist');
                end
                cond = obj.InitialHeading ~= obj.cMovingBicyclist.InitialHeading;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','InitialHeading','backscatterBicyclist');
                end
            end
        end
        
        function reflectMethodNontunablePropertiesCheck(obj)
            if obj.pIsInitialized
                cond = obj.PropagationSpeed ~= obj.cBackscatterBicyclistTarget.PropagationSpeed;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','PropagationSpeed','backscatterBicyclist');
                end
                cond = obj.OperatingFrequency ~= obj.cBackscatterBicyclistTarget.OperatingFrequency;
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','OperatingFrequency','backscatterBicyclist');
                end
                cond = any(obj.AzimuthAngles ~= obj.cBackscatterBicyclistTarget.AzimuthAngles);
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','AzimuthAngles','backscatterBicyclist');
                end
                cond = any(obj.ElevationAngles ~= obj.cBackscatterBicyclistTarget.ElevationAngles);
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','ElevationAngles','backscatterBicyclist');
                end
                cond = any(size(obj.RCSPattern(:)) ~= size(obj.cBackscatterBicyclistTarget.RCSPattern(:)));
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','RCSPattern','backscatterBicyclist');
                end
                cond = any(obj.RCSPattern(:) ~= obj.cBackscatterBicyclistTarget.RCSPattern(:));
                if cond
                    coder.internal.errorIf(cond,'phased:target:nontunableProperty','RCSPattern','backscatterBicyclist');
                end
            end
        end
        
        function s = saveobj(obj)
            s.Speed = obj.Speed;
            s.NumWheelSpokes = obj.NumWheelSpokes;
            s.GearTransmissionRatio = obj.GearTransmissionRatio;
            s.PropagationSpeed = obj.PropagationSpeed;
            s.OperatingFrequency = obj.OperatingFrequency;
            s.InitialPosition = obj.InitialPosition;
            s.InitialHeading = obj.InitialHeading;
            s.Coast = obj.Coast; 
            s.AzimuthAngles = obj.AzimuthAngles; 
            s.ElevationAngles = obj.ElevationAngles; 
            s.RCSPattern = obj.RCSPattern; 
            
            s.cMovingBicyclist = saveobj(obj.cMovingBicyclist);
            s.cBackscatterBicyclistTarget = saveobj(obj.cBackscatterBicyclistTarget);
            s.pIsInitialized = obj.pIsInitialized;
            s.pIsMotionInitialized = obj.pIsMotionInitialized;
            
            s.pBicyclistPosition = obj.pBicyclistPosition;
            s.pBicyclistVelocity = obj.pBicyclistVelocity;
            s.pBicyclistAxes = obj.pBicyclistAxes;
            s.pBicyclistInitialPosition = obj.pBicyclistInitialPosition;
            s.pNumScatterers = obj.pNumScatterers; 
           
            s.pFigInitialized = obj.pFigInitialized;
            s.pScatterObjInitialized = obj.pScatterObjInitialized;
            s.pFigHandle = obj.pFigHandle; 
            s.pScatterObjHandle = obj.pScatterObjHandle;
            s.pFigAxes = obj.pFigAxes;
            s.pFigAxesNV = obj.pFigAxesNV;
        end
        
        function loadPrivateProps(obj,s)
            obj.cMovingBicyclist = phased.internal.MovingBicyclist.loadobj(s.cMovingBicyclist);
            s = rmfield(s,'cMovingBicyclist');
            obj.cBackscatterBicyclistTarget = phased.internal.BackscatterBicyclistTarget.loadobj(s.cBackscatterBicyclistTarget);
            s = rmfield(s,'cBackscatterBicyclistTarget');
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function [speed,coast] = parseVarargin(obj,varargin)
            % Parse variable arguments
            if nargin == 3
                speed = varargin{1};
                coast = logical(varargin{2});
            elseif nargin == 2
                speed = varargin{1};
                coast = obj.Coast;
            else
                speed = obj.Speed;
                coast = obj.Coast;
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
            props = {'NumWheelSpokes','GearTransmissionRatio',...
                'PropagationSpeed','OperatingFrequency','InitialPosition',...
                'InitialHeading','AzimuthAngles','ElevationAngles',...
                'cMovingBicyclist','cBackscatterBicyclistTarget'};
        end
        function names = matlabCodegenOnceNames()
            names = {'pIsInitialized', 'callSetup', 'callSetup', 'setup', 'setup'};
        end
        % Load
        function obj = loadobj(s)
            obj = backscatterBicyclist;
            loadPrivateProps(obj,s);
        end
        % Default RCS pattern
        function rcs = defaultRCSPattern()
            rcs = phased.internal.BackscatterBicyclist.bicyclistrcs;
        end
    end
end