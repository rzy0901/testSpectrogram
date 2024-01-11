classdef (Sealed,StrictDefaults) ScenarioViewer < matlabshared.scopes.UnifiedSystemScope
%ScenarioViewer Display trajectories of radars and targets
%   SV = phased.ScenarioViewer returns a System object, SV, that can show
%   the trajectories of moving radars and targets over time.
%
%   SV = phased.ScenarioViewer('Name', Value, ...) returns a ScenarioViewer
%   System object, SV, with each specified property name set to the
%   specified value. You can specify name-value pair arguments in any order
%   as (Name 1, Value 1, ..., Name N, Value N).
%
%   Step method syntax:
%
%   step(sv,rpos,tpos) displays the trajectories of radars and targets
%   whose positions are specified in rpos and tpos respectively. This
%   syntax is applicable when the VelocityInputPort and the
%   OrientationInputPort properties are set to false. rpos and tpos are a
%   3xNr and 3xNt position matrices where each column is in the form of [x;
%   y; z] (in meters). Nr and Nt are the number of radars and targets to
%   track. Nr and Nt should be greater or equal to 1.
%
%   step(sv,rpos,rvel,tpos,tvel) specify the radars velocity vectors, rvel,
%   and the targets velocity vectors, tvel, when the VelocityInputPort
%   property is set to true. rvel and tvel have the same dimensions as rpos
%   and tpos respectively. Each column of rvel and tvel is in the form of
%   [Vx; Vy; Vz] (in meters/second). This is the default syntax.
%
%   step(sv,rpos,raxis,tpos,taxis) specify the radars orientation axis,
%   raxis, and the targets orientation axis, taxis, when the
%   OrientationInputPort property is set to true. raxis and taxis are in
%   the form of 3x3xNr and 3x3xNt matrices. Values along the first
%   dimension specify the direction of an axis in the form if [x; y; z]
%   coordinates. The index to the second dimension specifies the x, y or z
%   axis and the index to the third dimension specifies the corresponding
%   radar or target.
%
%   step(sv,rpos,rvel,raxis,tpos,tvel,taxis) specify both the velocity and
%   the orientation when the VelocityInputPort and OrientationInputPort
%   properties are set to true.
%
%   When the VelocityInputPort property is set to false, the velocity
%   vectors are estimated by computing the distance traveled between
%   consecutive positions during the elapsed time deduced from the
%   UpdateRate property.
%
%   When the OrientationInputPort property is set to false, the orientation
%   axis are aligned with the global coordinate axis.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ScenarioViewer methods:
%
%   step      - Plot trajectories in the ScenarioViewer figure (see above)
%   release   - Allow property value and input characteristics changes, and
%               release ScenarioViewer resources
%   clone     - Create ScenarioViewer object with same property values
%   isLocked  - Display locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>     - Clear ScenarioViewer figure
%   show      - Turn on visibility of ScenarioViewer figure
%   hide      - Turn off visibility of ScenarioViewer figure
%   isVisible - Return visibility of ScenarioViewer figure (logical)
%
%   ScenarioViewer properties:
%
%   Name                 - Caption to display on the ScenarioViewer window
%   ReferenceRadar       - Radar used as reference
%   ShowBeam             - Enable beam visualization
%   BeamWidth            - Beam width of reference radar
%   BeamRange            - Beam range of reference radar
%   BeamSteering         - Beam direction relative to reference radar
%   VelocityInputPort    - Enable velocity input
%   OrientationInputPort - Enable Orientation input
%   UpdateRate           - Viewer update rate
%   Title                - Display title
%   PlatformNames        - PlatformNames Names of radars and targets
%   TrailLength          - Length of visible trajectories
%   CameraPerspective    - Perspective of the camera
%   CameraPosition       - Position of the camera
%   CameraOrientation    - Orientation of the camera
%   CameraViewAngle      - View angle of the camera
%   ShowLegend           - Show legend
%   ShowGround           - Enable ground visualization
%   ShowName             - Annotate trajectories with names
%   ShowPosition         - Annotate trajectories with positions
%   ShowRange            - Annotate trajectories with range
%   ShowAltitude         - Annotate trajectories with altitude
%   ShowSpeed            - Annotate trajectories with speed
%   ShowRadialSpeed      - Annotate trajectories with radial speed
%   ShowAzEl             - Annotate trajectories with azimuth and elevation
%   Position             - Scope window position in pixels
%   ReducePlotRate       - Reduce plot rate to improve performance
%
%   % Example:
%   %   Visualize the trajectory of a radar and three targets.
%   radarPlatform = phased.Platform(...
%       'InitialPosition',[0;0;10], ...
%       'Velocity',[20;0;5]);
%   targetPlatforms = phased.Platform(...
%       'InitialPosition',[2000.66 3532.63 3845.04;0 500 0;1e3 1e3 10],...
%       'Velocity',[120 -120 -10; 10 -20 -10; 0 0 60]);
%   ur = .1;
%   sv = phased.ScenarioViewer('BeamRange',3e3,'UpdateRate',ur);
%   
%   for i = 1:1e2
%       [rPos,rVel] = radarPlatform(ur);
%       [tPos,tVel] = targetPlatforms(ur);
%       sv(rPos,rVel,tPos,tVel);
%       pause(.1);
%   end
%
%   See also phased.Platform, phased.IntensityScope


%   Copyright 2015-2018 The MathWorks, Inc.

properties(Nontunable)
    %Name Caption to display on the ScenarioViewer window
    %   Specify the caption to display on the ScenarioViewer window as any
    %   string. The default value of this property is 'Scenario Viewer'.
    Name = 'Scenario Viewer';
end

properties (Dependent)
    %ReferenceRadar Radar used as reference
    %   Specify the radar to be used as a reference. The property is an
    %   index corresponding to the radar at the specified column of the
    %   radar position matrix, rpos. rpos is passed to the input of the
    %   step method. The default value of this property is 1. This property
    %   is tunable.
    ReferenceRadar
    %ShowBeam Enable beam visualization
    %   Specify as 'None', 'ReferenceRadar', or 'All' which beams to
    %   visualize. When set to 'None' no beams are visualized, when set to
    %   'ReferenceRadar' the beam of the radar specified in the
    %   ReferenceRadar property is visualized, and when set to 'All' the
    %   beams on all radars are visualized. The default value of this
    %   property is 'ReferenceRadar'. This property is tunable.
    ShowBeam
    %BeamWidth Beam width of radars
    %   Specify as a scalar, an Nr element row vector, a 2 element column
    %   vector, or a 2-by-Nr matrix in degrees the beam width of the
    %   radars. The first row specifies the horizontal beam width and the
    %   second row specifies the vertical beam width. Each column
    %   corresponds to a radar. Nr is the number of radars. When one row is
    %   specified the vertical and horizontal beam width are assumed the
    %   same. When one column is specified then the beam widths of all
    %   radars are assumed the same. This property is only used for
    %   visualization. The angles are positive values less than 360
    %   degrees. Also this value is used to calculate the CameraViewAngle
    %   property when the CameraPerspective property is set to Radar. The
    %   default value of this property is 15. This property is tunable.
    BeamWidth
    %BeamRange Beam range of radars
    %   Specify as a positive scalar or an Nr row vector the beam range in
    %   meters of the radars. Nr is the number of radars. When specified as
    %   a vector, each column corresponds to a radar. When specified as a
    %   scalar, all radars are assumed to have the same beam range. This
    %   property is only used for visualization. The default value of this
    %   property is 1000. This property is tunable.
    BeamRange
    %BeamSteering Beam direction of radars
    %   Specify in degrees as a 2 element column vector or 2-by-Nr matrix
    %   in degrees the direction of the beam relative to the x-axis of the
    %   radars. The first row specifies the azimuth angle and the second
    %   row specifies the elevation angle. Nr is the number of radars. When
    %   one column is specified then the beam steering of all radars are
    %   assumed the same. The azimuth angle must be between -180 and 180
    %   degrees, and the elevation angle must be between -90 and 90
    %   degrees. The default value is this property is [0;0]. This property
    %   is tunable.
    BeamSteering
end
properties (Nontunable, Logical)
    %VelocityInputPort Enable velocity input
    %   Set this property to true to pass the velocity vectors of all
    %   radars and targets to the step method. Otherwise they are estimated
    %   by the viewer based on the difference in positions between the
    %   current and last input. The default value is true.
    VelocityInputPort = true
    %OrientationInputPort Enable Orientation input
    %   Set this property to true to pass the orientation axis of all
    %   radars and targets to the step method. Otherwise they are
    %   aligned with the global coordinate axis. The default value is
    %   false.
    OrientationInputPort = false
end
properties(Nontunable)
    %UpdateRate Viewer update rate (Hz)
    %   Specify the update rate (in Hz) as a positive scalar. The default
    %   value of this property is 1.
    UpdateRate = 1
end
properties (Dependent)
    %Title Display title
    %   Specify the display title as a string. The default value of this
    %   property is ''. This property is tunable.
    Title
end
properties (Nontunable, Dependent)
    %PlatformNames Names of radars and targets
    %   Specify as 1 by N cell array the names assigned to the radars and
    %   targets. The names are ordered by first the radars and then the
    %   targets. The names will appear on the legend and annotations. When
    %   set to 'Auto', the names would append sequential numbers starting
    %   from 'Radar 1' for the radars and 'Target 1' for the targets. N is
    %   the total number of radars and targets. The default value of this
    %   property is 'Auto'.
    PlatformNames
end
properties (Dependent)
    %TrailLength Length of visible trajectories
    %   Specify as a positive scalar or vector of N elements the length in
    %   points of the visible trajectories. One point is generated every
    %   call to the step method. When specified as a scalar, all
    %   trajectories will have the same length, otherwise each element of
    %   the vector specifies the length of the corresponding radar or
    %   target. The elements are ordered by first the radars and then the
    %   targets. N is the total number of radars and targets. The default
    %   value of this property is 500. This property is tunable.
    TrailLength
    %CameraPerspective Perspective of the camera
    %   Specify as 'Auto', 'Custom', or 'Radar' the perspective of the
    %   camera. When set to 'Auto', appropriate values are guessed for the
    %   camera's position, orientation and view angle in order to show all
    %   trajectories. When set to 'Custom', the camera's position,
    %   orientation and angles are enabled in order to be modified by the
    %   camera toolbar or the camera properties. When set to 'Radar', the
    %   camera's position, orientation and angles are determined by the
    %   radar position and its beam direction. The default value of this
    %   property is 'Auto'. This property is tunable.
    CameraPerspective
    %CameraPosition Position of the camera
    %   Specify as a 3 element vector of the form [x y x] in meters the
    %   position of the camera. This property is applicable when the
    %   CameraPerspective property is set to 'Custom'. The default value is
    %   system dependent. This property is tunable.
    CameraPosition
    %CameraOrientation Orientation of the camera
    %   Specify as a 3 element vector of the form [pan tilt roll] in
    %   degrees the orientation of the camera. Panning, tilting then
    %   rolling are performed in that order. pan and roll can have values
    %   between -180 and 180 degrees. tilt can have values between -90 and
    %   90 degrees. This property is applicable when the CameraPerspective
    %   property is set to 'Custom'. The default value is system dependent.
    %   This property is tunable.
    CameraOrientation
    %CameraViewAngle View angle of the camera
    %   Specify as a positive scalar in degrees the view angle of the
    %   camera. The CameraViewAngle property should be less than 360
    %   degrees. This property is applicable when the CameraPerspective
    %   property is set to 'Custom'. The default value is system dependent.
    %   This property is tunable.
    CameraViewAngle
    %ShowLegend Show legend
    %   Specify as true to show the legend. The default value of this
    %   property is false. This property is tunable.
    ShowLegend
    %ShowGround Enable ground visualization
    %   Specify as true to visualize the ground. The default value of this
    %   property is true. This property is tunable.
    ShowGround
    %ShowName Annotate trajectories with names
    %   Specify as true to annotate the trajectories with the corresponding
    %   radar or target name. The default value of this property is true.
    %   This property is tunable.
    ShowName
    %ShowPosition Annotate trajectories with positions
    %   Specify as true to annotate the trajectories with the corresponding
    %   radar or target position. The default value of this property is
    %   false. This property is tunable.
    ShowPosition
    %ShowRange Annotate trajectories with range
    %   Specify as true to annotate the trajectories with the corresponding
    %   range relative to the reference radar. The default value of this
    %   property is false. This property is tunable.
    ShowRange
    %ShowAltitude Annotate trajectories with altitude
    %   Specify as true to annotate the trajectories with the corresponding
    %   altitude. The default value of this property is false.
    %   This property is tunable.
    ShowAltitude
    %ShowSpeed Annotate trajectories with speed
    %   Specify as true to annotate the trajectories with the corresponding
    %   speed. The default value of this property is false.
    %   This property is tunable.
    ShowSpeed
    %ShowRadialSpeed Annotate trajectories with radial speed
    %   Specify as true to annotate the trajectories with the corresponding
    %   radial speed relative to the reference radar. The default value of
    %   this property is false. This property is tunable.
    ShowRadialSpeed
    %ShowAzEl Annotate trajectories with azimuth and elevation
    %   Specify as true to annotate the trajectories with the corresponding
    %   azimuth and elevation relative to the reference radar. The default
    %   value of this property is false. This property is tunable.
    ShowAzEl
    %ReducePlotRate Reduce plot rate to improve performance
    %   Set this property to false to update the viewer every time the step
    %   method is called. This will negatively impact performance. The
    %   default is true. This property is tunable.
    ReducePlotRate;
end
properties (Access = private)
    pNumRadars = 0;
    pNumTargets = 0;
end
properties(Constant, Hidden)
    CameraPerspectiveSet = matlab.system.StringSet({'Auto','Custom','Radar'});
    ShowBeamSet = matlab.system.StringSet({'None','ReferenceRadar','All'});
end
methods
    function this = ScenarioViewer(varargin)
        this@matlabshared.scopes.UnifiedSystemScope();
        this.Position = uiscopes.getDefaultPosition([800 450+30]);
        setProperties(this, nargin, varargin{:});
    end
    function set.Name(this, value)
        setScopeName(this, value);
        this.Name = value;
    end
    function set.ReferenceRadar(this,value)
        if this.pNumRadars
            extraChecks = {'scalar','<=',this.pNumRadars};
        else
            extraChecks = {'scalar'};
        end
        sigdatatypes.validateIndex(value,'','ReferenceRadar',extraChecks);
        this.getFramework.Visual.ReferenceRadar = value;
    end
    function value = get.ReferenceRadar(this)
        value = this.getFramework.Visual.ReferenceRadar;
    end
    function set.BeamWidth(this,value)
        sigdatatypes.validateAngle(value,'','BeamWidth',{'positive','<=',360});
        if isLocked(this)
            if size(value,1) > 2 || (size(value,2) ~= this.pNumRadars && size(value,2) ~= 1)
                if this.pNumRadars == 1
                    error(message('phased:scopes:expectedScalarVector','BeamWidth','column',2));
                else
                    error(message('phased:scopes:expectedScalarVectorMatrix','BeamWidth',2,this.pNumRadars));
                end
            end
        end
         this.getFramework.Visual.BeamWidth = value;
    end
    function value = get.BeamWidth(this)
            value = this.getFramework.Visual.BeamWidth;
    end
    function set.BeamRange(this,value)
        sigdatatypes.validateDistance(value,'','BeamRange');
        if isLocked(this)
            if size(value,1) ~= 1 || (size(value,2) ~= this.pNumRadars && size(value,2) ~= 1)
                if this.pNumRadars == 1
                    error(message('phased:scopes:expectedScalar','BeamRange'));
                else
                    error(message('phased:scopes:expectedScalarVector','BeamRange','row',this.pNumRadars));
                end
            end
        end
        this.getFramework.Visual.BeamRange = value;
    end
    function value = get.BeamRange(this)
        value = this.getFramework.Visual.BeamRange;
    end
    function set.BeamSteering(this,value)
        sigdatatypes.validateAzElAngle(value,'','BeamSteering');
        if isLocked(this)
            if size(value,2) ~= this.pNumRadars && size(value,2) ~= 1
                if this.pNumRadars == 1
                    error(message('phased:scopes:expectedOneColumn','BeamSteering'));
                else
                    error(message('phased:scopes:expectedOneOrNumColumns','BeamSteering',this.pNumRadars));
                end
            end
        end
        this.getFramework.Visual.BeamSteering = value;
    end
    function value = get.BeamSteering(this)
        value = this.getFramework.Visual.BeamSteering;
    end

    function set.VelocityInputPort(this, value)
        validateattributes(value,{'logical'}, {'scalar'},'','VelocityInputPort');
        this.VelocityInputPort = value;
        this.getFramework.Visual.VelocityInputPort = value;
    end
    function set.OrientationInputPort(this, value)
        validateattributes(value,{'logical'}, {'scalar'},'','VelocityInputPort');
        this.OrientationInputPort = value;
        this.getFramework.Visual.OrientationInputPort = value;
    end

    function set.Title(this,value)
        val = convertStringsToChars(value);
        if  ~ischar(val)
            error(message('phased:scopes:InvalidString','Title'));
        end
        this.getFramework.Visual.Title = val;
    end
    function value = get.Title(this)
        value = this.getFramework.Visual.Title;
    end

    function set.PlatformNames(this, value)
        if ~((iscellstr(value) || isstring(value)) && isvector(value)) && ~(ischar(value) && strcmp(value,'Auto'))
            error(message('phased:scopes:InvalidAutoCellString','PlatformNames'));
        end
        this.getFramework.Visual.PlatformNames = value;
    end
    function value = get.PlatformNames(this)
        value = this.getFramework.Visual.PlatformNames;
    end
    function set.TrailLength(this,value)
        sigdatatypes.validateIndex(value,'','TrailLength',{'vector'});
        numPlatf =  this.pNumRadars + this.pNumTargets;
        if numPlatf && ~isscalar(value) && numel(value) ~= numPlatf
            error(message('phased:scopes:expectedScalarVector','TrailLength','row',numPlatf));
        end
        this.getFramework.Visual.TrailLength = value;
    end
    function value = get.TrailLength(this)
        value = this.getFramework.Visual.TrailLength;
    end

    function set.CameraPerspective(this,value)
        %Redundant
        %validatestring(value,this.CameraPerspectiveSet.getAllowedValues,'','CameraPerspective',2)
        idx = this.CameraPerspectiveSet.getIndex(value);
        this.getFramework.Visual.CameraPerspective = idx;
    end
    function value = get.CameraPerspective(this)
        idx = this.getFramework.Visual.CameraPerspective;
        value = this.CameraPerspectiveSet.getValueFromIndex(idx);
    end
    function set.CameraPosition(this,value)
        validateattributes(value,{'double'}, ...
            {'finite','nonnan','nonempty','real','vector','numel',3}, ...
            '','CameraPosition');
        this.getFramework.Visual.CameraPosition = value;
    end
    function value = get.CameraPosition(this)
        value = this.getFramework.Visual.CameraPosition;
    end
    function set.CameraOrientation(this,value)
        sigdatatypes.validateAngle(value,'','CameraOrientation',{'vector','numel',3});
        validateattributes(value(1),{'double'},{'scalar','<=',180,'>=',-180},'',...
            'pan angle');

        validateattributes(value(2),{'double'},{'scalar','<=',90,'>=',-90},'',...
            'tilt angle');

        validateattributes(value(3),{'double'},{'scalar','<=',180,'>=',-180},'',...
            'roll angle');

        this.getFramework.Visual.CameraOrientation = value;
    end
    function value = get.CameraOrientation(this)
        value = this.getFramework.Visual.CameraOrientation;
    end
    function set.CameraViewAngle(this,value)
        sigdatatypes.validateAngle(value,'','CameraViewAngle',{'positive','scalar','<',360});
        this.getFramework.Visual.CameraViewAngle = value;
    end
    function value = get.CameraViewAngle(this)
        value = this.getFramework.Visual.CameraViewAngle;
    end

    function set.ShowBeam(this,value)
        idx = this.ShowBeamSet.getIndex(value);
        this.getFramework.Visual.ShowBeam = idx;
    end
    function value = get.ShowBeam(this)
        idx = this.getFramework.Visual.ShowBeam;
        value = this.ShowBeamSet.getValueFromIndex(idx);
    end
    function set.ShowLegend(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowLegend');
        this.getFramework.Visual.ShowLegend = value;
    end
    function value = get.ShowLegend(this)
        value = this.getFramework.Visual.ShowLegend;
    end
    function set.ShowGround(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowGround');
        this.getFramework.Visual.ShowGround = value;
    end
    function value = get.ShowGround(this)
        value = this.getFramework.Visual.ShowGround;
    end
    function set.ShowName(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowName');
        this.getFramework.Visual.ShowName = value;
    end
    function value = get.ShowName(this)
        value = this.getFramework.Visual.ShowName;
    end
    function set.ShowPosition(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowPosition');
        this.getFramework.Visual.ShowPosition = value;
    end
    function value = get.ShowPosition(this)
        value = this.getFramework.Visual.ShowPosition;
    end
    function set.ShowRange(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowRange');
        this.getFramework.Visual.ShowRange = value;
    end
    function value = get.ShowRange(this)
        value = this.getFramework.Visual.ShowRange;
    end
    function set.ShowAltitude(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowAltitude');
        this.getFramework.Visual.ShowAltitude = value;
    end
    function value = get.ShowAltitude(this)
        value = this.getFramework.Visual.ShowAltitude;
    end
    function set.ShowSpeed(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowSpeed');
        this.getFramework.Visual.ShowSpeed = value;
    end
    function value = get.ShowSpeed(this)
        value = this.getFramework.Visual.ShowSpeed;
    end
    function set.ShowRadialSpeed(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowRadialSpeed');
        this.getFramework.Visual.ShowRadialSpeed = value;
    end
    function value = get.ShowRadialSpeed(this)
        value = this.getFramework.Visual.ShowRadialSpeed;
    end
    function set.ShowAzEl(this,value)
        validateattributes(value,{'logical'}, {'scalar'},'','ShowAzEl');
        this.getFramework.Visual.ShowAzEl = value;
    end
    function value = get.ShowAzEl(this)
        value = this.getFramework.Visual.ShowAzEl;
    end
    function set.ReducePlotRate(obj, value)
        validateattributes(value,{'logical'}, {'scalar'},'','ReducePlotRate');
        setPropertyValue(obj.getSource, 'LockSynchronous', ~value);
    end
    function value = get.ReducePlotRate(obj)
        value = ~getPropertyValue(obj.getSource, 'LockSynchronous');
    end
end
methods (Access = protected)
    function h = getScopeCfg(~)
        h = phased.scopes.ScenarioViewerSpecification('AppName','Scenario Viewer');
    end
    function launchScope(this)
        launchScope@matlabshared.scopes.UnifiedSystemScope(this);
        this.getFramework.Handles.toolsMenu.Visible = 'off';
    end

    function num = getNumInputsImpl(this)
        num = 2;
        if this.VelocityInputPort
            num = num+2;
        end
        if this.OrientationInputPort
            num = num+2;
        end
        % update visual source with correct number
        % I expected this to be done automatically
        % by the uiscope infrastructure
        if isScopeLaunched(this)
            this.pSource.NumInputPorts = num;
        end
    end
    function setupImpl(this, varargin)
        this.pNumRadars = size(varargin{1},2);
        if this.VelocityInputPort && this.OrientationInputPort
            tposIdx = 4;
        elseif this.VelocityInputPort || this.OrientationInputPort
            tposIdx = 3;
        else
            tposIdx = 2;
        end
        %tposIdx tracks position of target position argument
        this.pNumTargets = size(varargin{tposIdx},2);
        setupImpl@matlabshared.scopes.UnifiedSystemScope(this, varargin{:});
    end
    
    function  releaseImpl(this)
        this.pNumRadars = 0;
        this.pNumTargets = 0;
        if isvalid(this.getFramework)
            release(this.getFramework.Visual);
            releaseImpl@matlabshared.scopes.UnifiedSystemScope(this);
        end
    end

    
    function validateInputsImpl(this,varargin)

        rArgs = {'rpos'};
        tArgs = {'tpos'};
        if this.VelocityInputPort
            rArgs{end+1} = 'rvel';
            tArgs{end+1} = 'tvel';
        end
        if this.OrientationInputPort
            rArgs{end+1} = 'rlaxis';
            tArgs{end+1} = 'tlaxis';
        end
        inArgs = [rArgs tArgs];


        for k = 1:numel(inArgs)
            value = varargin{k};
            argStr = inArgs{k};
            refArg = [argStr(1) 'pos'];
            expNumPlat = size(varargin{strcmp(inArgs,refArg)},2);
            switch argStr
                case {'rlaxis','tlaxis'}
                    validateattributes(value, ...
                        {'double'}, {'real','finite','nonnan','size',[3 3 expNumPlat]},'',argStr);
                otherwise
                    validateattributes(value, ...
                        {'double'}, {'real','finite','nonnan','size',[3 expNumPlat]},'',argStr);
            end
        end

        numRadars = size(varargin{1},2);
        numTargets = size(varargin{strcmp(inArgs,'tpos')},2);
        numPlatforms = numRadars + numTargets;
        extraChecks = {'scalar','<=',numRadars};
        sigdatatypes.validateIndex(this.ReferenceRadar,'','ReferenceRadar',extraChecks);
        if iscell(this.PlatformNames)
            validateattributes(this.PlatformNames,{'cell'}, {'numel',numPlatforms},'','PlatformNames');
        end
        if  ~isscalar(this.TrailLength) && numel(this.TrailLength) ~= numPlatforms
            error(message('phased:scopes:expectedScalarVector','TrailLength','row',numPlatforms));
        end
        value = this.BeamWidth;
        if size(value,1) > 2 || (size(value,2) ~= numRadars && size(value,2) ~= 1)
            if numRadars == 1
                error(message('phased:scopes:expectedScalarVector','BeamWidth','column',2));            
            else
                error(message('phased:scopes:expectedScalarVectorMatrix','BeamWidth',2,numRadars));
            end
        end
        value = this.BeamRange;
        if size(value,1) ~= 1 || (size(value,2) ~= numRadars && size(value,2) ~= 1)
            if numRadars == 1
                error(message('phased:scopes:expectedScalar','BeamRange'));
            else
                error(message('phased:scopes:expectedScalarVector','BeamRange','row',numRadars));
            end
        end
        value = this.BeamSteering;
        if size(value,2) ~= numRadars && size(value,2) ~= 1
            if numRadars == 1
                error(message('phased:scopes:expectedOneColumn','BeamSteering'));
            else
                error(message('phased:scopes:expectedOneOrNumColumns','BeamSteering',numRadars));
            end
        end

    end
    function loadObjectImpl(this, S, wasLocked)
        this.VelocityInputPort = S.VelocityInputPort;
        this.OrientationInputPort = S.OrientationInputPort;
        this.UpdateRate = S.UpdateRate;
        if wasLocked
            this.pNumRadars = S.pNumRadars;
            this.pNumTargets = S.pNumTargets;
        end
        this.BeamWidth = S.BeamWidth;
        this.ReferenceRadar = S.ReferenceRadar;
        this.PlatformNames = S.PlatformNames;
        this.ShowBeam = S.ShowBeam;
        this.BeamRange = S.BeamRange;
        this.BeamSteering = S.BeamSteering;
        this.TrailLength = S.TrailLength;
        loadObjectImpl@matlabshared.scopes.UnifiedSystemScope(this,S,wasLocked);
    end
    function S = saveObjectImpl(this)
        S = saveObjectImpl@matlabshared.scopes.UnifiedSystemScope(this);
        S.VelocityInputPort = this.VelocityInputPort;
        S.OrientationInputPort = this.OrientationInputPort;
        S.UpdateRate = this.UpdateRate;
        S.BeamWidth = this.BeamWidth;
        S.ReferenceRadar = this.ReferenceRadar;
        S.PlatformNames = this.PlatformNames;
        S.ShowBeam = this.ShowBeam;
        S.BeamRange = this.BeamRange;
        S.BeamSteering = this.BeamSteering;
        S.TrailLength = this.TrailLength;
        if isLocked(this)
            S.pNumRadars = this.pNumRadars;
            S.pNumTargets = this.pNumTargets;
        end        
    end
    function flag = isInactivePropertyImpl(this,prop)
        if any(strcmp(prop,{'CameraPosition','CameraOrientation','CameraViewAngle'})) ...
                && ~strcmp(this.CameraPerspective,'Custom')
            flag = true;
        else
            flag = false;
        end
    end
end
methods (Hidden)
    function ax = getAxes(this)
        ax = this.getFramework.Visual.Axes;
    end
    function st = getInputSampleTime(this)
        st = 1/this.UpdateRate;
    end
end

methods (Static, Hidden)    
    function flag = isAllowedInSystemBlock(~)
        flag = false;
    end
end
end
