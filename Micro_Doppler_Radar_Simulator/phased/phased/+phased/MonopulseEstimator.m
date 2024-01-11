classdef (Sealed,StrictDefaults) MonopulseEstimator < phased.internal.AbstractNarrowbandArrayProcessing & ...
                matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates
%MonopulseEstimator     Amplitude monopulse estimator
%   H = phased.MonopulseEstimator creates a System object, H, that
%   estimates the target direction using amplitude monopulse technique
%   based on the sum and difference signals.
%
%   H = phased.MonopulseEstimator(Name,Value) returns an amplitude
%   monopulse tracker object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   ANG = step(H,SUMCH,DAZCH,STEER) estimates the target direction, ANG (in
%   degrees), using the sum and azimuth difference channel signals
%   specified in SUMCH and DAZCH, respectively. Both SUMCH and DAZCH are
%   Nx1 vectors where N is the number of snapshots in the signal. STEER is
%   a scalar representing the azimuth steering angle (in degrees) of the
%   array when the sum and difference channels are formed. If you set the
%   OutputFormat property to 'Angle', ANG is a 1xN vector containing the
%   estimated azimuth angle (in degrees) for the corresponding signal
%   snapshots. If you set the OutputFormat property to 'Angle offset', ANG
%   is a 1xN vector containing the estimated azimuth angle offsets (in
%   degrees) against the steering angle for the corresponding signal
%   snapshots. This syntax is applicable when you set the Coverage property
%   to 'Azimuth'.
%
%   ANG = step(H,SUMCH,DAZCH,DELCH,STEER) estimates the target direction,
%   ANG (in degrees), using the sum, azimuth difference, and elevation
%   difference channel signals specified in SUMCH, DAZCH, and DELCH,
%   respectively. SUMCH, DAZCH, and DELCH are all Nx1 vectors where N is
%   the number of signal snapshots. STEER is a 2x1 column vector in the
%   form of [azimuth; elevation] representing the steer angle (in degrees)
%   of the array when sum and difference channels are formed. If you set
%   the OutputFormat to 'Angle', ANG is a 2xN matrix containing the
%   estimated azimuth and elevation angles (in degrees), with each column
%   in the form of [azimuth; elevation], of the corresponding signal
%   snapshots. If you set the OutputFormat to 'Angle offset', ANG is a 2xN
%   matrix containing the estimated azimuth and elevation angle offsets (in
%   degrees), with each column in the form of [azimuth; elevation], of the
%   corresponding signal snapshots. This syntax is applicable when you set
%   the Coverage property to '3D'.
%
%   [...,DRATIO] = step(...) also returns the sum and difference ratio in
%   DRATIO. When you set the Coverage property to 'Azimuth', DRATIO is a
%   1xN vector containing the ratios of azimuth difference channel and sum
%   channel for the corresponding signal snapshots. When you set the
%   Coverage property to '3D', DRATIO is a 2xN matrix. Its first row
%   contains the ratios of azimuth difference channel and sum channel for
%   the corresponding signal snapshots. Its second row contains the ratios
%   of elevation difference channel and sum channel for the corresponding
%   signal snapshots.
%   
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [ANG,DRATIO] = step(H,SUMCH,DAZCH,DELCH,STEER)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   MonopulseEstimator methods:
%   
%   step         - Estimate the target direction
%   release      - Allow property value and input characteristics changes
%   clone        - Create monopulse estimator object with same property
%                  values
%   isLocked     - Locked status (logical)
%
%   MonopulseEstimator properties:
%       
%   SensorArray                  - Sensor array
%   PropagationSpeed             - Signal propagation speed 
%   OperatingFrequency           - Operating frequency 
%   Coverage                     - Monopulse coverage
%   SquintAngle                  - Squint angle 
%   OutputFormat                 - Output format
%   SumDifferenceRatioOutputPort - Output sum difference ratio
%
%   % Example:
%   %   Determine the direction of a target using URA. The target echo is
%   %   first processed to identify the interested range cell before
%   %   applying monopulse processing.
% 
%   array = phased.URA('Size',4); 
%   collect = phased.Collector('Sensor',array);
%   feed = phased.MonopulseFeed('SensorArray',array,...
%           'Coverage','3D');
%   estimator = phased.MonopulseEstimator('SensorArray',array,...
%           'Coverage','3D');
%
%   % Source signal
%   x = sqrt(0.01/2)*(randn(100,1)+1i*randn(100,1));
%   x(20) = 1;
%   targetangle = [31;9];
%   steerangle = [30;10];
%
%   % Form echo and generate sum and difference channels
%   rx = collect(x,targetangle);
%   [sumch,azch,elch] = feed(rx,steerangle);
%
%   % Detection
%   [~,idx] = max(abs(sumch));
%
%   % Monopulse
%   est_dir = estimator(sumch(idx),azch(idx),elch(idx),steerangle)
%
%   See also phased, phased.MonopulseFeed,
%   phased.SumDifferenceMonopulseTracker,
%   phased.SumDifferenceMonopulseTracker2D.

%   Copyright 2018 The MathWorks, Inc.

%   References
%   [1] S. M. Sherman, Monopulse Principles and Techniques, Chapter in
%   Aspects of Modern Radar, Eli Brookner, Artech House, 1988

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    % Public, tunable properties
    properties (Nontunable)
        %Coverage   Monopulse coverage
        %   Specify the coverage of the corresponding monopulse feed as one
        %   of '3D' | 'Azimuth' where the default is '3D'. When you set
        %   this property to '3D', the inputs to the monopulse estimator
        %   contain the sum channel as well as both azimuth and elevation
        %   difference channel. When you set this property to 'Azimuth',
        %   the inputs to the monopulse estimator contain only the sum
        %   channel and the azimuth difference channel.
        Coverage = '3D'     % monopulse coordinates    
        %SquintAngle    Squint angle (degrees)
        %   Specify the squint angle (in degrees) as a scalar or a 2x1
        %   vector. The squint angle is the separation angle between the
        %   two beams along the azimuth and elevation directions. When you
        %   set the Coverage property to 'Azimuth', the SquintAngle
        %   property must be a scalar. When you set the Coverage property
        %   to '3D', if the SquintAngle property is a scalar, then the
        %   squint angle is the same along both the azimuth and elevation
        %   directions; if the SquintAngle property is a 2x1 vector, then
        %   its entries specify the squint angle along the azimuth and
        %   elevation directions, respectively. The default value of this
        %   property is 10.
        SquintAngle = 10    % squint angle [az; el] (degree)   
        %OutputFormat   Output format
        %   Specify the format of output as one of 'Angle' | 'Angle offset'
        %   where the default is 'Angle'. When you set this property to
        %   'Angle', the output is the direction of the target. When you
        %   set this property to 'Angle offset', the output is the angle
        %   offset from the array' steering direction.
        OutputFormat = 'Angle'
    end
    
    properties (Nontunable, Logical)
        %SumDifferenceRatioOutputPort  Output sum difference ratio
        %   Set this property to true to output the sum and difference
        %   ratios in azimuth and elevation directions. It is often used as
        %   an error control signal. Set this property to false to not
        %   output the sum and difference ratios.
        SumDifferenceRatioOutputPort = false
    end

    % Pre-computed constants
    properties(Access = private)
        cSteeringVector
        pTraverseGrid
        pElevationGrid
        pTraverseDictionary
        pElevationDictionary
        pSquintAngle
    end
    
    properties (Constant, Hidden)
        CoverageSet = matlab.system.StringSet(...
            {'Azimuth','3D'});
        OutputFormatSet = matlab.system.StringSet(...
            {'Angle','Angle offset'});
    end
    
    methods
       function obj = MonopulseEstimator(varargin)
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
       end               
    end

    methods
        function set.SquintAngle(obj,val)
            sigdatatypes.validateAngle(val,'','SquintAngle',{'positive','<=',90,'column'});
            obj.SquintAngle = val;
        end
    end
    
    methods(Access = protected)
        
        function privValidateSensorArray(~,val) 
        %privValidateSensorArray
        %   Certain array operation is limited to certain array geometries.
        %   Each operation can then overload this method to do its own
        %   validation. By default, any array geometry is ok.
        
            validateattributes( val, ...
                { 'phased.internal.AbstractArray', 'phased.internal.AbstractSubarray' },...
                { 'scalar' }, '', 'SensorArray');
        end
        
        function num = getNumInputsImpl(obj) 
            switch obj.Coverage
                case 'Azimuth'
                    num = 3;
                case '3D'
                    num = 4;
            end
        end
        
        function num = getNumOutputsImpl(obj) 
            switch obj.Coverage
                case 'Azimuth'
                    num = 1;
                    if obj.SumDifferenceRatioOutputPort
                        num = num+1;
                    end
                case '3D'
                    num = 1;
                    if obj.SumDifferenceRatioOutputPort
                        num = num+1;
                    end
            end
        end
        
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            if isLocked(obj)
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.pTraverseGrid = obj.pTraverseGrid;
                s.pElevationGrid = obj.pElevationGrid;
                s.pTraverseDictionary = obj.pTraverseDictionary;
                s.pElevationDictionary = obj.pElevationDictionary;
                s.pSquintAngle = obj.pSquintAngle;
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            s = loadSubObjects@phased.internal.AbstractNarrowbandArrayProcessing(obj,s);
            if wasLocked
                obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                s = rmfield(s,'cSteeringVector');
            end
            if isfield(s,'isLocked')
                s = rmfield(s,'isLocked');
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            
            s = loadSubObjects(obj,s,wasLocked);
            
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            resetImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            reset(obj.cSteeringVector);
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            releaseImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            release(obj.cSteeringVector);
        end

        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
            switch obj.Coverage
                case 'Azimuth'
                    cond = ~isscalar(obj.SquintAngle);
                    if cond
                        coder.internal.errorIf(cond,'MATLAB:system:inputMustBeScalar','SquintAngle');
                    end
                        
                case '3D'
                    cond = ~isscalar(obj.SquintAngle) && (length(obj.SquintAngle)~=2);
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:invalidColumnNumbers','SquintAngle',2);
                    end
            end
        end

        function validateInputsImpl(obj,sumch,varargin)
            % Validate inputs to the step method at initialization
            validateattributes(sumch,{'double'},{'column'},'','SIGMA');
            switch obj.Coverage
                case 'Azimuth'
                    deltaAz = varargin{1};
                    steerang = varargin{2};
                    validateattributes(deltaAz,{'double'},{'column'},'','DeltaAz');
                    cond = ~isequal(size(sumch),size(deltaAz));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:sizeMismatch','SIGMA','DeltaAz');
                    end
                    halfsqang = obj.SquintAngle/2;
                    angbound = 90-halfsqang;
                    sigdatatypes.validateAngle(steerang,'','STEER',{'scalar','>=',-angbound,'<=',angbound}); 
                case '3D'
                    deltaAz = varargin{1};
                    deltaEl = varargin{2};
                    steerang = varargin{3};
                    validateattributes(deltaAz,{'double'},{'column'},'','DeltaAz');
                    validateattributes(deltaEl,{'double'},{'column'},'','DeltaEl');
                    cond = ~isequal(size(sumch),size(deltaAz));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:sizeMismatch','SIGMA','DeltaAz');
                    end
                    cond = ~isequal(size(sumch),size(deltaEl));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:sizeMismatch','SIGMA','DeltaEl');
                    end
                    halfsqang = obj.SquintAngle/2;
                    if isscalar(halfsqang)
                        angbound = ones(2,1)*(90-halfsqang);
                    else
                        angbound = 90-halfsqang;
                    end
                    sigdatatypes.validateAngle(steerang,'','STEER',{'column','nrows',2,});
                    cond = (steerang(1)<-angbound(1)) || (steerang(1)>angbound(1)) || ...
                        (steerang(2)<-angbound(2)) || (steerang(2)>angbound(2));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:invalidAzElAngle',...
                            'Steer',num2str(-angbound(1)),num2str(angbound(1)),...
                            num2str(-angbound(2)),num2str(angbound(2)));
                    end
            end
            validateNumChannels(obj,sumch)
        end
        
        function varargout = trainingMonopulse(obj)
            
            fc = obj.OperatingFrequency;
            steervec = obj.cSteeringVector;
            squint = obj.pSquintAngle;
            lambda = obj.PropagationSpeed/fc;
            elempos = getElementPosition(obj.SensorArray);
            arrayspan = max(elempos,[],2)-min(elempos,[],2);
            gridN = 1024;
            switch obj.Coverage
                case 'Azimuth'
                    azbw = lambda/arrayspan(2);
                    tsinGrid = linspace(-1/2,1/2,gridN)*azbw;
                    tgrid = asind(tsinGrid);
                    nTr = numel(tgrid);
                    azSquint = squint/2;
                    
                    azstvmat = steervec(fc, [-azSquint azSquint tgrid; zeros(1,nTr+2)]);
                    sA = azstvmat(:,1);
                    sB = azstvmat(:,2);
                    patSum = sA+sB;
                    patDeltaTr = sA-sB;
                    sumbeam = patSum'*azstvmat(:,3:end);
                    deltaTr = patDeltaTr'*azstvmat(:,3:end);
                    tdict = real(deltaTr./sumbeam);
                    varargout{1} = tsinGrid;
                    varargout{2} = tdict;

                    
                case '3D'
                    azbw = lambda/arrayspan(2);
                    elbw = lambda/arrayspan(3);
                    azSquint = squint(1)/2;
                    elSquint = squint(2)/2;
                    
                    tsinGrid = linspace(-1/2,1/2,gridN)*azbw;
                    tgrid = asind(tsinGrid);
                    nTr = numel(tgrid);
                    azstvmat = steervec(fc, ...
                        [-azSquint azSquint -azSquint azSquint tgrid; ...
                        elSquint elSquint -elSquint -elSquint zeros(1,nTr)]);
                    sA = azstvmat(:,1);
                    sB = azstvmat(:,2);
                    sC = azstvmat(:,3);
                    sD = azstvmat(:,4);
                    patSum = sA+sB+sC+sD;
                    patDeltaTr = (sA+sC) - (sB+sD);
                    %patDeltaEl = (sA+sB - (sC+sD));
                    sumbeam = patSum'*azstvmat(:,5:end);
                    deltaTr = patDeltaTr'*azstvmat(:,5:end);
                    tdict = real(deltaTr./sumbeam);
                    varargout{1} = tsinGrid;
                    varargout{2} = tdict;

                    
                    esinGrid = linspace(-1/2,1/2,gridN)*elbw;
                    egrid = asind(esinGrid);
                    nEl = numel(egrid);
                    elstvmat = steervec(fc, ...
                        [-azSquint azSquint -azSquint azSquint zeros(1,nEl); ...
                        elSquint elSquint -elSquint -elSquint egrid]);
                    sA = elstvmat(:,1);
                    sB = elstvmat(:,2);
                    sC = elstvmat(:,3);
                    sD = elstvmat(:,4);
                    patSum = sA+sB+sC+sD;
                    %patDeltaTr = (sA+sC - (sB+sD));
                    patDeltaEl = (sA+sB - (sC+sD));
                    sumbeam = patSum'*elstvmat(:,5:end);
                    deltaEl = patDeltaEl'*elstvmat(:,5:end);
                    edict = real(deltaEl./sumbeam);
                    varargout{3} = esinGrid;
                    varargout{4} = edict;
            end
            
        end
        
        function setupImpl(obj,sumch,~,~,~)
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,'PropagationSpeed',obj.PropagationSpeed);
            if strcmp(obj.Coverage,'3D')
                if isscalar(obj.SquintAngle)
                    obj.pSquintAngle = [obj.SquintAngle;obj.SquintAngle];
                else
                    obj.pSquintAngle = obj.SquintAngle;
                end
            else
                obj.pSquintAngle = obj.SquintAngle;
            end
            
            switch obj.Coverage
                case 'Azimuth'
                    [tgrid,tdict] = trainingMonopulse(obj);
                    obj.pTraverseGrid = tgrid;
                    obj.pTraverseDictionary = tdict;
                case '3D'
                    [tgrid,tdict,egrid,edict] = trainingMonopulse(obj);
                    obj.pTraverseGrid = tgrid;
                    obj.pTraverseDictionary = tdict;
                    obj.pElevationGrid = egrid;
                    obj.pElevationDictionary = edict;
            end
            
            
            obj.pNumInputChannels = getNumChannels(obj,sumch);
            obj.pValidatedNumInputChannels = getNumChannels(obj,sumch);

        end

        function varargout = stepImpl(obj, sumch, dazch, delch, steering)
            % steering within bound

            switch obj.Coverage
                case 'Azimuth'
                    steerang = delch;
                    
                    azDev = real(dazch./sumch);
                    
                    tSinOffset = obj.getSinOffset(azDev,...
                        obj.pTraverseDictionary,obj.pTraverseGrid);
                    tAngEst = obj.addSinOffsetToRef(steerang,tSinOffset);
                    
                    if strcmp(obj.OutputFormat,'Angle')
                        varargout{1} = tAngEst.';
                    else
                        varargout{1} = tAngEst.'-steerang;
                    end
                    if obj.SumDifferenceRatioOutputPort
                        varargout{2} = azDev.';
                    end
                case '3D'
                    steerang = steering;
                    
                    azDev = real(dazch./sumch);
                    elDev = real(delch./sumch);
                    
                    elSinOffset = obj.getSinOffset(elDev,...
                        obj.pElevationDictionary,obj.pElevationGrid);
                    elAngEst = obj.addSinOffsetToRef(steerang(2),elSinOffset);
                                       
                    tSinOffset = obj.getSinOffset(azDev,...
                        obj.pTraverseDictionary,obj.pTraverseGrid);
                    azSinOffset = tSinOffset./cosd(elAngEst);
                    azAngEst = obj.addSinOffsetToRef(steerang(1),azSinOffset);
                    
                    if strcmp(obj.OutputFormat,'Angle')
                        varargout{1} = [azAngEst.';elAngEst.'];
                    else
                        varargout{1} = [azAngEst.';elAngEst.']-steerang;
                    end
                    if obj.SumDifferenceRatioOutputPort
                        varargout{2} = [azDev.';elDev.'];
                    end
            end

        end
    end
    
    methods (Access = protected)
        function flag = isInputComplexityLockedImpl(obj,index)  
            if strcmp(obj.Coverage,'Azimuth')
                if index == 3
                    flag = true;
                else
                    flag = false;
                end
            else % 3D
                if index == 4
                    flag = true;
                else % 
                    flag = false;
                end
            end
        end
        
        function icon = getIconImpl(obj) %#ok<MANU>
            % Define icon for System block
            icon = sprintf('Amplitude\nMonopulse\nEstimator');
        end

        function varargout = getInputNamesImpl(obj) 
            % Return input port names for System block
            switch obj.Coverage
                case 'Azimuth'
                    varargout = {'SIGMA','DeltaAz','STEER'};
                case '3D'
                    varargout = {'SIGMA','DeltaAz','DeltaEl','STEER'};
            end
        end

        function varargout = getOutputNamesImpl(obj)
            % Return output port names for System block
            switch obj.Coverage
                case 'Azimuth'
                    if strcmp(obj.OutputFormat,'Angle')
                        varargout{1} = 'Az';
                    else
                        varargout{1} = 'dAz';
                    end
                    if obj.SumDifferenceRatioOutputPort
                        varargout{2} = 'AzRatio';
                    end
                case '3D'
                    if strcmp(obj.OutputFormat,'Angle')
                        varargout{1} = 'AzEl';
                    else
                        varargout{1} = 'dAzEl';
                    end
                    if obj.SumDifferenceRatioOutputPort
                        varargout{2} = 'AzElRatio';
                    end
            end
        end

        function flag = isInputSizeLockedImpl(obj,ind) 
            switch obj.Coverage
                case 'Azimuth'
                    if ind == 3
                        flag = true;  % angle
                    else
                        flag = false;
                    end
                case '3D'
                    if ind == 4
                        flag = true;  % angle
                    else
                        flag = false;
                    end
            end
        end
        
        function varargout = getOutputSizeImpl(obj)
            % Return size for each output port
            sz = propagatedInputSize(obj,1);
            switch obj.Coverage
                case 'Azimuth'
                    varargout{1} = [1 sz(1)];
                    if obj.SumDifferenceRatioOutputPort
                        varargout{2} = [1 sz(1)];
                    end
                case '3D'
                    varargout{1} = [2 sz(1)];
                    if obj.SumDifferenceRatioOutputPort
                        varargout{2} = [2 sz(1)];
                    end
            end
            
        end

        function varargout = isOutputFixedSizeImpl(obj)
            flag = propagatedInputFixedSize(obj, 1);
            varargout{1} = flag;
            if obj.SumDifferenceRatioOutputPort
                varargout{2} = flag;
            end
        end
        
        function varargout = getOutputDataTypeImpl(obj)
            % Return data type for each output port
            varargout = cell(1,nargout);
            for k = 1:nargout
                varargout{k} = propagatedInputDataType(obj,1);
            end
        end

        function varargout = isOutputComplexImpl(obj)
            % Return true for each output port with complex data
            switch obj.Coverage
                case 'Azimuth'
                    varargout{1} = false;
                    if obj.SumDifferenceRatioOutputPort
                        varargout{2} = false;
                    end
                case '3D'
                    varargout{1} = false;
                    varargout{2} = false;
                    if obj.SumDifferenceRatioOutputPort
                        varargout{3} = false;
                        varargout{4} = false;
                    end
            end
        end
 

    end

    methods(Access = protected, Static)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractNarrowbandArrayProcessing('subarray');
            props = {...
                'Coverage',...
                'SquintAngle',...
                'OutputFormat',...
                'SumDifferenceRatioOutputPort'};
            groups(1).PropertyList = [groups(1).PropertyList props];
            
            action = matlab.system.display.Action(@phased.MonopulseEstimator.onActionCalled, ...
                'Label', 'Generate Monopulse Feed','Placement','last','Alignment','right');
            groups(1).Actions = action;
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:MonopulseEstimatorTitle')),...
                'Text',getString(message('phased:library:block:MonopulseEstimatorDesc')));
        end
        function onActionCalled(actionData, ~)
            SystemHandle = actionData.SystemHandle;
            sysname = get_param(SystemHandle,'Parent');
            b = add_block('phaseddoalib/Monopulse Feed',...
                [sysname '/Monopulse Feed'],'MakeNameUnique','on');
            set_param(b,'SensorArray',get_param(SystemHandle,'SensorArray'));
            set_param(b,'PropagationSpeed',get_param(SystemHandle,'PropagationSpeed'));
            set_param(b,'OperatingFrequency',get_param(SystemHandle,'OperatingFrequency'));
            set_param(b,'Coverage',get_param(SystemHandle,'Coverage'));
            set_param(b,'SquintAngle',get_param(SystemHandle,'SquintAngle'));
            estimator_name = get_param(SystemHandle,'Name');
            feed_name = get_param(b,'Name');
            add_line(sysname,sprintf('%s/1',feed_name),sprintf('%s/1',estimator_name),'AutoRouting','on');
            add_line(sysname,sprintf('%s/2',feed_name),sprintf('%s/2',estimator_name),'AutoRouting','on');
            add_line(sysname,sprintf('%s/3',feed_name),sprintf('%s/3',estimator_name),'AutoRouting','on');
        end
    end
    methods (Access = protected,Static)
        function sinOffset = getSinOffset(monopulseRatio,ratioGrid,sinOffsetGrid)
        %getSinOffset Obtain sine offset based on the monopulse ratio
        %   SINOFFSET = getSinOffset(R,RGRID,SINGRID) returns the offset
        %   SINOFFSET in the sin domain corresponding to the ratio R using the
        %   mapping between RGRID and SINGRID.
        %
        %   This is a static method.
            
            cond = any(abs(monopulseRatio) > abs(max(ratioGrid)));
            if cond
                coder.internal.errorIf(cond,'phased:SumDifferenceMonopulse:InvalidRatio','Steer');
            end
            sinOffset = interp1(ratioGrid,sinOffsetGrid,monopulseRatio);
        end

        function ang = addSinOffsetToRef(refAng,sinOffset)
        %addSinOffsetToRef Add sine offset to the reference angle
        %   ANG = addSinOffsetToRef(REFANG,SINOFFSET) returns the ANG (in
        %   degrees) that corresponding to the reference REFANG (in degrees)
        %   added with the offset SINOFFSET. Note that SINOFFSET is specified
        %   in the sine domain.
        %
        %   This is a static method.
            sinang = sind(refAng)*ones(1,size(sinOffset,2))+sinOffset;
            sinang(sinang>1) = 1;
            sinang(sinang<-1) = -1;
            ang = asind(sinang);
        end       
    end
end
