classdef (Sealed,StrictDefaults) MonopulseFeed < phased.internal.AbstractNarrowbandArrayProcessing & ...
                matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates
%MonopulseFeed     Amplitude monopulse feed
%   H = phased.MonopulseFeed creates a System object, H, that
%   simulates the feed system for the amplitude sum and difference
%   monopulse tracker. This object combines the received signal from an
%   array to sum and difference channels.
%
%   H = phased.MonopulseFeed(Name,Value) returns an amplitude
%   monopulse feed object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   [SUMCH,DAZCH] = step(H,X,STEER) forms both the sum channel, SUMCH, and
%   azimuth difference channel, DAZCH, from the received signal, X, at the
%   array elements. X is an N-column matrix where N represents the number
%   of elements of the receiving array or the number of subarrays if the
%   receiving array contains subarrays. STEER specifies the array's
%   steering direction as a 2x1 column vector in the form of [azimuth;
%   elevation] angles (in degrees). Both SUMCH and DAZCH are column vectors
%   and their number of rows are the same as the number of rows in X. This
%   syntax applies when the Coverage property is set to 'Azimuth'.
%
%   [SUMCH,DAZCH,DELCH] = step(H,X,STEER) also returns the elevation
%   difference channel, DELCH. DELCH is a column vector with the same
%   dimension as DAZCH. This syntax applies when the Coverage property is
%   set to '3D'.
%
%   [...,ANGEST] = step(...) returns the estimate of the target angle,
%   ANGEST, as a 2x1 vector in the form of [azimuth; elevation] angles (in
%   degree). This syntax applies when you set the AngleOutputPort property
%   to true.
%   
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [SUMCH,DAZCH,DELCH,ANGEST] = step(H,X,STEER)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   MonopulseFeed methods:
%   
%   step                  - Form sum and difference channel data
%   release               - Allow property value and input characteristics 
%                           changes
%   clone                 - Create monopulse feed object with same 
%                           property values
%   isLocked              - Locked status (logical)
%   getMonopulseEstimator - Create a monopulse estimator corresponding to 
%                           the feed
%
%   MonopulseFeed properties:
%       
%   SensorArray            - Sensor array
%   PropagationSpeed       - Signal propagation speed 
%   OperatingFrequency     - Operating frequency 
%   Coverage               - Monopulse coverage
%   SquintAngle            - Squint angle 
%   AngleOutputPort        - Output angle estimate
%
%   % Examples:
%
%   % Example 1:
%   %   Determine the direction of a target at around 60 degrees broadside
%   %   angle of a ULA
% 
%   array = phased.ULA('NumElements',16); 
%   steervector = phased.SteeringVector('SensorArray',array);
%   feed = phased.MonopulseFeed('SensorArray',array,...
%           'Coverage','Azimuth','AngleOutputPort',true);
%   x = steervector(feed.OperatingFrequency,60.1).';
%   [sumch,azch,est_dir] = feed(x,60)
%
%   % Example 2:
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
%   See also phased, phased.MonopulseEstimator,
%   phased.SumDifferenceMonopulseTracker,
%   phased.SumDifferenceMonopulseTracker2D.

%   Copyright 2018 The MathWorks, Inc.

%   References
%   [1] S. M. Sherman, Monopulse Principles and Techniques, Chapter in
%   Aspects of Modern Radar, Eli Brookner, Artech House, 1988


% Conceptually, the monopulse forms 4 closely-spacing beams; assume that
% the scan angle of the array, steering = [az, el], is placed at the origin
% of the azimuth-elevation plain, the center of the 4 beams, A,B,C,D, are
% located at [-squint(1), squint(2)], [squint(1), squint(2)], [-squint(1),
% -squint(2)], [squint(1), -squint(2)], respectively, as shown below.
% 
% |A|B|
% -----
% |C|D|  
%
% The beam patterns for the three channels, i.e., patSum, patDeltaAz,
% patDeltaEl are created by combining the 4 beams in such a way that patSum
% = sA+sB+sC+dD, patDeltaAz = sA+sC-(sB+sD), patDeltaEl = sA+sB-(sC+sD),
% where sA, sB, sC, sD are the steering vectors for the 4 beams,
% respectively.
% 
 

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    % Public, tunable properties
    properties (Nontunable)
        %Coverage   Monopulse coverage
        %   Specify the coverage of monopulse feed as one of '3D' |
        %   'Azimuth' where the default is '3D'. When you set this property
        %   to '3D', the monopulse feed forms the sum channel as well as
        %   both azimuth and elevation difference channel. When you set
        %   this property to 'Azimuth', the monopulse feed forms only the
        %   sum channel and the azimuth difference channel.
        Coverage = '3D'          % monopulse coordinates    
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
        SquintAngle = 10         % squint angle [az; el] (degree)   
    end
    
    properties (Nontunable, Logical)
        %AngleOutputPort    Output angle estimate
        %   Set this property to true to output the angle estimate in
        %   addition to sum and difference channels. Set this property to
        %   false to only output sum and difference channels. The default
        %   value of this property is false.
        AngleOutputPort = false
    end

    % Pre-computed constants
    properties(Access = private, Nontunable)
        cSteeringVector
        cMonopulseEstimator
        pSquintAngle
    end
    
    properties (Constant, Hidden)
        CoverageSet = matlab.system.StringSet(...
            {'Azimuth','3D'});
    end
    
    methods
        function obj = MonopulseFeed(varargin)
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
        end               
        
        function h = getMonopulseEstimator(obj)
        %getMonopulseEstimator Create a monopulse estimator corresponding to
        %the feed
        %   HE = getMonopulseEstimator(H) returns a amplitude monopulse
        %   estimator, HE, that corresponds to the amplitude monopulse
        %   feed, H.
        %
        %   % Example:
        %   %   Determine the direction of a target using URA. The target 
        %   %   echo is first processed to identify the interested range 
        %   %   cell before applying monopulse processing.
        % 
        %   array = phased.URA('Size',4); 
        %   collect = phased.Collector('Sensor',array);
        %   feed = phased.MonopulseFeed('SensorArray',array,...
        %           'Coverage','3D');
        %   estimator = getMonopulseEstimator(feed);
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
        %   See also phased.MonopulseFeed, phased.MonopulseEstimator.
        
            narginchk(1,1);
            validateProperties(obj);
            h = phased.MonopulseEstimator(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'OperatingFrequency',obj.OperatingFrequency,...
                'Coverage',obj.Coverage,...
                'SquintAngle',obj.SquintAngle);
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
        
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function num = getNumOutputsImpl(obj) 
            num = 1;
            switch obj.Coverage
                case 'Azimuth'
                    num = num+1;
                    if obj.AngleOutputPort
                        num = num+1;
                    end
                case '3D'
                    num = num+2;
                    if obj.AngleOutputPort
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
                s.cMonopulseEstimator = saveobj(obj.cMonopulseEstimator);
                s.pSquintAngle = obj.pSquintAngle;
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            s = loadSubObjects@phased.internal.AbstractNarrowbandArrayProcessing(obj,s);
            if wasLocked
                obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                s = rmfield(s,'cSteeringVector');
                obj.cMonopulseEstimator = phased.MonopulseEstimator.loadobj(s.cMonopulseEstimator);
                s = rmfield(s,'cMonopulseEstimator');
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
        
        
        function setupImpl(obj,x,~)
            setupImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
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
            
            if obj.AngleOutputPort 
                obj.cMonopulseEstimator = phased.MonopulseEstimator(...
                    'SensorArray',obj.SensorArray,'PropagationSpeed',obj.PropagationSpeed,...
                    'OperatingFrequency',obj.OperatingFrequency,'Coverage',obj.Coverage,...
                    'SquintAngle',obj.pSquintAngle,'OutputFormat','Angle');
            end
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            resetImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            reset(obj.cSteeringVector);
            if obj.AngleOutputPort
                reset(obj.cMonopulseEstimator);
            end
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            releaseImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            release(obj.cSteeringVector);
            if obj.AngleOutputPort
                release(obj.cMonopulseEstimator);
            end
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
                        coder.internal.errorIf(cond,'phased:phased:invalidRowNumbers','SquintAngle',2);
                    end
            end
        end

        function validateInputsImpl(obj,x,steerang)
            % Validate inputs to the step method at initialization
            validateattributes(x,{'double'},{'2d','ncols',getDOF(obj.SensorArray)},...
                '','X');
            switch obj.Coverage
                case 'Azimuth'
                    halfsqang = obj.SquintAngle/2;
                    angbound = 90-halfsqang;
                    sigdatatypes.validateAngle(steerang,'','STEER',{'scalar','>=',-angbound,'<=',angbound}); 
                case '3D'
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
            validateNumChannels(obj,x)
        end

        function  [sumCh, varargout] = stepImpl(obj, x, steerang)

            steervec = obj.cSteeringVector;
            fc = obj.OperatingFrequency;

            switch obj.Coverage
                case 'Azimuth'
                    azSteering = steerang;
                    azSquint = obj.pSquintAngle/2;
                    angA = azSteering-azSquint;
                    angB = azSteering+azSquint;
                    sA= steervec(fc, [angA; 0]);
                    sB= steervec(fc, [angB; 0]);
                    patSum = sA+sB;
                    patDeltaTr = sA-sB;

                    sumCh = x*conj(patSum);
                    deltaTr = x*conj(patDeltaTr);
                    
                    varargout{1} = deltaTr;
                    if obj.AngleOutputPort
                        az_est = step(obj.cMonopulseEstimator,sumCh,deltaTr,steerang);
                        varargout{2} = az_est;
                    end
                        
                case '3D'
                    azSteering = steerang(1);
                    elSteering = steerang(2);
                    azSquint = obj.pSquintAngle(1)/2;
                    elSquint = obj.pSquintAngle(2)/2;
                    angA = [azSteering-azSquint;elSteering+elSquint];
                    angB = [azSteering+azSquint;elSteering+elSquint];
                    angC = [azSteering-azSquint;elSteering-elSquint];
                    angD = [azSteering+azSquint;elSteering-elSquint];
                    
                    sA= steervec(fc, angA);
                    sB= steervec(fc, angB);
                    sC= steervec(fc, angC);
                    sD= steervec(fc, angD); 

                    patSum = sA+sB+sC+sD;
                    patDeltaTr = (sA+sC) - (sB+sD);
                    patDeltaEl = (sA+sB) - (sC+sD);

                    sumCh = x*conj(patSum);
                    deltaTr = x*conj(patDeltaTr);
                    deltaEl = x*conj(patDeltaEl);     
                    
                    varargout{1} = deltaTr;
                    varargout{2} = deltaEl;
                    if obj.AngleOutputPort
                        ang_est = step(obj.cMonopulseEstimator,sumCh,deltaTr,deltaEl,steerang);
                        varargout{3} = ang_est;
                    end
            end
            
              
        end
        
    end
    
    methods (Access = protected)
        function flag = isInputComplexityLockedImpl(obj,index)  %#ok<INUSL>
            if index == 1
                flag = false;
            else % index == 2
                flag = true;
            end
        end
        
        function icon = getIconImpl(obj) %#ok<MANU>
            % Define icon for System block
            icon = sprintf('Amplitude\nMonopulse\nFeed');
        end

        function varargout = getInputNamesImpl(obj) %#ok<MANU>
            % Return input port names for System block
            varargout = {'X','STEER'};
        end

        function varargout = getOutputNamesImpl(obj)
            % Return output port names for System block
            switch obj.Coverage
                case 'Azimuth'
                    varargout = {'SIGMA','DeltaAz'};
                    if obj.AngleOutputPort
                        varargout{3} = 'ANG';
                    end
                case '3D'
                    varargout = {'SIGMA','DeltaAz','DeltaEl'};
                    if obj.AngleOutputPort
                        varargout{4} = 'ANG';
                    end
            end
        end

        function flag = isInputSizeLockedImpl(obj,ind) %#ok<INUSL>
            if ind == 1
                flag = false;
            else
                flag = true;  % angle
            end
        end
        
        function varargout = getOutputSizeImpl(obj)
            % Return size for each output port
            varargout = cell(1,nargout);
            sz = propagatedInputSize(obj,1);
            switch obj.Coverage
                case 'Azimuth'
                    varargout{1} = [sz(1) 1];
                    varargout{2} = [sz(1) 1];
                    if obj.AngleOutputPort
                        varargout{3} = [1 sz(1)];
                    end
                case '3D'
                    varargout{1} = [sz(1) 1];
                    varargout{2} = [sz(1) 1];
                    varargout{3} = [sz(1) 1];
                    if obj.AngleOutputPort
                        varargout{4} = [2 sz(1)];
                    end
            end
            
        end

        function varargout = isOutputFixedSizeImpl(obj)
            flag = propagatedInputFixedSize(obj, 1);
            switch obj.Coverage
                case 'Azimuth'
                    varargout{1} = flag;
                    varargout{2} = flag;
                    if obj.AngleOutputPort
                        varargout{3} = flag;
                    end
                case '3D'
                    varargout{1} = flag;
                    varargout{2} = flag;
                    varargout{3} = flag;
                    if obj.AngleOutputPort
                        varargout{4} = flag;
                    end
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
            flag = propagatedInputComplexity(obj,1);

            switch obj.Coverage
                case 'Azimuth'
                    varargout{1} = flag;
                    varargout{2} = flag;
                    if obj.AngleOutputPort
                        varargout{3} = false;
                    end
                case '3D'
                    varargout{1} = flag;
                    varargout{2} = flag;
                    varargout{3} = flag;
                    if obj.AngleOutputPort
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
                'AngleOutputPort'};
            groups(1).PropertyList = [groups(1).PropertyList props];
            
            action = matlab.system.display.Action(@phased.MonopulseFeed.onActionCalled, ...
                'Label', 'Generate Monopulse Tracker','Placement','last','Alignment','right');
            groups(1).Actions = action;
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:MonopulseFeedTitle')),...
                'Text',getString(message('phased:library:block:MonopulseFeedDesc')));
        end
        function onActionCalled(actionData, ~)
            SystemHandle = actionData.SystemHandle;
            sysname = get_param(SystemHandle,'Parent');
            b = add_block('phaseddoalib/Monopulse Estimator',...
                [sysname '/Monopulse Estimator'],'MakeNameUnique','on');
            set_param(b,'SensorArray',get_param(SystemHandle,'SensorArray'));
            set_param(b,'PropagationSpeed',get_param(SystemHandle,'PropagationSpeed'));
            set_param(b,'OperatingFrequency',get_param(SystemHandle,'OperatingFrequency'));
            set_param(b,'Coverage',get_param(SystemHandle,'Coverage'));
            set_param(b,'SquintAngle',get_param(SystemHandle,'SquintAngle'));
            feed_name = get_param(SystemHandle,'Name');
            estimator_name = get_param(b,'Name');
            add_line(sysname,sprintf('%s/1',feed_name),sprintf('%s/1',estimator_name),'AutoRouting','on');
            add_line(sysname,sprintf('%s/2',feed_name),sprintf('%s/2',estimator_name),'AutoRouting','on');
            add_line(sysname,sprintf('%s/3',feed_name),sprintf('%s/3',estimator_name),'AutoRouting','on');
        end
    end
    
end
