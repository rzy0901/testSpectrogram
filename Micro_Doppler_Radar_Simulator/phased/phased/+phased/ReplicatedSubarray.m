classdef (Sealed,StrictDefaults) ReplicatedSubarray < phased.internal.AbstractSubarray
%ReplicatedSubarray   Phased array formed by replicated subarrays
%   H = phased.ReplicatedSubarray creates a replicated subarray System
%   object, H. This object models an array formed by replicated subarrays.
%   The default array is a 4-element ULA formed by two 2-element ULAs.
%
%   H = phased.ReplicatedSubarray(Name,Value) creates a replicated subarray
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   RESP = step(H,FREQ,ANGLE,V) returns the subarrays' responses, RESP,
%   given the array's operating frequency FREQ (in Hz) and the directions
%   specified in ANGLE (in degrees). FREQ is a row vector of length L and
%   ANGLE can be either a row vector of length M or a 2xM matrix. V is a
%   scalar representing the propagation speed (in m/s). 
%
%   If the subarray is not capable of simulating polarization, RESP has a
%   dimension of NxMxL where N is the number of subarrays in the phased
%   array. Each column of RESP contains the responses of the array elements
%   for the corresponding direction specified in ANGLE. Each page of RESP
%   contains the responses of the array for the given frequency specified
%   in FREQ. If the subarray is capable of simulating polarization, RESP is
%   a structure containing two fields, H and V. H represents the array's
%   response in horizontal polarization and V represents the array's
%   response in vertical polarization. Each field is an NxMxL array whose
%   columns contains the responses of the array elements for the
%   corresponding direction and frequency.
%
%   When ANGLE is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in the [azimuth; elevation] form. The azimuth
%   angle should be between [-180 180] degrees and the elevation angle
%   should be between [-90 90] degrees. If ANGLE is a length-M row vector,
%   then each element specifies a direction's azimuth angle, and the
%   corresponding elevation angle is assumed to be 0.
%
%   Note that the elements within each subarray are connected to the
%   subarray phase center using an equal-path feed.
%   
%   RESP = step(H,FREQ,ANGLE,V,STEERANGLE) uses STEERANGLE as the
%   subarray's steering direction. STEERANGLE can be either a scalar or a
%   2-element column vector. This syntax is only applicable if you set the
%   SubarraySteering property to either 'Phase' or 'Time'.
%
%   If STEERANGLE is a column vector, it specifies the steering direction
%   in the [azimuth; elevation] form. The azimuth angle should be between
%   [-180 180] degrees and the elevation angle should be between [-90 90]
%   degrees. If STEERANGLE is a scalar, it specifies the direction's
%   azimuth angle, and the corresponding elevation angle is assumed to be
%   0.
%   
%   RESP = step(H,FREQ,ANGLE,V,WS) uses WS as the weights applied to each
%   element in the subarray. WS must be an NSExN matrix where NSE is the
%   number of elements in each individual subarray and N is the number of
%   subarrays. Each column in WS specifies the weights for the elements in
%   the corresponding subarray. This syntax is only applicable when you set
%   the SubarraySteering property to 'Custom'.
%   
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   RESP = step(H,FREQ,ANGLE,V,STEERANGLE)
%
%   or
%
%   RESP = step(H,FREQ,ANGLE,V,WS)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ReplicatedSubarray methods:
%
%   step                  - Output responses of the subarrays (see above)
%   release               - Allow property value and input characteristics
%                           changes
%   clone                 - Create a replicated subarray object with same 
%                           property values
%   collectPlaneWave      - Simulate received plane waves
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the subarray is capable of 
%                           simulating polarization
%   getNumElements        - Number of elements in array
%   getNumSubarrays       - Number of subarrays in array
%   getElementPosition    - Element positions in array
%   getSubarrayPosition   - Subarray positions in array
%   directivity           - Compute array directivity
%   pattern               - Plot subarray response pattern
%   patternAzimuth        - Plot azimuth pattern
%   patternElevation      - Plot elevation pattern
%   viewArray             - View array geometry
%
%   ReplicatedSubarray properties:
%
%   Subarray              - Subarray to replicate
%   Layout                - Subarrays layout
%   GridSize              - Grid size
%   GridSpacing           - Grid spacing
%   SubarrayPosition      - Subarray positions
%   SubarrayNormal        - Subarray normals
%   SubarraySteering      - Subarray steering method
%   PhaseShifterFrequency - Subarray phase shifter frequency
%   NumPhaseShifterBits   - Number of bits in phase shifters
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a 4-element ULA formed by two 2-element ULAs. Plot its 
%   %   azimuth response. Assume the operating frequency is 1 GHz and the 
%   %   wave propagation speed is 3e8 m/s.
%
%   array = phased.ULA(2,0.5);
%   subarray = phased.ReplicatedSubarray('Subarray',array,...
%               'Layout','Rectangular','GridSize',[2 1],...
%               'GridSpacing','Auto');
%   fc = 1e9; c = 3e8;
%   pattern(subarray,fc,-180:180,0,'CoordinateSystem','polar',...
%       'PropagationSpeed',c);
%
%   % Example 2:
%   %   Find the response of each subarray in the above array at the 
%   %   boresight.
%
%   array = phased.ULA(2,0.5);
%   subarray = phased.ReplicatedSubarray('Subarray',array,...
%               'Layout','Rectangular','GridSize',[2 1],...
%               'GridSpacing','Auto');
%   fc = 1e9; c = 3e8; ang = [0;0];
%   resp = subarray(fc,ang,c)
%
%   See also phased.PartitionedArray, phased.ULA, phased.URA,
%   phased.ConformalArray.

%   Copyright 2011-2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)

        %Subarray   Subarray to replicate
        %   Specify the subarray used to form the array. Subarray must be
        %   of type phased.ULA, phased.URA, or phased.ConformalArray. The
        %   default value of this property is a 2-element ULA.
        Subarray
        %Layout     Subarrays layout
        %   Specify the layout of the replicated subarrays as one of
        %   'Rectangular' | 'Custom', where the default is 'Rectangular'.
        Layout = 'Rectangular'
        %GridSize   Grid size
        %   Specify the size of the rectangular grid as a length-2 row
        %   vector in the form of [NumberOfRows NumberOfColumns] where rows
        %   are along y-axis and columns are along z-axis. There is one
        %   subarray in each grid. If GridSize is a scalar, the array has
        %   the same number of subarrays in each row and column. The
        %   default value of this property is [1 2]. This property is only
        %   applicable when you set the Layout property to 'Rectangular'.
        %
        %   Note that the row is along the local y axis and the column is
        %   along the local z axis.
        GridSize = [1 2]
        %GridSpacing   Grid spacing (m)
        %   Specify the rectangular grid spacing (in meters) of the array
        %   as 'Auto', a 1x2 vector or a scalar. The default value of this
        %   property is 'Auto'. This property is only applicable when you
        %   set the Layout property to 'Rectangular'.
        %
        %   If GridSpacing is 'Auto', the replication preserves the element
        %   spacing in both row and column. This option is only applicable
        %   if you use either phased.ULA or phased.URA as the subarray.
        %
        %   If GridSpacing is a 1x2 vector, it is in the form of
        %   [SpacingBetweenRows SpacingBetweenColumn]. If GridSpacing is a
        %   scalar, the spacing along the row and the spacing along the
        %   column are the same.
        GridSpacing = 'Auto'
        %SubarrayPosition   Subarray positions (m)
        %   Specify the positions of the subarrays in the custom grid.
        %   SubarrayPosition must be a 3xN matrix where N indicates the
        %   number of subarrays in the array. Each column of
        %   SubarrayPosition represents the position, in the form of [x; y;
        %   z] (in meters), of a single subarray in the array's local
        %   coordinate system. This property is only applicable if you set
        %   the Layout property to 'Custom'. The default value of this
        %   property is [0 0;-0.5 0.5;0 0].
        SubarrayPosition = [0 0;-0.5 0.5;0 0];
        %SubarrayNormal   Subarray normals (deg)
        %   Specify the normal directions of the subarrays in the array.
        %   SubarrayNormal must be a 2xN matrix where N indicates the
        %   number of subarrays in the array. Each column of SubarrayNormal
        %   specifies the normal direction of the corresponding subarray in
        %   the form of [azimuth; elevation] (in degrees) defined in the
        %   local coordinate system. This property is only applicable if
        %   you set the Layout property to 'Custom'. The default value of
        %   this property is [0 0;0 0].
        SubarrayNormal = [0 0;0 0];
        
    end
    
    properties (Constant,Hidden)
        LayoutSet = matlab.system.StringSet(...
            {'Rectangular','Custom'});       
    end
    
    properties (Access = private, Logical)
        pRectangularLayoutFlag = false
    end
    
    methods
        function set.Subarray(obj,value)
            validateattributes(value,{'phased.internal.AbstractArray'},...
                {'scalar'},'','Subarray');
            obj.Subarray = value;
        end
        
        function set.GridSize(obj,valueArg)
            if isscalar(valueArg)
                value = [valueArg valueArg];
            else
                value = valueArg;
            end
            sigdatatypes.validateIndex(value,'','GridSize',{'size',[1 2]});
            cond = all(value==1);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:invalidSubarrayGridSize','GridSize');
            end
            obj.GridSize = value;
        end
        
        function set.GridSpacing(obj,value)
            validateattributes(value,{'double','char'},{},'','GridSpacing');
            if ischar(value)
                validatestring(value,{'Auto'},'','GridSpacing');
                valueOut = value;
            else
                if isscalar(value)
                    valueOut = [value value];
                else
                    valueOut = value;
                end
                sigdatatypes.validateDistance(valueOut,'','GridSpacing',...
                    {'size',[1 2],'positive'});
            end
            obj.GridSpacing = valueOut;
        end
        
        function set.SubarrayPosition(obj,value)
            sigdatatypes.validate3DCartCoord(value,'','SubarrayPosition');
            cond = size(value,2)<2;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:expectMoreColumns','SubarrayPosition',2);
            end
            obj.SubarrayPosition = value;
        end
        
        function set.SubarrayNormal(obj,value)
            sigdatatypes.validateAzElAngle(value,'','SubarrayNormal');
            cond = size(value,2)<2;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:expectMoreColumns','SubarrayNormal',2);
            end
            obj.SubarrayNormal = value;
        end
    end
    
    methods

        function obj = ReplicatedSubarray(varargin)
            obj@phased.internal.AbstractSubarray(varargin{:});
            if isempty(coder.target)
                if isempty(obj.Subarray)
                    obj.Subarray = phased.ULA;
                end
            else
                if ~coder.internal.is_defined(obj.Subarray)
                    obj.Subarray = phased.ULA;
                end
            end

        end
        
    end
    methods(Hidden)
        function cl = clonecg(obj)
            cl = phased.ReplicatedSubarray(...
                'Subarray',clonecg(obj.Subarray), ...
                'Layout',obj.Layout, ...    
                'SubarraySteering',obj.SubarraySteering);

            if strcmp(obj.Layout,'Rectangular' )
                cl.GridSize = obj.GridSize;                
                cl.GridSpacing = obj.GridSpacing;
            else
                cl.SubarrayPosition = obj.SubarrayPosition;
                cl.SubarrayNormal = obj.SubarrayNormal;                
            end

            if strcmp(obj.SubarraySteering,'Phase')
                cl.PhaseShifterFrequency = obj.PhaseShifterFrequency; 
                cl.NumPhaseShifterBits = obj.NumPhaseShifterBits;
            end
        end
    end
    methods (Access = protected)
        function num = calcNumSubarrays(obj)
            if strcmp(obj.Layout,'Rectangular')
                num = prod(obj.GridSize);
            else
                num = size(obj.SubarrayPosition,2);
            end
        end
        
        function num = calcNumElements(obj,subarrayidx)
            NumSubarray = getNumSubarrays(obj);
            NumPerSubarray = getNumElements(obj.Subarray);
            if nargin > 1
                sigdatatypes.validateIndex(subarrayidx,'getNumElements',...
                    'SubarrayIndex',{'row','<=',NumSubarray});
                num = NumPerSubarray*ones(1,numel(subarrayidx));
            else
                num = NumPerSubarray*NumSubarray;
            end
        end
        
        function pos = calcSubarrayPosition(obj)
            if strcmp(obj.Layout,'Rectangular')
                if ischar(obj.GridSpacing)
                    if isa(obj.Subarray,'phased.ULA') || ...
                            isa(obj.Subarray,'phased.HeterogeneousULA')
                        spacing = [1 getNumElements(obj.Subarray)]*...
                            obj.Subarray.ElementSpacing;
                    elseif isa(obj.Subarray,'phased.URA')
                        spacing = obj.Subarray.Size.*obj.Subarray.ElementSpacing;
                    else % HeterogeneousURA
                        spacing = size(obj.Subarray.ElementIndices).*obj.Subarray.ElementSpacing;
                    end
                else
                    spacing = obj.GridSpacing;
                end
                coder.internal.const(spacing);
                N = coder.internal.const(getNumSubarrays(obj));       
                pos = coder.internal.const(calcSubarrayPositionCG(N,spacing,obj.GridSize));

            else
                pos = obj.SubarrayPosition;
            end
        end
         
        function pos = calcElementPosition(obj)
            pos_element = getElementPosition(obj.Subarray);
            pos_subarray = getSubarrayPosition(obj);
            NPerSubarray = getNumElements(obj.Subarray);
            NSubarray = getNumSubarrays(obj);
            pos = kron(pos_subarray,ones(1,NPerSubarray));
            if strncmpi(obj.Layout,'Rectangular',1)
                pos = pos + kron(ones(1,NSubarray),pos_element);
            else
                normal_subarray = obj.SubarrayNormal;
                for m = 1:NSubarray
                    idx = (m-1)*NPerSubarray+1:m*NPerSubarray;
                    pos(:,idx) = pos(:,idx) + ...
                        phased.internal.rotazel(pos_element,normal_subarray(:,m));
                end
           end
        end
    end
    
    methods (Access = protected)
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSubarray(obj);
            s.Subarray = saveobj(obj.Subarray);
            if isLocked(obj)
                s.pRectangularLayoutFlag = obj.pRectangularLayoutFlag;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSubarray(obj,s);
            obj.Subarray = eval(...
                sprintf('%s.loadobj(s.Subarray)',s.Subarray.ClassNameForLoadTimeEval));
            s = rmfield(s,'Subarray');
            if isfield(s,'isLocked')
                s = rmfield(s,'isLocked');
            end
            if isfield(s,'pPropagationSpeed')  % obsolete property
                s = rmfield(s,'pPropagationSpeed');
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function validatePropertiesImpl(obj)
            if strncmpi(obj.Layout,'Custom',1)
                cond = size(obj.SubarrayPosition,2) ~= size(obj.SubarrayNormal,2);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:NumColumnsMismatch',...
                        'SubarrayPosition','SubarrayNormal');
                end
            end
            cond = strncmpi(obj.Layout,'Rectangular',1) && ...
                    strncmpi(obj.GridSpacing,'auto',1) && ...
                    ~(isa(obj.Subarray,'phased.ULA') || ...
                    isa(obj.Subarray,'phased.URA') || ...
                    isa(obj.Subarray,'phased.HeterogeneousULA') || ...
                    isa(obj.Subarray,'phased.HeterogeneousURA'));
            if cond
                coder.internal.errorIf(cond,'phased:system:array:invalidSubarrayForAutoGrid',...
                    'Subarray','phased.ULA','phased.URA','GridSpacing','Auto');
            end
        end

        function setupImpl(obj,freq,ang,varargin)

            if isempty(coder.target)
                obj.cArray = clone(obj.Subarray);
                release(obj.cArray);
            else
                obj.cArray = clonecg(obj.Subarray);
            end

            setupImpl@phased.internal.AbstractSubarray(obj,freq,ang,varargin{:});
            if isPolarizationEnabled(obj)
                obj.cModuleResponse = phased.ArrayResponse(...
                    'SensorArray',obj.cArray,...
                    'WeightsInputPort',true,...
                    'EnablePolarization',true);
                obj.cModuleResponse.PropagationSpeedSource = 'Input port';
            else
                obj.cModuleResponse = phased.ArrayResponse(...
                    'SensorArray',obj.cArray,...
                    'WeightsInputPort',true);
                obj.cModuleResponse.PropagationSpeedSource = 'Input port';
            end
            obj.cModuleSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.cArray,...
                'IncludeElementResponse',false);
            if obj.pPhaseSteeringFlag
                obj.cModuleSteeringVector.NumPhaseShifterBits = ...
                    obj.NumPhaseShifterBits;
            end
            obj.cModuleSteeringVector.PropagationSpeedSource = 'Input port';
            if strncmp(obj.Layout,'Rectangular',1)
                obj.pRectangularLayoutFlag = true;
            else
                obj.pRectangularLayoutFlag = false;
            end
            
        end

        function resp = stepImpl(obj, freq, angArg, c, stang)

            sigdatatypes.validateSpeed(c,'','V',{'positive'});
            
            if size(angArg,1) == 1
                ang = [angArg;zeros(size(angArg))];
            else
                ang = angArg;
            end
            
            if obj.pTimeSteeringFlag || obj.pPhaseSteeringFlag
                [f,theta] = getSteeringConfiguration(obj,freq,stang);
            elseif ~obj.pNoSteeringFlag  % element weights
                cond = ~all(isfinite(stang(:)));
                if cond        
                    coder.internal.errorIf(cond,'phased:step:expectedFinite','W');
                end
            end
            
            [Na,N,M,L] = calcRespDimensions(obj);
            
            if obj.pRectangularLayoutFlag
                if obj.pNoSteeringFlag || (~obj.pNoSteeringFlag && ...
                        (obj.pPhaseSteeringFlag || obj.pTimeSteeringFlag))
                    if obj.pNoSteeringFlag
                        w = ones(Na,1);
                    else % phase or time steering
                        w = squeeze(step(obj.cModuleSteeringVector,f,theta,c));
                    end
                    if isPolarizationEnabled(obj)
                        subarrayresp = step(obj.cModuleResponse,freq,ang,w,c);
                        subarrayresp_h = subarrayresp.H;
                        subarrayresp_v = subarrayresp.V;
                        resp_h = reshape(bsxfun(@times,ones(N,M*L),...
                            subarrayresp_h(:).'),N,M,[]);
                        resp_v = reshape(bsxfun(@times,ones(N,M*L),...
                            subarrayresp_v(:).'),N,M,[]);
                        resp.H = resp_h;
                        resp.V = resp_v;
                    else
                        subarrayresp = step(obj.cModuleResponse,freq,ang,w,c);
                        resp = reshape(bsxfun(@times,ones(N,M*L),...
                            subarrayresp(:).'),N,M,[]);
                    end
                else  % custom weights
                    w = conj(stang);    % this is like a taper, needs conj
                    if isPolarizationEnabled(obj)
                        resp_h = complex(zeros(N,M,L));
                        resp_v = complex(zeros(N,M,L));
                        for m = N:-1:1
                            wele = getElementWeightsForSubarray(obj,w,m);
                            subarrayresp = step(obj.cModuleResponse,freq,ang,wele,c);
                            resp_h(m,:,:) = subarrayresp.H;
                            resp_v(m,:,:) = subarrayresp.V;
                        end
                        resp.H = resp_h;
                        resp.V = resp_v;
                    else
                        resp = complex(zeros(N,M,L));
                        for m = N:-1:1
                            wele = getElementWeightsForSubarray(obj,w,m);
                            resp(m,:,:) = step(obj.cModuleResponse,freq,ang,wele,c);
                        end
                    end
                end
            else
                % calculate the corresponding az/el for the incoming
                % direction at each subarray

                % Note that in conformal case, no steering is not the
                % same as pointing to [0;0]

                subarray_incang = zeros(2,M,N);
                subarray_steerang = zeros(2,N);
                subarray_normal = obj.SubarrayNormal;
                for n = N:-1:1
                    tempaxes = phased.internal.rotazel(eye(3),...
                        subarray_normal(:,n));
                    subarray_incang(:,:,n) = ...
                        phased.internal.incident2azel(ang,tempaxes);
                    if ~obj.pNoSteeringFlag && ...
                            (obj.pPhaseSteeringFlag || obj.pTimeSteeringFlag)
                        subarray_steerang(:,n) = ...
                            phased.internal.incident2azel(theta,tempaxes);
                    end
                end

                % calculate array response toward that az/el given the
                % weights steered to the steering direction
                if isPolarizationEnabled(obj)
                    subarrayresp_h = complex(zeros(N,M,L));
                    subarrayresp_v = complex(zeros(N,M,L));
                    [subarray_phi_hat, subarray_theta_hat] = ...
                        phased.internal.azel2vec(ang);
                    for n = N:-1:1
                        if obj.pNoSteeringFlag
                            w = ones(Na,1);
                        elseif obj.pPhaseSteeringFlag || obj.pTimeSteeringFlag
                            w = squeeze(step(obj.cModuleSteeringVector,...
                                f,subarray_steerang(:,n),c));
                        else % custom weights
                            % this is like a taper, needs conj
                            w = conj(getElementWeightsForSubarray(obj,stang,n));
                        end
                        tempaxes = phased.internal.rotazel(eye(3),...
                            subarray_normal(:,n));
                        subarrayresp_l = ...
                            step(obj.cModuleResponse,freq,...
                            subarray_incang(:,:,n),w,c);
                        [phi_hat_l, theta_hat_l] = ...
                            phased.internal.azel2vec(...
                            subarray_incang(:,:,n));
                        phi_hat_g = phased.internal.local2globalvec(...
                            phi_hat_l,tempaxes);
                        theta_hat_g = phased.internal.local2globalvec(...
                            theta_hat_l,tempaxes);
                        for m = 1:M
                            map_matrix = [...
                                phi_hat_g(:,m)'*subarray_phi_hat(:,m) ...
                                theta_hat_g(:,m)'*subarray_phi_hat(:,m); ...
                                phi_hat_g(:,m)'*subarray_theta_hat(:,m) ...
                                theta_hat_g(:,m)'*subarray_theta_hat(:,m)];
                            temp = map_matrix*[subarrayresp_l.H(m,:);...
                                subarrayresp_l.V(m,:)];
                            subarrayresp_h(n,m,:) = temp(1,:);
                            subarrayresp_v(n,m,:) = temp(2,:);
                        end
                    end

                    resp.H = subarrayresp_h;
                    resp.V = subarrayresp_v;
                else
                    subarrayresp = complex(zeros(N,M,L));
                    for n = N:-1:1
                        if obj.pNoSteeringFlag
                            w = ones(Na,1);
                        elseif obj.pPhaseSteeringFlag || obj.pTimeSteeringFlag
                            w = squeeze(step(obj.cModuleSteeringVector,...
                                f,subarray_steerang(:,n),c));
                        else % custom weights
                            % this is like a taper, needs conj
                            w = conj(getElementWeightsForSubarray(obj,stang,n));
                        end
                        subarrayresp(n,:,:) = step(obj.cModuleResponse,freq,...
                            subarray_incang(:,:,n),w,c);
                    end

                    resp = subarrayresp;
                end

            end
        end
        
        function validateElementWeights(obj,ws)
            sz_ws = size(ws);
            cond =  ~ismatrix(ws) || isempty(ws);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','WS');
            end
            Nse = getNumElements(obj.Subarray);
            Ns = getNumSubarrays(obj);
            cond = ~isequal(sz_ws,[Nse Ns]);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:phased:expectedMatrixSize','WS',Nse,Ns);
            end
            cond = ~isa(ws,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','WS','double');
            end
            
        end

        function wele = getElementWeightsForSubarray(obj,w,subidx) %#ok<INUSL>

            % subidx is the index of the subarray
            wele = w(:,subidx);
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractSubarray(obj, prop);
            if strncmp(obj.Layout,'Custom',1) && ...
                    (strcmp(prop,'GridSpacing') || ...
                    strcmp(prop,'GridSize'))
                flag = true;
            end
            if strncmp(obj.Layout,'Rectangular',1) && ...
                    (strcmp(prop,'SubarrayPosition') || ...
                    strcmp(prop,'SubarrayNormal'))
                flag = true;
            end
        end
        
    end
    
    methods (Hidden, Access = {?phased.gpu.internal.AbstractClutterSimulator, ?phased.internal.AbstractSubarray})
        function [pos, az, el] = getSubarrayPosAzEl(obj, subArrIdx, azin, elin)
            %Transform to new pos, az, and el for subarray at index
            %subArrIdx based on the SubarrayNormal. 
            pos =getElementPosition(obj.Subarray); 
            subarray_normal = obj.SubarrayNormal;
            if ~strcmpi(obj.Layout, 'Rectangular')
                azAx = phased.internal.deg2rad(subarray_normal(1, subArrIdx));
                elAx = phased.internal.deg2rad(subarray_normal(2, subArrIdx));
                [az, el] = phased.gpu.internal.angleAxesConversion(azin, elin, ...
                                                           azAx, elAx); 
            else
                az = azin;
                el = elin;
            end
          
            
        end
        function sv = getSubArraySteeringVec(obj, subArrIdx, freq, c, stang )
            if strncmp(obj.SubarraySteering,'Phase',1) || ...
                    strncmp(obj.SubarraySteering,'Time',1)
                %Get the steering vector for this subarray at steering angle
                %stang.
                if strncmp(obj.SubarraySteering,'Phase',1)
                    f = obj.PhaseShifterFrequency;
                    theta = stang;
                else %Time
                    f = freq;
                    theta = stang;
                end
                
                %Do this on the CPU, it's small.
                if ~strcmpi(obj.Layout, 'Rectangular')
                    subarray_normal = obj.SubarrayNormal(:, subArrIdx);
                    tempaxes = phased.internal.rotazel(eye(3),...
                        subarray_normal);
                    subarr_steerang = phased.internal.incident2azel(theta,tempaxes);
                else
                    subarr_steerang = theta;
                end
                posSub =getElementPosition(obj.Subarray);
                
                sv = phased.internal.steeringvec(posSub, f,c, subarr_steerang);
                sv = conj(sv);
            else % custom weights
                sv = getElementWeightsForSubarray(obj,stang,subArrIdx);
            end
        
        end
        
    end
    
    methods (Access = {?phased.ArrayGain, ?phased.gpu.internal.AbstractClutterSimulator})
        function Rn = getNoiseGainMatrix(obj,freq,c,stang)
            if strcmp(obj.SubarraySteering,'None')
                w = ones(getNumElements(obj.Subarray),1);
                Rn = real(w'*w)*eye(getDOF(obj)); 
            elseif strcmp(obj.SubarraySteering,'Phase')
                hstv = phased.SteeringVector(...
                    'SensorArray',obj.Subarray,...
                    'PropagationSpeed',c,...
                    'NumPhaseShifterBits',obj.NumPhaseShifterBits);
                w = step(hstv,obj.PhaseShifterFrequency,stang);
                Rn = real(w'*w)*eye(getDOF(obj)); 
            elseif strcmp(obj.SubarraySteering,'Time')
                hstv = phased.SteeringVector(...
                    'SensorArray',obj.Subarray,...
                    'PropagationSpeed',c);
                w = step(hstv,freq,stang);
                Rn = real(w'*w)*eye(getDOF(obj)); 
            else % custom weights
                Rn = eye(getDOF(obj));
                for m = 1:getDOF(obj)
                    w = stang(:,m);
                    Rn(m,m) = real(w'*w);
                end
            end
        end
    end
    
    
    methods 
        function flag = isPolarizationCapable(obj)
        %isPolarizationCapable Indicate if the subarray is capable of 
        %simulating polarization
        %   F = isPolarizationCapable(H) returns the flag, F, which
        %   indicates whether the subarray H is capable of simulating
        %   polarization.
        %
        %   % Example: 
        %   %   Determine whether a subarray is capable of simulating 
        %   %   polarization.
        %   
        %   h = phased.ReplicatedSubarray;
        %   f = isPolarizationCapable(h)
        %
        %   See also phased.
            flag = isPolarizationCapable(obj.Subarray);
        end
        
    end
    
    methods (Hidden)
        function flag = isPolarizationEnabled(obj)
        %isPolarizationEnabled Indicate if the subarray is enabled to 
        %simulate polarization
        %   F = isPolarizationEnabled(H) returns the flag, F, which
        %   indicates whether the subarray H is enabled to simulate
        %   polarization.
        %
        %   % Example: 
        %   %   Determine whether a subarray is enabled to simulate 
        %   %   polarization.
        %   
        %   h = phased.ReplicatedSubarray;
        %   f = isPolarizationEnabled(h)
        %
        %   See also phased.
            flag = isPolarizationEnabled(obj.Subarray);
        end
        
    end
    
    methods (Hidden)
        function h = getElementHandle(obj)
            if isa(obj.Subarray,'phased.internal.AbstractHeterogeneousArray')
                h = obj.Subarray.ElementSet;
            else
                h = obj.Subarray.Element;
            end
        end
        function newObj = cloneSensor(obj)
            if isElementFromAntenna(obj)
                newObj = clone(obj); 
                release(newObj); 
                newObj.Subarray = cloneSensor(obj.Subarray);
            else
                newObj = clone(obj);
                release(newObj);
            end
        end
        function flag = isElementFromAntenna(obj) 
            flag = isElementFromAntenna(obj.Subarray);
        end
    end
    methods(Static, Hidden, Access=protected)
        function groups = getPropertyGroupsImpl
            sensorClassSet = matlab.system.display.internal.ClassStringSet(...
                {
                    'phased.ULA',...
                    'phased.URA',...
                    'phased.UCA',...
                    'phased.ConformalArray'...
                }, ...
                'PropertiesTitle', '', ...
                'Labels', { 
                    'ULA',...
                    'URA',...
                    'UCA',...                    
                    'Conformal Array'...
                          });
            sensorProp = matlab.system.display.internal.Property('Subarray', ...
                'IsObjectDisplayOnly', true, ...
                'Description', 'Sensor array', ...
                'ClassStringSet', sensorClassSet);  
            props = {sensorProp,...
                'Layout',...
                'GridSize',...
                phased.internal.GridSpacingProperty('GridSpacing'),...
                'SubarrayPosition',...
                'SubarrayNormal',...
                'SubarraySteering',...
                'PhaseShifterFrequency'...
                'NumPhaseShifterBits'};
            groups = matlab.system.display.Section('Title', 'Parameters', ...
                'PropertyList', props);
        end
    end
end

function pos = calcSubarrayPositionCG(N,spacing,gridSize)
                NPerRow = gridSize(2);     % number of elements in each row
                NPerCol = gridSize(1);     % number of elements in each column
                pos = zeros(3,N);
                RowEnd = (NPerRow-1)/2;
                ColEnd = (NPerCol-1)/2;
                [row_temp, col_temp] = meshgrid(...
                    (-RowEnd:RowEnd)*spacing(2),...
                    (ColEnd:-1:-ColEnd)*spacing(1));
                pos(2,:) = row_temp(:).';
                pos(3,:) = col_temp(:).';
end


% [EOF]
