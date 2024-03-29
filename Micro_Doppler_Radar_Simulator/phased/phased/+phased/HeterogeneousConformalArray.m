classdef (Sealed,StrictDefaults) HeterogeneousConformalArray < ...
phased.internal.AbstractHeterogeneousArray
%HeterogeneousConformalArray   Heterogeneous conformal array
%   H = phased.HeterogeneousConformalArray creates a heterogeneous
%   conformal array System object, H. This object models a heterogeneous
%   conformal array formed with sensor elements whose pattern may vary. The
%   default array is a single isotropic antenna element.
%
%   H = phased.HeterogeneousConformalArray(Name,Value) creates a
%   heterogeneous conformal array object, H, with the specified property
%   Name set to the specified Value. You can specify additional name-value
%   pair arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   HeterogeneousConformalArray can be used to model an array with
%   arbitrary geometry whose elements point in different directions.
%
%   Step method syntax:
%
%   RESP = step(H,FREQ,ANGLE) returns the array elements' responses, RESP,
%   given the array's operating frequency FREQ (in Hz) and the directions
%   specified in ANGLE (in degrees). FREQ is a row vector of length L and
%   ANGLE can be either a row vector of length M or a 2xM matrix. 
%
%   If the array is not capable of simulating polarization, RESP has a
%   dimension of NxMxL where N is the number of elements in the phased
%   array. Each column of RESP contains the responses of the array elements
%   for the corresponding direction specified in ANGLE. Each page of RESP
%   contains the responses of the array for the given frequency specified
%   in FREQ. If the array is capable of simulating polarization, RESP is a
%   structure containing two fields, H and V. H represents the array's
%   response in horizontal polarization and V represents the array's
%   response in vertical polarization. Each field is an NxMxL array whose
%   columns contains the responses of the array elements for the
%   corresponding direction and frequency.
%
%   When ANGLE is a 2xM matrix, each column of the matrix specifies the
%   direction in space in [azimuth; elevation] form. The azimuth angle
%   should be between [-180 180] degrees and the elevation angle should be
%   between [-90 90] degrees. If ANGLE is a length M row vector, then each
%   element specifies a direction's azimuth angle, and the corresponding
%   elevation angle is assumed to be 0.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   HeterogeneousConformalArray methods:
%
%   step                  - Output responses of the array elements
%   release               - Allow property name and input characteristics
%                           changes
%   clone                 - Create a conformal array object with same 
%                           property values
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the array is capable of 
%                           simulating polarization
%   getNumElements        - Get number of elements in the array
%   getElementPosition    - Get positions of array elements
%   getElementNormal      - Get normal directions of array elements
%   getTaper              - Get array element tapers
%   directivity           - Compute array directivity
%   pattern               - Plot array response pattern
%   patternAzimuth        - Plot azimuth pattern
%   patternElevation      - Plot elevation pattern
%   collectPlaneWave      - Simulate received plane waves
%   viewArray             - View array geometry
%
%   HeterogeneousConformalArray properties:
%
%   ElementSet            - Set of elements used in the array
%   ElementIndices        - Assign elements to their locations
%   ElementPosition       - Element positions 
%   ElementNormal         - Element normal directions 
%   Taper                 - Taper on elements
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a 8-element heterogeneous uniform circular array (UCA). 
%   %   Four of the elements have a cosine pattern with a power of 1.6 
%   %   while the rest of elements have a cosine pattern with a power of 
%   %   1.5. Plot its azimuth response. Assume the operating frequency is 
%   %   1 GHz and the wave propagation speed is 3e8 m/s.
%
%   sElement1 = phased.CosineAntennaElement('CosinePower',1.6);
%   sElement2 = phased.CosineAntennaElement('CosinePower',1.5);
%   sArray = phased.HeterogeneousConformalArray(...
%       'ElementSet',{sElement1,sElement2},...
%       'ElementIndices',[1 1 1 1 2 2 2 2]);
%   N = 8; azang = (0:N-1)*360/N-180;
%   sArray.ElementPosition = [cosd(azang);sind(azang);zeros(1,N)];
%   sArray.ElementNormal = [azang;zeros(1,N)];
%   fc = 1e9; c = 3e8;
%   pattern(sArray,fc,-180:180,0,'CoordinateSystem','polar',...
%       'PropagationSpeed',c);
%
%   % Example 2:
%   %   Find the response of each element in the above array at the 
%   %   boresight.
%
%   sElement1 = phased.CosineAntennaElement('CosinePower',1.6);
%   sElement2 = phased.CosineAntennaElement('CosinePower',1.5);
%   sArray = phased.HeterogeneousConformalArray(...
%       'ElementSet',{sElement1,sElement2},...
%       'ElementIndices',[1 1 1 1 2 2 2 2]);
%   N = 8; azang = (0:N-1)*360/N-180;
%   sArray.ElementPosition = [cosd(azang);sind(azang);zeros(1,N)];
%   sArray.ElementNormal = [azang;zeros(1,N)];
%   fc = 1e9; ang = [0;0];
%   resp = sArray(fc,ang)
%
%   See also phased, phased.ULA, phased.URA, phased.ConformalArray.

%   Copyright 2012-2016 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] Josefsson and Persson, Conformal Array Antenna Theory and Design,
%       Wiley, 2006
%   [3] Harold Mott, Antennas for Radar and Communications, A Polarimetric
%       Approach, John Wiley & Sons, 1992

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable) 
        %ElementIndices  Assign elements to their locations
        %   This property specifies the mapping of elements in the array.
        %   The property assigns elements to their locations in the array
        %   using the indices into the ElementSet property. ElementIndices
        %   must be a length-N row vector. N is the number of elements in
        %   the sensor array. The values in ElementIndices should be less
        %   than or equal to the number of entries in the ElementSet
        %   property. The default value of this property is 1.
        ElementIndices = 1
        %ElementPosition   Element positions
        %   ElementPosition specifies the positions of the elements in the
        %   conformal array. ElementPosition must be a 3xN matrix where N
        %   indicates the number of elements in the conformal array. Each
        %   column of ElementPosition represents the position, in the form
        %   of [x; y; z] (in meters), of a single element in the array's
        %   local coordinate system. The default value of this property is
        %   [0; 0; 0], i.e., a single element.
        %
        %   For details about the local coordinate system of a
        %   heterogeneous conformal array, type
        %   phased.ConformalArray.coordinateSystemInfo.
        ElementPosition = [0; 0; 0]
        %ElementNormal   Element normal directions
        %   ElementNormal specifies the normal directions of the elements
        %   in the conformal array. ElementNormal must be either a 2xN
        %   matrix where N indicates the number of elements in the array,
        %   or a 2x1 column vector. If ElementNormal is a matrix, each
        %   column of ElementNormal specifies the normal direction of the
        %   corresponding element in the form of [azimuth; elevation] (in
        %   degrees) defined in the local coordinate system. If
        %   ElementNormal is a vector, it specifies the pointing direction
        %   of all elements in the array in the form of [azimuth;
        %   elevation] (in degrees). The default value of this property is
        %   [0; 0].
        %
        %   For details about the local coordinate system of a
        %   heterogeneous conformal array, type
        %   phased.ConformalArray.coordinateSystemInfo.
        ElementNormal = [0; 0]
        %Taper  Taper on elements
        %    Specify the Taper property as a scalar or a length-N vector of
        %    complex weights applied to each element in the sensor array. N
        %    is the number of elements in the array. If Taper is a scalar
        %    the same weights will be applied to each element. If Taper is
        %    a vector, each weight will be applied to the corresponding
        %    element. The default value of Taper is 1.
        Taper = 1;
    end
    
    properties (Access = private, Nontunable)
        % private property to hold number of elements. Gets set at the
        % compiler time so that we don't have to constantly validate the
        % ElementPosition vs. ElementNormal.
        pNumElements
        % private property to hold element axes for each element. Once the
        % normals are set, the element axes are fixed. So we don't have to
        % compute it at each simulation step.
        pElementAxes
    end
    
    methods
        function set.ElementPosition(obj,val)
            sigdatatypes.validate3DCartCoord(val,...
                '','ElementPosition');
            obj.ElementPosition = val;
        end
        
        function set.ElementNormal(obj,val)
            sigdatatypes.validateAzElAngle(val,...
                '','ElementNormal');
            obj.ElementNormal = val;
        end
        
        function set.ElementIndices(obj,val)
            sigdatatypes.validateIndex(val,'','ElementIndices',...
                {'row'});
            obj.ElementIndices = val;
        end
        
        function set.Taper(obj,val)
            validateattributes(val,{'double'},...
                               {'nonnan','nonempty','finite','vector'},...
                               '','Taper');
            obj.Taper = val;
        end
        
    end
      
    methods

        function obj = HeterogeneousConformalArray(varargin)

            obj@phased.internal.AbstractHeterogeneousArray(nargin,...
                varargin{:});
            

        end

        function N = getNumElements(obj)
        %getNumElements   Number of elements in the array
        %   N = getNumElements(H) returns the number of elements, N, in
        %   the heterogeneous conformal array object H. 
        %
        %   % Example:
        %   %   Construct a default heterogeneous conformal array and 
        %   %   obtain its number of elements.
        %
        %   ha = phased.HeterogeneousConformalArray;
        %   N = getNumElements(ha)
            
            if isempty(coder.target)
                if ~isLocked(obj)
                    setNumElements(obj);
                end
            else
                if ~coder.internal.is_defined(obj.pNumElements)
                    setNumElements(obj);
                end
            end
            N = obj.pNumElements;
        end
               
        function pos = getElementPosition(obj,EleIdx)
        %getElementPosition Element positions of the array
        %   POS = getElementPosition(H) returns the element positions of
        %   the heterogeneous conformal array H. POS is a 3xN matrix where
        %   N is the number of elements in H. Each column of POS defines
        %   the position, in the form of [x; y; z] (in meters), of an
        %   element in the local coordinate system.
        %
        %   For details regarding the local coordinate system of a
        %   heterogeneous conformal array, type
        %   phased.ConformalArray.coordinateSystemInfo.
        %
        %   POS = getElementPosition(H,ELEIDX) returns the positions of the
        %   elements that are specified in the element index vector ELEIDX.
        %
        %   % Example:
        %   %   Construct a default heterogeneous conformal array and 
        %   %   obtain its element positions.
        %
        %   ha = phased.HeterogeneousConformalArray;
        %   pos = getElementPosition(ha)
            
            narginchk(1,2);
            validateattributes(obj,{'phased.HeterogeneousConformalArray'},...
                {'scalar'},'getElementPosition','H');
            N = getNumElements(obj);
            if nargin < 2
                EleIdx = 1:N;
            end
            sigdatatypes.validateIndex(EleIdx,'getElementPosition',...
                'ELEIDX',{'vector','<=',N});
            pos = obj.ElementPosition(:,EleIdx);
        end
        
        function nv = getElementNormal(obj,EleIdx)
        %getElementNormal Get normal directions of array elements
        %   NV = getElementNormal(H) returns the element normal directions
        %   of the heterogeneous conformal array H. NV is a 2xN matrix
        %   where N is the number of elements of H. Each column of NV
        %   defines the normal direction, in the form of [azimuth;
        %   elevation] (in degrees), of an element in the local coordinate
        %   system.
        %
        %   For details regarding the local coordinate system of a
        %   heterogeneous conformal array, type
        %   phased.ConformalArray.coordinateSystemInfo.
        %
        %   NV = getElementNormal(H,ELEIDX) returns the normal directions
        %   of the elements that are specified in the element index vector
        %   ELEIDX.
        %
        %   % Example:
        %   %   Construct a default heterogeneous conformal array and 
        %   %   obtain its element normal directions.
        %
        %   ha = phased.HeterogeneousConformalArray;
        %   nv = getElementNormal(ha)

            N = obj.pNumElements;
            if nargin < 2
                EleIdx = 1:N;
            end
            if iscolumn(obj.ElementNormal)
                nv = repmat(obj.ElementNormal,1,numel(EleIdx));
            else
                nv = obj.ElementNormal(:,EleIdx);
            end
        end
        
        function w = getTaper(obj)
        %getTaper Get array element tapers
        %   W = getTaper(H) returns the element taper of the heterogeneous
        %   conformal array H. W is a length-N column vector where N is the
        %   number of elements in H.
        %
        %   % Example: 
        %   %   Construct a default heterogeneous conformal array and 
        %   %   obtain its element taper.
        %
        %   ha = phased.HeterogeneousConformalArray;
        %   pos = getTaper(ha)
            
            wArg = obj.Taper(:);
            if isscalar(wArg)
                w = wArg*ones(getNumElements(obj),1);
            else
                w = wArg;
            end
        end
        
    end
    
    methods(Access = private)
        function setNumElements(obj)
            N = size(obj.ElementPosition,2);
            cond = ~iscolumn(obj.ElementNormal) && (N ~= size(obj.ElementNormal,2));
            if cond
                coder.internal.errorIf(cond,'phased:system:DimensionMismatch',...
                    'ElementPosition','ElementNormal');
            end
            cond =  (N ~= numel(obj.ElementIndices));
            if cond
                coder.internal.errorIf(cond,'phased:system:DimensionMismatch',...
                    'ElementPosition','ElementIndices');
            end
            obj.pNumElements = N;
        end
    end
            
    methods(Access = protected)
        
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractHeterogeneousArray(obj);
            sz = getNumElements(obj);
            
            cond = ~isscalar(obj.Taper) ...
                && any(size(obj.Taper) ~= [1 sz]) ...
                && any(size(obj.Taper) ~= [sz 1]);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidTaper','Taper',...
                    sprintf('%d',sz));
            end
        end
        
        function setupImpl(obj, varargin)
            setupImpl@phased.internal.AbstractHeterogeneousArray(obj, varargin{:});
            setNumElements(obj);
            nv = getElementNormal(obj);
            eleaxes = computeElementAxes(nv,obj.pNumElements);
            obj.pElementAxes = eleaxes;
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractHeterogeneousArray(obj);
            if isLocked(obj)
                s.pNumElements = obj.pNumElements;
                s.pElementAxes = obj.pElementAxes;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function az_el = convertIncidentToAzEl(obj,inc_ang)
        %convertIncidentToAzEl Converts the incident angle to corresponding
        %normal directions for an element in the array.
        %   AZ_EL = convertIncidentToAzEl(H,INC_ANG) converts the incident
        %   directions specified in INC_ANG (in degrees), in the form of
        %   [azimuth; elevation] to the local azimuth and elevation angles
        %   (in degrees) of sensor elements in the array H. AZ_EL is a
        %   2x(N*M) matrix where M is the number of rows in INC_ANG and N
        %   is the number of elements in the array. Each column of the
        %   AZ_EL represents the translated azimuth and elevation angle, in
        %   the form of [azimuth; elevation], of the corresponding incident
        %   angle in INC_ANG. The first N columns of AZ_EL is the [az; el]
        %   for the elements corresponding to the first angle specified in
        %   INC_ANG. The second N columns of AZ_EL is the [az; el] for the
        %   elements corresponding to the second angle specified in INC_ANG,
        %   and so on.
        
            % obtain the normal direction
            M = size(inc_ang,2);
            N = getNumElements(obj);
            az_el = zeros(2,N*M);
            for m = 1:N
                az_el(:,(0:M-1)*N+m) = incidentAngleToAzEl(obj,inc_ang,m);
            end
        end

        function arPl = getArrayPlaneOrAxis(obj) %#ok<MANU>
           arPl = 'Not Specified'; 
        end
        
        function resp = getPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles)
            N = getNumElements(obj);
            [array_phi_hat, array_theta_hat] = phased.internal.azel2vec(ang);
            for m = N:-1:1
                ElementAxes = squeeze(obj.pElementAxes(:,:,m));
                inc_angle = phased.internal.incident2azel(ang,ElementAxes);

                lElementTypeMap = obj.ElementIndices;
                if isscalar(lElementTypeMap)
                    resp_l = ...
                        step(obj.cElement{lElementTypeMap},freq,inc_angle);
                else
                    resp_l = ...
                        step(obj.cElement{lElementTypeMap(m)},freq,inc_angle);
                end
                    
                [phi_hat_l, theta_hat_l] = phased.internal.azel2vec(inc_angle);
                phi_hat_g = phased.internal.local2globalvec(...
                    phi_hat_l,ElementAxes);
                theta_hat_g = phased.internal.local2globalvec(...
                    theta_hat_l,ElementAxes);
                for n = num_angles:-1:1
                    map_matrix = [phi_hat_g(:,n)'*array_phi_hat(:,n) ...
                        theta_hat_g(:,n)'*array_phi_hat(:,n); ...
                        phi_hat_g(:,n)'*array_theta_hat(:,n) ...
                        theta_hat_g(:,n)'*array_theta_hat(:,n)];
                    temp = map_matrix*[resp_l.H(n,:);resp_l.V(n,:)];
                    respidx = (n-1)*N+m;
                    resp.H(respidx,:) = temp(1,:);
                    resp.V(respidx,:) = temp(2,:);
                end
            end
            resp.H = resp.H.*obj.pTaper;
            resp.H = reshape(resp.H,N,num_angles,[]);
            resp.V = resp.V.*obj.pTaper;
            resp.V = reshape(resp.V,N,num_angles,[]);
        end
    end
    
    methods (Hidden)
        function cl = clonecg(obj)
            Ne = numel(obj.ElementSet);
            newElementSet = cell(1,Ne);
            for m = 1:Ne
                newElementSet{m} = clonecg(obj.ElementSet{m});
            end
            cl = phased.HeterogeneousConformalArray(...
                    'ElementSet',newElementSet, ...
                    'ElementIndices',obj.ElementIndices, ...                
                    'ElementPosition',obj.ElementPosition, ...
                    'ElementNormal', obj.ElementNormal, ...
                    'Taper',obj.Taper);
        end
        
        function az_el = incidentAngleToAzEl(obj,inc_ang,m)
        %incidentAngleToAzEl Convert incident angle to local azimuth and
        %elevation angles
        %   AZ_EL = incidentAngleToAzEl(H,INC_ANG,N) converts the incident
        %   angle in the form of [azimuth; elevation] to the local
        %   [azimuth; elevation] pair, AZ_EL, for the Nth element in the
        %   conformal array H.
        %
        %   % Example:
        %   %   Convert direction of [180; 80] to the local azimuth and
        %   %   elevation angle when the heterogeneous conformal array's 
        %   %   element normal points to [0; 90].
        %   
        %   ha = phased.HeterogeneousConformalArray(...
        %       'ElementNormal',[0; 90]);
        %   az_el = incidentAngleToAzEl(ha,[180; 80],1);
            
            % construct the local coordinate system of the element that
            % pointing to that normal direction. 
            ElementAxes = squeeze(obj.pElementAxes(:,:,m));
            % translate incoming direction to the local coordinate system
            az_el = phased.internal.incident2azel(inc_ang,ElementAxes);
        end
    end
            
    methods(Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractHeterogeneousArray;
            props = {...
              'ElementIndices',...
              'ElementPosition',...
              'ElementNormal',...
              'Taper'};
            groups.PropertyList = [groups.PropertyList props];
        end
    end      
    
    methods (Hidden, Access = {?phased.internal.AbstractArray, ?phased.gpu.internal.AbstractClutterSimulator})
        %Methods used by the GPU ConstantGammaClutter model.
        function xscaled = scaleByElemResponse(obj, azin, elin, freq, x, idx)
            %scale x - the steering vector, by the element pattern for the
            %elements selected by the indices in idx
            
            %Convert the ElementNormal to radians and transform the input
            %azimuth and elevation angles to the ElementNormal axes.
            
            nv = getElementNormal(obj);
            normAz = phased.internal.deg2rad(nv(1,idx));
            normEl = phased.internal.deg2rad(nv(2,idx));
            
            [az, el] = phased.gpu.internal.angleAxesConversion(azin, elin, normAz, normEl); 
           
            xscaled = scaleByElemResponse@phased.internal.AbstractHeterogeneousArray(obj, az, el, freq, x, idx);
        end     
       
    end
end

function eleaxes = computeElementAxes(nv,Ne)
    eleaxes = zeros(3,3,Ne);
    for m = 1:Ne
        eleaxes(:,:,m) = phased.internal.rotazel(...
            eye(3),nv(:,m));
    end
end

% [EOF]
