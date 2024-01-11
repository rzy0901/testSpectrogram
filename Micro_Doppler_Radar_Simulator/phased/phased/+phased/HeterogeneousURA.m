classdef (Sealed,StrictDefaults) HeterogeneousURA < phased.internal.AbstractHeterogeneousArray
%HeterogeneousURA   Heterogeneous uniform rectangular array
%   H = phased.HeterogeneousURA creates a heterogeneous uniform rectangular
%   array (URA) System object, H. This object models a heterogeneous URA
%   formed with sensor elements whose pattern may vary. The default array
%   is a 2x2 URA of isotropic antenna elements.
%
%   H = phased.HeterogeneousURA(Name,Value) creates a heterogeneous URA
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   An MxN heterogeneous URA has M rows and N columns.
%
%   Step method syntax:
%
%   RESP = step(H,FREQ,ANGLE) returns the array elements' responses, RESP,
%   given the array's operating frequency FREQ (in Hz) and the directions
%   specified by ANGLE (in degrees). FREQ is a row vector of length L and
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
%   HeterogeneousURA methods:
%
%   step                  - Output responses of the array elements
%   release               - Allow property value and input characteristics
%                           changes
%   clone                 - Create a HeterogeneousURA object with same 
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
%   HeterogeneousURA properties:
%
%   ElementSet            - Set of elements used in the array
%   ElementIndices        - Assign elements to their locations
%   ElementSpacing        - Element spacing of the array
%   Lattice               - Element lattice
%   ArrayNormal           - Array normal
%   Taper                 - Taper on elements
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a 3x4 heterogeneous URA. The elements along the edge 
%   %   have a cosine pattern with a power of 1.5 while the rest elements
%   %   have a cosine pattern with a power of 1.6. Plot its azimuth 
%   %   response. Assume the operating frequency is 1 GHz and the wave 
%   %   propagation speed is 3e8 m/s.
%
%   sElement1 = phased.CosineAntennaElement('CosinePower',1.6);
%   sElement2 = phased.CosineAntennaElement('CosinePower',1.5);
%   sArray = phased.HeterogeneousURA(...
%       'ElementSet',{sElement1,sElement2},...
%       'ElementIndices',[2 2 2 2;2 1 1 2;2 2 2 2]);
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
%   sArray = phased.HeterogeneousURA(...
%       'ElementSet',{sElement1,sElement2},...
%       'ElementIndices',[2 2 2 2;2 1 1 2;2 2 2 2]);
%   fc = 1e9; ang = [0;0];
%   resp = sArray(fc,ang)
%
%   See also phased, phased.ULA, phased.URA, phased.ConformalArray.

%   Copyright 2012-2016 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] Eli Brookner, Radar Technology, Lex Book, 1996
%   [3] Eli Brookner, Practical Phased Array Antenna Systems, Artech House,
%       1991
%   [4] Robert Mailloux, Phased Array Theory and Technology, Proceedings of
%       The IEEE, Vol. 70, No. 3, 1982
%   [5] Harold Mott, Antennas for Radar and Communications, A Polarimetric
%       Approach, John Wiley & Sons, 1992

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)

        %ElementIndices  Assign elements to their locations
        %   This property specifies the mapping of elements in the array.
        %   The property assigns elements to their locations in the array
        %   using the indices into the ElementSet property. ElementIndices
        %   must be an MxN matrix where M is the number of rows and N is
        %   the number of columns. Rows are along y-axis and columns are
        %   along z-axis. The values in ElementIndices should be less than
        %   or equal to the number of entries in the Element property. The
        %   default value of this property is [1 1; 1 1].
        ElementIndices = ones([2 2])
        %ElementSpacing   Element spacing of the array
        %   Specify the element spacing (in meters) of the array as a 1x2
        %   vector or a scalar. If ElementSpacing is a 1x2 vector, it is in
        %   the form of [SpacingBetweenRows SpacingBetweenColumns]. Rows
        %   are along y-axis and columns are along z-axis. If
        %   ElementSpacing is a scalar, the spacing along both axes are the
        %   same. The default value of this property is [0.5 0.5].
        ElementSpacing = [0.5 0.5];
        %Lattice    Element lattice
        %   Specify the element lattice as one of 'Rectangular' |
        %   'Triangular', where the default is 'Rectangular'. When you set
        %   the Lattice property to 'Rectangular', all elements are aligned
        %   in both row and column directions. When you set the Lattice
        %   property to 'Triangular', the elements in even rows are shifted
        %   toward the positive row axis direction by a distance of half
        %   element spacing.
        Lattice = 'Rectangular';
        %ArrayNormal Array normal
        %   Specify the normal of the plane on which the elements are
        %   placed as one of 'x'|'y'|'z', where the default is 'x'.
        ArrayNormal = 'x';
        %Taper  Taper
        %    Specify the Taper property as a scalar, a P element vector, or
        %    an MxN matrix of complex weights applied to each element in
        %    the sensor array. P is the total number of elements, M is the
        %    number of elements along z-axis and N is the number of
        %    elements along y-axis. If Taper is a scalar the same weights
        %    will be applied to each element. If Taper is a matrix, each
        %    weight will be applied to the corresponding element. The
        %    default value of Taper is 1.
        Taper = 1;
        

    end

    properties (Access = private, Nontunable)
        pNumElements;     % private property to hold number of elements
    end

    properties(Constant, Hidden)
        LatticeSet = matlab.system.StringSet({'Rectangular','Triangular'});
        ArrayNormalSet = matlab.system.StringSet({'x','y','z'});
    end
    
    methods
        function set.ElementSpacing(obj,valueArg)
            if isscalar(valueArg)
                value = [valueArg valueArg];
            else
                value = valueArg;
            end
            sigdatatypes.validateDistance(value,...
                '','ElementSpacing',...
                {'size',[1 2],'positive','finite'});
            obj.ElementSpacing = value;
        end

        function set.ElementIndices(obj,val)
            sigdatatypes.validateIndex(val,'','ElementIndices',...
                {'2d'});
            cond = (any(size(val) < 2));
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidURAElementMap','ElementIndices');
            end
            obj.ElementIndices = val;
        end
        
        function set.Taper(obj,val)
            validateattributes(val,{'double'},...
                               {'nonnan','nonempty','finite','2d'},...
                               '','Taper');
            obj.Taper = val;
        end
        
    end
    
    methods (Access = protected)
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractHeterogeneousArray(obj);
            sz = getURASize(obj);
            cond =  ~isscalar(obj.Taper) && any(size(obj.Taper) ~= sz) ...
                    && any(size(obj.Taper) ~= [1 sz(1)*sz(2)]) ...
                    && any(size(obj.Taper) ~= [sz(1)*sz(2) 1]);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidURATaper','Taper',...
                    sprintf('%d',sz(1)*sz(2)), ...
                    sprintf('%d %d',sz(1),sz(2)));
            end    
        end
        
        function setupImpl(obj, varargin)
            setupImpl@phased.internal.AbstractHeterogeneousArray(obj, varargin{:});
            setNumElements(obj);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractHeterogeneousArray(obj);
            if isLocked(obj)
                s.pNumElements = obj.pNumElements;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function arNr = getArrayNormal(obj)
            arNr = obj.ArrayNormal;
        end
        
        function arPl = getArrayPlaneOrAxis(obj)
            arNr = obj.ArrayNormal;
            if strcmpi(arNr,'x')
                arPl = 'YZ';
            elseif strcmpi(arNr,'y')
                arPl = 'XZ';
            else
                arPl = 'XY';
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
        
            Ns = numel(obj.cElement);
            Na = size(inc_ang,2);
            az_el_set = zeros(2,Ns*Na);
            
            inc_ang_rot = incidentAngleToAzEl(obj,inc_ang);
            
            for m = 1:Ns
                % check if the element normal and array normal are aligned. 
                isAligned = isElementNormalArrayNormalAligned(obj.cElement{m});    
                if(isAligned)
                    az_el_set(:,(1:Na)+(m-1)*Na) = inc_ang_rot;
                else
                    az_el_set(:,(1:Na)+(m-1)*Na) = inc_ang;
                end
            end
            
            % obtain the normal direction
            
            lMap = bsxfun(@plus,(obj.ElementIndices(:)-1)*Na,1:Na);
            
            az_el = az_el_set(:,lMap(:).');
        end
    end
    
    methods (Access = private)
        function setNumElements(obj)
            obj.pNumElements = numel(obj.ElementIndices);
        end
        
        function sz = getURASize(obj)
            sz = size(obj.ElementIndices);
        end
        
        function az_el = incidentAngleToAzEl(obj,inc_ang)
        %incidentAngleToAzEl Convert incident angle to local azimuth and
        %elevation angles
        %   AZ_EL = incidentAngleToAzEl(H,INC_ANG,N) converts the incident
        %   angle in the form of [azimuth; elevation] to the local
        %   [azimuth; elevation] pair, AZ_EL, for the Nth element in the
        %   heterogeneous uniform rectangular array H.
            
            % construct the local coordinate system of the element that
            % pointing to that normal direction. 
            arNr = obj.getArrayNormal;
            if (strcmpi(arNr,'x'))
               az_el = inc_ang;
            elseif(strcmpi(arNr,'y'))
                ElementAxes = [0 -1 0;1 0 0;0 0 1];
                % translate incoming direction to the local coordinate system
                az_el = phased.internal.incident2azel(inc_ang,ElementAxes);
            else
                ElementAxes = [0 0 -1;0 1 0;1 0 0];
                % translate incoming direction to the local coordinate system
                az_el = phased.internal.incident2azel(inc_ang,ElementAxes);
            end
            
        end
    end
    
    methods

        function obj = HeterogeneousURA(varargin)
        %HeterogeneousURA   Construct the HeterogeneousURA class.
            obj@phased.internal.AbstractHeterogeneousArray(nargin,...
                varargin{:});

        end

        function N = getNumElements(obj)
        %getNumElements   Number of elements in the array
        %   N = getNumElements(H) returns the number of elements, N, in
        %   the heterogeneous URA object H. 
        %
        %   % Example:
        %   %   Construct a default heterogeneous URA and obtain its number
        %   %   of elements.
        %
        %   ha = phased.HeterogeneousURA;
        %   N = getNumElements(ha)

            N = numel(obj.ElementIndices);
        end
        
        function pos = getElementPosition(obj,EleIdx)
        %getElementPosition Element positions of the array
        %   POS = getElementPosition(H) returns the element positions of
        %   the heterogeneous URA H. POS is a 3xN matrix where N is the
        %   number of elements in H. Each column of POS defines the
        %   position, in the form of [x; y; z] (in meters), of an element
        %   in the local coordinate system.
        %
        %   For details regarding the local coordinate system of a
        %   heterogeneous URA, type phased.URA.coordinateSystemInfo.
        %
        %   POS = getElementPosition(H,ELEIDX) returns the positions of the
        %   elements that are specified in the element index vector ELEIDX.
        %   The index of a heterogeneous URA runs through each column, one
        %   after another. For example, in a heterogeneous URA with 4
        %   elements along z-axis (column) and 3 elements along y-axis
        %   (row), the element in the 2nd row and 3rd column has the index
        %   value of 10.
        %
        %   % Example:
        %   %   Construct a default heterogeneous URA and obtain its 
        %   %   element positions.
        %
        %   ha = phased.HeterogeneousURA;
        %   pos = getElementPosition(ha)
            
            narginchk(1,2);
            validateattributes(obj,{'phased.HeterogeneousURA'},{'scalar'},...
                'getElementPosition','H');
            N = getNumElements(obj);
            if nargin < 2
                EleIdx = 1:N;
            end
            sigdatatypes.validateIndex(EleIdx,...
                'getElementPosition','ELEIDX',{'vector','<=',N});
            EleIdx = EleIdx(:);
            array_sz = getURASize(obj);
            NPerRow = array_sz(2);     % number of elements in each row
            NPerCol = array_sz(1);     % number of elements in each column
            deltaRows = (NPerRow-1)/2+1;   % delta along the row
            deltaCols = (NPerCol-1)/2+1;   % delta along the column
            [IdxInCol, IdxInRow] = ind2sub(array_sz,EleIdx);
            pos = zeros(3,numel(EleIdx));
            pos2 = (IdxInRow-deltaRows)*obj.ElementSpacing(2);
            if strcmp(obj.Lattice,'Triangular')
                evenrowidx = ~rem(IdxInCol,2);
                pos2(evenrowidx) = pos2(evenrowidx) + ...
                    obj.ElementSpacing(2)/2;
            end
            arNr = obj.ArrayNormal;
            if strcmpi(arNr,'x')
               pos(2,:) = pos2;
               pos(3,:) = (deltaCols-IdxInCol)*obj.ElementSpacing(1); 
            elseif strcmpi(arNr,'y')
               pos(1,:) = pos2;
               pos(3,:) = (deltaCols-IdxInCol)*obj.ElementSpacing(1); 
            else
               pos(1,:) = pos2;
               pos(2,:) = (deltaCols-IdxInCol)*obj.ElementSpacing(1);
            end
        end
        
        function nv = getElementNormal(obj, EleIdx)
        %getElementNormal Get normal directions of array elements
        %   NV = getElementNormal(H) returns the element normals of the
        %   heterogeneous URA H. NORM is a 2xN matrix where N is the number 
        %   of elements in H. each column of NORM specifies the normal 
        %   direction of the corresponding element in the form of 
        %   [azimuth; elevation] (in degrees) defined in the local 
        %   coordinate system.
        %
        %   For details regarding the local coordinate system of a
        %   heterogeneous URA, type phased.URA.coordinateSystemInfo.
        %
        %   NV = getElementNormal(H,ELEIDX) returns the normals of the
        %   elements that are specified in the element index vector ELEIDX.
        %
        %   % Example:
        %   %   Construct a default heterogeneous URA and obtain its 
        %   %   element normals.
        %
        %   ha = phased.HeterogeneousURA;
        %   nv = getElementNormal(ha)    
            phased.internal.narginchk(1,2,nargin);
            arNr = obj.ArrayNormal;
            N = getNumElements(obj);
            
            if nargin < 2
                EleIdx = 1:N;
            end
            sigdatatypes.validateIndex(EleIdx,...
                'getElementPosition','ELEIDX',{'vector','<=',N});

            nv_temp = zeros(2,N);
            
            if ~strcmp(arNr,'x')
            
                for m = 1:numel(obj.ElementSet)
                    % check if the element normal and array normal are aligned
                    isAligned = ~isElementFromAntenna(obj.ElementSet{m}) && ...
                        isElementNormalArrayNormalAligned(obj.ElementSet{m});

                    if isAligned
                        rotidx = (obj.ElementIndices(:)==m);
                        if strcmp(arNr,'y')
                            nv_temp(:,rotidx) = repmat([90;0],[1 sum(rotidx)]);
                        else % z
                            nv_temp(:,rotidx) = repmat([0;90],[1 sum(rotidx)]);   
                        end
                    end
                end

            end

            nv = nv_temp(:,EleIdx);
        end
        
        function w = getTaper(obj)
        %getTaper Get array element tapers
        %   W = getTaper(H) returns the element taper of the heterogeneous
        %   URA H. W is a length-N column vector where N is the number of
        %   elements in H.
        %
        %   % Example:
        %   %   Construct a default heterogeneous URA and obtain its 
        %   %   element taper.
        %
        %   ha = phased.HeterogeneousURA;
        %   pos = getTaper(ha)
            
            wArg = obj.Taper(:);
            if isscalar(wArg)
                w = wArg*ones(getNumElements(obj),1);
            else
                w = wArg;
            end
        end
    end
    
    methods(Hidden)
        function cl = clonecg(obj)
            Ne = numel(obj.ElementSet);
            newElementSet = cell(1,Ne);
            for m = 1:Ne
                newElementSet{m} = clonecg(obj.ElementSet{m});
            end
            cl = phased.HeterogeneousURA(...
                    'ElementSet',newElementSet, ...
                    'ElementIndices',obj.ElementIndices, ...                
                    'ElementSpacing',obj.ElementSpacing, ...
                    'Lattice',obj.Lattice, ...
                    'ArrayNormal', obj.ArrayNormal, ...
                    'Taper',obj.Taper);
        end
        
    end
    
    methods(Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractHeterogeneousArray;
            props = {...
                'ElementIndices',...
                'ElementSpacing',...
                'Lattice',...
                'ArrayNormal',...
                'Taper'};
            groups.PropertyList = [groups.PropertyList props];
        end
          
    end
    
end

% [EOF]
