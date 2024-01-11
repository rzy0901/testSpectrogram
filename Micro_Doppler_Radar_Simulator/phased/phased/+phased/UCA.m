classdef (Sealed,StrictDefaults) UCA < phased.internal.AbstractHomogeneousArray
%UCA   Uniform circular array
%   H = phased.UCA creates a uniform circular array (UCA) System object, H.
%   This object models a UCA formed with identical sensor elements. The
%   default array is a 5-element UCA of five isotropic antenna elements.
%
%   H = phased.UCA(Name,Value) creates a UCA object, H, with the specified
%   property Name set to the specified Value. You can specify additional
%   name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   H = phased.UCA(N,R,Name,Value) creates a UCA object, H, with the
%   NumElements property set to N, Radius property set to R and other
%   specified property Names set to the specified Values. N and R are
%   value-only arguments. To specify a value-only argument, you must also
%   specify all preceding value-only arguments. You can specify name-value
%   pair arguments in any order.
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
%   UCA methods:
%
%   step                   - Output responses of the array elements
%   release                - Allow property value and input characteristics
%                            changes
%   clone                  - Create a URA object with same property values
%   isLocked               - Locked status (logical)
%   isPolarizationCapable  - Indicate if the array is capable of 
%                            simulating polarization
%   getNumElements         - Get number of elements in the array
%   getElementPosition     - Get positions of array elements
%   getElementNormal       - Get normal directions of array elements
%   getElementSpacing      - Get distance between array elements    
%   getTaper               - Get array element tapers
%   directivity            - Compute array directivity
%   pattern                - Plot array response pattern
%   patternAzimuth         - Plot azimuth pattern
%   patternElevation       - Plot elevation pattern
%   collectPlaneWave       - Simulate received plane waves
%   viewArray              - View array geometry
%                          
%   UCA properties:        
%                          
%   Element                - Element of the array
%   NumElements            - Number of elements
%   Radius                 - Radius of UCA
%   ArrayNormal            - Array normal
%   Taper                  - Taper
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a 15 element UCA with a radius of 35 cm and plot its  
%   %   pattern. Assume the operating frequency is 1 GHz and the wave
%   %   propagation speed is 3e8 m/s.
%
%   array = phased.UCA(15,0.35);
%   fc = 1e9; c = 3e8;
%   pattern(array,fc,'CoordinateSystem','polar','PropagationSpeed',c);
%
%   % Example 2:
%   %   Find the response of each element in the above array at the 
%   %   boresight.
%
%   array = phased.UCA(15,0.35);
%   fc = 1e9; ang = [0;0];
%   resp = array(fc,ang)
%
%   See also phased, phased.ULA, phased.URA,
%   phased.ConformalArray.

%   Copyright 2014-2016 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)

        %NumElements   Number of elements
        %   An integer containing the number of elements in the array. The
        %   default value of this property is 5. 
        NumElements = 5
        %Radius   Radius of UCA (m)
        %   A scalar specifying the radius(in meters) of the UCA. The
        %   default value of this property is 0.5
        Radius = 0.5
        %ArrayNormal Array normal
        %   Specify the normal of the plane on which the elements are
        %   placed as one of 'x'|'y'|'z', where the default is 'z'.
        ArrayNormal = 'z';
        %Taper  Taper
        %    Specify the Taper property as a scalar or a length-N vector of
        %    complex weights applied to each element in the sensor array. N
        %    is the number of elements in the array. If Taper is a scalar
        %    the same weights will be applied to each element. If Taper is
        %    a vector, each weight will be applied to the corresponding
        %    element. The default value of Taper is 1.
        Taper = 1;

    end
    
    properties (Constant, Hidden)
        ArrayNormalSet = matlab.system.StringSet({'x','y','z'});
    end
    
    properties (Access = private, Nontunable)
        pElementNormal
    end
    
    properties (Access = private)
        % private property to hold element axes for each element. Once the
        % normals are set, the element axes are fixed. So we don't have to
        % compute it at each simulation step.
        pElementAxes
    end
    
    methods
        function set.NumElements(obj,value)
            sigdatatypes.validateIndex(value,...
                '','NumElements',...
                {'scalar','>=',2});
            obj.NumElements = value;
        end
        
        function set.Radius(obj,value)
            sigdatatypes.validateDistance(value,...
                '',...
                'Radius',{'scalar','positive','finite'});
            obj.Radius = value;
        end
        
        function set.Taper(obj,val)
            validateattributes(val,{'double'},...
                               {'nonnan','nonempty','finite','vector'},...
                               '','Taper');
            obj.Taper = val;
        end
        
    end
    
    methods

        function obj = UCA(varargin)
        %UCA   Construct the UCA class.
            obj@phased.internal.AbstractHomogeneousArray(nargin, ...
               varargin{:},'NumElements','Radius');
        end
    end
    methods(Hidden)
        function cl = clonecg(obj)
            cl = phased.UCA(...
                    'Element',clonecg(obj.Element), ...
                    'NumElements',obj.NumElements, ...                
                    'Radius',obj.Radius, ...
                    'ArrayNormal', obj.ArrayNormal,...            
                    'Taper',obj.Taper);
        end
        
        function az_el = incidentAngleToAzEl(obj,inc_ang,m)
        %incidentAngleToAzEl Convert incident angle to local azimuth and
        %elevation angles
        %   AZ_EL = incidentAngleToAzEl(H,INC_ANG,N) converts the incident
        %   angle in the form of [azimuth; elevation] to the local
        %   [azimuth; elevation] pair, AZ_EL, for the Nth element in the
        %   uniform circular array H.
          
            % construct the local coordinate system of the element that
            % pointing to that normal direction. 
            ElementAxes = squeeze(obj.pElementAxes(:,:,m));
            % translate incoming direction to the local coordinate system
            az_el = phased.internal.incident2azel(inc_ang,ElementAxes);
        end
    end
    
    methods
        function N = getNumElements(obj)
        %getNumElements   Number of elements in the array
        %   N = getNumElements(H) returns the number of elements, N, in
        %   the UCA object H. 
        %
        %   % Example:
        %   %   Construct a default UCA and obtain its number of elements.
        %
        %   ha = phased.UCA;
        %   N = getNumElements(ha)
            N = obj.NumElements;
        end
               
        function pos  = getElementPosition(obj,varargin)
        %getElementPosition Element positions of the array
        %   POS = getElementPosition(H) returns the element positions of
        %   the UCA H. POS is a 3xN matrix where N is the number of
        %   elements in H. Each column of POS defines the position, in the
        %   form of [x; y; z] (in meters), of an element in the local
        %   coordinate system.
        %
        %   For details regarding the local coordinate system of a UCA,
        %   type phased.UCA.coordinateSystemInfo.
        %
        %   POS = getElementPosition(H,ELEIDX) returns the positions of the
        %   elements that are specified in the element index vector ELEIDX.
        %
        %   % Example:
        %   %   Construct a default UCA and obtain its element positions.
        %
        %   ha = phased.UCA;
        %   pos = getElementPosition(ha)
        
            phased.internal.narginchk(1,2,nargin);
            pos  = getElementPositionAndNormal(obj,varargin{:});
            
        end
        
        function nv  = getElementNormal(obj,varargin)
        %getElementNormal Get normal directions of array elements
        %   NV = getElementNormal(H) returns the element normals of the
        %   UCA H. NORM is a 2xN matrix where N is the number of elements
        %   in H. each column of NORM specifies the normal direction of the
        %   corresponding element in the form of [azimuth; elevation] (in
        %   degrees) defined in the local coordinate system.
        %
        %   For details regarding the local coordinate system of a UCA,
        %   type phased.UCA.coordinateSystemInfo.
        %
        %   NV = getElementNormal(H,ELEIDX) returns the normals of the
        %   elements that are specified in the element index vector ELEIDX.
        %
        %   % Example:
        %   %   Construct a default UCA and obtain its element normals.
        %
        %   ha = phased.UCA;
        %   nv = getElementNormal(ha)
        
            phased.internal.narginchk(1,2,nargin);
            [~, nv]  = getElementPositionAndNormal(obj,varargin{:});   
        end
        
        function dist  = getElementSpacing(obj,type)
        %getElementSpacing Distance between array elements    
        %   DIST = getElementSpacing(H) returns the arc length, DIST, in
        %   meters between two adjacent elements of the UCA, H.  
        %
        %   Note that the syntax getElementSpacing(H) is equivalent to the
        %   syntax getElementSpacing(H,'arc').
        %
        %   DIST = getElementSpacing(H,'chord') returns the chord length 
        %   in meters of two adjacent elements. 
        %
        %   % Example:
        %   %   Construct a default UCA and obtain the distance along the 
        %   %   circumference between its element positions.
        %
        %   ha = phased.UCA;
        %   dist = getElementSpacing(ha,'arc')
        
            phased.internal.narginchk(1,2,nargin);
            if  nargin == 2
                type = validatestring(type,{'arc','chord'},'getElementSpacing','',2);
            else
                type = 'arc';
            end
            if strcmp(type,'arc')
                dist = 2*pi*obj.Radius/obj.NumElements;
            else
                dist = 2*obj.Radius*sin(pi/obj.NumElements);
            end
        end
        
        function w = getTaper(obj)
        %getTaper Get array element tapers
        %   W = getTaper(H) returns the element taper of the UCA H. W is a
        %   length-N column vector where N is the number of elements in H.
        %
        %   % Example:
        %   %   Construct a default UCA and obtain its element taper.
        %
        %   ha = phased.UCA;
        %   pos = getTaper(ha)
            
            ws = obj.Taper(:);
            if isscalar(ws)
                w = ws*ones(getNumElements(obj),1);
            else
                w = ws;
            end
        end
        
    end
    
    methods(Access = private)
     function [pos, nv] = getElementPositionAndNormal(obj,EleIdx)
         coder.extrinsic('phased.UCA.calcElementPosition');  
            phased.internal.narginchk(1,2,nargin);
            N = getNumElements(obj);
            arNr = getArrayNormal(obj);
            if nargin < 2
                EleIdx = 1:N;
            end
            sigdatatypes.validateIndex(EleIdx,...
                'getElementPosition','ELEIDX',{'vector','<=',N});
            
            % compute normal directions along circle
            EleIdx = EleIdx(:);
            allNormals = (-(N-1)/2:(N-1)/2)*360/N;
            nv_temp = allNormals(EleIdx);
            numElIDX = numel(EleIdx);
            
            % compute positions
            pos = ...
                coder.const(obj.calcElementPosition(numElIDX,nv_temp,obj.Radius,arNr));
            [az,el] =cart2sph(pos(1,:),pos(2,:),pos(3,:));
            nv = [rad2deg(az);rad2deg(el)];
         
        end
    end
    
    methods (Access = protected)
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractHomogeneousArray(obj);
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
            setupImpl@phased.internal.AbstractHomogeneousArray(obj, varargin{:});
            [~, nv] = getElementPositionAndNormal(obj);
            obj.pElementNormal = nv;
            obj.pElementAxes = getElementAxes(getNumElements(obj),nv);
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractHomogeneousArray(obj);
            if isLocked(obj)
                s.pElementNormal = obj.pElementNormal;
                s.pElementAxes = obj.pElementAxes;
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
            if strcmp(obj.ArrayNormal,'z')
 
                nv = obj.pElementNormal; 
                elemNormal = nv(1,:);
                for m = 1:M 
                    frameNum = N*(m-1); 
                    az = inc_ang(1,m); 
                    el = inc_ang(2,m); 
                    for n = 1:N 
                        relativeAz = az - elemNormal(n); 
                        if(relativeAz < -180) 
                           az_el(1,n+frameNum) = relativeAz + 360;   
                        elseif(relativeAz > 180) 
                           az_el(1,n+frameNum) = relativeAz - 360; 
                        else 
                           az_el(1,n+frameNum) = relativeAz;  
                        end 
                        az_el(2,n+frameNum) = el; 
                    end 
                end
            else
                for m = 1:N
                    az_el(:,(0:N:N*(M-1))+m) = incidentAngleToAzEl(obj,inc_ang,m);
                end
            end
         end 
        
        function arNr = getArrayNormal(obj)
        %getArrayNormal Get array normal
        %   arNr = getArrayNormal(H) returns the normal of the UCA H. arNr
        %   is a char containing the normal to the plane within which the 
        %   array is positioned
        %
        %   % Example:
        %   % Construct a default UCA and obtain its normal.
        %   
        %   ha = phased.UCA;
        %   arNr = getArrayNormal(ha);
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
        
    end
    
    methods (Static, Hidden)
        function coordinateSystemInfo
            %coordinateSystemInfo Displays local coordinate system
            %information
            
            fprintf([...
            '    For a UCA, the origin of the local coordinate system, i.e.,\n',...
            '    the phase center of the array, is defined at the center of\n',...
            '    the array. In local coordinate system, the default array\n',...
            '    normal is along z axis so the elements are located on the\n', ... 
            '    xy plane. The element normals are in the xy plane and \n', ...
            '    are pointing radially outwards from the center of the UCA.\n\n', ...
            '    If the array normal is set to x axis, then the elements are\n', ... 
            '    located on the yz plane. The element normals are in the yz \n', ...
            '    plane and are pointing radially outwards from the center of\n',...
            '    the UCA. If the array normal is set to y axis, then the \n', ... 
            '    elements are located on the xz plane. The element normals \n', ...
            '    are in the xz plane and are pointing radially outwards from\n',...
            '    the center of the UCA.\n\n',... 
            '    The azimuth angle is defined as the angle from x axis\n',...
            '    toward y axis. The elevation angle is defined as the angle\n',...
            '    from xy plane toward z axis.\n']);
        end
        function pos = calcElementPosition(numElIDX,nv,radius, arNr)    
            
            if strcmpi(arNr,'x')
            pos = ... 
                [zeros(1,numElIDX);...
                 radius*cosd(nv); ...
                 radius*sind(nv)];
            
            elseif strcmpi(arNr,'y')
            pos = ... 
                [radius*cosd(nv); ...
                 zeros(1,numElIDX);...
                 radius*sind(nv)];
            else
            pos = ... 
                [radius*cosd(nv); ...
                 radius*sind(nv); ...
                 zeros(1,numElIDX)];
            end
         end
    end
    
    methods(Static,Hidden,Access=protected)        
        function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractHomogeneousArray;
        props = {'NumElements',...
                 'Radius',...
                 'ArrayNormal',...
                 'Taper'};
        groups.PropertyList = [groups.PropertyList props];
        end
    end
end

function eleAxes = getElementAxes(numElements,elementNormal)
    eleAxes = zeros(3,3,numElements);
    for m = 1:numElements
        eleAxes(:,:,m) = phased.internal.rotazel(...
            eye(3),elementNormal(:,m));
    end
end

