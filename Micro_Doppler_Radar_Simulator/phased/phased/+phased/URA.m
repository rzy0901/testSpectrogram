classdef (Sealed,StrictDefaults) URA < phased.internal.AbstractHomogeneousArray
%URA   Uniform rectangular array
%   H = phased.URA creates a uniform rectangular array (URA) System object,
%   H. This object models a URA formed with identical sensor elements. The
%   default array is a 2x2 URA of isotropic antenna elements.
%
%   H = phased.URA(Name,Value) creates a URA object, H, with the specified
%   property Name set to the specified Value. You can specify additional
%   name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   H = phased.URA(SZ,D,Name,Value) creates a URA object, H, with the Size
%   property set to SZ, ElementSpacing property set to D and other
%   specified property Names set to the specified Values. SZ and D are
%   value-only arguments. To specify a value-only argument, you must also
%   specify all preceding value-only arguments. You can specify name-value
%   pair arguments in any order.
%
%   An MxN URA has M rows and N columns.
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
%   URA methods:
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
%   getTaper               - Get array element tapers
%   directivity            - Compute array directivity
%   pattern                - Plot array response pattern
%   patternAzimuth         - Plot azimuth pattern
%   patternElevation       - Plot elevation pattern
%   plotGratingLobeDiagram - Plot the grating lobe diagram of the array
%   collectPlaneWave       - Simulate received plane waves
%   viewArray              - View array geometry
%                          
%   URA properties:        
%                          
%   Element                - Element of the array
%   Size                   - Array size
%   ElementSpacing         - Element spacing
%   Lattice                - Element lattice
%   ArrayNormal            - Array normal
%   Taper                  - Taper
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a 2x3 URA and plot its azimuth response. Assume the 
%   %   operating frequency is 1 GHz and the wave propagation speed is 3e8 
%   %   m/s.
%
%   array = phased.URA([2 3]);
%   fc = 1e9; c = 3e8;
%   pattern(array,fc,-180:180,0,'CoordinateSystem','polar',...
%       'PropagationSpeed',c);
%
%   % Example 2:
%   %   Find the response of each element in the above array at the 
%   %   boresight.
%
%   array = phased.URA([2 3]);
%   fc = 1e9; ang = [0;0];
%   resp = array(fc,ang)
%
%   See also phased, phased.ULA, phased.ConformalArray, 
%   phased.HeterogeneousURA.

%   Copyright 2009-2016 The MathWorks, Inc.

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

        %Size   Array size
        %   Specify the size of the array as a 1x2 integer vector or an
        %   integer. If Size is a 1x2 vector, it is in the form of
        %   [NumberOfRows NumberOfColumns] where rows are along y-axis and
        %   columns are along z-axis. If Size is a scalar, the array has
        %   the same number of elements along both axes. The default value
        %   of this property is [2 2]. 
        Size = [2 2]
        %ElementSpacing   Element spacing (m)
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
        function set.Size(obj,valueArg)
            if isscalar(valueArg)
                value = [valueArg valueArg];
            else
                value = valueArg;
            end
            sigdatatypes.validateIndex(value,...
                '','Size',...
                {'size',[1 2],'>=',2});
            obj.Size = value;
        end

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

        function set.Taper(obj,val)
            validateattributes(val,{'double'},...
                               {'nonnan','nonempty','finite','2d'},...
                               '','Taper');
            obj.Taper = val;
        end
        
    end
    
    methods (Access = protected)
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractHomogeneousArray(obj);
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
            setupImpl@phased.internal.AbstractHomogeneousArray(obj,varargin{:});
            setNumElements(obj);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractHomogeneousArray(obj);
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
        
            % check if the element normal and array normal are aligned. 
            isAligned = isElementNormalArrayNormalAligned(obj.cElement);
            
            % obtain the normal direction
            N = getNumElements(obj);
            if (isAligned) 
                az_el_temp = incidentAngleToAzEl(obj,inc_ang);
            else
                az_el_temp = inc_ang;
            end
            az_el = repmat(az_el_temp,N,1);
            az_el = reshape(az_el,2,[]);

        end
    end
    
    methods (Access = private)
        function setNumElements(obj)
            obj.pNumElements = prod(obj.Size);
        end
        
        function sz = getURASize(obj)
            sz = obj.Size;
        end
    end
    
    methods

        function obj = URA(varargin)
        %URA   Construct the URA class.
            obj@phased.internal.AbstractHomogeneousArray(nargin, ...
                 varargin{:},'Size', 'ElementSpacing');

        end
        
        function N = getNumElements(obj)
        %getNumElements   Number of elements in the array
        %   N = getNumElements(H) returns the number of elements, N, in
        %   the URA object H.
        %
        %   % Example:
        %   %   Construct a default URA and obtain its number of elements.
        %
        %   ha = phased.URA;
        %   N = getNumElements(ha)
            
            %if ~isempty(coder.target)|| ~isLocked(obj)
            %    setNumElements(obj);
            %end
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
        %   the URA H. POS is a 3xN matrix where N is the number of
        %   elements in H. Each column of POS defines the position, in the
        %   form of [x; y; z] (in meters), of an element in the local
        %   coordinate system.
        %
        %   For details regarding the local coordinate system of a URA,
        %   type phased.URA.coordinateSystemInfo.
        %
        %   POS = getElementPosition(H,ELEIDX) returns the positions of the
        %   elements that are specified in the element index vector ELEIDX.
        %   The index of a URA runs through each column, one after another.
        %   For example, in a URA with 4 elements along z-axis (column) and
        %   3 elements along y-axis (row), the element in the 2nd row and
        %   3rd column has the index value of 10.
        %
        %   % Example:
        %   %   Construct a default URA and obtain its element positions.
        %
        %   ha = phased.URA;
        %   pos = getElementPosition(ha)
            
            coder.extrinsic('phased.URA.calcElementPosition');
            phased.internal.narginchk(1,2,nargin);
            validateattributes(obj,{'phased.URA'},{'scalar'},...
                'getElementPosition','H');
            N = getNumElements(obj);
            arNr = getArrayNormal(obj);
            if nargin < 2
                EleIdx = 1:N;
            end
            sigdatatypes.validateIndex(EleIdx,...
                'getElementPosition','ELEIDX',{'vector','<=',N});
            EleIdx = EleIdx(:);
            array_sz = getURASize(obj);
            pos = ...
                coder.internal.const(obj.calcElementPosition(EleIdx,array_sz,obj.ElementSpacing,obj.Lattice, arNr));
        end 
        
        function nv = getElementNormal(obj, EleIdx)
        %getElementNormal Get normal directions of array elements
        %   NV = getElementNormal(H) returns the element normals of the
        %   URA H. NORM is a 2xN matrix where N is the number of elements
        %   in H. each column of NORM specifies the normal direction of the
        %   corresponding element in the form of [azimuth; elevation] (in
        %   degrees) defined in the local coordinate system.
        %
        %   For details regarding the local coordinate system of a ULA,
        %   type phased.URA.coordinateSystemInfo.
        %
        %   NV = getElementNormal(H,ELEIDX) returns the normals of the
        %   elements that are specified in the element index vector ELEIDX.
        %
        %   % Example:
        %   %   Construct a default URA and obtain its element normals.
        %
        %   ha = phased.URA;
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
            
            % check if the element normal and array normal are aligned
            if ~strcmp(arNr,'x') && ...
                    ~isElementFromAntenna(obj.Element) && ...
                    isElementNormalArrayNormalAligned(obj.Element); 
                                        
                if strcmp(arNr,'y')
                    nv_temp = repmat([90;0],[1 N]);
                else % 'z'
                    nv_temp = repmat([0;90],[1 N]);
                end
            end

            nv = nv_temp(:,EleIdx);
            
        end
        
        function w = getTaper(obj)
        %getTaper Get array element tapers
        %   W = getTaper(H) returns the element taper of the URA H. W is a
        %   length-N column vector where N is the number of elements in H.
        %
        %   % Example:
        %   %   Construct a default URA and obtain its element taper.
        %
        %   ha = phased.URA;
        %   pos = getTaper(ha)

            ws = obj.Taper(:);
            if isscalar(ws)
                w = ws*ones(getNumElements(obj),1);
            else
                w = ws;            
            end
        end
        

        
        function varargout = plotGratingLobeDiagram(obj,f,ang,c,f0)
            %plotGratingLobeDiagram  Plot grating lobe diagram in sine
            %   space (u v coordinates)
            %
            %   plotGratingLobeDiagram(H,FREQ) plots the grating lobe
            %   diagram of the array, H, at the specified frequency, FREQ,
            %   in Hz and at 0 degrees azimuth and elevation.
            %
            %   plotGratingLobeDiagram(H,FREQ,ANGLE) plots the grating lobe
            %   diagram at a specified steering angle, ANGLE. ANGLE can be
            %   either a 2x1 vector or a scalar. If ANGLE is a vector, the
            %   form is [azimuth; elevation] (in degrees) where the azimuth
            %   angle must be between [-180 180] and the elevation angle
            %   must be between [-90 90]. If ANGLE is a scalar, ANGLE
            %   specifies the azimuth angle and the corresponding elevation
            %   angle is assumed to be 0 degrees. 
            %
            %   plotGratingLobeDiagram(H,FREQ,ANGLE,C) plots the grating
            %   lobe diagram at a specified propagation speed, C, in m/s. 
            %   The default propagation speed is the speed of light.
            %
            %   plotGratingLobeDiagram(H,FREQ,ANGLE,C,F0) plots the grating
            %   lobe diagram at a specified phased-shift frequency, F0, in
            %   Hz. The default value of F0 is FREQ.
            %   
            %   h = plotGratingLobeDiagram(H,FREQ,ANGLE,C,F0) plots the
            %   grating lobe diagram and returns its axes handle.
            %
            %   % Examples:
            %
            %   % Example 1: 
            %   %   Construct a rectangular lattice URA with elements 
            %   %   spaced at 0.8 wavelength at 5 GHz operating frequency.
            %   %   The array is steered to 45 degrees azimuth and 20 
            %   %   degrees elevation. Plot the grating lobe diagram.
            %
            %   c = 3e8;f = 5e9;lambda = c/f;
            %   ha = phased.URA('ElementSpacing',[lambda*0.8 lambda*0.8]);
            %   plotGratingLobeDiagram(ha,f,[45;20])
            %
            %   % Example 2: 
            %   %   Construct a 0.5-wavelength by 0.9-wavelength element 
            %   %   spacing triangular lattice URA operating at 3 GHz and 
            %   %   steered to [30; 0] degrees. The propagation speed is 
            %   %   2.8e8 m/s and the phase shifter frequency is 3.1 GHz. 
            %   %   Plot the grating lobe diagram.
            %
            %   c = 2.8e8; f = 3e9; f0 = 3.1e9; lambda = c/f;
            %   ha = phased.URA('Lattice','triangular','ElementSpacing',...
            %       [lambda*0.5 lambda*0.9]);
            %   plotGratingLobeDiagram(ha,f,[30;0],c,f0);
            
            if ~isempty(coder.target)
                coder.internal.errorIf(true, ...
                    'phased:Waveform:CodegenNotSupported','plotGratingLobeDiagram');
            end            
            narginchk(2,5);
            % Ensure all input arguments are defined
            if (nargin < 5) 
                f0 = f;
            end
            if (nargin < 4) 
                c = physconst('LightSpeed');
            end
            if (nargin < 3) 
                ang = [0;0];
            end
            
            % Validate inputs
            sigdatatypes.validateFrequency(f,'plotGratingLobeDiagram',...
                'FREQ',{'scalar'});
            sigdatatypes.validateFrequency(f0,'plotGratingLobeDiagram',...
                'F0',{'scalar'});
            sigdatatypes.validateSpeed(c,'plotGratingLobeDiagram','C',{'scalar'});
            if isscalar(ang)
                ang = [ang;0];
            end
            sigdatatypes.validateAzElAngle(ang,'plotGratingLobeDiagram',...
                'ANGLE',{'vector'}); 
            az = ang(1);
            el = ang(2);
            
            lambda = c/f;
            u0 = cosd(el)*sind(az);
            u_center = f0/f*u0;
            v0 = sind(el);
            v_center = f0/f*v0;
            
            gl_u_spacing = lambda/obj.ElementSpacing(2);
            if strcmp(obj.Lattice,'Rectangular')
                gl_v_spacing = lambda/obj.ElementSpacing(1);
            else
                gl_v_spacing = lambda/(2*obj.ElementSpacing(1));
            end
            
            nu = 10;
            nl = -10;
            u_grat_grid = (nl:nu)*gl_u_spacing;%
            v_grat_grid = (nl:nu)*gl_v_spacing;
            
            cla reset;
            hold on;
            theta = 0:5:360;
            c_x = cosd(theta);
            c_y = sind(theta);
            
            % plot visible region patch
            hgreen = patch(c_x,c_y,'g','Tag','green_patch');
    
            % plot grating lobe region
            alphaVal = 0.8;
            for m = 1:(nu-nl+1)
                xref = u_grat_grid(m);
                for n = 1:(nu-nl+1)
                    yref = v_grat_grid(n);
                    if strcmp(obj.Lattice,'Rectangular')
                        if xref ~= 0 || yref ~= 0
                            tag_patch = sprintf('red_patch_%d_%d',m,n);
                            hred = patch(xref+c_x,yref+c_y,[1 0.54 0.54],...
                                'EdgeColor',[1 0.64 0.64],'FaceAlpha',alphaVal,'Tag',tag_patch);
                        end
                    else
                        if (xref ~= 0 || yref ~= 0) && ~rem(2*nl+m+n-2,2)
                            tag_patch = sprintf('red_patch_%d_%d',m,n);
                            hred = patch(xref+c_x,yref+c_y,[1 0.54 0.54],...
                                'EdgeColor',[1 0.64 0.64],'FaceAlpha',alphaVal,'Tag',tag_patch);
                        end
                    end
                end
            end

            % plot grating lobe point
            [u_plot,v_plot] = meshgrid(u_grat_grid,v_grat_grid);
            if ~strcmp(obj.Lattice,'Rectangular')
                u_plot(2:2:end) = nan;
                v_plot(2:2:end) = nan;
            end
            h = plot(u_center+u_plot(:),v_center+v_plot(:),'ko','MarkerSize',10,'Tag','GratingLobe_plot');
            haxes = get(h,'Parent');
            
            % plot main lobe point
            hmainlobe = plot(u_center,v_center,'ko', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','k','MarkerSize',10,'Tag','MainLobe_plot');
            
            % plot a unit circle edge
            plot(c_x,c_y,'k','Tag','Edge_plot');
                     
            % information of scan area
            if strcmp(obj.Lattice,'Rectangular')
                if (gl_u_spacing < 1)||(gl_v_spacing < 1)
                    arrayspan_str = sprintf('%s\n%s',...
                        getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                        getString(message('phased:system:array:AlwaysGratingLobe')));
                elseif (gl_u_spacing >= 2)&&(gl_v_spacing >= 2)
                    arrayspan_str = sprintf('%s\n%s',...
                        getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                        getString(message('phased:system:array:NoGratingLobe')));
                else
                    apu = min((gl_u_spacing-1),1);
                    apv = min((gl_v_spacing-1),1);
                    arrayspan_str = sprintf('%s\nU: [%.2f %.2f] (Az: [%.1f %.1f] deg)\nV: [%.2f %.2f] (El: [%.1f %.1f] deg)',...
                        getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                        -apu,apu, -asind(apu),asind(apu),-apv,apv, -asind(apv),asind(apv));
                end
            else
                % triangular lattice
                if ((gl_u_spacing >= 2)&&(gl_v_spacing >= 1))||...
                        ((gl_v_spacing >= 2)&&(gl_u_spacing >= 1))
                    arrayspan_str = sprintf('%s\n%s',...
                        getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                        getString(message('phased:system:array:NoGratingLobe')));
                elseif (gl_v_spacing < 0.5)||...
                       (gl_u_spacing < 0.5)||...
                       ((gl_u_spacing < 1)&&(gl_v_spacing < 1)&&((gl_u_spacing^2+gl_v_spacing^2)<1))
                    % (gl_u_spacing^2+gl_v_spacing^2) is the distance from the center
                    % to the grating lobe. when it is <=1 the four diagonal corner grating 
                    % lobes enter visible area (unit circle)
                    % The grating lobe located at the top and bottom enter are at a distance
                    % of 2*gl_v_spacing. When gl_v_spacing < 0.5 those lobes enter the visible 
                    % area. 
                    % Similarly when gl_u_spacing <0.5 the  left and right grating lobes enter 
                    % visible area.
                    arrayspan_str = sprintf('%s\n%s',...
                        getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                        getString(message('phased:system:array:AlwaysGratingLobe')));
                else 
                    if (gl_u_spacing < 1)&&(gl_v_spacing < 1)&&~(((gl_u_spacing^2+gl_v_spacing^2)<1))
                        apu = min((2*gl_u_spacing-1),(gl_u_spacing-sqrt(1-gl_v_spacing^2))); 
                        apv = min((2*gl_v_spacing-1),(gl_v_spacing-sqrt(1-gl_u_spacing^2)));
                        % Here we are trying to determine which edge is closer to the center of
                        % the visible area. (2*gl_u_spacing-1) calculates the distance from the
                        % center to the edge of the left (or right lobes).
                        % (gl_u_spacing-sqrt(1-gl_v_spacing^2)) calculates the distance from the
                        % center to the intersection of the upper left and lower right lobes.
                        % min((2*gl_u_spacing-1),(gl_u_spacing-sqrt(1-gl_v_spacing^2))) will
                        % give us the closer distance
                        % try [0.52 1.6] lambda spacing and then try [0.52 1.2] lambda spacing
                        % in sensorArrayAnalyzer to visualize those lobes.
                    elseif ((gl_u_spacing >= 2)&&(gl_v_spacing < 1)&&(gl_v_spacing >= 0.5))||...
                            ((gl_v_spacing >= 2)&&(gl_u_spacing < 1)&&(gl_u_spacing >= 0.5))||...
                            ((gl_u_spacing < 2)&&(gl_u_spacing >= 1)&&(gl_v_spacing < 2)&&(gl_v_spacing >= 1))
                        apu = min((2*gl_u_spacing-1),1);
                        apv = min((2*gl_v_spacing-1),1);  
                    elseif (gl_v_spacing < 1)&&(gl_u_spacing < 2)&&(gl_u_spacing >= 1)
                        apu = min((gl_u_spacing-sqrt(1-gl_v_spacing^2)),1);
                        apv = 2*gl_v_spacing-1;                
                    elseif (gl_u_spacing < 1)&&(gl_v_spacing < 2)&&(gl_v_spacing >= 1)
                        apu = 2*gl_u_spacing-1;
                        apv = min((gl_v_spacing-sqrt(1-gl_u_spacing^2)),1);  
                    end
                    arrayspan_str = sprintf('%s\nU: [%.2f %.2f] (Az: [%.1f %.1f] deg)\nV: [%.2f %.2f] (El: [%.1f %.1f] deg)',...
                        getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                        -apu,apu, -asind(apu),asind(apu),-apv,apv, -asind(apv),asind(apv));
                end
            end
            xlabel('U','Tag','ULabel');
            ylabel('V','Tag','VLabel')
            title(getString(message('phased:system:array:GratingLobeDiagram','U-V')),'Tag','GratingLobeTitle');
            set(haxes,'Color',get(gcf,'Color'));
            hleg = legend([hmainlobe,h,hgreen,hred],...
                {getString(message('phased:system:array:MainLobeLegend')),...
                getString(message('phased:system:array:GratingLobeLegend')),...
                getString(message('phased:system:array:GratingLobeFreeAreaLegend')),...
                getString(message('phased:system:array:GratingLobeAreaLegend'))},...
                'Location','NorthEastOutside','Tag','GratingLobeLegend','AutoUpdate','off');
            set(hleg,'FontSize',9);
            insetdist = get(haxes,'TightInset');
            text('Unit','Normalized','Position',[0.52 -3*insetdist(2)],'String',arrayspan_str,...
                'Color',[.501 .501 .501],'FontSize',9,'Tag','GratingLobeAnnotation');
            hold off;
            axis equal
            axis([-3 3 -3 3]);
            set(haxes,'OuterPosition',[0.1 0.25 0.8 0.7]);
            zoom reset;
            set(hleg,'Position',[0.2 0.1 0.28 0.15]);
            if nargout
                varargout{1} = haxes;
            end   
        end
    end
    
    methods(Hidden)
        function cl = clonecg(obj)
            cl = phased.URA(...
                    'Element',clonecg(obj.Element), ...
                    'Size',obj.Size, ...                
                    'ElementSpacing',obj.ElementSpacing, ...
                    'Lattice',obj.Lattice, ...
                    'ArrayNormal', obj.ArrayNormal,...
                    'Taper',obj.Taper);
        end
        
        function az_el = incidentAngleToAzEl(obj,inc_ang)
        %incidentAngleToAzEl Convert incident angle to local azimuth and
        %elevation angles
        %   AZ_EL = incidentAngleToAzEl(H,INC_ANG,N) converts the incident
        %   angle in the form of [azimuth; elevation] to the local
        %   [azimuth; elevation] pair, AZ_EL, for the Nth element in the
        %   uniform rectangular array H.

            
            % get array normal
            arNr = getArrayNormal(obj);
            
            % construct the local coordinate system of the element that
            % pointing to that normal direction.
            
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
    
    methods (Static,Hidden)
        function coordinateSystemInfo
            %coordinateSystemInfo Displays local coordinate system
            %information
            
            fprintf([...
            '    For a URA, the origin of the local coordinate system, i.e.,\n',...
            '    the phase center of the array, is defined at the center of\n',...
            '    the array. In local coordinate system, x axis is the default\n',...
            '    array normal of the URA, y axis defines the row direction of \n',...
            '    the URA and the z axis defines the column direction of the URA.\n',...
            '    Therefore, the elements are located on the yz plane.\n\n',...
            '    If the array normal is set to y axis, then the elements are\n',...
            '    located on the xz plane. Similarly, if the array normal is set\n',...
            '    to z axis, the elements are located on the xy plane.\n\n',...
            '    The azimuth angle is defined as the angle from x axis\n',...
            '    toward y axis. The elevation angle is defined as the angle\n',...
            '    from xy plane toward z axis.\n']);
        end
        
        function pos = calcElementPosition(EleIdx,array_sz,elementSpacing,lattice,arNr)
            NPerRow = array_sz(2);     % number of elements in each row
            NPerCol = array_sz(1);     % number of elements in each column
            deltaRows = (NPerRow-1)/2+1;   % delta along the row
            deltaCols = (NPerCol-1)/2+1;   % delta along the column
            [IdxInCol, IdxInRow] = ind2sub(array_sz,EleIdx);
            pos = zeros(3,numel(EleIdx));
            pos2 = (IdxInRow-deltaRows)*elementSpacing(2);
            
            if strcmp(lattice,'Triangular')
                evenrowidx = ~rem(IdxInCol,2);
                pos2(evenrowidx) = pos2(evenrowidx) + ...
                    elementSpacing(2)/2;
            end
            if strcmpi(arNr,'x')
                pos(2,:) = pos2;
                pos(3,:) = (deltaCols-IdxInCol)*elementSpacing(1);
               
            elseif strcmpi(arNr,'y')
               pos(1,:) = pos2;
               pos(3,:) = (deltaCols-IdxInCol)*elementSpacing(1); 
            else
               pos(1,:) = pos2;
               pos(2,:) = (deltaCols-IdxInCol)*elementSpacing(1);
                
            end
        end
    end
    
    methods(Static,Hidden,Access=protected)        
        function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractHomogeneousArray;
            props = {'Size',...
                     'ElementSpacing',...
                     'Lattice',...
                     'ArrayNormal',...
                     'Taper'};
        groups.PropertyList = [groups.PropertyList props];
        end
    end
    
end




