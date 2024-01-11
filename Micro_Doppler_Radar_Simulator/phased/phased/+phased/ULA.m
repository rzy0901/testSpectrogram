classdef (Sealed,StrictDefaults) ULA < phased.internal.AbstractHomogeneousArray
%ULA   Uniform linear array
%   H = phased.ULA creates a uniform linear array (ULA) System object, H.
%   This object models a ULA formed with identical sensor elements. The
%   default array is a 2-element ULA of two isotropic antenna elements.
%
%   H = phased.ULA(Name,Value) creates a ULA object, H, with the specified
%   property Name set to the specified Value. You can specify additional
%   name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   H = phased.ULA(N,D,Name,Value) creates a ULA object, H, with the
%   NumElements property set to N, ElementSpacing property set to D and
%   other specified property Names set to the specified Values. N and D are
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
%   ULA methods:
%
%   step                   - Output responses of the array elements (see 
%                            above)
%   release                - Allow property value and input characteristics
%                            changes
%   clone                  - Create a ULA object with same property values
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
%   ULA properties:
%
%   Element                - Element of the array
%   NumElements            - Number of elements
%   ElementSpacing         - Element spacing
%   ArrayAxis              - Array axis
%   Taper                  - Taper
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a 4-element ULA and plot its azimuth response. Assume the
%   %   operating frequency is 1 GHz and the wave propagation speed is 3e8 
%   %   m/s.
%
%   array = phased.ULA(4);
%   fc = 1e9; c = 3e8;
%   pattern(array,fc,-180:180,0,'CoordinateSystem','polar',...
%       'PropagationSpeed',c);
%
%   % Example 2:
%   %   Find the response of each element in the above array at the 
%   %   boresight.
%
%   array = phased.ULA(4);
%   fc = 1e9; ang = [0;0];
%   resp = array(fc,ang)
%
%   See also phased, phased.URA, phased.ConformalArray,
%   phased.HeterogeneousULA.

%   Copyright 2009-2016 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] Brookner, Radar Technology, Lex Book, 1996
%   [3] Harold Mott, Antennas for Radar and Communications, A Polarimetric
%       Approach, John Wiley & Sons, 1992


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)

        %NumElements   Number of elements
        %   An integer containing the number of elements in the array. The
        %   default value of this property is 2. 
        NumElements = 2
        %ElementSpacing   Element spacing (m)
        %   A scalar containing the spacing (in meters) between two
        %   adjacent elements in the array. The default value of this
        %   property is 0.5.
        ElementSpacing = 0.5
        %ArrayAxis Array axis
        %   Specify the axis on which the elements are placed as one of
        %   'x'|'y'|'z', where the default is 'y'.
        ArrayAxis = 'y';
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
        ArrayAxisSet = matlab.system.StringSet({'x','y','z'});
    end
    
    methods
        function set.NumElements(obj,value)
            sigdatatypes.validateIndex(value,...
                '','NumElements',...
                {'scalar','>=',2});
            obj.NumElements = value;
        end
        
        function set.ElementSpacing(obj,value)
            sigdatatypes.validateDistance(value,...
                '',...
                'ElementSpacing',{'scalar','positive','finite'});
            obj.ElementSpacing = value;
        end
        
        function set.Taper(obj,val)
            validateattributes(val,{'double'},...
                               {'nonnan','nonempty','finite','vector'},...
                               '','Taper');
            obj.Taper = val;
        end
        
    end
    
    methods

        function obj = ULA(varargin)
        %ULA   Construct the ULA class.
            obj@phased.internal.AbstractHomogeneousArray(nargin, ...
               varargin{:},'NumElements','ElementSpacing');
        end
    end
    methods(Hidden)
        function cl = clonecg(obj)
            cl = phased.ULA(...
                    'Element',clonecg(obj.Element), ...
                    'NumElements',obj.NumElements, ...                
                    'ElementSpacing',obj.ElementSpacing, ...
                    'ArrayAxis', obj.ArrayAxis,...
                    'Taper',obj.Taper);
        end
        
        function az_el = incidentAngleToAzEl(obj,inc_ang)
        %incidentAngleToAzEl Convert incident angle to local azimuth and
        %elevation angles
        %   AZ_EL = incidentAngleToAzEl(H,INC_ANG,N) converts the incident
        %   angle in the form of [azimuth; elevation] to the local
        %   [azimuth; elevation] pair, AZ_EL, for the Nth element in the
        %   uniform linear array H.

            % get array axis
            arAx = getArrayPlaneOrAxis(obj);    
        
            % construct the local coordinate system of the element that
            % pointing to that normal direction. 
            if (strcmpi(arAx,'x'))
                ElementAxes = [0 -1 0;1 0 0;0 0 1];
                % translate incoming direction to the local coordinate system
                az_el = phased.internal.incident2azel(inc_ang,ElementAxes);
            else
                az_el = inc_ang;
            end
            
        end
    end
    
    methods
        function N = getNumElements(obj)
        %getNumElements   Number of elements in the array
        %   N = getNumElements(H) returns the number of elements, N, in
        %   the ULA object H. 
        %
        %   % Example:
        %   %   Construct a default ULA and obtain its number of elements.
        %
        %   ha = phased.ULA;
        %   N = getNumElements(ha)
            N = obj.NumElements;
        end
               
        function pos = getElementPosition(obj,EleIdx)
        %getElementPosition Element positions of the array
        %   POS = getElementPosition(H) returns the element positions of
        %   the ULA H. POS is a 3xN matrix where N is the number of
        %   elements in H. Each column of POS defines the position, in the
        %   form of [x; y; z] (in meters), of an element in the local
        %   coordinate system.
        %
        %   For details regarding the local coordinate system of a ULA,
        %   type phased.ULA.coordinateSystemInfo.
        %
        %   POS = getElementPosition(H,ELEIDX) returns the positions of the
        %   elements that are specified in the element index vector ELEIDX.
        %
        %   % Example:
        %   %   Construct a default ULA and obtain its element positions.
        %
        %   ha = phased.ULA;
        %   pos = getElementPosition(ha)
        
            coder.extrinsic('phased.ULA.calcElementPosition');        
            phased.internal.narginchk(1,2,nargin);
            validateattributes(obj,{'phased.ULA'},{'scalar'},...
                'getElementPosition','H');
            N = getNumElements(obj);
            elSp = obj.ElementSpacing;
            arAx = obj.ArrayAxis;
            if nargin < 2
                EleIdx = 1:N;
            end
            sigdatatypes.validateIndex(EleIdx,...
                'getElementPosition','ELEIDX',{'vector','<=',N});
            pos = coder.internal.const(obj.calcElementPosition(N,elSp,EleIdx,arAx));

        end
        
        function nv = getElementNormal(obj, EleIdx)
        %getElementNormal Get normal directions of array elements
        %   NV = getElementNormal(H) returns the element normals of the
        %   ULA H. NORM is a 2xN matrix where N is the number of elements
        %   in H. each column of NORM specifies the normal direction of the
        %   corresponding element in the form of [azimuth; elevation] (in
        %   degrees) defined in the local coordinate system.
        %
        %   For details regarding the local coordinate system of a ULA,
        %   type phased.ULA.coordinateSystemInfo.
        %
        %   NV = getElementNormal(H,ELEIDX) returns the normals of the
        %   elements that are specified in the element index vector ELEIDX.
        %
        %   % Example:
        %   %   Construct a default ULA and obtain its element normals.
        %
        %   ha = phased.ULA;
        %   nv = getElementNormal(ha)
        
            phased.internal.narginchk(1,2,nargin);
            arAx = obj.ArrayAxis;
            N = getNumElements(obj);
            
            if nargin < 2
                EleIdx = 1:N;
            end
            sigdatatypes.validateIndex(EleIdx,...
                'getElementPosition','ELEIDX',{'vector','<=',N});
            
            % check if the element normal and array normal are aligned
            if strcmpi(arAx,'x') && ...
                    ~isElementFromAntenna(obj.Element) && ...
                    isElementNormalArrayNormalAligned(obj.Element)
                nv_temp = repmat([90;0],[1 N]);
            else
                nv_temp = zeros(2,N);
            end    
            nv = nv_temp(:,EleIdx);
        end
        
        function w = getTaper(obj)
        %getTaper Get array element tapers
        %   W = getTaper(H) returns the element taper of the ULA H. W is a
        %   length-N column vector where N is the number of elements in H.
        %
        %   % Example:
        %   %   Construct a default ULA and obtain its element taper.
        %
        %   ha = phased.ULA;
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
            %   %   Construct a ULA with elements spaced at 0.8 wavelength  
            %   %   at 5 GHz operating frequency. The array is steered to  
            %   %   45 degrees azimuth and 20 degrees elevation. Plot the 
            %   %   grating lobe diagram.
            %
            %   c = 3e8; f = 5e9; lambda = c/f;
            %   ha = phased.ULA('ElementSpacing',lambda*0.8);
            %   plotGratingLobeDiagram(ha,f,[45;20])
            %
            %   % Example 2: 
            %   %   Construct a ULA with elements spacing of 0.5 wavelength 
            %   %   at 3 GHz operating frequency. The array is steered to 
            %   %   30 degrees azimuth and 0 degrees elevation. The 
            %   %   propagation speed is 2.8e8 m/s and the phase shifter 
            %   %   frequency is 3.1 GHz. Plot the grating lobe diagram.  
            %
            %   c = 2.8e8; f = 3e9; f0 = 3.1e9; lambda = c/f;
            %   ha = phased.ULA('ElementSpacing',lambda*0.5 );
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
            u0 = sind(az2broadside(az,el));
            u_center = f0/f*u0;
            gl_spacing = lambda/obj.ElementSpacing;
            nu = 10;
            nl = -10;
            u_grat_grid = (nl:nu)*gl_spacing;  
            cla reset;
            hold on;
            hgreen = patch([-1 -1 1 1],[-0.5 0.5 0.5 -0.5],'g','Tag','green_patch');
            alphaVal = 0.8;
            
            for m = 1:(nu-nl+1)
                xref = u_grat_grid(m);
                if xref ~= 0
                    tag_patch = sprintf('red_patch_%d',m);
                    hred =  patch(xref+[-1 -1 1 1],[-0.5 0.5 0.5 -0.5],[1 0.54 0.54],...
                        'EdgeColor',[1 0.64 0.64],'FaceAlpha',alphaVal,'Tag',tag_patch);
                end
            end
            
            h = plot(u_grat_grid+u_center,zeros(size(u_grat_grid)),'ko','MarkerSize',10,'Tag','GratingLobe_plot');
            haxes = get(h,'Parent');
            %Turn off panning (g1159168)
            setAllowAxesPan(pan,haxes,false);
            hmainlobe = plot(u_center,0,'ko', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','k','MarkerSize',10,'Tag','MainLobe_plot');
            plot([-1 -1 1 1],[-0.5 0.5 0.5 -0.5],'k','Tag','Edge_plot');
            % information of scan area
            if gl_spacing >= 2
                arrayspan_str = sprintf('%s\n%s',...
                    getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                    getString(message('phased:system:array:NoGratingLobe')));
            elseif gl_spacing <= 1
                arrayspan_str = sprintf('%s\n%s',...
                    getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                    getString(message('phased:system:array:AlwaysGratingLobe')));
            else
                apu = gl_spacing-1;
                arrayspan_str = sprintf('%s\nU: [%.2f %.2f] (Az: [%.1f %.1f] deg)',...
                    getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
                    -apu,apu, -asind(apu),asind(apu));
            end
            title(getString(message('phased:system:array:GratingLobeDiagram','U')),'Tag','GratingLobeTitle');
            set(haxes,'Color',get(gcf,'Color'));
            hleg = legend([hmainlobe,h,hgreen,hred],...
                {getString(message('phased:system:array:MainLobeLegend')),...
                getString(message('phased:system:array:GratingLobeLegend')),...
                getString(message('phased:system:array:GratingLobeFreeAreaLegend')),...
                getString(message('phased:system:array:GratingLobeAreaLegend'))},...
                'Location','SouthWest','Tag','GratingLobeLegend','AutoUpdate','off');
            set(hleg,'FontSize',9);
            line([-3 3],[-1 -1],'Color','k','Tag','UAxis');
            for m = -3:3
                line([m m],[-1 -0.9],'Color','k');
                text(m-0.05-(m<0)*0.05,-1.2,num2str(m));
            end
            text(-0.05,-1.5,'U','Tag','ULabel');
            text(0,-2,arrayspan_str,'Color',[.501 .501 .501],'FontSize',9,'Tag','GratingLobeAnnotation');
            hold off;
            set(haxes,'OuterPosition',[0.1 0.1 0.8 0.8]);
            axis equal
            axis([-3 3 -2.5 1]);
            axis off
            zoom reset;
            set(hleg,'Position',[0.2 0.25 0.28 0.15]);
            if nargout
                varargout{1} = haxes;
            end
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
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function arAx = getArrayPlaneOrAxis(obj)        
            arAx = obj.ArrayAxis;
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
        
            % check if the element normal and array normal are aligned
            isAligned = isElementNormalArrayNormalAligned(obj.cElement); 
            
            % obtain the normal direction
            N = getNumElements(obj);
            if(isAligned)
                az_el_temp = incidentAngleToAzEl(obj,inc_ang);
            else
                az_el_temp = inc_ang;
            end
            az_el = repmat(az_el_temp,N,1);
            az_el = reshape(az_el,2,[]);

        end

        
    end
    
    methods (Static, Hidden)
        function coordinateSystemInfo
            %coordinateSystemInfo Displays local coordinate system
            %information
            
            fprintf([...
            '    For a ULA, the origin of the local coordinate system, i.e.,\n',...
            '    the phase center of the array, is defined at the center of\n',...
            '    the array. In local coordinate system, x axis is the default\n',...
            '    array normal of the ULA, y axis is the default array axis of \n',...
            '    the ULA and the z axis is the cross product of x and y axes. \n',...
            '    Therefore,the elements are located on the y axis.\n\n',...
            '    If the array axis is set to along x axis, then the array normal\n',...
            '    is along y axis. If the array axis is set to along z axis, then\n',...
            '    the array normal is along x axis.\n\n',...
            '    The azimuth angle is defined as the angle from x axis\n',...
            '    toward y axis. The elevation angle is defined as the angle\n',...
            '    from xy plane toward z axis.\n\n',...
            '    Many processing algorithms defined on a ULA use broadside \n',...
            '    angle instead of azimuth and elevation pair. The broadside\n',...
            '    angle is defined as the angle measured from the normal \n',...
            '    direction, which is the array normal direction projected \n',...
            '    onto the plane determined by the signal incident direction\n',...
            '    and the array axis, to the signal incident direction. The \n',...
            '    broadside angle is measured positive toward the positive \n',...
            '    direction of the array axis.\n']);
        end
        function pos = calcElementPosition(N,elSp,EleIdx,arAx)
            EleIdx = EleIdx(:);
            delta = (N-1)/2+1;
            numElIDX = numel(EleIdx);
            if strcmpi(arAx,'x')
                pos = [(EleIdx'-delta)*elSp;...
                       zeros(1,numElIDX);...
                       zeros(1,numElIDX)];
            elseif strcmpi(arAx,'y')
                pos = [zeros(1,numElIDX);...
                      (EleIdx'-delta)*elSp;...
                       zeros(1,numElIDX)];
            else
                pos = [zeros(1,numElIDX);...
                      zeros(1,numElIDX);...
                      (EleIdx'-delta)*elSp];
            end
        end
    end
    
    methods(Static,Hidden,Access=protected)        
        function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractHomogeneousArray;
        props = {'NumElements',...
                 'ElementSpacing',...
                 'ArrayAxis',...
                 'Taper'};
        groups.PropertyList = [groups.PropertyList props];
        end
    end
    
end

