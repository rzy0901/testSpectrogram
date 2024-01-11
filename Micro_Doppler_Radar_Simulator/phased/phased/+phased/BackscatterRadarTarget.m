classdef (Sealed,StrictDefaults) BackscatterRadarTarget < phased.internal.AbstractRadarTarget
%BackscatterRadarTarget Backscatter point radar target
%   H = phased.BackscatterRadarTarget creates a backscatter point target
%   System object, H, that computes the reflected signal from a target. A
%   backscatter target is designed to be used in a monostatic setting where
%   the incident and reflect angles of the signal at the target are the
%   same.
%
%   H = phased.BackscatterRadarTarget(Name,Value) creates a radar target
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax when the EnablePolarization property is false
%
%   Y = step(H,X,ANG) returns the reflected signal Y due to the incident
%   signal X when you set the Model property to 'Nonfluctuating'. In this
%   case, the values specified in the RCSPattern property are used to
%   compute the RCS values at given directions.
%
%   X is a column vector or a matrix representing the incident signal. If X
%   is a matrix, then each column of X represents an independent signal.
%   When X represents multiple signals, the RCSPattern property contains
%   either one pattern or M patterns where M is the number of columns in X.
%   If there is only one pattern, then the multiple signals are reflected
%   off the same target. If there are M patterns, then the multiple signals
%   are reflected off the corresponding pattern.
%
%   ANG is a 2-row matrix representing the signal's incident direction.
%   Each column of ANG specifies the incident direction of the
%   corresponding signal in the form of an [AzimuthAngle; ElevationAngle]
%   pair (in degrees). The number of columns in ANG must be the same as the
%   number of entries in X.
%
%   The reflected signal is calculated as
%
%   Y = sqrt(G)*X where G = 4*pi*RCS/(LAMBDA^2).
%
%   G is the target gain factor corresponding to the target's radar cross
%   section RCS and LAMBDA is the wavelength of the incoming signal. Note
%   that G is dimensionless.
%
%   Y = step(H,X,ANG,UPDATE) uses UPDATE as the indicator of whether to
%   update the RCS values when you set the Model property to 'Swerling1',
%   'Swerling2', 'Swerling3' or 'Swerling4'. If UPDATE is true, a new RCS
%   value is generated. If UPDATE is false, the previous RCS value is used.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANG,UPDATE)
%
%   Step method syntax when the EnablePolarization property is true
%
%   Y = step(H,X,ANG,LAXES) returns the reflected signal Y due to the
%   incident signal X when you set the Model property to 'Nonfluctuating'.
%   In this case, the values in the ShhPattern, SvvPattern, and ShvPattern
%   properties are used to compute the scattering matrices at given
%   directions.
%
%   X is a row struct array. Each struct represents a separate incoming
%   signal and contains three fields: X, Y, and Z. Each field contains a
%   column vector representing the X, Y, and Z component of the polarized
%   input signal, also measured in the global coordinate system,
%   respectively. The signals in X, Y, and Z fields must have the same
%   dimension.
%
%   When X represents multiple signals, the ShhPattern, SvvPattern, and
%   ShvPattern properties contain either one pattern or M patterns where M
%   is the number of signals in X. If there is only one pattern in each of
%   these properties, then the multiple signals are reflected off the same
%   target. If there are M patterns in each of these properties, then the
%   multiple signals are reflected off the corresponding pattern. Note that
%   the size of the patterns specified in each of these properties must be
%   the same.
%
%   ANG is a 2-row matrix representing the signal's incoming direction.
%   Each column of ANG specifies the incident direction of the
%   corresponding signal in the form of an [AzimuthAngle; ElevationAngle]
%   pair (in degrees). The number of columns in ANG must be the same as the
%   number of entries in X.
%
%   LAXES is a 3x3 matrix or a 3x3xM array. Each page of the array defines
%   a local coordinate system and the columns within each page specify the
%   local coordinate system's orthonormal x, y, and z axes, respectively.
%   Each axis is specified in [x;y;z] form measured in the global
%   coordinate system. If there is only a single signal in X, LAXES is a
%   3x3 matrix. If there are multiple signals in X, LAXES can either be a
%   3x3 matrix, in which case all targets share the same local coordinate
%   systems; or a 3x3xM matrix where each page defines the local coordinate
%   system for the corresponding target.
%
%   Y is also a row struct array. Each struct represents the reflected
%   signal of the corresponding incoming signal in X and contains three
%   fields: X, Y, and Z. Each field contains the X, Y, and Z component,
%   also measured in the global coordinate system, of the polarized
%   reflected signal, respectively.
%
%   The reflected signal is calculated as
%
%   Y_HV = sqrt(4*pi/(LAMBDA^2))*SMAT*X_HV 
%
%   where SMAT is the target's scattering matrix and LAMBDA is the
%   wavelength of the incoming signal. X_HV and Y_HV are X and Y in the
%   appropriate HV coordinates. The target's radar cross section can be
%   considered as the magnitude squared of the entries in the scattering
%   matrix.
%
%   Y = step(H,X,ANG,LAXES,UPDATE) uses UPDATE as the indicator of whether
%   to update the scattering matrix value when you set the Model property
%   to 'Swerling1', 'Swerling2', 'Swerling3' or 'Swerling4'. If UPDATE is
%   true, a new scattering matrix is generated. If UPDATE is false, the
%   previous scattering matrix is used.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANG,LAXES,UPDATE)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   BackscatterRadarTarget methods:
%
%   step     - Reflect the incoming signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create backscatter target object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of the radar target object
%   
%   BackscatterRadarTarget properties:
%
%   EnablePolarization  - Enable polarization
%   AzimuthAngles       - Azimuth angles 
%   ElevationAngles     - Elevation angles 
%   RCSPattern          - Radar cross section pattern 
%   ShhPattern          - Radar scattering pattern for HH 
%   SvvPattern          - Radar scattering pattern for VV 
%   ShvPattern          - Radar scattering pattern for HV 
%   Model               - Fluctuation model
%   PropagationSpeed    - Propagation speed
%   OperatingFrequency  - Operating frequency 
%   SeedSource          - Source of seed for random number generator
%   Seed                - Seed for random number generator
%
%   % Examples:
%   
%   % Example 1:
%   %   Calculate the reflected signal from a nonfluctuating point target
%   %   with RCS 10. Assume the incident angle is 10 degrees azimuth.
%
%   x = ones(10,1);
%   target = phased.BackscatterRadarTarget('Model','Nonfluctuating');
%   y = target(x,[10;0]);
%
%   % Example 2:
%   %   Calculate the reflected signal from a Swerling 1 point target. The
%   %   target's scattering matrix is given by the measured patterns. 
%   %   Assume the target is facing -x direction and the incident angle is
%   %   10 degrees azimuth and 20 degrees elevation respective to the 
%   %   target's orientation.
%
%   tgtaxes = azelaxes(180,0);
%   ang = [10;20];
%   shhpat = 6.3*ones(181,361);
%   shvpat = ones(181,361);
%   svvpat = 3*ones(181,361);
%   x = struct('X',ones(10,1),'Y',ones(10,1),'Z',ones(10,1));
%   target = phased.BackscatterRadarTarget('EnablePolarization',true,...
%       'Model','Swerling1','ShhPattern',shhpat,'ShvPattern',shvpat,...
%       'SvvPattern',svvpat);
%   y = target(x,[10;0],tgtaxes,true);
%
%   See also phased, phased.FreeSpace, phased.Platform,
%   phased.WidebandBackscatterRadarTarget.

%   Copyright 2015-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %AzimuthAngles   Azimuth angles (deg)
        %   Specify the azimuth angles (in degrees) as a length P vector.
        %   These are the azimuth angles where the custom pattern is
        %   evaluated. P must be greater than 2. The default value of this
        %   property is -180:180.
        AzimuthAngles = -180:180;
        %ElevationAngles   Elevation angles (deg)
        %   Specify the elevation angles (in degrees) as a length Q vector.
        %   These are the elevation angles where the custom pattern is
        %   evaluated. The default value of this property is -90:90.
        ElevationAngles = -90:90;
    end

    properties (Access = private)
        %pRCSPattern - private property which holds the pattern info
        pInitialRCSPattern 
        pInitialShhPattern 
        pInitialSvvPattern 
        pInitialShvPattern 
        pRCSPattern 
        pShhPattern 
        pSvvPattern 
        pShvPattern 
    end
    
    properties (Access = private, Nontunable, Logical)
        pIsPatternSingleElevation
        pIsSinglePattern
    end
    
    properties (Nontunable)
        %RCSPattern   Radar cross section pattern (square meters)
        %   Specify the radar cross section pattern (in square meters) as a
        %   QxP matrix or a QxPxL array, where Q is the number of elements
        %   presented in the ElevationAngles property, P is the number of
        %   elements presented in the AzimuthAngles property, and L is the
        %   number of patterns. Each page of the array defines an RCS
        %   pattern across the region defined by the azimuth and elevation
        %   angles. The default value of this property is a 181x361 matrix
        %   with all elements equal to 1.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xP
        %   vector or an LxP matrix. In this case, each row is a separate
        %   pattern.
        %
        %   It is useful to specify multiple pattern in this property when
        %   you want to use one object to represent multiple targets in the
        %   simulation.
        %
        %   This property applies when the EnablePolarization property is
        %   false.
        RCSPattern = ones(181,361)
        %ShhPattern   Radar scattering pattern for HH (meters)
        %   Specify the scattering pattern for horizontal transmission and
        %   horizontal reception (in meters) as a QxP matrix or a QxPxL
        %   array, where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and L is the number of patterns.
        %   Each page of the array defines a scattering pattern across the
        %   region defined by the azimuth and elevation angles. The default
        %   value of this property is a 181x361 matrix with all elements
        %   equal to 1.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xP
        %   vector or an LxP matrix. In this case, each row is a separate
        %   pattern.
        %
        %   It is useful to specify multiple pattern in this property when
        %   you want to use one object to represent multiple targets in the
        %   simulation.
        %
        %   This property applies when the EnablePolarization property is
        %   true.
        ShhPattern = ones(181,361)
        %ShvPattern   Radar scattering pattern for HV (meters)
        %   Specify the scattering pattern for vertical transmission and
        %   horizontal reception (in meters) as a QxP matrix or a QxPxL
        %   array, where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and L is the number of patterns.
        %   Each page of the array defines a scattering pattern across the
        %   region defined by the azimuth and elevation angles. The default
        %   value of this property is a 181x361 matrix with all elements
        %   equal to 0. Note that for a backscatterer, the pattern for
        %   horizontal transmission and vertical reception is the same as
        %   the pattern for vertical transmission and horizontal reception.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xP
        %   vector or an LxP matrix. In this case, each row is a separate
        %   pattern.
        %
        %   It is useful to specify multiple pattern in this property when
        %   you want to use one object to represent multiple targets in the
        %   simulation.
        %
        %   This property applies when the EnablePolarization property is
        %   true.
        ShvPattern = zeros(181,361)
        %SvvPattern   Radar scattering pattern for VV (meters)
        %   Specify the scattering pattern for vertical transmission and
        %   vertical reception (in meters) as a QxP matrix or a QxPxL
        %   array, where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and L is the number of patterns.
        %   Each page of the array defines a scattering pattern across the
        %   region defined by the azimuth and elevation angles. The default
        %   value of this property is a 181x361 matrix with all elements
        %   equal to 1.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xP
        %   vector or an LxP matrix. In this case, each row is a separate
        %   pattern.
        %
        %   It is useful to specify multiple pattern in this property when
        %   you want to use one object to represent multiple targets in the
        %   simulation.
        %
        %   This property applies when the EnablePolarization property is
        %   true.
        SvvPattern = ones(181,361)
    end
    
    methods
        function obj = BackscatterRadarTarget(varargin)
            obj@phased.internal.AbstractRadarTarget(varargin{:});
        end
    end

    methods
        
        function set.AzimuthAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.BackscatterRadarTarget','AzimuthAngles',...
                {'vector','>=',-180,'<=',180});
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'AzimuthAngles');
            end
            obj.AzimuthAngles = value;
        end

        function set.ElevationAngles(obj,value)
            % allow single elevation cut, i.e., azimuth only pattern
            sigdatatypes.validateAngle(value,...
                'phased.BackscatterRadarTarget','ElevationAngles',...
                {'vector','>=',-90,'<=',90});
            obj.ElevationAngles = value;
        end
        
        function set.RCSPattern(obj,value)
            validateattributes(value,{'double'},...
                {'real','nonnegative','nonempty','3d'},...
                'phased.BackscatterRadarTarget','RCSPattern');

            obj.RCSPattern  = value;
        end
        
        function set.ShhPattern(obj,value)
            validateattributes(value,{'double'},...
                {'finite','nonempty','3d'},...
                'phased.BackscatterRadarTarget','ShhPattern');

            obj.ShhPattern  = value;
        end
        
        function set.SvvPattern(obj,value)
            validateattributes(value,{'double'},...
                {'finite','nonempty','3d'},...
                'phased.BackscatterRadarTarget','SvvPattern');

            obj.SvvPattern  = value;
        end
        
        function set.ShvPattern(obj,value)
            validateattributes(value,{'double'},...
                {'finite','nonempty','3d'},...
                'phased.BackscatterRadarTarget','ShvPattern');

            obj.ShvPattern  = value;
        end
    end
    
    methods(Access = protected)
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractRadarTarget(obj);
            pattern_size = [numel(obj.ElevationAngles)...
                numel(obj.AzimuthAngles)];
            if obj.EnablePolarization
                if numel(obj.ElevationAngles) == 1
                    if numel(size(obj.ShhPattern))==3
                        validateattributes(obj.ShhPattern,{'double'},...
                            {'nrows',pattern_size(1),'ncols',pattern_size(2)},...
                            'phased.BackscatterRadarTarget','ShhPattern');
                    else
                        validateattributes(obj.ShhPattern,{'double'},...
                            {'ncols',pattern_size(2)},...
                            'phased.BackscatterRadarTarget','ShhPattern');
                    end
                else
                    validateattributes(obj.ShhPattern,{'double'},...
                        {'nrows',pattern_size(1),'ncols',pattern_size(2)},...
                        'phased.BackscatterRadarTarget','ShhPattern');
                end
                validateattributes(obj.SvvPattern,{'double'},...
                    {'size',size(obj.ShhPattern)},...
                    'phased.BackscatterRadarTarget','SvvPattern');
                validateattributes(obj.ShvPattern,{'double'},...
                    {'size',size(obj.ShhPattern)},...
                    'phased.BackscatterRadarTarget','ShvPattern');
            else
                if numel(obj.ElevationAngles) == 1
                    if numel(size(obj.RCSPattern))==3
                        validateattributes(obj.RCSPattern,{'double'},...
                            {'nrows',pattern_size(1),'ncols',pattern_size(2)},...
                            'phased.BackscatterRadarTarget','RCSPattern');
                    else
                        validateattributes(obj.RCSPattern,{'double'},...
                            {'ncols',pattern_size(2)},...
                            'phased.BackscatterRadarTarget','RCSPattern');
                    end
                else
                    validateattributes(obj.RCSPattern,{'double'},...
                        {'nrows',pattern_size(1),'ncols',pattern_size(2)},...
                        'phased.BackscatterRadarTarget','RCSPattern');
                end
            end
        end

        function validateInputsImpl(obj,x,ang,lclaxes,rcsupdate)
            if obj.EnablePolarization
                cond = ~isa(x,'struct');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','X','struct');
                end
                flag_hasXYZ = isfield(x(1),'X') && isfield(x(1),'Y') && isfield(x(1),'Z');
                cond = ~flag_hasXYZ;
                if cond
                    coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationXYZStruct');
                end
                
                xsize = size(x);
                for m = 1:xsize(2)
                    x_x = x(m).X;
                    x_y = x(m).Y;
                    x_z = x(m).Z;
                    cond =  ~isa(x_x,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType',sprintf('X(%d).X',m),'double');
                    end
                    cond =  ~iscolumn(x_x) || isempty(x_x);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector',sprintf('X(%d).X',m));
                    end
                    cond =  ~isa(x_y,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType',sprintf('X(%d).Y',m),'double');
                    end
                    cond =  ~iscolumn(x_y) || isempty(x_y);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Y',m));
                    end
                    cond =  ~isa(x_z,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType',sprintf('X(%d).Z',m),'double');
                    end
                    cond =  ~iscolumn(x_z) || isempty(x_z);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Z',m));
                    end
                    cond = numel(x_x)~=numel(x_y) || numel(x_x)~=numel(x_z);
                    if cond
                        coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                            'X,Y,Z',sprintf('X(%d)',m));
                    end
                end
                
                if xsize(2)==1 || ismatrix(lclaxes)
                    sigdatatypes.validate3DCartCoord(lclaxes,'','LAxes',...
                        {'size',[3 3]});
                else
                    validateattributes(lclaxes,{'double'},{'finite','nonnan','nonempty','real',...
                        'size',[3 3 xsize(2)]},'','LAxes');
                end
                
            else
                xsize = size(x);
                cond =  ~isa(x,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','X','double');
                end
                cond =  ~ismatrix(x) || isempty(x);
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:inputMustBeMatrix','X');
                end
            end
            
            if obj.EnablePolarization
                if numel(obj.ElevationAngles)==1 
                    if ismatrix(obj.ShhPattern)
                        npat = size(obj.ShhPattern,1);
                    else
                        npat = size(obj.ShhPattern,3);
                    end
                else
                    npat = size(obj.ShhPattern,3);
                end
            else
                if numel(obj.ElevationAngles)==1 
                    if ismatrix(obj.RCSPattern)
                        npat = size(obj.RCSPattern,1);
                    else
                        npat = size(obj.RCSPattern,3);
                    end
                else
                    npat = size(obj.RCSPattern,3);
                end
            end
            cond = (npat ~= 1) && (xsize(2) ~= npat);
            if cond
                coder.internal.errorIf(cond,'phased:target:signalPatternMismatch');
            end
                
            validateNumChannels(obj,x);
            
            angsize = size(ang);
            cond = xsize(2) ~= angsize(2);
            if cond
                coder.internal.errorIf(cond,'phased:phased:collector:AngleDimensionMismatch');
            end
            sigdatatypes.validateAzElAngle(ang,'','Ang');

            fluctuateFlag = (obj.Model(1) == 'S');%Swerling
            

            if fluctuateFlag
                if ~obj.EnablePolarization
                    rcsupdate = lclaxes;
                end
                cond =   ~islogical(rcsupdate) && ~isa(rcsupdate,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                          'MATLAB:system:invalidInputDataType','Update','logical');
                end
                cond =   ~isscalar(rcsupdate);
                if cond
                    coder.internal.errorIf(cond, ...
                          'MATLAB:system:inputMustBeScalar','Update');
                end
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractRadarTarget(obj);
            if isLocked(obj)
                s.pLambda = obj.pLambda;
                s.pRCSPattern = obj.pRCSPattern;
                s.pShhPattern = obj.pShhPattern;
                s.pSvvPattern = obj.pSvvPattern;
                s.pShvPattern = obj.pShvPattern;
                s.pInitialRCSPattern = obj.pInitialRCSPattern;
                s.pInitialShhPattern = obj.pInitialShhPattern;
                s.pInitialSvvPattern = obj.pInitialSvvPattern;
                s.pInitialShvPattern = obj.pInitialShvPattern;
                s.pIsPatternSingleElevation = obj.pIsPatternSingleElevation;
                s.pIsSinglePattern = obj.pIsSinglePattern;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractRadarTarget(obj,prop);
            if obj.EnablePolarization
                if strcmp(prop,'RCSPattern') 
                    flag = true;
                end
            else
                if strcmp(prop,'ShhPattern') || ...
                        strcmp(prop,'SvvPattern') || ...
                        strcmp(prop,'ShvPattern')
                    flag = true;
                end
            end
        end
        
        function setupImpl(obj,x,~,~,~)
            setupImpl@phased.internal.AbstractRadarTarget(obj,x);
            obj.pLambda = obj.PropagationSpeed/obj.OperatingFrequency;
            obj.pIsPatternSingleElevation = (numel(obj.ElevationAngles)==1);
            if obj.EnablePolarization
                if obj.pIsPatternSingleElevation && ismatrix(obj.ShhPattern)
                    obj.pInitialShhPattern = permute(obj.ShhPattern,[3 2 1]);
                    obj.pInitialSvvPattern = permute(obj.SvvPattern,[3 2 1]);
                    obj.pInitialShvPattern = permute(obj.ShvPattern,[3 2 1]);
                else
                    obj.pInitialShhPattern = obj.ShhPattern;
                    obj.pInitialSvvPattern = obj.SvvPattern;
                    obj.pInitialShvPattern = obj.ShvPattern;
                end
                obj.pIsSinglePattern = (size(obj.pInitialShhPattern,3)==1);
            else
                if obj.pIsPatternSingleElevation && ismatrix(obj.RCSPattern)
                    obj.pInitialRCSPattern = permute(obj.RCSPattern,[3 2 1]);
                else
                    obj.pInitialRCSPattern = obj.RCSPattern;
                end
                obj.pIsSinglePattern = (size(obj.pInitialRCSPattern,3)==1);
            end
        end
        
        function resetImpl(obj)
            %setup random stream
            resetImpl@phased.internal.AbstractRadarTarget(obj);
            if obj.EnablePolarization
                obj.pShhPattern = obj.pInitialShhPattern;
                obj.pShvPattern = obj.pInitialShvPattern;
                obj.pSvvPattern = obj.pInitialSvvPattern;
            else
                obj.pRCSPattern = obj.pInitialRCSPattern;
            end

        end

        function y = stepImpl(obj,x,ang,lclaxesArg,updateRCSArg)
            az = obj.AzimuthAngles;
            el = obj.ElevationAngles;
            isSingleEl = obj.pIsPatternSingleElevation;
            if obj.EnablePolarization
                lclaxes_in = lclaxesArg;

                if obj.pFluctuate
                    updateRCS = updateRCSArg;
                    % update the entire pattern at each trigger because the
                    % incident angle is unknown
                    if obj.pNeedRCSInit
                        obj.pShhPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShhPattern,false);
                        obj.pSvvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialSvvPattern,false);
                        obj.pShvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShvPattern,false);
                        obj.pNeedRCSInit = false;
                    elseif updateRCS
                        obj.pShhPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShhPattern,false);
                        obj.pSvvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialSvvPattern,false);
                        obj.pShvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShvPattern,false);
                    end
                end
                
                pat_hh = obj.pShhPattern;
                pat_vv = obj.pSvvPattern;
                pat_hv = obj.pShvPattern;
                
                % x should be all same size
                % ang should match x column numbers
                num_ang = size(ang,2);
                M = size(x,2);
                if ismatrix(lclaxes_in) && M > 1
                    lclaxes = repmat(lclaxes_in,1,1,M);
                else
                    lclaxes = lclaxes_in;
                end
                x_x = complex(zeros(size(x(1).X,1),M));
                x_y = x_x; x_z = x_x;
                for m = 1:M
                    x_x(:,m) = x(m).X;
                    x_y(:,m) = x(m).Y;
                    x_z(:,m) = x(m).Z;
                end
                if obj.pIsCodeGen
                    initYfields = complex(zeros(size(x(1).X,1),1));
                    y = repmat(struct('X',initYfields, 'Y', initYfields,'Z', initYfields),1, M);
                else
                    %keep extra fields in MATLAB
                    y = x;
                end
                for m = num_ang:-1:1
                    
                    % rcs could be complex, use direct computation.
                
                    if obj.pIsSinglePattern
                        % convert incoming signal to local HV
                        lcl_x = phased.internal.global2localvec(...
                            [x_x(:,m) x_y(:,m) x_z(:,m)].',lclaxes);
                        lcl_x = cart2sphvec(lcl_x,ang(1,m),ang(2,m));
                        lcl_x_sph = lcl_x(1:2,:);
                        
                        tempShh = interpolatePattern(...
                            az,el,pat_hh,ang(1,m), ang(2,m),isSingleEl);
                        tempSvv = interpolatePattern(...
                            az,el,pat_vv,ang(1,m), ang(2,m),isSingleEl);
                        tempShv = interpolatePattern(...
                            az,el,pat_hv,ang(1,m), ang(2,m),isSingleEl);
                    else
                        % convert incoming signal to local HV
                        lcl_x = phased.internal.global2localvec(...
                            [x_x(:,m) x_y(:,m) x_z(:,m)].',lclaxes(:,:,m));
                        lcl_x = cart2sphvec(lcl_x,ang(1,m),ang(2,m));
                        lcl_x_sph = lcl_x(1:2,:);
                        
                        tempShh = interpolatePattern(...
                            az,el,pat_hh(:,:,m),ang(1,m), ang(2,m),isSingleEl);
                        tempSvv = interpolatePattern(...
                            az,el,pat_vv(:,:,m),ang(1,m), ang(2,m),isSingleEl);
                        tempShv = interpolatePattern(...
                            az,el,pat_hv(:,:,m),ang(1,m), ang(2,m),isSingleEl);
                    end
                    
                    % for backscatter, Shv = Svh
                    tempSmat = [tempShh tempShv;tempShv tempSvv]; 
                    
                    g = sqrt(4*pi/obj.pLambda^2)*tempSmat;
                    
                    lcl_y_sph = g*lcl_x_sph;
                    
                    lcl_y = sph2cartvec(...
                        [lcl_y_sph;zeros(1,size(lcl_y_sph,2))],...
                        ang(1,m),ang(2,m));
                    
                    if obj.pIsSinglePattern
                        g_y = phased.internal.local2globalvec(...
                            lcl_y,lclaxes);
                    else
                        g_y = phased.internal.local2globalvec(...
                            lcl_y,lclaxes(:,:,m));
                    end
                    y(m).X = g_y(1,:).';
                    y(m).Y = g_y(2,:).';
                    y(m).Z = g_y(3,:).';
                end
                
            else  % no polarization
                if obj.pFluctuate
                    updateRCS = lclaxesArg;  % 3rd input is updateRCS
                end

                if obj.pFluctuate
                    if obj.pNeedRCSInit
                        obj.pRCSPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialRCSPattern);
                        obj.pNeedRCSInit = false;
                    elseif updateRCS
                        obj.pRCSPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialRCSPattern);
                    end
                end  
                
                pat_rcs = obj.pRCSPattern;
                
                num_ang = size(ang,2);
                g = zeros(1,num_ang);
                for m = 1:num_ang
                    if obj.pIsSinglePattern
                        tempRCS = interpolatePattern(...
                            az,el,pat_rcs,ang(1,m), ang(2,m),isSingleEl);
                    else
                        tempRCS = interpolatePattern(...
                            az,el,pat_rcs(:,:,m),ang(1,m), ang(2,m),isSingleEl);
                    end
                    g(m) = tempRCS;
                end
                
                g = sqrt(db2pow(aperture2gain(g,obj.pLambda)));
                
                %one rcs per channel
                y = bsxfun(@times, x, g);
            end

        end

        function num = getNumInputsImpl(obj)
            num  = 2;
            if (obj.Model(1) == 'S');%Swerling
                num = num+1;
            end
            if obj.EnablePolarization
                num = num+1;
            end
        end
        
    end
    
    methods (Static, Hidden, Access = protected)        
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:BackScatterRadarTargetTitle')),...
              'Text',getString(message('phased:library:block:BackScatterRadarTargetDesc')));
        end
        
        function groups = getPropertyGroupsImpl
            dEnablePolarization = matlab.system.display.internal.Property(...
                'EnablePolarization', 'IsGraphical', false);
            dRCSPattern = matlab.system.display.internal.Property(...
                'RCSPattern','Description', 'RCS pattern (m^2)', ...
                'UseClassDefault',false,'Default', 'ones(181,361)');
            dShhPattern = matlab.system.display.internal.Property(...
                'ShhPattern','Description', 'Shh pattern (m)', ...
                'UseClassDefault',false,'Default', 'ones(181,361)');
            dSvvPattern = matlab.system.display.internal.Property(...
                'SvvPattern','Description', 'Svv pattern (m)', ...
                'UseClassDefault',false,'Default', 'ones(181,361)');
            dShvPattern = matlab.system.display.internal.Property(...
                'ShvPattern','Description', 'Shv pattern (m)', ...
                'UseClassDefault',false,'Default', 'ones(181,361)');

            dSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false,'UseClassDefault',false,...
                'Default','Property');
            dSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false,'UseClassDefault',false,...
                'Default','randi(65535,1)');

            groups = matlab.system.display.Section(...
                            'PropertyList',...
                            {...
                            dEnablePolarization,...
                            'AzimuthAngles',...
                            'ElevationAngles',...
                            dRCSPattern,...
                            dShhPattern,...
                            dShvPattern,...
                            dSvvPattern,...
                            'Model', ...
                            'PropagationSpeed', ...
                            'OperatingFrequency', ...
                            dSeedSource, ...
                            dSeed});   
                       
        end
    end
    
    methods (Access = protected) %for Simulink
        function varargout = getInputNamesImpl(obj)   
            if obj.EnablePolarization
                if (obj.Model(1) == 'S') %Swerling
                    varargout = {'X','Ang','LAxes','Update'};
                else
                    varargout = {'X','Ang','LAxes'};
                end
            else
                if (obj.Model(1) == 'S')
                    varargout = {'X','Ang','Update'};
                else
                    varargout = {'X','Ang'};
                end     
            end
        end

        function varargout = getOutputNamesImpl(obj)  %#ok<MANU>
            varargout = {''};
        end
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Backscatter\n Target');
        end
        function varargout = getOutputSizeImpl(obj)
           varargout{1} = propagatedInputSize(obj,1);
        end
        function varargout = getOutputDataTypeImpl(obj)
            varargout{1} = propagatedInputDataType(obj,1);
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout{1} = true;
        end    
    end
end

function pat = interpolatePattern(az,el,pattern,az_q,el_q,interpIn1D)
    if interpIn1D
        pat = interpolatePatternRadians1D(phased.internal.deg2rad(az), ...
            pattern, phased.internal.deg2rad(az_q(:)));
    else
        pat = interpolatePatternRadians2D(phased.internal.deg2rad(az), phased.internal.deg2rad(el), ...
            pattern, phased.internal.deg2rad(az_q(:)), phased.internal.deg2rad(el_q(:)));
    end
end

function pat = interpolatePatternRadians2D(azr, elr, pattern, az_qr, el_qr)
 pat = interp2(azr,elr,pattern,...
        az_qr, el_qr,'nearest',0);
end

function pat = interpolatePatternRadians1D(azr, pattern, az_qr)
 pat = interp1(azr,pattern,...
        az_qr,'nearest',0);
end

