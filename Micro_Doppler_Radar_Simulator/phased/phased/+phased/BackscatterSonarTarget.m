classdef (Sealed,StrictDefaults) BackscatterSonarTarget < phased.internal.AbstractSonarTarget
%BackscatterSonarTarget Backscatter sonar point target
%   H = phased.BackscatterSonarTarget creates a backscatter sonar point
%   target System object, H, that computes the reflected signal from a
%   target. A backscatter target is designed to be used in a monostatic
%   setting where the incident and reflected angles of the signal at the
%   target are the same.
%
%   H = phased.BackscatterSonarTarget(Name,Value) creates a sonar target
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Y = step(H,X,ANG) returns the reflected signal Y due to the incident
%   signal X when you set the Model property to 'Nonfluctuating'. In this
%   case, the values specified in the TSPattern property are used to
%   compute the target strength values at given directions.
%
%   X is a column vector or a matrix representing the incident signal. If X
%   is a matrix, then each column of X represents an independent signal.
%   When X represents multiple signals, the TSPattern property contains
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
%   Y = sqrt(G)*X where G = 10^(Ts/10).
%
%   G is the target gain factor corresponding to the target's sonar target
%   strength. Note that G is dimensionless.
%
%   Y = step(H,X,ANG,UPDATE) uses UPDATE as the indicator of whether to
%   update the TS values when you set the Model property to 'Swerling1',
%   'Swerling2', 'Swerling3' or 'Swerling4'. If UPDATE is true, a new TS
%   value is generated. If UPDATE is false, the previous TS value is used.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANG,UPDATE)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   BackscatterSonarTarget methods:
%
%   step     - Reflect the incoming signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create backscatter target object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of the sonar target object
%   
%   BackscatterSonarTarget properties:
%
%   AzimuthAngles       - Azimuth angles 
%   ElevationAngles     - Elevation angles 
%   TSPattern           - Target strength pattern
%   Model               - Fluctuation model
%   SeedSource          - Source of seed for random number generator
%   Seed                - Seed for random number generator
%
%   % Examples:
% 
%   % Example 1:
%   %   Calculate the reflected signal from a nonfluctuating point target
%   %   with target strength of -3 dB. Assume the incident angle is 10 
%   %   degrees azimuth.
% 
%   x = ones(10,1);
%   TSPattern = -3*ones(181,361);
%   target = phased.BackscatterSonarTarget('TSPattern',TSPattern,...
%       'Model','Nonfluctuating');
%   y = target(x,[10;0]);
% 
%   % Example 2:
%   %   Calculate the reflected signal from a Swerling 1 point target. The
%   %   target's target strength is -10 dB. The incident angle is 10 
%   %   degrees azimuth and 20 degrees elevation.
% 
%   ang = [10;20];
%   x = ones(10,1);
%   TSPattern = -10*ones(181,361);
%   target = phased.BackscatterSonarTarget('Model','Swerling1',...
%       'TSPattern',TSPattern);
%   y = target(x,ang,true);
%
%   See also phased.IsotropicHydrophone, phased.IsotropicProjector,
%   phased.IsoSpeedUnderwaterPaths, phased.MultipathChannel

%   Copyright 2016 The MathWorks, Inc.

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
        %TSPattern   Target strength pattern (dB)
        %   Specify the target strength (TS) pattern in dB as a QxP matrix
        %   or a QxPxL array, where Q is the number of elements presented
        %   in the ElevationAngles property, P is the number of elements
        %   presented in the AzimuthAngles property, and L is the number of
        %   patterns. Each page of the array defines a TS pattern across
        %   the region defined by the azimuth and elevation angles. The
        %   default value of this property is a 181x361 matrix with all
        %   elements equal to 0.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xP
        %   vector or an LxP matrix. In this case, each row is a separate
        %   pattern.
        %
        %   It is useful to specify multiple patterns in this property when
        %   you want to use one object to represent multiple targets in the
        %   simulation.
        %
        TSPattern = zeros(181,361)
    end
    
    properties (Access = private)
        %pTSPattern - private property which holds the pattern info
        pInitialTSPattern 
        pTSPattern 
    end
    
    properties (Access = private, Nontunable, Logical)
        pIsPatternSingleElevation
        pIsSinglePattern
    end
    
    methods
        function obj = BackscatterSonarTarget(varargin)
            obj@phased.internal.AbstractSonarTarget(varargin{:});
        end
    end

    methods
        
        function set.AzimuthAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.BackscatterSonarTarget','AzimuthAngles',...
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
                'phased.BackscatterSonarTarget','ElevationAngles',...
                {'vector','>=',-90,'<=',90});
            obj.ElevationAngles = value;
        end
        
        function set.TSPattern(obj,value)
            validateattributes(value,{'double'},...
                {'real','nonempty','3d'},...
                'phased.BackscatterSonarTarget','TSPattern');
            obj.TSPattern  = value;
        end
    end
    
    methods(Access = protected)
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractSonarTarget(obj);
            pattern_size = [numel(obj.ElevationAngles)...
            numel(obj.AzimuthAngles)];
            if numel(obj.ElevationAngles) == 1
                if numel(size(obj.TSPattern))==3
                    validateattributes(obj.TSPattern,{'double'},...
                        {'nrows',pattern_size(1),'ncols',pattern_size(2)},...
                        'phased.BackscatterSonarTarget','TSPattern');
                else
                    validateattributes(obj.TSPattern,{'double'},...
                        {'ncols',pattern_size(2)},...
                        'phased.BackscatterSonarTarget','TSPattern');
                end
            else
                validateattributes(obj.TSPattern,{'double'},...
                    {'nrows',pattern_size(1),'ncols',pattern_size(2)},...
                    'phased.BackscatterSonarTarget','TSPattern');
            end
        end

        function validateInputsImpl(obj,x,ang,tsupdate)
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
            
            if numel(obj.ElevationAngles)==1 
                if ismatrix(obj.TSPattern)
                    npat = size(obj.TSPattern,1);
                else
                    npat = size(obj.TSPattern,3);
                end
            else
                npat = size(obj.TSPattern,3);
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
                cond =   ~islogical(tsupdate) && ~isa(tsupdate,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                          'MATLAB:system:invalidInputDataType','Update','logical');
                end
                cond =   ~isscalar(tsupdate);
                if cond
                    coder.internal.errorIf(cond, ...
                          'MATLAB:system:inputMustBeScalar','Update');
                end
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSonarTarget(obj);
            if isLocked(obj)
                s.pTSPattern = obj.pTSPattern;
                s.pInitialTSPattern = obj.pInitialTSPattern;
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
        
        function setupImpl(obj,x,~,~,~)
            setupImpl@phased.internal.AbstractSonarTarget(obj,x);
            obj.pIsPatternSingleElevation = (numel(obj.ElevationAngles)==1);

            if obj.pIsPatternSingleElevation && ismatrix(obj.TSPattern)
                obj.pInitialTSPattern = permute(obj.TSPattern,[3 2 1]);
            else
                obj.pInitialTSPattern = obj.TSPattern;
            end
            obj.pIsSinglePattern = (size(obj.pInitialTSPattern,3)==1);
        end
        
        function resetImpl(obj)
            %setup random stream
            resetImpl@phased.internal.AbstractSonarTarget(obj);
                obj.pTSPattern = obj.pInitialTSPattern;
        end

        function y = stepImpl(obj,x,ang,updateTS)
            az = obj.AzimuthAngles;
            el = obj.ElevationAngles;
            isSingleEl = obj.pIsPatternSingleElevation;
          
            if obj.pFluctuate
                if obj.pNeedTSInit
                    obj.pTSPattern = 10*log10(obj.pFluctuateFunc(obj.cNoiseSource,10.^(obj.pInitialTSPattern/10)));
                    obj.pNeedTSInit = false;
                elseif updateTS
                    obj.pTSPattern = 10*log10(obj.pFluctuateFunc(obj.cNoiseSource,10.^(obj.pInitialTSPattern/10)));
                end
            end  

            pat_rcs = obj.pTSPattern;

            num_ang = size(ang,2);
            g = zeros(1,num_ang);
            for m = 1:num_ang
                if obj.pIsSinglePattern
                    tempTS = interpolatePattern(...
                        az,el,pat_rcs,ang(1,m), ang(2,m),isSingleEl);
                else
                    tempTS = interpolatePattern(...
                        az,el,pat_rcs(:,:,m),ang(1,m), ang(2,m),isSingleEl);
                end
                g(m) = tempTS;
            end

            g = sqrt(db2pow(g));

            %one rcs per channel
            y = bsxfun(@times, x, g);

        end

        function num = getNumInputsImpl(obj)
            num  = 2;
            if (obj.Model(1) == 'S');%Swerling
                num = num+1;
            end
        end
        
    end
    
    methods (Static, Hidden, Access = protected)        

        function groups = getPropertyGroupsImpl
            dTSPattern = matlab.system.display.internal.Property(...
                'TSPattern','Description', 'TS pattern (dB)', ...
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
                            'AzimuthAngles',...
                            'ElevationAngles',...
                            dTSPattern,...
                            'Model', ...
                            dSeedSource, ...
                            dSeed});   
        end
    end
    
    methods (Access = protected) %for Simulink
        function varargout = getInputNamesImpl(obj)   
            if (obj.Model(1) == 'S')
                varargout = {'X','Ang','Update'};
            else
                varargout = {'X','Ang'};
            end     
        end

        function varargout = getOutputNamesImpl(obj)  %#ok<MANU>
            varargout = {''};
        end
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Backscatter\n Sonar Target');
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
    
    methods (Hidden, Static)
        function flag = isAllowedInSystemBlock(obj) %#ok<INUSD>
            flag = false;
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

