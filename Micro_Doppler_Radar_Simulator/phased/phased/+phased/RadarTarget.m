classdef (Sealed,StrictDefaults) RadarTarget < phased.internal.AbstractRadarTarget
%RadarTarget Radar target
%   H = phased.RadarTarget creates a radar target System object, H, that
%   computes the reflected signal from a target.
%
%   H = phased.RadarTarget(Name,Value) creates a radar target object, H,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax when the EnablePolarization property is false
%
%   Y = step(H,X) returns the reflected signal Y due to the incident signal
%   X when you set the Model property to 'Nonfluctuating'. In this case,
%   the value of MeanRCS property is used as the RCS value.
%
%   The reflected signal is calculated as
%
%   Y = sqrt(G)*X where G = 4*pi*RCS/(LAMBDA^2).
%
%   G is the target gain factor corresponding to the target's radar cross
%   section RCS and LAMBDA is the wavelength of the incoming signal. Note
%   that G is dimensionless.
%
%   Y = step(H,X,RCS) uses RCS as the mean radar cross-section values, when
%   you set the MeanRCSSource property to 'Input port'. RCS must be a 1xN
%   vector of positive values. N is the number of targets to model.
%   
%   Y = step(H,X,UPDATE) uses UPDATE as the indicator of whether to update
%   the RCS values when you set the Model property to 'Swerling1',
%   'Swerling2', 'Swerling3' or 'Swerling4'. If UPDATE is true, a new RCS
%   value is generated. If UPDATE is false, the previous RCS value is used.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,RCS,UPDATE)
%
%   Step method syntax when the EnablePolarization property is true
%
%   Y = step(H,X,ANGIN,LAXES) returns the reflected signal Y due to the
%   incident signal X when you set the Model property to 'Nonfluctuating'
%   and the Mode property to 'Monostatic'. In this case, the value of
%   ScatteringMatrix property is used as the scattering matrix value.
%
%   X is a row struct array. Each struct represents a separate incoming
%   signal and contains three fields: X, Y, and Z. Each field contains a
%   column vector representing the X, Y, and Z component of the polarized
%   input signal, also measured in the global coordinate system,
%   respectively. The signals in X, Y, and Z fields must have the same
%   dimension.
%
%   ANGIN is a 2-row matrix representing the signal's incoming direction.
%   Each column of ANGIN specifies the incident direction of the
%   corresponding signal in the form of an [AzimuthAngle; ElevationAngle]
%   pair (in degrees). The number of columns in ANGIN must be the same as
%   the number of entries in X.
%
%   LAXES is a 3x3 matrix whose columns specify the local coordinate
%   system's orthonormal x, y, and z axes, respectively. Each axis is
%   specified in [x;y;z] form measured in the global coordinate system.
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
%   considered as the magnitude square of the entries in the scattering
%   matrix.
%
%   Y = step(H,X,ANGIN,ANGOUT,LAXES) specifies the reflection angle, ANGOUT
%   (in degrees), when you set the Mode property to 'Bistatic'. ANGOUT is a
%   2-row matrix representing the signal's reflecting direction. Each
%   column of ANGOUT specifies the incident direction of the corresponding
%   signal in the form of an [AzimuthAngle; ElevationAngle] pair (in
%   degrees). The number of columns in ANGOUT must be the same as the
%   number of entries in X.
%
%   Y = step(H,X,ANGIN,LAXES,SMAT) uses SMAT as the scattering matrix,
%   when you set the ScatteringMatrixSource property to 'Input port'. SMAT
%   must be a 2x2 matrix.
%   
%   Y = step(H,X,ANGIN,LAXES,UPDATE) uses UPDATE as the indicator of
%   whether to update the scattering matrix value when you set the Model
%   property to 'Swerling1', 'Swerling2', 'Swerling3' or 'Swerling4'. If
%   UPDATE is true, a new scattering matrix is generated. If UPDATE is
%   false, the previous scattering matrix is used.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANGIN,ANGOUT,LAXES,SMAT,UPDATE)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   RadarTarget methods:
%
%   step     - Reflect the incoming signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create radar target object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of the radar target object
%   
%   RadarTarget properties:
%
%   EnablePolarization       - Enable polarization
%   Mode                     - Scattering mode
%   MeanRCSSource            - Source of mean radar cross section
%   ScatteringMatrixSource   - Source of mean scattering matrix
%   ScatteringMatrix         - Mean scattering matrix
%   MeanRCS                  - Mean radar cross section 
%   Model                    - Fluctuation model
%   PropagationSpeed         - Propagation speed
%   OperatingFrequency       - Operating frequency 
%   SeedSource               - Source of seed for random number generator
%   Seed                     - Seed for random number generator
%
%   % Example:
%   %   Calculate the reflected signal from a nonfluctuating point target
%   %   with RCS 10.
%
%   x = ones(10,1);
%   target = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',10);
%   y = target(x);
%
%   See also phased, phased.BackscatterRadarTarget, phased.FreeSpace,
%   phased.Platform.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [2] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001
%   [3] Harold Mott, Antenna for Radar and Communications: A polarimetric
%       Approach, Wiley-Interscience, 1992
%   [4] Eugene Knott, et al. Radar Cross Section, Scitech Publishing, 2004


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Mode   Scattering mode
        %   Specify the target scattering mode as one of 'Monostatic' |
        %   'Bistatic', where the default is 'Monostatic'. If you set this
        %   property to 'Monostatic', the signal's reflection direction is
        %   the same as its incoming direction. If you set this property to
        %   'Bistatic', the signal's reflection direction may be different
        %   to its incoming direction. This property applies when you set
        %   the EnablePolarization property to true.
        Mode = 'Monostatic';
    end

    properties (Nontunable)
        %MeanRCSSource Source of mean radar cross section
        %   Specify how to determine the target's mean RCS value as one of
        %   'Property' | 'Input port', where the default is 'Property'. If
        %   you set this property to 'Property', the mean RCS value is
        %   determined by the value of the MeanRCS property. If you set
        %   this property to 'Input port', the mean RCS value is determined
        %   by the input argument. This property applies when you set the
        %   EnablePolarization property to false
        MeanRCSSource = 'Property'
        %ScatteringMatrixSource Source of mean scattering matrix
        %   Specify how to determine the target's mean scattering matrix as
        %   one of 'Property' | 'Input port', where the default is
        %   'Property'. If you set this property to 'Property', the mean
        %   scattering matrix is determined by the value of the
        %   ScatteringMatrix property. If you set this property to 'Input
        %   port', the mean RCS value is determined by the input argument.
        %   This property applies when you set the EnablePolarization
        %   property to true.
        ScatteringMatrixSource = 'Property'
    end

    properties
        %MeanRCS    Mean radar cross section (m^2)
        %   Specify one or several mean values of the targets' radar cross
        %   section (in square meters) as a 1xN vector of nonnegative
        %   values. N is the number of targets to model. This property
        %   applies when you set the MeanRCSSource property to 'Property'
        %   and the EnablePolarization property to false. The default value
        %   of this property is one. This property is tunable.
        MeanRCS = 1
        %ScatteringMatrix  Mean scattering matrix (m)
        %   Specify the mean value of the target's radar cross section (in
        %   square meters) matrix as a 2x2 matrix. The matrix is in the
        %   form of [s_hh s_hv;s_vh s_vv], where s_hv specifies the
        %   scattering when the input signal is vertically polarized and
        %   the reflected signal is horizontally polarized. This property
        %   applies when you set the ScatteringMatrixSource property to
        %   'Property' and the EnablePolarization property to true. The
        %   default value of this property is [1 0;0 1]. This property is
        %   tunable.
        ScatteringMatrix = eye(2)
    end

    properties(Constant, Hidden)
        MeanRCSSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        ScatteringMatrixSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        ModeSet = matlab.system.StringSet({'Monostatic','Bistatic'});
    end

    properties(Access = private)
        % Private RCS value
        pRCS;
    end

    properties(Access = private, Nontunable)
        % Private flag about whether to specify mean rcs through input
        pInputRCS;
    end

    methods
        function set.MeanRCS(obj, value)
            validateattributes(value,{'double'},{'row','finite',...
                'nonnegative'},'phased.RadarTarget','MeanRCS');
            obj.MeanRCS = value;
        end
        function set.ScatteringMatrix(obj, value)
            validateattributes(value,{'double'},{'finite',...
                'size',[2 2]},'phased.RadarTarget','ScatteringMatrix');
            obj.ScatteringMatrix = value;
        end
    end

    methods
        function obj = RadarTarget(varargin)
            obj@phased.internal.AbstractRadarTarget(varargin{:});
        end
    end

    methods (Access = protected)
        
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractRadarTarget(obj);
            cond = obj.EnablePolarization && (obj.Mode(1) == 'M') && ... %Monostatic
                    (obj.ScatteringMatrixSource(1) == 'P') && ... %Property
                    obj.ScatteringMatrix(2) ~= obj.ScatteringMatrix(3) ;
            if cond
                coder.internal.errorIf(cond,'phased:target:ExpectedSymmetric',...
                    'ScatteringMatrix','Mode','Monostatic');
            end
        end

        function validateInputsImpl(obj,x,ang_in,ang_out,lclaxes,meanrcsArg,rcsupdateArg)
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

                angsize = size(ang_in);
                cond = xsize(2) ~= angsize(2);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:collector:AngleDimensionMismatch');
                end
                sigdatatypes.validateAzElAngle(ang_in,'','AngIn');
                if (obj.Mode(1) == 'B') %Bistatic
                    angsize = size(ang_out);
                    cond = xsize(2) ~= angsize(2);
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:collector:AngleDimensionMismatch');
                    end
                    sigdatatypes.validateAzElAngle(ang_out,'','AngOut');
                end
                
                if (obj.Mode(1) == 'B') %Bistatic
                    laxes = lclaxes;
                else
                    laxes = ang_out;
                end
                sigdatatypes.validate3DCartCoord(laxes,'','LAxes',...
                    {'size',[3 3]});
                
            else
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
                
            validateNumChannels(obj,x);

            fluctuateFlag = (obj.Model(1) == 'S');%Swerling
            
            if (obj.EnablePolarization && obj.ScatteringMatrixSource(1) == 'I') || ...
                    (~obj.EnablePolarization && obj.MeanRCSSource(1) == 'I')
                RCSViaInputPort = true;
            else
                RCSViaInputPort = false;
            end
            
            
            if RCSViaInputPort
                if obj.EnablePolarization
                    if (obj.Mode(1) == 'M') %Monostatic
                        if fluctuateFlag
                            rcsupdate = meanrcsArg;
                        end
                        meanrcs = lclaxes;
                    else
                        if fluctuateFlag
                            rcsupdate = rcsupdateArg;
                        end
                         meanrcs = meanrcsArg;
                    end
                    validateattributes(meanrcs,{'double'},{'finite',...
                        'size',[2 2]},'','SMat');
                else
                    meanrcs = ang_in;
                    if fluctuateFlag
                        rcsupdate = ang_out;
                    end
                    cond =   ~isa(meanrcs,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','RCS','double');
                    end
                    cond =   ~isscalar(meanrcs) && ~isequal(size(meanrcs),[1 size(x,2)]);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:target:mustBeScalarOrSameColumns','RCS','X');
                    end
                    cond =   ~isreal(meanrcs);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:NeedReal', 'RCS');
                    end
                end

            else
                if fluctuateFlag
                    if obj.EnablePolarization
                        if (obj.Mode(1) == 'M') %Monostatic
                            rcsupdate = lclaxes;
                        else
                            rcsupdate = meanrcsArg;
                        end
                    else
                        rcsupdate = ang_in;
                    end
                end
                cond =   ~obj.EnablePolarization && ...
                         ~isscalar(obj.MeanRCS) && ...
                         ~isequal(size(obj.MeanRCS),[1 size(x,2)]);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:target:mustBeScalarOrSameColumns','MeanRCS','X');
                end

            end

            if fluctuateFlag
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

        function setupImpl(obj,x,ang_in,ang_out,lclaxes,meanrcs,updateRCS) %#ok<INUSD>
            setupImpl@phased.internal.AbstractRadarTarget(obj,x);
            obj.pLambda = obj.PropagationSpeed/obj.OperatingFrequency;
            if obj.EnablePolarization
                obj.pInputRCS = (obj.ScatteringMatrixSource(1) == 'I'); %Input port
                obj.pRCS = complex(zeros(2,2)); %needed to codegen real value will be set in reset
            else
                obj.pInputRCS = (obj.MeanRCSSource(1) == 'I'); %Input port
                if obj.pInputRCS
                    obj.pRCS = ang_in;
                else
                    obj.pRCS = obj.MeanRCS; %needed to codegen real value will be set in reset                     
                end
            end
        end

        function resetImpl(obj)
            %setup random stream
            resetImpl@phased.internal.AbstractRadarTarget(obj);
            if ~obj.pFluctuate && ~obj.pInputRCS
                if obj.EnablePolarization
                    obj.pRCS = obj.ScatteringMatrix;
                else
                    obj.pRCS = obj.MeanRCS;
                end
            end

        end

        function processTunedPropertiesImpl(obj)
            if ~obj.pInputRCS && ~obj.pFluctuate
                if obj.EnablePolarization
                    obj.pRCS = obj.ScatteringMatrix;
                else
                    obj.pRCS = obj.MeanRCS;
                end
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractRadarTarget(obj);
            if isLocked(obj)
                s.pInputRCS = obj.pInputRCS;
                s.pRCS = obj.pRCS;
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s,wasLocked);
            if isfield(s,'pIsPrivateRCSEmpty')
                obj.pNeedRCSInit = s.pIsPrivateRCSEmpty;
                s = rmfield(s,'pIsPrivateRCSEmpty');
            end
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function y = stepImpl(obj,x,ang_in,ang_out,lclaxesArg,meanrcsArg,updateRCSArg)
            if obj.EnablePolarization
                if obj.pInputRCS
                    if (obj.Mode(1) == 'M') %Monostatic
                        if obj.pFluctuate
                            updateRCS = meanrcsArg;
                        end
                        meanrcs = lclaxesArg;
                        cond = meanrcs(2) ~= meanrcs(3); %symmetric
                        if cond
                            coder.internal.errorIf(cond,'phased:target:ExpectedSymmetric',...
                                'SMat','Mode','Monostatic');
                        end
                    else
                        if obj.pFluctuate
                           updateRCS = updateRCSArg;
                        end
                        meanrcs = meanrcsArg;
                    end
                    if obj.pFluctuate
                        sigma = meanrcs;
                    else
                        obj.pRCS = meanrcs;
                    end
                else
                    if obj.pFluctuate
                        if (obj.Mode(1) == 'M') %Monostatic
                            updateRCS = lclaxesArg;  
                        else
                            updateRCS = meanrcsArg;
                        end
                    end
                    sigma = obj.ScatteringMatrix;
                end

                if obj.pFluctuate
                    if obj.pNeedRCSInit
                        if (obj.Mode(1) == 'M') %Monostatic
                            obj.pRCS = obj.pFluctuateFunc(obj.cNoiseSource,sigma,true);
                        else
                            obj.pRCS = obj.pFluctuateFunc(obj.cNoiseSource,sigma,false);
                        end
                        obj.pNeedRCSInit = false;
                    elseif updateRCS
                        if (obj.Mode(1) == 'M') %Monostatic
                            obj.pRCS = obj.pFluctuateFunc(obj.cNoiseSource,sigma,true);
                        else
                            obj.pRCS = obj.pFluctuateFunc(obj.cNoiseSource,sigma,false);
                        end
                    end
                end
                
                if (obj.Mode(1) == 'M') %Monostatic
                    lclaxes = ang_out;
                else
                    lclaxes = lclaxesArg;
                end
                
                % x should be all same size
                % ang_in and ang_out should match x column numbers
                num_ang = size(ang_in,2);
                M = size(x,2);
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
                    % convert incoming signal to local HV
                    lcl_x = phased.internal.global2localvec(...
                        [x_x(:,m) x_y(:,m) x_z(:,m)].',lclaxes);
                    lcl_x = cart2sphvec(lcl_x,ang_in(1,m),ang_in(2,m));
                    lcl_x_sph = lcl_x(1:2,:);
                    % rcs could be complex, use direct computation.
                    g = reshape(sqrt(4*pi/obj.pLambda^2)*obj.pRCS(:),2,2);
                    lcl_y_sph = g*lcl_x_sph;
                    if (obj.Mode(1) == 'M') %Monostatic
                        lcl_y = sph2cartvec(...
                            [lcl_y_sph;zeros(1,size(lcl_y_sph,2))],...
                            ang_in(1,m),ang_in(2,m));
                    else
                        lcl_y = sph2cartvec(...
                            [lcl_y_sph;zeros(1,size(lcl_y_sph,2))],...
                            ang_out(1,m),ang_out(2,m));
                    end
                    g_y = phased.internal.local2globalvec(...
                        lcl_y,lclaxes);
                    y(m).X = g_y(1,:).';
                    y(m).Y = g_y(2,:).';
                    y(m).Z = g_y(3,:).';
                end
                
            else  % no polarization
                if obj.pInputRCS
                    meanrcs = ang_in;

                    cond =  any(meanrcs < 0);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:step:expectedNonnegative', 'RCS');                    
                    end
                    
                    if obj.pFluctuate
                        updateRCS = ang_out;
                        sigma = meanrcs;
                    else
                        obj.pRCS = meanrcs;
                    end
                else
                    if obj.pFluctuate
                        updateRCS = ang_in;  % 2nd input is updateRCS
                    end
                    sigma = obj.MeanRCS;
                end

                if obj.pFluctuate
                    if obj.pNeedRCSInit
                        obj.pRCS = obj.pFluctuateFunc(obj.cNoiseSource,sigma);
                        obj.pNeedRCSInit = false;
                    elseif updateRCS
                        obj.pRCS = obj.pFluctuateFunc(obj.cNoiseSource,sigma);
                    end
                end              
                g = sqrt(db2pow(aperture2gain(obj.pRCS,obj.pLambda)));
                if isscalar(g)
                    %one rcs, several channels
                    %This is when we have signals coming from several directions 
                    y = g*x;
                else
                    %one rcs per channel
                    y = bsxfun(@times, x, g);
                end
            end

        end

        function num = getNumInputsImpl(obj)
            num  = 1;
            if (obj.Model(1) == 'S');%Swerling
                num = num+1;
            end
            if obj.EnablePolarization
                num = num+2;
                if (obj.Mode(1) == 'B') %Bistatic
                    num = num+1;
                end
                if (obj.ScatteringMatrixSource(1) == 'I'); %Input port
                    num = num+1;
                end
            else
                if (obj.MeanRCSSource(1) == 'I'); %Input port
                    num = num+1;
                end
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractRadarTarget(obj,prop);
            if obj.EnablePolarization
                if strcmp(prop,'MeanRCS') || ...
                        strcmp(prop,'MeanRCSSource')
                    flag = true;
                end
                if  (obj.ScatteringMatrixSource(1) == 'I') && ...
                        strcmp(prop,'ScatteringMatrix')
                    flag = true;
                end
            else
                if strcmp(prop,'ScatteringMatrix') || ...
                        strcmp(prop,'Mode') || ...
                        strcmp(prop,'ScatteringMatrixSource')
                    flag = true;
                end
                if  (obj.MeanRCSSource(1) == 'I') && ... %Input port
                        strcmp(prop,'MeanRCS')
                    flag = true;
                end
            end
        end
    end
    
    methods (Access = protected, Static, Hidden)
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:RadarTargetTitle')),...
              'Text',getString(message('phased:library:block:RadarTargetDesc')));
        end
        
        function groups = getPropertyGroupsImpl
            dEnablePolarization = ...
                matlab.system.display.internal.Property('EnablePolarization', ...
                                                        'IsGraphical', false);

            dSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false,'UseClassDefault',false,...
                'Default','Property');
            dSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false,'UseClassDefault',false,...
                'Default','randi(65535,1)');
            props = {...
                dEnablePolarization, ...
                'Mode', ...
                'MeanRCSSource', ...
                'ScatteringMatrixSource', ...
                'MeanRCS', ...
                'ScatteringMatrix', ...
                'Model', ...
                'PropagationSpeed', ...
                'OperatingFrequency', ...
                dSeedSource, ...
                dSeed};
            groups = matlab.system.display.Section(...
                'phased.RadarTarget','PropertyList',props);
       end
    end
    
    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)   
            if obj.EnablePolarization
                if (obj.ScatteringMatrixSource(1) == 'I') %Input port
                    if (obj.Model(1) == 'S') %Swerling
                        if (obj.Mode(1) == 'B')
                            varargout = {'X','AngIn','AngOut','LAxes',...
                                'SMat','Update'};
                        else
                            varargout = {'X','AngIn','LAxes',...
                                'SMat','Update'};
                        end
                    else
                        if (obj.Mode(1) == 'B')
                            varargout = {'X','AngIn','AngOut','LAxes',...
                                'SMat'};
                        else
                            varargout = {'X','AngIn','LAxes','SMat'};
                        end
                    end
                else
                     if (obj.Model(1) == 'S') %Swerling
                        if (obj.Mode(1) == 'B')
                            varargout = {'X','AngIn','AngOut','LAxes',...
                                'Update'};
                        else
                            varargout = {'X','AngIn','LAxes',...
                                'Update'};
                        end
                    else
                        if (obj.Mode(1) == 'B')
                            varargout = {'X','AngIn','AngOut','LAxes'};
                        else
                            varargout = {'X','AngIn','LAxes'};
                        end
                    end
                end
            else
                if (obj.MeanRCSSource(1) == 'I') %Input port
                    if (obj.Model(1) == 'S')
                        varargout = {'X','RCS','Update'};
                    else
                        varargout = {'X','RCS'};
                    end
                elseif (obj.Model(1) == 'S')
                    varargout = {'X','Update'};
                else
                    varargout = {'X'};
                end     
            end
        end

        function varargout = getOutputNamesImpl(obj)  %#ok<MANU>
            varargout = {''};
        end
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Target');
        end
        function varargout = getOutputSizeImpl(obj)
           varargout{1} = propagatedInputSize(obj,1);
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj, 1);
        end
        function varargout = getOutputDataTypeImpl(obj)
            varargout{1} = propagatedInputDataType(obj,1);
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout{1} = true;
        end    
    end
end

