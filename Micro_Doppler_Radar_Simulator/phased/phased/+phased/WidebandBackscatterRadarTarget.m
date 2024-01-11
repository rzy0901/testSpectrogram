classdef (Sealed,StrictDefaults) WidebandBackscatterRadarTarget < phased.internal.AbstractRadarTarget
%WidebandBackscatterRadarTarget Wideband backscatter point radar target
%   H = phased.WidebandBackscatterRadarTarget creates a wideband
%   backscatter point target System object, H, that computes the reflected
%   wideband signal from a target. A backscatter target is designed to be
%   used in a monostatic setting where the incident and reflect angles of
%   the signal at the target are the same. The target radar cross section
%   (RCS) pattern is both angle and frequency dependent. The signal is
%   divided into subbands so the appropriate reflected signal in each
%   subband can be computed and then combined to form the total reflected
%   signal.
%
%   H = phased.WidebandBackscatterRadarTarget(Name,Value) creates a
%   wideband radar target object, H, with the specified property Name set
%   to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
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
%   The reflected signal in each subband is calculated as
%
%   Y = sqrt(G)*X where G = 4*pi*RCS/(LAMBDA^2).
%
%   G is the target gain factor corresponding to the target's radar cross
%   section RCS and LAMBDA is the wavelength of the subband. Note that G is
%   dimensionless.
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
%   The reflected signal in each subband is calculated as
%
%   Y_HV = sqrt(4*pi/(LAMBDA^2))*SMAT*X_HV 
%
%   where SMAT is the target's scattering matrix and LAMBDA is the
%   wavelength of the subband. X_HV and Y_HV are X and Y in the appropriate
%   HV coordinates. The target's radar cross section can be considered as
%   the magnitude squared of the entries in the scattering matrix.
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
%   WidebandBackscatterRadarTarget methods:
%
%   step     - Reflect the incoming signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create backscatter target object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of the radar target object
%   
%   WidebandBackscatterRadarTarget properties:
%
%   EnablePolarization  - Enable polarization
%   FrequencyVector     - Backscatter pattern frequency vector
%   AzimuthAngles       - Azimuth angles 
%   ElevationAngles     - Elevation angles 
%   RCSPattern          - Radar cross section pattern 
%   ShhPattern          - Radar scattering pattern for HH 
%   SvvPattern          - Radar scattering pattern for VV 
%   ShvPattern          - Radar scattering pattern for HV 
%   Model               - Fluctuation model
%   PropagationSpeed    - Propagation speed
%   OperatingFrequency  - Operating frequency 
%   SampleRate          - Sample rate 
%   NumSubbands         - Number of subbands
%   SeedSource          - Source of seed for random number generator
%   Seed                - Seed for random number generator
%
%   % Examples:
%   
%   % Example 1:
%   %   Calculate the reflected signal from a nonfluctuating point target
%   %   with RCS 10. Assume the incident angle is 10 degrees azimuth.
%
%   x = randn(10,1);
%   tgt = phased.WidebandBackscatterRadarTarget('Model','Nonfluctuating');
%   y = tgt(x,[10;0]);
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
%   x = struct('X',randn(10,1),'Y',randn(10,1),'Z',randn(10,1));
%   tgt = phased.WidebandBackscatterRadarTarget(...
%       'EnablePolarization',true,'Model','Swerling1',...
%       'ShhPattern',shhpat,'ShvPattern',shvpat,'SvvPattern',svvpat);
%   y = tgt(x,[10;0],tgtaxes,true);
%
%   See also phased, phased.BackscatterRadarTarget, phased.FreeSpace,
%   phased.Platform.

%   Copyright 2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %FrequencyVector    Backscatter pattern frequency vector (Hz)
        %   Specify the frequencies (in Hz) where the frequency responses
        %   of element are measured as a row vector. The elements of the
        %   vector must be increasing. The default of this property is [0
        %   1e20]. The antenna element has no response outside the
        %   specified frequency range.
        FrequencyVector = [0 1e20];
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
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a positive scalar. The
        %   default value is 1e6 (1 MHz).
        SampleRate = 1e6
    end

    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from Simulink time engine. Set SampleRateFromInputCheckbox to
        %   false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true
    end
    
    properties (Nontunable, PositiveInteger)
        %NumSubbands    Number of subbands
        %   Specify the number of subbands used in the subband processing
        %   as a positive integer. The default value of this property is
        %   64.
        NumSubbands = 64
    end
    
    properties (Access = private)
        %pRCSPattern - private property which holds the pattern info
        pInitialRCSPattern 
        pInitialShhPattern 
        pInitialSvvPattern 
        pInitialShvPattern 
    end
    
    properties (Access = private)
        pRCSPattern 
        pShhPattern 
        pSvvPattern 
        pShvPattern 
    end
    
    properties (Access = private, Nontunable, Logical)
        pIsPatternSingleElevation
        pIsSinglePattern
    end
    
    properties (Access = private, Nontunable)
        pSampleRate
        pSubbandFreqs
        pSubbandFreqVecIndex
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end

    properties (Nontunable)
        %RCSPattern   Radar cross section pattern (square meters)
        %   Specify the radar cross section pattern (in square meters).
        %
        %   For a given target, the pattern is specified as a QxPxK array,
        %   where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and K is the number of
        %   frequencies specified in the FrequencyVector property. Each
        %   page of the array defines an RCS pattern across the region
        %   defined by the azimuth and elevation angles at a given
        %   frequency. If K is 1, the same pattern is used across all
        %   frequencies. The default value of this property is a 181x361
        %   matrix with all elements equal to 1. 
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xPxK
        %   array or a KxP matrix. In this case, each row is a pattern for
        %   a given frequency
        %
        %   If you want to use the object to represent L targets in the
        %   simulation, put L patterns in a cell array and specify it to
        %   the property. All patterns should use the same format.
        %
        %   This property applies when the EnablePolarization property is
        %   false.
        RCSPattern = ones(181,361)
        %ShhPattern   Radar scattering pattern for HH (meters)
        %   Specify the scattering pattern for horizontal transmission and
        %   horizontal reception (in meters).
        %
        %   For a given target, the pattern is specified as a QxPxK array,
        %   where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and K is the number of
        %   frequencies specified in the FrequencyVector property. Each
        %   page of the array defines an RCS pattern across the region
        %   defined by the azimuth and elevation angles at a given
        %   frequency. If K is 1, the same pattern is used across all
        %   frequencies. The default value of this property is a 181x361
        %   matrix with all elements equal to 1. 
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xPxK
        %   array or a KxP matrix. In this case, each row is a pattern for
        %   a given frequency
        %
        %   If you want to use the object to represent L targets in the
        %   simulation, put L patterns in a cell array and specify it to
        %   the property. All patterns should use the same format.
        %
        %   This property applies when the EnablePolarization property is
        %   true.
        ShhPattern = ones(181,361)
        %ShvPattern   Radar scattering pattern for HV (meters)
        %   Specify the scattering pattern for vertical transmission and
        %   horizontal reception (in meters).
        %
        %   For a given target, the pattern is specified as a QxPxK array,
        %   where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and K is the number of
        %   frequencies specified in the FrequencyVector property. Each
        %   page of the array defines an RCS pattern across the region
        %   defined by the azimuth and elevation angles at a given
        %   frequency. If K is 1, the same pattern is used across all
        %   frequencies. The default value of this property is a 181x361
        %   matrix with all elements equal to 1. 
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xPxK
        %   array or a KxP matrix. In this case, each row is a pattern for
        %   a given frequency
        %
        %   If you want to use the object to represent L targets in the
        %   simulation, put L patterns in a cell array and specify it to
        %   the property. All patterns should use the same format.
        %
        %   This property applies when the EnablePolarization property is
        %   true.
        ShvPattern = zeros(181,361)
        %SvvPattern   Radar scattering pattern for VV (meters)
        %   Specify the scattering pattern for vertical transmission and
        %   vertical reception (in meters).
        %
        %   For a given target, the pattern is specified as a QxPxK array,
        %   where Q is the number of elements presented in the
        %   ElevationAngles property, P is the number of elements presented
        %   in the AzimuthAngles property, and K is the number of
        %   frequencies specified in the FrequencyVector property. Each
        %   page of the array defines an RCS pattern across the region
        %   defined by the azimuth and elevation angles at a given
        %   frequency. If K is 1, the same pattern is used across all
        %   frequencies. The default value of this property is a 181x361
        %   matrix with all elements equal to 1. 
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can also be given as either a 1xPxK
        %   array or a KxP matrix. In this case, each row is a pattern for
        %   a given frequency
        %
        %   If you want to use the object to represent L targets in the
        %   simulation, put L patterns in a cell array and specify it to
        %   the property. All patterns should use the same format.
        %
        %   This property applies when the EnablePolarization property is
        %   true.
        SvvPattern = ones(181,361)
    end
    
    properties (Access = private)
        cSubbandDivider
        cSubbandCombiner
    end
    
    methods
        function obj = WidebandBackscatterRadarTarget(varargin)
            obj@phased.internal.AbstractRadarTarget(varargin{:});
        end
    end

    methods
        
        function set.FrequencyVector(obj,value)
            validateattributes( value, { 'double' }, ...
                {'nonempty','finite','row','nonnegative','nondecreasing'}, '', 'FrequencyVector');
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'FrequencyVector');
            end
            obj.FrequencyVector = value;
        end
        
        function set.AzimuthAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.WidebandBackscatterRadarTarget','AzimuthAngles',...
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
                'phased.WidebandBackscatterRadarTarget','ElevationAngles',...
                {'vector','>=',-90,'<=',90});
            obj.ElevationAngles = value;
        end
        
        function set.RCSPattern(obj,value)
            validateattributes(value,{'double','cell'},{},...
                'phased.WidebandBackscatterRadarTarget','RCSPattern');
            if iscell(value)
                for m = 1:numel(value)
                    validateattributes(value{m},{'double'},...
                        {'real','nonnegative','nonempty','3d'},...
                        'phased.WidebandBackscatterRadarTarget',...
                        feval( 'sprintf' , 'RCSPattern{%d}',m));
                end
            else
                validateattributes(value,{'double'},...
                    {'real','nonnegative','nonempty','3d'},...
                    'phased.WidebandBackscatterRadarTarget',...
                    'RCSPattern');
            end

            obj.RCSPattern  = value;
        end
        
        function set.ShhPattern(obj,value)
            validateattributes(value,{'double','cell'},{},...
                'phased.WidebandBackscatterRadarTarget','ShhPattern');
            if iscell(value)
                for m = 1:numel(value)
                    validateattributes(value{m},{'double'},...
                        {'finite','nonempty','3d'},...
                        'phased.WidebandBackscatterRadarTarget',...
                        feval( 'sprintf' , 'ShhPattern{%d}',m));
                end
            else
                validateattributes(value,{'double'},...
                    {'finite','nonempty','3d'},...
                    'phased.WidebandBackscatterRadarTarget',...
                    'ShhPattern');
            end

            obj.ShhPattern  = value;
        end
        
        function set.SvvPattern(obj,value)
            validateattributes(value,{'double','cell'},{},...
                'phased.WidebandBackscatterRadarTarget','SvvPattern');
            if iscell(value)
                for m = 1:numel(value)
                    validateattributes(value{m},{'double'},...
                        {'finite','nonempty','3d'},...
                        'phased.WidebandBackscatterRadarTarget',...
                        feval( 'sprintf' , 'SvvPattern{%d}',m));
                end
            else
                validateattributes(value,{'double'},...
                    {'finite','nonempty','3d'},...
                    'phased.WidebandBackscatterRadarTarget',...
                    'SvvPattern');
            end

            obj.SvvPattern  = value;
        end
        
        function set.ShvPattern(obj,value)
            validateattributes(value,{'double','cell'},{},...
                'phased.WidebandBackscatterRadarTarget','ShvPattern');
            if iscell(value)
                for m = 1:numel(value)
                    validateattributes(value{m},{'double'},...
                        {'finite','nonempty','3d'},...
                        'phased.WidebandBackscatterRadarTarget',...
                        feval( 'sprintf' , 'ShvPattern{%d}',m));
                end
            else
                validateattributes(value,{'double'},...
                    {'finite','nonempty','3d'},...
                    'phased.WidebandBackscatterRadarTarget',...
                    'ShvPattern');
            end

            obj.ShvPattern  = value;
        end
        
        function set.SampleRate(obj,value)
            validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'SampleRate');
            obj.SampleRate = value;
        end
    end
    
    methods(Access = protected)
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractRadarTarget(obj);
            pattern_size = [numel(obj.ElevationAngles)...
                numel(obj.AzimuthAngles)];
            num_freq = numel(obj.FrequencyVector);
            if obj.EnablePolarization
                % polarization
                if numel(obj.ElevationAngles) == 1
                    % azimuth pattern only
                    if iscell(obj.ShhPattern)
                        % multiple targets
                        validateattributes(obj.SvvPattern,{'cell'},...
                            {'size',size(obj.ShhPattern)},...
                            'phased.WidebandBackscatterRadarTarget','SvvPattern');
                        validateattributes(obj.ShvPattern,{'cell'},...
                            {'size',size(obj.ShhPattern)},...
                            'phased.WidebandBackscatterRadarTarget','ShvPattern');
                        if numel(size(obj.ShhPattern{1})) == 3
                            % multiple 3D azimuth wideband pattern
                            % {1xAzxFreq}
                            for m = 1:numel(obj.ShhPattern)
                                validateattributes(obj.ShhPattern{m},{'double'},...
                                    {'size',[pattern_size num_freq]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShhPattern(%d)',m));
                                validateattributes(obj.SvvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'SvvPattern(%d)',m));
                                validateattributes(obj.ShvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShvPattern(%d)',m));
                            end
                        elseif size(obj.ShhPattern{1},1) == 1
                            % multiple azimuth wideband pattern, same
                            % across frequencies, {1xAz}
                            for m = 1:numel(obj.ShhPattern)
                                validateattributes(obj.ShhPattern{m},{'double'},...
                                    {'size',[1 pattern_size(2)]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShhPattern(%d)',m));
                                validateattributes(obj.SvvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'SvvPattern(%d)',m));
                                validateattributes(obj.ShvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShvPattern(%d)',m));
                            end
                        else
                            % multiple azimuth wideband pattern
                            % {FreqxAz}
                            for m = 1:numel(obj.ShhPattern)
                                validateattributes(obj.ShhPattern{m},{'double'},...
                                    {'size',[num_freq pattern_size(2)]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShhPattern(%d)',m));
                                validateattributes(obj.SvvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'SvvPattern(%d)',m));
                                validateattributes(obj.ShvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShvPattern(%d)',m));
                            end
                        end
                    else
                        % single target
                        if numel(size(obj.ShhPattern)) == 3
                            % single 3D wideband pattern
                            % 1xAzxFreq
                            validateattributes(obj.ShhPattern,{'double'},...
                                {'size',[pattern_size num_freq]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'ShhPattern');
                        elseif size(obj.ShhPattern,1) == 1
                            % single azimuth wideband pattern, same
                            % across frequencies, 1xAz
                            validateattributes(obj.ShhPattern,{'double'},...
                                {'size',[1 pattern_size(2)]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'ShhPattern');
                        else
                            % single azimuth wideband pattern
                            % FreqxAz
                            validateattributes(obj.ShhPattern,{'double'},...
                                {'size',[num_freq pattern_size(2)]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'ShhPattern');
                        end
                        validateattributes(obj.SvvPattern,{'double'},...
                            {'size',size(obj.ShhPattern)},...
                            'phased.BackscatterRadarTarget','SvvPattern');
                        validateattributes(obj.ShvPattern,{'double'},...
                            {'size',size(obj.ShhPattern)},...
                            'phased.BackscatterRadarTarget','ShvPattern');
                        
                    end
                else
                    % 3D pattern
                    if iscell(obj.ShhPattern)
                        % multiple targets
                        validateattributes(obj.SvvPattern,{'cell'},...
                            {'size',size(obj.ShhPattern)},...
                            'phased.WidebandBackscatterRadarTarget','SvvPattern');
                        validateattributes(obj.ShvPattern,{'cell'},...
                            {'size',size(obj.ShhPattern)},...
                            'phased.WidebandBackscatterRadarTarget','ShvPattern');
                        if ismatrix(obj.ShhPattern{1})
                            % multiple 3D wideband pattern, share across
                            % frequencies {ElxAz}
                            for m = 1:numel(obj.ShhPattern)
                                validateattributes(obj.ShhPattern{m},{'double'},...
                                    {'size',pattern_size},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShhPattern(%d)',m));
                                validateattributes(obj.SvvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'SvvPattern(%d)',m));
                                validateattributes(obj.ShvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShvPattern(%d)',m));
                            end
                        else
                            % multiple 3D wideband pattern
                            % {ElxAzxFreq}
                            for m = 1:numel(obj.ShhPattern)
                                validateattributes(obj.ShhPattern{m},{'double'},...
                                    {'size',[pattern_size num_freq]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShhPattern(%d)',m));
                                validateattributes(obj.SvvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'SvvPattern(%d)',m));
                                validateattributes(obj.ShvPattern{m},{'double'},...
                                    {'size',size(obj.ShhPattern{m})},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'ShvPattern(%d)',m));
                            end
                        end
                    else
                        if ismatrix(obj.ShhPattern)
                            % single 3D wideband pattern, share across
                            % frequencies ElxAz
                            validateattributes(obj.ShhPattern,{'double'},...
                                {'size',pattern_size},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'ShhPattern');
                            validateattributes(obj.SvvPattern,{'double'},...
                                {'size',size(obj.ShhPattern)},...
                                'phased.BackscatterRadarTarget','SvvPattern');
                            validateattributes(obj.ShvPattern,{'double'},...
                                {'size',size(obj.ShhPattern)},...
                                'phased.BackscatterRadarTarget','ShvPattern');
                        else
                            % single 3D wideband pattern
                            % ElxAzxFreq
                            validateattributes(obj.ShhPattern,{'double'},...
                                {'size',[pattern_size num_freq]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'ShhPattern');
                            validateattributes(obj.SvvPattern,{'double'},...
                                {'size',size(obj.ShhPattern)},...
                                'phased.BackscatterRadarTarget','SvvPattern');
                            validateattributes(obj.ShvPattern,{'double'},...
                                {'size',size(obj.ShhPattern)},...
                                'phased.BackscatterRadarTarget','ShvPattern');
                        end
                        
                    end
                end
            else
                % no polarization
                if numel(obj.ElevationAngles) == 1
                    % azimuth pattern only
                    if iscell(obj.RCSPattern)
                        % multiple target
                        if numel(size(obj.RCSPattern{1}))==3
                            for m = 1:numel(obj.RCSPattern)
                                % multiple 3D azimuth wideband pattern
                                validateattributes(obj.RCSPattern{m},{'double'},...
                                    {'size',[pattern_size num_freq]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'RCSPattern(%d)',m));
                            end
                        elseif size(obj.RCSPattern{1},1)==1
                            for m = 1:numel(obj.RCSPattern)
                                % multiple 3D azimuth wideband pattern,
                                % share across frequencies
                                validateattributes(obj.RCSPattern{m},{'double'},...
                                    {'size',[1 pattern_size(2)]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'RCSPattern(%d)',m));
                            end
                        else
                            for m = 1:numel(obj.RCSPattern)
                                % multiple azimuth wideband pattern
                                validateattributes(obj.RCSPattern{m},{'double'},...
                                    {'size',[num_freq pattern_size(2)]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'RCSPattern(%d)',m));
                            end
                        end
                    else
                        % single target
                        if numel(size(obj.RCSPattern))==3
                            % multiple 3D azimuth wideband pattern
                            validateattributes(obj.RCSPattern,{'double'},...
                                {'size',[pattern_size num_freq]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'RCSPattern');
                        elseif size(obj.RCSPattern,1)==1
                            % multiple 3D azimuth wideband pattern,
                            % share across frequencies
                            validateattributes(obj.RCSPattern,{'double'},...
                                {'size',[1 pattern_size(2)]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'RCSPattern(%d)');
                        else
                            % multiple azimuth wideband pattern
                            validateattributes(obj.RCSPattern,{'double'},...
                                {'size',[num_freq pattern_size(2)]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'RCSPattern');
                        end
                    end
                else
                    % 3D pattern
                    if iscell(obj.RCSPattern)
                        if numel(size(obj.RCSPattern{1}))==3
                            % multiple target
                            for m = 1:numel(obj.RCSPattern)
                                % multiple 3D wideband pattern
                                validateattributes(obj.RCSPattern{m},{'double'},...
                                    {'size',[pattern_size num_freq]},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'RCSPattern(%d)',m));
                            end
                        else
                             % multiple target, share pattern across
                             % frequencies
                            for m = 1:numel(obj.RCSPattern)
                                % multiple 3D wideband pattern
                                validateattributes(obj.RCSPattern{m},{'double'},...
                                    {'size',pattern_size},...
                                    'phased.WidebandBackscatterRadarTarget',...
                                    feval( 'sprintf' , 'RCSPattern(%d)',m));
                            end
                       end
                    else
                        if numel(size(obj.RCSPattern))==3
                            % single 3D wideband pattern
                            validateattributes(obj.RCSPattern,{'double'},...
                                {'size',[pattern_size num_freq]},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'RCSPattern');
                        else
                            % single 3D wideband pattern across frequencies
                            validateattributes(obj.RCSPattern,{'double'},...
                                {'size',pattern_size},...
                                'phased.WidebandBackscatterRadarTarget',...
                                'RCSPattern');
                        end
                    end
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
                             'MATLAB:system:invalidInputDataType',feval( 'sprintf' , 'X(%d).X',m),'double');
                    end
                    cond =  ~iscolumn(x_x) || isempty(x_x);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector',feval( 'sprintf' , 'X(%d).X',m));
                    end
                    cond =  ~isa(x_y,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType',feval( 'sprintf' , 'X(%d).Y',m),'double');
                    end
                    cond =  ~iscolumn(x_y) || isempty(x_y);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector',feval( 'sprintf' , 'X(%d).Y',m));
                    end
                    cond =  ~isa(x_z,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType',feval( 'sprintf' , 'X(%d).Z',m),'double');
                    end
                    cond =  ~iscolumn(x_z) || isempty(x_z);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector',feval( 'sprintf' , 'X(%d).Z',m));
                    end
                    cond = numel(x_x)~=numel(x_y) || numel(x_x)~=numel(x_z);
                    if cond
                        coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                            'X,Y,Z',feval( 'sprintf' , 'X(%d)',m));
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
                if iscell(obj.ShhPattern)
                    npat = numel(obj.ShhPattern);
                else
                    npat = 1;
                end
            else
                if iscell(obj.RCSPattern)
                    npat = numel(obj.RCSPattern);
                else
                    npat = 1;
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
                s.cSubbandDivider = saveobj(obj.cSubbandDivider);
                s.cSubbandCombiner = saveobj(obj.cSubbandCombiner);
                s.pSubbandFreqs = obj.pSubbandFreqs;
                s.pSubbandFreqVecIndex = obj.pSubbandFreqVecIndex;
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            s = loadSubObjects@phased.internal.AbstractRadarTarget(obj,s,wasLocked);
            if wasLocked
                obj.cSubbandDivider = phased.internal.SubbandDivider.loadobj(s.cSubbandDivider);
                s = rmfield(s,'cSubbandDivider');
                obj.cSubbandCombiner = phased.internal.SubbandCombiner.loadobj(s.cSubbandCombiner);
                s = rmfield(s,'cSubbandCombiner');
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
            subbandfreq = phased.internal.subbandCenterFrequency(...
                obj.OperatingFrequency,obj.SampleRate,obj.NumSubbands);
            obj.pLambda = obj.PropagationSpeed./subbandfreq;
            obj.pIsPatternSingleElevation = (numel(obj.ElevationAngles)==1);
            
            % convert all pattern to the form of ElxAzxFreq
            % If there is only one pattern, the result is a 3D array
            % if there are multiple patterns, the result is a cell array of
            % 3D arrays
            num_freq = numel(obj.FrequencyVector);
            if obj.EnablePolarization
                % polarized pattern
                if iscell(obj.ShhPattern)
                    if numel(obj.ShhPattern) == 1
                        % single pattern
                        if obj.pIsPatternSingleElevation && ismatrix(obj.ShhPattern{1})
                            if size(obj.ShhPattern{1},1)==1
                                obj.pInitialShhPattern = repmat(obj.ShhPattern{1},1,1,num_freq);
                                obj.pInitialSvvPattern = repmat(obj.SvvPattern{1},1,1,num_freq);
                                obj.pInitialShvPattern = repmat(obj.ShvPattern{1},1,1,num_freq);
                            else
                                obj.pInitialShhPattern = permute(obj.ShhPattern{1},[3 2 1]);
                                obj.pInitialSvvPattern = permute(obj.SvvPattern{1},[3 2 1]);
                                obj.pInitialShvPattern = permute(obj.ShvPattern{1},[3 2 1]);
                            end
                        elseif ismatrix(obj.ShhPattern{1})
                            obj.pInitialShhPattern = repmat(obj.ShhPattern{1},1,1,num_freq);
                            obj.pInitialSvvPattern = repmat(obj.SvvPattern{1},1,1,num_freq);
                            obj.pInitialShvPattern = repmat(obj.ShvPattern{1},1,1,num_freq);
                        else
                            obj.pInitialShhPattern = obj.ShhPattern{1};
                            obj.pInitialSvvPattern = obj.SvvPattern{1};
                            obj.pInitialShvPattern = obj.ShvPattern{1};
                        end
                        obj.pIsSinglePattern = true;
                    else
                        % multiple pattern
                        obj.pInitialShhPattern = cell(1,numel(obj.ShhPattern));
                        obj.pInitialSvvPattern = cell(1,numel(obj.SvvPattern));
                        obj.pInitialShvPattern = cell(1,numel(obj.ShvPattern));
                        if obj.pIsPatternSingleElevation && ismatrix(obj.ShhPattern{1})
                            if size(obj.ShhPattern{1},1)==1
                                for m = numel(obj.ShhPattern):-1:1
                                    obj.pInitialShhPattern{m} = repmat(obj.ShhPattern{m},1,1,num_freq);
                                    obj.pInitialSvvPattern{m} = repmat(obj.SvvPattern{m},1,1,num_freq);
                                    obj.pInitialShvPattern{m} = repmat(obj.ShvPattern{m},1,1,num_freq);
                                end
                            else
                                for m = numel(obj.ShhPattern):-1:1
                                    obj.pInitialShhPattern{m} = permute(obj.ShhPattern{m},[3 2 1]);
                                    obj.pInitialSvvPattern{m} = permute(obj.SvvPattern{m},[3 2 1]);
                                    obj.pInitialShvPattern{m} = permute(obj.ShvPattern{m},[3 2 1]);
                                end
                            end
                        elseif ismatrix(obj.ShhPattern{1})
                            for m = numel(obj.ShhPattern):-1:1
                                obj.pInitialShhPattern{m} = repmat(obj.ShhPattern{m},1,1,num_freq);
                                obj.pInitialSvvPattern{m} = repmat(obj.SvvPattern{m},1,1,num_freq);
                                obj.pInitialShvPattern{m} = repmat(obj.ShvPattern{m},1,1,num_freq);
                            end
                        else
                            obj.pInitialShhPattern = obj.ShhPattern;
                            obj.pInitialSvvPattern = obj.SvvPattern;
                            obj.pInitialShvPattern = obj.ShvPattern;
                        end
                        obj.pIsSinglePattern = false;
                    end
                else
                    % single pattern
                    if obj.pIsPatternSingleElevation && ismatrix(obj.ShhPattern)
                        if size(obj.ShhPattern,1) == 1
                            obj.pInitialShhPattern = repmat(obj.ShhPattern,1,1,num_freq);
                            obj.pInitialSvvPattern = repmat(obj.SvvPattern,1,1,num_freq);
                            obj.pInitialShvPattern = repmat(obj.ShvPattern,1,1,num_freq);
                        else
                            obj.pInitialShhPattern = permute(obj.ShhPattern,[3 2 1]);
                            obj.pInitialSvvPattern = permute(obj.SvvPattern,[3 2 1]);
                            obj.pInitialShvPattern = permute(obj.ShvPattern,[3 2 1]);
                        end
                    elseif ismatrix(obj.ShhPattern)
                        obj.pInitialShhPattern = repmat(obj.ShhPattern,1,1,num_freq);
                        obj.pInitialSvvPattern = repmat(obj.SvvPattern,1,1,num_freq);
                        obj.pInitialShvPattern = repmat(obj.ShvPattern,1,1,num_freq);
                    else
                        obj.pInitialShhPattern = obj.ShhPattern;
                        obj.pInitialSvvPattern = obj.SvvPattern;
                        obj.pInitialShvPattern = obj.ShvPattern;
                    end
                    obj.pIsSinglePattern = true;
                end
                obj.pShhPattern = obj.pInitialShhPattern;
                obj.pSvvPattern = obj.pInitialSvvPattern;
                obj.pShvPattern = obj.pInitialShvPattern;
            else
                % non-polarized pattern
                if iscell(obj.RCSPattern)
                    if numel(obj.RCSPattern) == 1
                        % single pattern
                        if obj.pIsPatternSingleElevation && ismatrix(obj.RCSPattern{1})
                            if size(obj.RCSPattern{1},1)==1
                                obj.pInitialRCSPattern = repmat(obj.RCSPattern{1},1,1,num_freq);
                            else
                                obj.pInitialRCSPattern = permute(obj.RCSPattern{1},[3 2 1]);
                            end
                        elseif ismatrix(obj.RCSPattern{1})
                            obj.pInitialRCSPattern = repmat(obj.RCSPattern{1},1,1,num_freq);
                        else
                            obj.pInitialRCSPattern = obj.RCSPattern{1};
                        end
                        obj.pIsSinglePattern = true;
                    else
                        % multiple pattern
                        obj.pInitialRCSPattern = cell(1,numel(obj.RCSPattern));
                        if obj.pIsPatternSingleElevation && ismatrix(obj.RCSPattern{1})
                            if size(obj.RCSPattern{1},1)==1
                                for m = numel(obj.RCSPattern):-1:1
                                    obj.pInitialRCSPattern{m} = repmat(obj.RCSPattern{m},1,1,num_freq);
                                end
                            else
                                for m = numel(obj.RCSPattern):-1:1
                                    obj.pInitialRCSPattern{m} = permute(obj.RCSPattern{m},[3 2 1]);
                                end
                            end
                        elseif ismatrix(obj.RCSPattern{1})
                            for m = numel(obj.RCSPattern):-1:1
                                obj.pInitialRCSPattern{m} = repmat(obj.RCSPattern{m},1,1,num_freq);
                            end
                        else
                            obj.pInitialRCSPattern = obj.RCSPattern;
                        end
                        obj.pIsSinglePattern = false;
                    end
                else
                    % single pattern
                    if obj.pIsPatternSingleElevation && ismatrix(obj.RCSPattern)
                        if size(obj.RCSPattern,1)==1
                            obj.pInitialRCSPattern = repmat(obj.RCSPattern,1,1,num_freq);
                        else
                            obj.pInitialRCSPattern = permute(obj.RCSPattern,[3 2 1]);
                        end
                    elseif ismatrix(obj.RCSPattern)
                        obj.pInitialRCSPattern = repmat(obj.RCSPattern,1,1,num_freq);
                    else
                        obj.pInitialRCSPattern = obj.RCSPattern;
                    end
                    obj.pIsSinglePattern = true;
                end
                obj.pRCSPattern = obj.pInitialRCSPattern;
            end
            
            % obj.pSampleRate = getSampleRate(obj,xsize(1),1,obj.SampleRate);
            fs = obj.SampleRate; % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            obj.pSampleRate = fs;

            nfft = obj.NumSubbands;
            obj.pSubbandFreqs = phased.internal.subbandCenterFrequency(...
                obj.OperatingFrequency,obj.pSampleRate,nfft).';
            
            obj.pSubbandFreqVecIndex = interp1(obj.FrequencyVector,...
                1:numel(obj.FrequencyVector),obj.pSubbandFreqs,'nearest',0);
            
           
            cond = obj.pSampleRate > 2*obj.OperatingFrequency;
            if cond
                coder.internal.errorIf(cond,'phased:phased:WidebandCollector:FsTooHigh');
            end
            
            obj.cSubbandDivider = phased.internal.SubbandDivider(...
                'OperatingFrequency',obj.OperatingFrequency,...
                'SampleRate',obj.pSampleRate,...
                'NumSubbands',nfft,'EnableWarning',false);
            obj.cSubbandCombiner = phased.internal.SubbandCombiner(...
                'NumSubbands',nfft,'TimeSignalLengthSource','Inherit',...
                'EnableWarning',false);
        end
        
        function resetImpl(obj)
            %setup random stream
            resetImpl@phased.internal.AbstractRadarTarget(obj);
            reset(obj.cSubbandDivider);
            reset(obj.cSubbandCombiner);
            if obj.EnablePolarization
                obj.pShhPattern = obj.pInitialShhPattern;
                obj.pShvPattern = obj.pInitialShvPattern;
                obj.pSvvPattern = obj.pInitialSvvPattern;
            else
                obj.pRCSPattern = obj.pInitialRCSPattern;
            end

        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractRadarTarget(obj);
            release(obj.cSubbandDivider);
            release(obj.cSubbandCombiner);
        end

        function y = stepImpl(obj,x,ang,lclaxesArg,updateRCSArg)
            az = obj.AzimuthAngles;
            el = obj.ElevationAngles;
            subbandidx = obj.pSubbandFreqVecIndex;
            isSingleEl = obj.pIsPatternSingleElevation;
            if obj.EnablePolarization
                % polarization
                lclaxes_in = lclaxesArg;

                if obj.pFluctuate
                    updateRCS = updateRCSArg;
                    % update the entire pattern at each trigger because the
                    % incident angle is unknown
                    if obj.pNeedRCSInit
                        if iscell(obj.pInitialShhPattern)
                            for m = numel(obj.pInitialShhPattern):-1:1
                                obj.pShhPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShhPattern{m},false);
                                obj.pSvvPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialSvvPattern{m},false);
                                obj.pShvPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShvPattern{m},false);
                            end
                        else
                            obj.pShhPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShhPattern,false);
                            obj.pSvvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialSvvPattern,false);
                            obj.pShvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShvPattern,false);
                        end
                        obj.pNeedRCSInit = false;
                    elseif updateRCS
                        if iscell(obj.pInitialShhPattern)
                            for m = numel(obj.pInitialShhPattern):-1:1
                                obj.pShhPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShhPattern{m},false);
                                obj.pSvvPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialSvvPattern{m},false);
                                obj.pShvPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShvPattern{m},false);
                            end
                        else
                            obj.pShhPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShhPattern,false);
                            obj.pSvvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialSvvPattern,false);
                            obj.pShvPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialShvPattern,false);
                        end
                    end
                end
                
                pat_hh = obj.pShhPattern;
                pat_vv = obj.pSvvPattern;
                pat_hv = obj.pShvPattern;
                
                % x should be all same size
                % ang should match x column numbers
                num_ang = size(ang,2);
                M = obj.pValidatedNumInputChannels;
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
                
                % convert to frequency domain nfft x chan x freq
                xsubband = step(obj.cSubbandDivider,[x_x,x_y,x_z]);
                xsubband_x = xsubband(:,1:M,:);
                xsubband_y = xsubband(:,M+1:2*M,:);
                xsubband_z = xsubband(:,2*M+1:3*M,:);
                
                ysubband_x = complex(zeros(size(xsubband_x),'like',xsubband_x));
                ysubband_y = complex(zeros(size(xsubband_y),'like',xsubband_y));
                ysubband_z = complex(zeros(size(xsubband_z),'like',xsubband_z));
                
                for m = num_ang:-1:1
                    
                    % rcs could be complex, use direct computation.
                
                    for n = obj.NumSubbands:-1:1
                        nidx = subbandidx(n);
                        if obj.pIsSinglePattern
                        
                            % convert incoming signal to local HV
                            lcl_x = phased.internal.global2localvec(...
                                [xsubband_x(:,m,n) xsubband_y(:,m,n) xsubband_z(:,m,n)].',lclaxes);
                            lcl_x = cart2sphvec(lcl_x,ang(1,m),ang(2,m));
                            lcl_x_sph = lcl_x(1:2,:);

                            if nidx==0
                                tempShh = 0;
                                tempSvv = 0;
                                tempShv = 0;
                            else
                                tempShh = interpolatePattern(...
                                    az,el,pat_hh(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                                tempSvv = interpolatePattern(...
                                    az,el,pat_vv(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                                tempShv = interpolatePattern(...
                                    az,el,pat_hv(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                            end
                        
                        else
                            % convert incoming signal to local HV
                            lcl_x = phased.internal.global2localvec(...
                                [xsubband_x(:,m,n) xsubband_y(:,m,n) xsubband_z(:,m,n)].',lclaxes(:,:,m));
                            lcl_x = cart2sphvec(lcl_x,ang(1,m),ang(2,m));
                            lcl_x_sph = lcl_x(1:2,:);

                            if nidx==0
                                tempShh = 0;
                                tempSvv = 0;
                                tempShv = 0;
                            else
                                tempShh = interpolatePattern(...
                                    az,el,pat_hh{m}(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                                tempSvv = interpolatePattern(...
                                    az,el,pat_vv{m}(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                                tempShv = interpolatePattern(...
                                    az,el,pat_hv{m}(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                            end
                        end
                    
                        % for backscatter, Shv = Svh
                        tempSmat = [tempShh tempShv;tempShv tempSvv]; 

                        g = sqrt(4*pi/obj.pLambda(n)^2)*tempSmat;

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
                        
                        ysubband_x(:,m,n) = g_y(1,:).';
                        ysubband_y(:,m,n) = g_y(2,:).';
                        ysubband_z(:,m,n) = g_y(3,:).';
                        
                    end
                    
                end
                ysubband = zeros([size(xsubband_y,1),3*M,size(ysubband_x,3)],'like',ysubband_x);
                ysubband(:,1:M,:) = ysubband_x;
                ysubband(:,M+1:2*M,:) = ysubband_y;
                ysubband(:,2*M+1:3*M,:) = ysubband_z;
                yt = step(obj.cSubbandCombiner,ysubband,complex(x(1).X));
                ylen = size(y(1).X,1);
                for m = M:-1:1
                    y(m).X = yt(1:ylen,m);
                    y(m).Y = yt(1:ylen,M+m);
                    y(m).Z = yt(1:ylen,2*M+m);
                end
                
            else  % no polarization
                if obj.pFluctuate
                    updateRCS = lclaxesArg;  % 3rd input is updateRCS
                end

                if obj.pFluctuate
                    if obj.pNeedRCSInit
                        if iscell(obj.pInitialRCSPattern)
                            for m = numel(obj.pInitialRCSPattern):-1:1
                                obj.pRCSPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialRCSPattern{m});
                            end
                        else
                            obj.pRCSPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialRCSPattern);
                        end
                        obj.pNeedRCSInit = false;
                    elseif updateRCS
                        if iscell(obj.pInitialRCSPattern)
                            for m = numel(obj.pInitialRCSPattern):-1:1
                                obj.pRCSPattern{m} = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialRCSPattern{m});
                            end
                        else
                            obj.pRCSPattern = obj.pFluctuateFunc(obj.cNoiseSource,obj.pInitialRCSPattern);
                        end
                    end
                end  
                
                pat_rcs = obj.pRCSPattern;
                
                num_ang = size(ang,2);
                
                % convert to frequency domain nfft x chan x freq
                xsubband = step(obj.cSubbandDivider,x);
                
                g = zeros(1,num_ang,obj.NumSubbands);
                for m = 1:num_ang
                    for n = 1:obj.NumSubbands
                        nidx = subbandidx(n);
                        if obj.pIsSinglePattern
                            if nidx==0
                                tempRCS = 0;
                            else
                                tempRCS = interpolatePattern(...
                                    az,el,pat_rcs(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                            end
                        else
                            if nidx==0
                                tempRCS = 0;
                            else
                                tempRCS = interpolatePattern(...
                                    az,el,pat_rcs{m}(:,:,nidx),ang(1,m), ang(2,m),isSingleEl);
                            end
                        end
                        g(1,m,n) = tempRCS;
                    end
                end
                
                g = sqrt(bsxfun(@rdivide,4*pi*g,permute(obj.pLambda.^2,[3 2 1])));
                
                %one rcs per channel
                ysubband = bsxfun(@times, xsubband, g);
                
                y = step(obj.cSubbandCombiner,ysubband,complex(x));
            end

        end

        function num = getNumInputsImpl(obj)
            num  = 2;
            if (obj.Model(1) == 'S') %Swerling
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
              'Title',getString(message('phased:library:block:WidebandBackScatterRadarTargetTitle')),...
              'Text',getString(message('phased:library:block:WidebandBackScatterRadarTargetDesc')));
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

            % dSampleRate = matlab.system.display.internal.Property(...
            %     'SampleRate','IsObjectDisplayOnly',true);

            groups = matlab.system.display.Section(...
                            'PropertyList',...
                            {...
                            dEnablePolarization,...
                            'FrequencyVector',...
                            'AzimuthAngles',...
                            'ElevationAngles',...
                            dRCSPattern,...
                            dShhPattern,...
                            dShvPattern,...
                            dSvvPattern,...
                            'Model', ...
                            'PropagationSpeed', ...
                            'OperatingFrequency', ...
                            'SampleRateFromInputCheckbox', ...
                            'SampleRate', ...
                            'NumSubbands', ...
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
            str = feval( 'sprintf' , 'Wideband\n Backscatter\n Target');
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

