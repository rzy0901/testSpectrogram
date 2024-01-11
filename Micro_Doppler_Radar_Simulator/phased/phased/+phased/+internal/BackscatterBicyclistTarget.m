classdef (Sealed) BackscatterBicyclistTarget < phased.internal.AbstractSampleRateEngine & matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.

%BackscatterBicyclistTarget    Signal reflection off bicyclist
%   H = phased.internal.BackscatterBicyclistTarget creates a backscatter
%   bicyclist target System object, H, that computes the reflected signal
%   from a bicyclist. A backscatter target is designed to be used in a
%   monostatic setting where the incident and reflect angles of the signal
%   at the target are the same.
%
%   H = phased.internal.BackscatterBicyclistTarget(Name,Value) creates a
%   radar target object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The bicyclist is derived from a multi-scatterer model developed for a
%   77 GHz radar system. It is composed of 5 primary parts: bicycle frame
%   and upper-body of rider, pedals, legs, front wheel, and rear wheel.
%   Each part is composed of individual point scatterers. The total number
%   of point scatterers in the bicycle target model is dependent on the
%   number of spokes defined. The bicycle model is the aggregate of all of
%   the point scatterers.
%   
%   Step method syntax:
%
%   Y = step(H,X,ANG) returns the reflected signal Y off the bicyclist
%   target due to the input signal X.
%
%   X is an MxN matrix where M is the number of samples in the signal
%   incident to each point scatterer, and N is the number of point
%   scatterers.
%
%   ANG is a 2xN matrix representing the signal's incident direction to
%   each point scatterer. Each column of ANG specifies the incident
%   direction of the corresponding signal in the form of an [AzimuthAngle;
%   ElevationAngle] pair (in degrees).
%
%   Y is an M-length vector containing the combined reflected signals from
%   all point scatterers. The number of rows in Y is the same as the number
%   of rows in X.
%
%   RCS values depend on the viewing angle. The RCS values for the
%   individual point scatterers are computed from the angle of view of the
%   whole bicyclist. The default value of this property is a 1x361 matrix
%   with values that were derived from measurements taken at 77 GHz.
%
%   Internal occlusions within the bicyclist are not modeled.
%
%   System objects may be called directly like a function instead of
%   using the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   BackscatterBicyclistTarget methods:
%
%   step     - Signal reflection off bicyclist 
%   release  - Allow property value and input characteristics changes
%   clone    - Create a backscatter bicyclist target object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset the backscatter bicyclist target to its initial state
%
%   BackscatterBicyclistTarget properties:
%
%   OperatingFrequency     - Signal carrier frequency (Hz)
%   PropagationSpeed       - Propagation speed (m/s)
%   RCSPatternSource       - Source of RCS pattern 
%   AzimuthAngles          - Azimuth angles (deg)
%   ElevationAngles        - Elevation angles (deg)
%   RCSPattern             - Radar cross section pattern (square meters)
%
%   % Example:
%   %   Compute the reflected signal from a bicyclist moving along the
%   %   x axis. Assume the radar works at 24 GHz and the signal has a 1 GHz
%   %   bandwidth. The signal is captured at the moment the bicyclist
%   %   starts to move and 1 second into the movement. Assume the radar is
%   %   at the origin.
%             
%   % Define parameters
%   c = 3e8;
%   bw = 3e8;
%   fs = bw;
%   fc = 24e9;
%             
%   % Initialize bicyclist object
%   bicycleMove = phased.internal.MovingBicyclist(...
%       'InitialPosition',[5;0;0]);
%   N = bicycleMove.getNumScatterers;
%   bicycleReflect = phased.internal.BackscatterBicyclistTarget(...
%       'OperatingFrequency',fc);
%             
%   % Initialize radar and propagation channel
%   rpos = [0;0;0];
%   wav = phased.LinearFMWaveform('SampleRate',fs,'SweepBandwidth',bw);
%   x = wav();
%   chan = phased.FreeSpace('OperatingFrequency',fc,'SampleRate',fs,...
%       'TwoWayPropagation',true);
%             
%   % Time 0
%   [bbpos,bbvel,bbax] = step(bicycleMove,1);
%   xp = chan(repmat(x,1,N),rpos,bbpos,[0;0;0],bbvel);
%   [~,ang] = rangeangle(rpos,bbpos,bbax);
%   y0 = step(bicycleReflect,xp,ang);
%             
%   % Time 1
%   [bbpos,bbvel,bbax] = step(bicycleMove,1);
%   xp = chan(repmat(x,1,N),rpos,bbpos,[0;0;0],bbvel);
%   [~,ang] = rangeangle(rpos,bbpos,bbax);
%   y1 = step(bicycleReflect,xp,ang);
%             
%   % Plot
%   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(wav));
%   ymf = mf([y0 y1]);
%   t = (0:size(ymf,1)-1)/fs;
%   plot(t,mag2db(abs(ymf)));
%   ylim([-200 0])
%   xlabel('Time (s)');
%   ylabel('Magnitude (dB)');
%
%   See also phased, backscatterPedestrian, phased.RadarTarget.
    
%   Copyright 2019 The MathWorks, Inc.

%   Reference
%   [1] Stolz, M. et al. "Multi-Target Reflection Point Model of Cyclists
%   for Automotive Radar." 2017 European Radar Conference (EURAD),
%   Nuremberg, 2017, pp. 94-97.
%
%   [2] Chen, V., D. Tahmoush, and W. J. Miceli. Radar Micro-Doppler
%   Signatures: Processing and Applications. The Institution of Engineering
%   and Technology: London, 2014.
%
%   [3] Belgiovane, D., and C. C. Chen. "Bicycles and  Human Riders
%   Backscattering at 77 GHz for Automotive Radar." 2016 10th European
%   Conference on Antennas and Propagation (EuCAP), Davos, 2016, pp. 1-5.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %OperatingFrequency     Signal carrier frequency (Hz)
        %   Specify the carrier frequency (in Hz) of the narrowband signal
        %   as a scalar. The default value of this property is 77e9 (77
        %   GHz).
        OperatingFrequency = 77e9
        %PropagationSpeed       Propagation speed (m/s)
        %   Specify the wave propagation speed (in m/s) in free space as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('lightspeed')
    end
    
    properties (Nontunable)
        %RCSPatternSource       Source of RCS pattern
        % Specify the source of the RCS pattern as one of 'Auto' ||
        % 'Property'. When you set the RCSPatternSource property to 'Auto',
        % the properties AzimuthAngles, ElevationAngles, and RCSPattern are
        % inactive, and default values are used for these properties. If
        % RCSPatternSource is set to 'Property', the properties
        % AzimuthAngles, ElevationAngles, and RCSPattern are active, and a
        % custom RCS Pattern can be defined. Defaults to 'Auto'.
        RCSPatternSource = 'Auto'
    end
    
    properties (Nontunable)
        %AzimuthAngles          Azimuth angles (deg)
        %   Specify the azimuth angles (in degrees) as a length P vector.
        %   These are the azimuth angles where the custom pattern is
        %   evaluated. P must be greater than 2. The default value of this
        %   property is -180:180.
        AzimuthAngles = -180:180
        %ElevationAngles        Elevation angles (deg)
        %   Specify the elevation angles (in degrees) as a length Q vector.
        %   These are the elevation angles where the custom pattern is
        %   evaluated. The default value of this property is 0.
        ElevationAngles = 0
    end
    
    properties
        %RCSPattern             Radar cross section pattern (square meters)
        %   Specify the radar cross section pattern of the whole bicyclist
        %   (in square meters) as a QxP matrix, where Q is the number of
        %   elements presented in the ElevationAngles property, and P is
        %   the number of elements presented in the AzimuthAngles property.
        %
        %   Alternatively, if the pattern is only given at a specific
        %   elevation, then the pattern can be given as a 1xP vector.
        %
        %   The individual point scatterers get an RCS value based on the
        %   angle of view of the whole bicyclist. The default value of this
        %   property is a 1x361 matrix with values that were derived from
        %   measurements taken at 77 GHz.
        RCSPattern = backscatterBicyclist.defaultRCSPattern
    end
    
    properties(Constant, Hidden)
        RCSPatternSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
   
    properties (Hidden)
        NumScatterers = -1
    end
    
    properties (Access = protected)
        pLambda
        pRCSPattern
        pNumScatterers = -1 
    end
    
    methods
        function obj = BackscatterBicyclistTarget(varargin)
            setProperties(obj, nargin, varargin{:});
            cond = obj.NumScatterers~=-1;  
            if cond
                obj.pNumScatterers = obj.NumScatterers; 
            end
        end
    end
    
    methods
        function set.OperatingFrequency(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'nonnegative','real'},'phased.internal.BackscatterBicyclistTarget','OperatingFrequency');
            obj.OperatingFrequency = value;
        end
        function set.PropagationSpeed(obj,value)
            sigdatatypes.validateSpeed(value,'phased.internal.BackscatterBicyclistTarget',...
                'PropagationSpeed',{'scalar','positive','real'});
            obj.PropagationSpeed = value;
        end
        function set.AzimuthAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.internal.BackscatterBicyclistTarget','AzimuthAngles',...
                {'vector','>=',-180,'<=',180});
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'AzimuthAngles');
            end
            obj.AzimuthAngles = value;
        end
        function set.ElevationAngles(obj,value)
            % Allow single elevation cut, i.e., azimuth only pattern
            sigdatatypes.validateAngle(value,...
                'phased.internal.BackscatterBicyclistTarget','ElevationAngles',...
                {'vector','>=',-90,'<=',90});
            obj.ElevationAngles = value;
        end
        function set.RCSPattern(obj,value)
            validateattributes(value,{'double'},...
                {'real','nonnegative','finite','nonempty','2d'},...
                'phased.internal.BackscatterBicyclistTarget','RCSPattern');
            obj.RCSPattern  = value;
        end
    end

    methods (Access = protected)
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and system block dialog
            flag = false; 
            cond = strcmp(obj.RCSPatternSource,'Auto');
            if cond
                propChk = strcmpi(prop,'AzimuthAngles') || ...
                    strcmpi(prop,'ElevationAngles') || ...
                    strcmpi(prop,'RCSPattern');
                if propChk
                    flag = true;
                end
            end
        end
        
        function setupImpl(obj,x)
            setupImpl@phased.internal.AbstractSampleRateEngine(obj);
            obj.pLambda = obj.PropagationSpeed/obj.OperatingFrequency;
            obj.pRCSPattern = obj.RCSPattern;  % Used to check for changes to a "nontunable" property; this is done to support codegen
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        end
        
        function y = stepImpl(obj,x,ang)    
            reflectMethodNontunablePropertiesCheck(obj); % Check to make sure RCS pattern has not changed

            % Set number of scatterers 
            cond = obj.pNumScatterers == -1; 
            if cond
                obj.pNumScatterers = size(x,2); 
            end
            
            % Determine RCS for overall bicycle 
            isSingleEl = (numel(obj.ElevationAngles)==1); 
            pat_rcs = obj.RCSPattern;
            thisAng = median(ang,2);
            az = obj.AzimuthAngles; 
            el = obj.ElevationAngles;
            thisRCS = interpolatePattern(...
                az,el,pat_rcs,thisAng(1),thisAng(2),isSingleEl);
            thisRCS = thisRCS/size(x,2); 
            
            % Sum reflections
            g = sqrt(db2pow(aperture2gain(thisRCS,obj.pLambda)));
            y = g*x;
            y = sum(y,2);
        end
        
        function validateInputsImpl(obj,x,ang)
            % Verify that x is a double
            cond =  ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','X','double');
            end
            
            % Check to see if x is a matrix and nonempty
            cond =  ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeMatrix','X');
            end
            
            % Verify ang and x dim 2 matches 
            angsize = size(ang);
            xsize = size(x);
            cond = xsize(2) ~= angsize(2);
            if cond
                coder.internal.errorIf(cond,'phased:phased:collector:AngleDimensionMismatch');
            end
            sigdatatypes.validateAzElAngle(ang,'','Ang');
            
            % Verify size of x
            validateNumChannels(obj,x); 
            
            % Check angle values
            cond = any(ang(1,:)<-180) || any(ang(1,:)>180) || ...
                any(ang(2,:)<-90) || any(ang(2,:)>90);
            if cond
                coder.internal.errorIf(cond,'phased:phased:invalidAzElAngle','ANG',...
                    '-180','180','-90','90');
            end
        end

        function flag = isInputSizeMutableImpl(obj,index) %#ok<INUSL>
            % Return false if input size cannot change
            % between calls to the System object
            if index == 1
                flag = true;
            else
                flag = false;
            end
        end

        function flag = isInputComplexityMutableImpl(obj,index) %#ok<INUSL>
            % Return false if input complexity cannot change
            % between calls to the System object
            if index == 1
                flag = true;
            else
                flag = false;
            end
        end

        function flag = isInputDataTypeMutableImpl(obj,index) %#ok<INUSD>
            % Return false if input data type cannot change
            % between calls to the System object
            flag = false;
        end

        function num = getNumInputsImpl(obj) %#ok<MANU>
            % Define total number of inputs for system with optional inputs
            num = 2;
        end
    
        function loadObjectImpl(obj,s,~)
            % Set properties in object obj to values in structure s
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            
            if isLocked(obj)
                s.pNumScatterers = obj.pNumScatterers; 
                s.pLambda = obj.pLambda; 
                s.pRCSPattern = obj.pRCSPattern; 
            end

        end

        function icon = getIconImpl(obj) %#ok<MANU>
            % Define icon for System block
            icon = "Backscatter Bicyclist Target"; 
        end

        function [name,name2] = getInputNamesImpl(obj) %#ok<MANU>
            % Return input port names for System block
            name = 'X';
            name2 = 'ang';
        end

        function name = getOutputNamesImpl(obj) %#ok<MANU>
            % Return output port names for System block
            name = 'Y';
        end

        function out = getOutputSizeImpl(obj)
            % Return size for each output port
            out = propagatedInputSize(obj,1);
            out(2) = 1; 
        end

        function out = getOutputDataTypeImpl(obj) %#ok<MANU>
            % Return data type for each output port
            out = "double";
        end

        function out = isOutputComplexImpl(obj)
            % Return true for each output port with complex data
            out = propagatedInputComplexity(obj,1);
        end

        function out = isOutputFixedSizeImpl(obj)
            % Return true for each output port with fixed size
            out = propagatedInputFixedSize(obj,1);
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractSampleRateEngine(obj);
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
        end
    end
    methods (Hidden)
        function reflectMethodNontunablePropertiesCheck(obj)
            cond = any(size(obj.RCSPattern(:)) ~= size(obj.pRCSPattern(:)));
            if cond
                coder.internal.errorIf(cond,'phased:target:nontunableProperty','RCSPattern','BackscatterBicyclistTarget');
            end
            
            cond = any(obj.RCSPattern(:) ~= obj.pRCSPattern(:));
            if cond
                coder.internal.errorIf(cond,'phased:target:nontunableProperty','RCSPattern','BackscatterBicyclistTarget');
            end
        end
    end
end

% Interpolate pattern
function pat = interpolatePattern(az,el,pattern,az_q,el_q,interpIn1D)
if interpIn1D
    pat = interpolatePatternRadians1D(phased.internal.deg2rad(az), ...
        pattern, phased.internal.deg2rad(az_q(:)));
else
    pat = interpolatePatternRadians2D(phased.internal.deg2rad(az), phased.internal.deg2rad(el), ...
        pattern, phased.internal.deg2rad(az_q(:)), phased.internal.deg2rad(el_q(:)));
end
end

% Interpolate pattern 2D
function pat = interpolatePatternRadians2D(azr, elr, pattern, az_qr, el_qr)
pat = interp2(azr,elr,pattern,...
    az_qr, el_qr,'nearest',0);
end

% Interpolate pattern 1D
function pat = interpolatePatternRadians1D(azr, pattern, az_qr)
pat = interp1(azr,pattern,...
    az_qr,'nearest',0);
end