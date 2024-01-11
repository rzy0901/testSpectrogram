classdef (Sealed,StrictDefaults) TimeDelayLCMVBeamformer < phased.internal.AbstractTimeDomainSMIBeamformer
%TimeDelayLCMVBeamformer    Time delay LCMV beamformer
%   H = phased.TimeDelayLCMVBeamformer returns a time delay linear
%   constraint minimum variance (LCMV) beamformer System object, H. This
%   object performs time delay LCMV beamforming on the received signal.
%
%   H = phased.TimeDelayLCMVBeamformer(Name,Value) creates a time delay
%   LCMV beamformer object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The beamforming algorithm is the time domain counterpart of the
%   narrowband linear constraint minimum variance (LCMV) beamformer. The
%   received samples at each element in the sensor array are fed into an
%   FIR filter. The algorithm first steers the array to the beamforming
%   direction and then applies the FIR filter to the output of each sensor
%   to achieve the specified constraints.
%
%   Step method syntax:
%
%   Y = step(H,X) performs time delay LCMV beamforming on the input X, and
%   returns the beamformed output in Y. X is an MxN matrix where N is the
%   number of elements of the sensor array. Y is a length-M column vector.
%   M must be larger than the FIR filter length specified in the
%   FilterLength property.
%
%   [Y,W] = step(H,X) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is a
%   length-L column vector where L is the degrees of freedom of the
%   beamformer. For a time delay LCMV beamformer, L is given by the product
%   of the number of elements in the sensor array specified by the
%   SensorArray property and the FIR filter length specified by the
%   FilterLength property.
%
%   Y = step(H,X,XT) uses XT as the training samples to calculate the
%   beamforming weights when you set the TrainingInputPort property to
%   true. XT is an MxN matrix where N is the number of elements of the
%   sensor array. M must be larger than the FIR filter length specified in
%   the FilterLength property.
%
%   Y = step(H,X,ANG) uses ANG as the beamforming direction, when you
%   set the DirectionSource property to 'Input port'. ANG is a length-2
%   column vector in the form of [AzimuthAngle; ElevationAngle] (in
%   degrees). The azimuth angle must be between [-180 180] and the
%   elevation angle must be between [-90 90].
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [Y,W] = step(H,X,XT,ANG)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   TimeDelayLCMVBeamformer methods:
%
%   step     - Perform time delay LCMV beamforming (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create phase shift beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   TimeDelayLCMVBeamformer properties:
%
%   SensorArray           - Sensor array
%   PropagationSpeed      - Signal propagation speed
%   SampleRate            - Sample rate
%   FilterLength          - FIR filter length
%   Constraint            - Constraint matrix
%   DesiredResponse       - Desired response vector
%   DiagonalLoadingFactor - Diagonal loading factor
%   TrainingInputPort     - Enable training data input
%   DirectionSource       - Source of beamforming direction
%   Direction             - Beamforming direction
%   WeightsOutputPort     - Enable weights output
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Apply a time delay LCMV beamformer to an 11-element array. The 
%   %   incident angle of the signal is -50 degrees in azimuth and 30 
%   %   degrees in elevation.
%
%   % signal simulation
%   array = phased.ULA('NumElements',11,'ElementSpacing',0.04);
%   array.Element.FrequencyRange = [20 20000];
%   fs = 8e3; t = 0:1/fs:0.3;
%   x = chirp(t,0,1,500);
%   c = 340; % Wave propagation speed (m/s)
%   sigcol = phased.WidebandCollector('Sensor',array,...
%            'PropagationSpeed',c,'SampleRate',fs,'ModulatedInput',false);
%   incidentAngle = [-50; 30];
%   x = sigcol(x.',incidentAngle);
%   noise = 0.2*randn(size(x));
%   rx = x+noise;
%
%   % beamforming
%   beamformer = phased.TimeDelayLCMVBeamformer('SensorArray',array,...
%         'PropagationSpeed',c,'SampleRate',fs,'FilterLength',5,...
%         'Direction',incidentAngle);
%   beamformer.Constraint = kron(eye(5),ones(11,1));
%   beamformer.DesiredResponse = eye(5, 1);
%   y = beamformer(rx);
%   plot(t,rx(:,6),'r:',t,y);
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
%
%   See also phased, phased.TimeDelayBeamformer, phased.FrostBeamformer,
%   phased.LCMVBeamformer.

%   Copyright 2009-2016 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] Otis Frost, An Algorithm For Linearly Constrained Adaptive Array
%       Processing, Proceedings of the IEEE, Vol 60, Issue 8, August 1972


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %Constraint   Constraint matrix
        %   Specify the constraint matrix used for time delay LCMV
        %   beamformer as an MxK matrix. Each column of the matrix is a
        %   constraint and M is the degrees of freedom of the beamformer.
        %   For a time delay LCMV beamformer, the degrees of freedom is the
        %   product of the number of elements in the sensor array and the
        %   order of the FIR filter specified in the FilterLength property.
        %   The default value of this property is [1;1].
        Constraint = [1;1];
        %DesiredResponse    Desired response vector
        %   Specify the desired response used for time delay LCMV
        %   beamformer as a length-K column vector where K is the number of
        %   constraints in the Constraint property. Each element in the
        %   vector defines the desired response of the constraint specified
        %   in the corresponding column of the Constraint property. The
        %   default value of this property is 1, which is equivalent to a
        %   distortionless response.
        DesiredResponse = 1;
    end
    
    methods
        function set.Constraint(obj,val)
            validateattributes( val, { 'double','single' }, { '2d', 'finite', 'nonempty' }, '', 'Constraint');
            cond = any(all(val==0));
            if cond
                coder.internal.errorIf(cond,'phased:beamformer:SMI:expectedNonZero');
            end
            obj.Constraint = val;
        end
        function set.DesiredResponse(obj,val)
            validateattributes( val, { 'double','single' }, { 'column', 'finite', 'nonempty' }, '', 'DesiredResponse');
            obj.DesiredResponse = val;
        end
    end
    
    methods
        function obj = TimeDelayLCMVBeamformer(varargin)
            obj = obj@phased.internal.AbstractTimeDomainSMIBeamformer(varargin{:});
        end
    end

    methods (Access = 'protected')
        
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractTimeDomainSMIBeamformer(obj);
            dim = getNumElements(obj.SensorArray)*obj.FilterLength;
            cond = size(obj.Constraint,1) ~= dim;
            if cond
                coder.internal.errorIf(cond,'phased:TimeDomainBeamformer:InvalidConstraint', dim);
            end
        end
            

        function setupImpl(obj, varargin)
            setupImpl@phased.internal.AbstractTimeDomainSMIBeamformer(obj, varargin{:});
            obj.cWeightsEstimator.DesiredResponse = obj.DesiredResponse;
            obj.cWeightsEstimator.Constraint = obj.Constraint;
        end

        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
    end

    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractTimeDomainBeamformer('array');
            props = {...
                     'SampleRateFromInputCheckbox',...
                     'SampleRate',...
                     'FilterLength',...
                     'Constraint',...
                     'DesiredResponse',...
                     'DiagonalLoadingFactor',...
                     'TrainingInputPort',...
                     'DirectionSource',...
                     'Direction',...
                     'WeightsOutputPort'};
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:TimeDelayLCMVBeamformerTitle')),...
                'Text',getString(message('phased:library:block:TimeDelayLCMVBeamformerDesc')));
        end        
    end
    methods (Access = protected) %for Simulink    
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Time Delay LCMV\nBeamformer');
        end        
        function flag = isInputSizeLockedImpl(obj,index) 
            if index == 1
                flag = false;
            elseif index == 2 && obj.TrainingInputPort
                flag = false;
            else
                flag = true;
            end
        end
    end        
end

