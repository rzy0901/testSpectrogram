classdef (Sealed,StrictDefaults) FrostBeamformer < phased.internal.AbstractTimeDomainSMIBeamformer
%FrostBeamformer    Frost beamformer
%   H = phased.FrostBeamformer returns a Frost beamformer System object, H.
%   This object performs Frost beamforming on the received signal. 
%
%   H = phased.FrostBeamformer(Name,Value) creates a Frost beamformer
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   The beamforming algorithm is proposed by Frost. It can be considered
%   the time domain counterpart of the minimum variance distortionless
%   response (MVDR) beamformer. The received samples at each element in the
%   sensor array are fed into an FIR filter. The algorithm first steers the
%   array to the beamforming direction and then applies the FIR filter to
%   the output of each sensor to achieve the distortionless response
%   constraint.
%
%   Step method syntax:
%
%   Y = step(H,X) performs Frost beamforming on the input X, and returns
%   the beamformed output in Y. X is an MxN matrix where N is the number of
%   elements of the sensor array. Y is a length-M column vector. M must be
%   larger than the FIR filter length specified in the FilterLength
%   property.
%
%   [Y,W] = step(H,X) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is a
%   length-L column vector where L is the degrees of freedom of the
%   beamformer. For a Frost beamformer, L is given by the product of the
%   number of elements in the sensor array specified by the SensorArray
%   property and the FIR filter length specified by the FilterLength
%   property.
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
%   FrostBeamformer methods:
%
%   step     - Perform Frost beamforming (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create phase shift beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   FrostBeamformer properties:
%
%   SensorArray           - Sensor array
%   PropagationSpeed      - Signal propagation speed
%   SampleRate            - Sample rate
%   FilterLength          - FIR filter length
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
%   %   Apply a Frost beamformer to an 11-element array. The incident 
%   %   angle of the signal is -50 degrees in azimuth and 30 degrees in
%   %   elevation.
%
%   % signal simulation
%   array = phased.ULA('NumElements',11,'ElementSpacing',0.04);
%   array.Element.FrequencyRange = [20 20000];
%   fs = 8e3; t = 0:1/fs:0.3;
%   x = chirp(t,0,1,500);
%   c = 340; % Wave propagation speed (m/s)
%   collector = phased.WidebandCollector('Sensor',array,...
%               'PropagationSpeed',c,'SampleRate',fs,...
%               'ModulatedInput',false);
%   incidentAngle = [-50; 30];
%   x = collector(x.',incidentAngle);
%   noise = 0.2*randn(size(x));
%   rx = x+noise;
%
%   % beamforming
%   frost = phased.FrostBeamformer('SensorArray',array,...
%         'PropagationSpeed',c,'SampleRate',fs,...
%         'Direction',incidentAngle,'FilterLength',5);
%   y = frost(rx);
%   plot(t,rx(:,6),'r:',t,y);
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
%
%   See also phased, phased.TimeDelayBeamformer, phased.MVDRBeamformer,
%   phased.TimeDelayLCMVBeamformer, phased.SubbandPhaseShiftBeamformer.

%   Copyright 2009-2016 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] Otis Frost, An Algorithm For Linearly Constrained Adaptive Array
%       Processing, Proceedings of the IEEE, Vol 60, Issue 8, August 1972


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Hidden)
        hasSampleRateSet = true
    end
    
    methods
        function obj = FrostBeamformer(varargin)
            obj = obj@phased.internal.AbstractTimeDomainSMIBeamformer(varargin{:});
        end

    end

    methods (Access = 'protected')
        
        function setupImpl(obj, x, varargin)
            coder.extrinsic('phased.FrostBeamformer.initConstr');
            setupImpl@phased.internal.AbstractTimeDomainSMIBeamformer(obj,x,varargin{:});
            obj.cWeightsEstimator.DesiredResponse = eye(obj.FilterLength, 1);
            obj.cWeightsEstimator.Constraint = ...
                coder.const(obj.initConstr(obj.FilterLength,obj.pDOF));
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
                'DiagonalLoadingFactor',...
                'TrainingInputPort',...
                'DirectionSource',...
                'Direction',...
                'WeightsOutputPort'};
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:FrostBeamformerTitle')),...
                'Text',getString(message('phased:library:block:FrostBeamformerDesc')));
        end    
    end
    methods (Static,Hidden)
      function constr = initConstr(filtLen,DOF)
          constr = kron(eye(filtLen),ones(DOF,1));
      end
    end
    methods (Access = protected) %for Simulink
        function flag = isInputSizeLockedImpl(obj,index) 
            if index == 1
                flag = false;
            elseif index == 2 && obj.TrainingInputPort
                flag = false;
            else
                flag = true;
            end
        end
        function varargout = isOutputComplexImpl(obj)
            varargout = { propagatedInputComplexity(obj,1), propagatedInputComplexity(obj,1)};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Frost\nBeamformer');
        end
    end
end

