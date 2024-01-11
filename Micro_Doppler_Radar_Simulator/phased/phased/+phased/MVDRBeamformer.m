classdef (Sealed,StrictDefaults) MVDRBeamformer < phased.internal.AbstractNarrowbandBeamformer
%MVDRBeamformer Narrowband MVDR beamformer
%   H = phased.MVDRBeamformer creates a minimum variance distortionless
%   response (MVDR) beamformer System object, H. This object performs MVDR
%   beamforming on the received signal.
%
%   H = phased.MVDRBeamformer(Name,Value) creates an MVDR beamformer
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) performs MVDR beamforming on the input X, and returns the
%   beamformed output in Y. X is an MxN matrix where N is the number of
%   subarrays if SensorArray contains subarrays, or the number of elements
%   otherwise. Y is an MxL matrix where L is the number of beamforming
%   directions.
%
%   If you set the TrainingInputPort to false, then X is used to do the
%   training and M must be larger than N. If you set the TrainingInputPort
%   to true, then M can be any positive integer.
%
%   [Y,W] = step(H,X) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is an
%   NxL matrix.
%
%   Y = step(H,X,XT) uses XT as the training samples to calculate the
%   beamforming weights when you set the TrainingInputPort property to
%   true. XT is a PxN matrix where P must be larger than N.
%
%   Y = step(H,X,ANG) uses ANG as the beamforming direction, when you
%   set the DirectionSource property to 'Input port'. ANG is a 2-row
%   matrix whose columns are in the form of [AzimuthAngle; ElevationAngle]
%   (in degrees). The azimuth angle must be between [-180 180] and the
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
%   MVDRBeamformer methods:
%
%   step     - Perform MVDR beamforming (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create phase shift beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   MVDRBeamformer properties:
%
%   SensorArray           - Sensor array
%   PropagationSpeed      - Signal propagation speed
%   OperatingFrequency    - Operating frequency
%   DiagonalLoadingFactor - Diagonal loading factor
%   TrainingInputPort     - Enable training data input
%   DirectionSource       - Source of beamforming direction
%   Direction             - Beamforming direction
%   NumPhaseShifterBits   - Number of bits in phase shifters
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
%   %    Apply an MVDR beamformer to a 5-element ULA. The incident angle 
%   %    of the signal is 45 degrees in azimuth and 0 degree in elevation.
%
%   % signal simulation
%   t = (0:1000)';
%   x = sin(2*pi*0.01*t);
%   c = 3e8; Fc = 3e8;
%   incidentAngle = [45; 0];
%   array = phased.ULA('NumElements',5);
%   x = collectPlaneWave(array,x,incidentAngle,Fc,c);
%   noise = 0.1*(randn(size(x)) + 1j*randn(size(x)));
%   rx = x+noise;
%
%   % beamforming
%   beamformer = phased.MVDRBeamformer('SensorArray',array,...
%           'OperatingFrequency',Fc,'Direction',incidentAngle,...
%           'WeightsOutputPort',true,'PropagationSpeed',c);
%   [y,w] = beamformer(rx);
%   plot(t,real(rx(:,3)),'r:',t,real(y));
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
%
%   % plot response pattern
%   figure; pattern(array,Fc,-180:180,0,'PropagationSpeed',c,'Weights',w);
%
%   See also phased, phased.PhaseShiftBeamformer, phased.LCMVBeamformer,
%   phased.FrostBeamformer.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties
        %DiagonalLoadingFactor  Diagonal loading factor
        %   Specify the diagonal loading factor as a positive scalar. The
        %   default value of this property is 0. Diagonal loading is a
        %   technique used to achieve robust beamforming performance,
        %   especially when the sample support is small. This property is
        %   tunable.
        DiagonalLoadingFactor = 0;
    end

    properties (Nontunable, Logical) 
        %TrainingInputPort  Enable training data input
        %   Set this property to true to add input to specify the
        %   additional training data. Set this property to false to not add
        %   input to specify the training data. The default value of this
        %   property is false. When this property is false, the input
        %   signal itself is used as the training data.
        TrainingInputPort = false;
    end

    methods
        function set.DiagonalLoadingFactor(obj,val)
            validateattributes( val, { 'double','single' }, { 'scalar', 'nonnegative', 'real', 'finite', 'nonempty' }, '', 'DiagonalLoadingFactor');
            obj.DiagonalLoadingFactor = val;
        end
    end

    methods
        function obj = MVDRBeamformer(varargin)
            obj@phased.internal.AbstractNarrowbandBeamformer(varargin{:});
        end
    end

    methods (Access = protected)

        function flag = isMultipleInputAnglesAllowed(obj) %#ok<MANU>
            flag = true;
        end
        
        function flag = isSubarraySupported(obj) %#ok<MANU>
            flag = true;
        end
        
        function setupImpl(obj,x,varargin)
            setupImpl@phased.internal.AbstractNarrowbandBeamformer(obj,x);
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
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = isInputComplexityLockedImpl@phased.internal.AbstractNarrowbandBeamformer(obj,index);
            if obj.TrainingInputPort
                if index == 2
                    flag = false;
                end
                if strncmpi(obj.DirectionSource,'i',1) && (index == 3)
                    flag = true;
                end
            else
                if strncmpi(obj.DirectionSource,'i',1) && (index == 2)
                    flag = true;
                end
            end
        end

        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractNarrowbandBeamformer(obj);
        end

        function resetImpl(obj)
            resetImpl@phased.internal.AbstractNarrowbandBeamformer(obj);
        end
        
        function validateInputsImpl(obj,x,varargin)
            validateInputsImpl@phased.internal.AbstractNarrowbandBeamformer(obj,x);

            if obj.TrainingInputPort
                xt = varargin{1};                
                validateInputSignal(obj,xt,'XT');
                
                if strncmpi(obj.DirectionSource,'i',1)
                   ang = varargin{2};
                   validateInputAngle(obj,ang);
                end
                
            else
                if strncmpi(obj.DirectionSource,'i',1)
                    ang = varargin{1};
                    validateInputAngle(obj,ang);
                end
                xt = x;
            end

            sz_xt = size(xt);
            if sz_xt(1) < 2*sz_xt(2) && isempty(coder.target)
                warning(message('phased:beamformer:SMI:InsufficientTrainingSample'));
            end

        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractNarrowbandBeamformer(obj);
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractNarrowbandBeamformer(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked && isfield(s,'cWeightsEstimator')
                    s = rmfield(s,'cWeightsEstimator');
                end
                s = rmfield(s,'isLocked');
            end
        end
        
        function loadObjectImpl(obj,s,~)
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function [y, w] = stepImpl(obj,x,varargin)
            
            classtouse = class(x);
            anginputflag = (obj.DirectionSource(1) == 'I'); %Input
            %strncmpi(obj.DirectionSource,'i',1);
            
            if obj.TrainingInputPort
                xt = cast(varargin{1},classtouse);
                if anginputflag
                    ang = varargin{2};
                else
                    ang = obj.Direction;
                end
            else
                xt = x;
                if anginputflag
                    ang = varargin{1};
                else
                    ang = obj.Direction;
                end
            end
            w = complex(zeros(obj.pDOF,size(ang,2),classtouse));
            for m = size(ang,2):-1:1
                Constraint = ...
                    step(obj.cSteeringVector,cast(obj.OperatingFrequency,'double'),cast(ang(:,m),'double'));
                w(:,m) = phased.internal.lcmvweights(xt,Constraint,...
                    1, obj.DiagonalLoadingFactor);
            end
            y = x*conj(w);

        end

        function num = getNumInputsImpl(obj)
            num = getNumInputsImpl@phased.internal.AbstractNarrowbandBeamformer(obj);
            if obj.TrainingInputPort
                num = num + 1;
            end
        end
        
    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractNarrowbandBeamformer('subarray');
        props = {...
          'OperatingFrequency',...
          'DiagonalLoadingFactor',...
          'TrainingInputPort',...
          'DirectionSource',...
          'Direction',...
          'NumPhaseShifterBits',...
          'WeightsOutputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:MVDRBeamformerTitle')),...
              'Text',getString(message('phased:library:block:MVDRBeamformerDesc')));
      end            
    end
    methods (Access = protected) %for Simulink
        function varargout = getInputNamesImpl(obj)
            %Insert XT
            if obj.TrainingInputPort
                [varargout{1:nargout-1}] = getInputNamesImpl@phased.internal.AbstractNarrowbandBeamformer(obj);
                varargout = {varargout{1} 'XT' varargout{2:end}};
            else
                [varargout{1:nargout}] = getInputNamesImpl@phased.internal.AbstractNarrowbandBeamformer(obj);
            end
        end
         function aIdx = getAngInputIdx(obj)
            if obj.TrainingInputPort
                aIdx = 3;
            else
                aIdx = 2;
            end
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('MVDR\nBeamformer');
        end        
    end    
end

