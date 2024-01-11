classdef (Sealed,StrictDefaults) LCMVBeamformer < phased.internal.AbstractVarSizeEngine  & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%LCMVBeamformer     Narrowband LCMV beamformer
%   H = phased.LCMVBeamformer returns a linear constraint minimum variance
%   (LCMV) beamformer System object, H. This object performs narrowband
%   LCMV beamforming on the received signal.
%
%   H = phased.LCMVBeamformer(Name,Value) creates an LCMV beamformer
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) performs LCMV beamforming on the input X, and returns the
%   beamformed output in Y. X is an MxN matrix where N is the number of
%   elements of the sensor array. Y is a length-M column vector.
%
%   [Y,W] = step(H,X) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is a
%   length-N column vector where N is the number of elements in the sensor
%   array.
%
%   Y = step(H,X,XT) uses XT as the training samples to calculate the
%   beamforming weights when you set the TrainingInputPort property to
%   true. XT is a PxN matrix where N is the number of elements of the
%   sensor array. P must be greater than N.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [Y,W] = step(H,X,XT)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   LCMVBeamformer methods:
%
%   step     - Perform LCMV beamforming (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create phase shift beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   LCMVBeamformer properties:
%
%   Constraint            - Constraint matrix
%   DesiredResponse       - Desired response vector
%   DiagonalLoadingFactor - Diagonal loading factor
%   TrainingInputPort     - Enable training data input
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
%   %   Apply an LCMV beamformer to a 5-element ULA, preserving the signal
%   %   from the desired direction.
%
%   % simulate signal
%   t = (0:1000)';
%   x = sin(2*pi*0.01*t);
%   c = 3e8; Fc = 3e8;
%   incidentAngle = [45; 0];
%   array = phased.ULA('NumElements',5);
%   x = collectPlaneWave(array,x,incidentAngle,Fc,c);
%   noise = 0.1*(randn(size(x)) + 1j*randn(size(x)));
%   rx = x + noise;
%
%   % beamforming
%   steer = phased.SteeringVector('SensorArray',array,...
%       'PropagationSpeed',c);
%   lcmv = phased.LCMVBeamformer;
%   lcmv.Constraint = steer(Fc,incidentAngle);
%   lcmv.DesiredResponse = 1;
%   y = lcmv(rx);
%   plot(t,real(rx(:,3)),'r:',t,real(y));
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
%
%   See also phased, phased.PhaseShiftBeamformer, phased.MVDRBeamformer,
%   phased.TimeDelayLCMVBeamformer.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
   
    properties (Nontunable)
        %Constraint   Constraint matrix
        %   Specify the constraint matrix used for LCMV beamforming as an
        %   NxK matrix. Each column of the matrix is a constraint and N is
        %   the number of elements in the sensor array or the number of
        %   subarrays if the sensor array contains subarrays. The number of
        %   constraints, K, should be less than or equal to N. The default
        %   value of this property is [1;1].
        Constraint = complex([1;1]);
        %DesiredResponse    Desired response vector
        %   Specify the desired response used for LCMV beamforming as a
        %   length-K column vector where K is the number of constraints in
        %   the Constraint property. Each element in the vector defines the
        %   desired response of the constraint specified in the
        %   corresponding column of the Constraint property. The default
        %   value of this property is 1, i.e., distortionless response.
        DesiredResponse = 1;
    end

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
        %WeightsOutputPort    Enable weights output
        %   Set this property to true to output the weights used in the
        %   beamformer. Set this property to false to not output the
        %   weights. The default value of this property is false.
        WeightsOutputPort = false;
    end
    
    properties (Access = private, Nontunable)
        cWeightsEstimator;
        pSensorArrayNumElements;
    end
    
    methods
        function set.Constraint(obj,val)
            validateattributes( val, { 'double','single' }, { '2d', 'finite', 'nonempty' }, '', 'Constraint');
            cond =  any(all(val==0));
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:beamformer:SMI:expectedNonZero');
            end
            cond = size(val,2) > size(val,1);
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:LCMVBeamformer:TooManyConstraints');
            end
            obj.Constraint = val;
        end
        function set.DesiredResponse(obj,val)
            validateattributes( val, { 'double','single' }, { 'column', 'finite', 'nonempty' }, '', 'DesiredResponse');
            obj.DesiredResponse = val;
        end
        function set.DiagonalLoadingFactor(obj,val)
            validateattributes( val, { 'double','single' }, { 'scalar', 'nonnegative', 'real', 'finite', 'nonempty' }, '', 'DiagonalLoadingFactor');
            obj.DiagonalLoadingFactor = val;
        end
    end

    methods
        function obj = LCMVBeamformer(varargin)
            setProperties(obj, nargin, varargin{:});
        end

    end

    methods (Access = protected)

        function validatePropertiesImpl(obj)
            cond =  size(obj.Constraint,2) ~= size(obj.DesiredResponse,1);
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:LCMVBeamformer:ConstraintMismatch');
            end
        end

        function setupImpl(obj,x,varargin) 
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
            obj.cWeightsEstimator = phased.internal.BlockSMIWeightsEstimator(...
                'DesiredResponse',obj.DesiredResponse,...
                'Constraint',obj.Constraint);
            obj.pSensorArrayNumElements = size(obj.Constraint,1);
            processTunedPropertiesImpl(obj);
            
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = false;
            if index == 1
                flag = false;
            else
                if obj.TrainingInputPort && (index == 2)
                    flag = false;
                end
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~) %#ok
            flag = false;
        end
        
        function releaseImpl(obj)
            release(obj.cWeightsEstimator);
        end
        
        function resetImpl(obj)
            reset(obj.cWeightsEstimator);
        end
        
        function validateInputsImpl(obj,x,varargin)            
            validateInputSignal(obj,x,'X');
            if obj.TrainingInputPort
                xt = varargin{1};
                validateInputSignal(obj,xt,'XT');
            end
        end
        
        function processTunedPropertiesImpl(obj)
            set(obj.cWeightsEstimator,'DiagonalLoadingFactor', ...
                obj.DiagonalLoadingFactor);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.cWeightsEstimator = saveobj(obj.cWeightsEstimator);
                s.pSensorArrayNumElements = obj.pSensorArrayNumElements;
            end
        end
        
        function s = loadSubObjects(obj,s)
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cWeightsEstimator = phased.internal.BlockSMIWeightsEstimator.loadobj(s.cWeightsEstimator);                
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
        
        function flag = isInputSizeLockedImpl(obj,index) 
            if index == 1
                flag = false;
            elseif index == 2 && obj.TrainingInputPort
                flag = false;
            else
                flag = true;
            end
        end
        
        function [y, w] = stepImpl(obj,x,varargin)
            if obj.TrainingInputPort
                xt = cast(varargin{1},class(x));
            else
                xt = x;
            end
            w = step(obj.cWeightsEstimator,xt);
            y = x*conj(w);

        end

        function num = getNumOutputsImpl(obj) 
            num = 1;
            if obj.WeightsOutputPort
                num = num+1;
            end
        end
        
        function num = getNumInputsImpl(obj)
            num = 1;
            if obj.TrainingInputPort
                num = num+1;
            end
        end        
    end
    
    methods (Access=private)
        function validateInputSignal(obj,x,str_x)
            cond =  ~(isa(x,'float'));
            if cond
                coder.internal.errorIf(cond, ...
                                       'MATLAB:system:invalidInputDataType',str_x,'float');
            end
            cond =  ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                                       'MATLAB:system:inputMustBeMatrix',str_x);
            end
            sz_x = size(x);
            cond =  sz_x(2) ~= size(obj.Constraint,1);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:beamformer:InvalidDataDimension', str_x, size( obj.Constraint, 1 ));
            end
            validateNumChannels(obj,x)
        end
    end
    methods (Access = protected) %for Simulink
        function varargout = getOutputNamesImpl(~)
            varargout = {'Y','W'};
        end
        function varargout = getInputNamesImpl(~)
            varargout = {'X','XT'};
        end
        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            varargout{1} = [szX(1) 1];
            if obj.WeightsOutputPort
               varargout{2} = [szX(2) 1];
            end
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj, 1);
            varargout{2} = true;
        end
        function varargout = getOutputDataTypeImpl(obj)
            dt = propagatedInputDataType(obj,1);
            varargout = {dt dt};
        end
        function varargout = isOutputComplexImpl(obj) %#ok<MANU>
            varargout = {true, true};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('LCMV\nBeamformer');
        end
    end
    methods (Static,Hidden,Access=protected)
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:LCMVBeamformerTitle')),...
              'Text',getString(message('phased:library:block:LCMVBeamformerDesc')));
      end        
    end
    
end


