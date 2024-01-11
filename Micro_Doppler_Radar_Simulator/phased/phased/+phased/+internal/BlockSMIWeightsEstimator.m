classdef (Hidden, Sealed, StrictDefaults) BlockSMIWeightsEstimator < phased.internal.AbstractVarSizeEngine
%This class is for internal use only. It may be removed in the future.

%BlockSMIWeightsEstimator Estimate weights using sample matrix inversion
%                          (SMI)
%   H = phased.internal.BlockSMIWeightsEstimator creates a block SMI
%   weights estimator System object, H. The object estimates the weights
%   using the SMI algorithm. In addition, the weights are estimated using
%   each block of the training signal.
%
%   H = phased.internal.BlockSMIWeightsEstimator(Name,Value) creates a
%   block SMI weights estimator object, H, with the specified property Name
%   set to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   W = step(H,X) estimates the weights W using the SMI algorithm with the
%   training data specified in X. X must be a N-column matrix where N is
%   the degrees of freedom defined in the Constraint property. W is a
%   length-N column vector.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   BlockSMIWeightsEstimator methods:
%
%   step     - Calculate SMI weights (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create SMI weights estimator object with same property
%              values
%   isLocked - Locked status (logical)
%
%   BlockSMIWeightsEstimator properties:
%
%   Constraint            - Constraint matrix
%   DesiredResponse       - Desired response vector
%   DiagonalLoadingFactor - Diagonal loading factor
%
%   % Example:
%   %   Estimate the MVDR weights of a 5-element ULA. Assume that the
%   %   beamforming direction is 10 degrees of azimuth and 0 degrees of
%   %   elevation.
%   ha = phased.ULA('NumElements',5);
%   hstv = phased.SteeringVector('SensorArray',ha);
%   fc = 3e8; ang = [10 0];
%   hw = phased.internal.BlockSMIWeightsEstimator(...
%           'Constraint',step(hstv,fc,ang));
%   x = phased.internal.cwgn(1,100,5);
%   w = step(hw,x)
%
%   See also phased.

%   Copyright 2010-2016 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Constraint   Constraint matrix
        %   Specify the constraint matrix used for SMI weights calculation
        %   as an MxN matrix. Each column of the matrix is a constraint and
        %   M is the degrees of freedom of the weights, i.e., the length of
        %   the weights. The default value of this property is 1. This
        %   property is tunable.
        Constraint = complex(1);
        %DesiredResponse    Desired response vector
        %   Specify the desired response used for SMI weights calculation
        %   as a length-N column vector where N is the number of
        %   constraints in the Constraint property. Each element in the
        %   vector defines the desired response of the constraint specified
        %   in the corresponding column of the Constraint property. The
        %   default value of this property is 1, i.e., distortionless
        %   response. This property is tunable.
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

    methods

        function set.DiagonalLoadingFactor(obj,val)
            validateattributes( val, { 'double','single'}, { 'scalar', 'nonnegative', 'real', 'finite', 'nonempty' }, '', 'DiagonalLoadingFactor');
            obj.DiagonalLoadingFactor = val;
        end

        function set.DesiredResponse(obj,val)
            validateattributes(val,{'double','single'},{'column','finite','nonempty'},...
                'phased.beamformer.SMI','DesiredResponse');
            obj.DesiredResponse = val;
        end

        function set.Constraint(obj,val)
            validateattributes(val,{'double','single'},{'2d','finite','nonempty'},...
                'phased.beamformer.SMI','Constraint');
            cond = any(all(val==0));
            if cond
                coder.internal.errorIf(cond,...
                     'phased:beamformer:SMI:expectedNonZero');            
            end
            obj.Constraint = val;
        end

    end

    methods

        function obj = BlockSMIWeightsEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

    end

    methods (Access = 'protected')

        function validatePropertiesImpl(obj)
            cond = size(obj.Constraint,1) < size(obj.Constraint,2);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:beamformer:SMI:TooManyConstraints');
            end
            cond = size(obj.Constraint,2) ~= size(obj.DesiredResponse,1);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:beamformer:SMI:ConstraintMismatch');
            end
        end

        function validateInputsImpl(obj,x)
            cond =  ~isa(x,'float');
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:invalidInputDataType','X','float');
            end
            sz_x = size(x);
            cond = ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeMatrix','X');
            end
            cond =  sz_x(2) ~= size(obj.Constraint,1);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:beamformer:SMI:TrainingDimensionMismatch', size( obj.Constraint, 1 ));
            end
            if sz_x(1) <= 2*sz_x(2) && isempty(coder.target)
                warning(message('phased:beamformer:SMI:InsufficientTrainingSample'));
            end
            validateNumChannels(obj,x)
           
        end
        
        function flag = isInputComplexityLockedImpl(obj,~)  %#ok
            flag = false;
        end

        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok
            flag = false;
        end
        
        function setupImpl(obj,x)
            sz_x = size(x);
            obj.pNumInputChannels = sz_x(2);
            % obj.pValidatedNumInputChannels = sz_x(2);
        end
        
        function w = stepImpl(obj,x)
            classtouse=class(x);
            w = phased.internal.lcmvweights(x,cast(obj.Constraint,classtouse),...
                cast(obj.DesiredResponse,classtouse),cast(obj.DiagonalLoadingFactor,classtouse));

        end

        function flag = isInputSizeLockedImpl(~,~)
            flag = false;
        end

    end
    
end
