classdef (Sealed, Hidden, StrictDefaults)STAPTraining < phased.internal.AbstractTrainingInterface
%This class is for internal use only. It may be removed in the future.

%STAPTraining   STAP training data retriever
%   H = phased.internal.STAPTraining creates a STAP training data retriever
%   System object, H. The object extracts the training data for STAP
%   algorithm from input data.
%
%   H = phased.internal.STAPTraining(Name,Value) creates a training data
%   retriever object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN)
%
%   Step method syntax:
%
%   Y = step(H, X, CUTIDX) retrieves the training data from the input X. X
%   must be a matrix whose columns are space time snapshots. Each column of
%   X is one cell. The training data is used to train the STAP processor
%   for the cell specified by the index CUTIDX.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   STAPTraining methods:
%
%   step     - Retrieve training data (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a training data object with same property values
%   isLocked - Locked status (logical)
%
%   STAPTraining properties:
%
%   NumGuardCells    - Number of guarding cells
%   NumTrainingCells - Number of training cells
%
%   % Example:
%   %   Retrieve training data toward 10th cell of the input data.
%   Ht = phased.internal.STAPTraining; x = rand(10,20);
%   tdata = step(Ht,x,10);

%   Copyright 2009-2016 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    methods
        function obj = STAPTraining(varargin)
            obj@phased.internal.AbstractTrainingInterface(varargin{:});
        end
    end

    methods (Access = protected)
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function validateInputsImpl(~,x,cutidx)
            cond = ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','double');
            end
            cond = ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeMatrix','X');
            end
            cond = ~isa(cutidx,'double');
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','CUTIDX','double');
            end
            cond = ~isscalar(cutidx);
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeScalar','CUTIDX');
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index) %#ok<MANU>
            flag = false;  % index == 1
            if (index == 2)
                flag  =true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~) %#ok<MANU>
            flag = false;  % (index == 1) 
        end
        
        function trnData = stepImpl(obj,datavec,CUTIdx)
           
            % Assume input data are valid as this is a protected method
            
            NumCells = size(datavec,2);
            
            [FrontTrainingCellIdx, RearTrainingCellIdx] = ...
                get1DTrainingIdx(obj,NumCells,CUTIdx);
                           
            % Retrieve data
            trnData = datavec(:,...
                [FrontTrainingCellIdx,RearTrainingCellIdx]);
        end
    end

end




