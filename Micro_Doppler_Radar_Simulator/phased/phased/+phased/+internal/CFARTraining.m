classdef (Sealed, Hidden, StrictDefaults)CFARTraining < phased.internal.AbstractTrainingInterface
%This class is for internal use only. It may be removed in the future.

%CFARTraining   CFAR training data retriever
%   H = phased.internal.CFARTraining creates a CFAR training data retriever
%   System object, H. The object extracts training data for CFAR algorithm
%   from input data.
%
%   H = phased.internal.CFARTraining(Name,Value) creates a training data
%   retriever object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN)
%
%   Step method syntax:
%
%   Y = step(H, X, CUTIDX) retrieves the training data from the input X. X
%   must be a matrix whose columns are space time snapshots. Each column of
%   X is one cell. The training data, Y, contains data from all training
%   cells and is used to train the CFAR processor for the cell specified by
%   the index CUTIDX. This syntax applies when you set the
%   CombineTrainingData property to true.
%
%   [YL, YT] = step(H, X, CUTIDX) retrieves the training data from the
%   input X. X must be a matrix whose columns are space time snapshots.
%   Each column of X is one cell. YL contains the data from leading
%   training cells and YT contains the data from the trailing training
%   cells for the cell specified by the index CUTIDX. This syntax applies
%   when you set the CombineTrainingData property to false.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   CFARTraining methods:
%
%   step     - Retrieve training data (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a training data object with same property values
%   isLocked - Locked status (logical)
%
%   CFARTraining properties:
%
%   NumGuardCells       - Number of guarding cells
%   NumTrainingCells    - Number of training cells
%   CombineTrainingData - Output data from all training cells
%
%   % Example:
%   %   Retrieve training data toward 10th cell of the input data.
%   Ht = phased.internal.CFARTraining; x = rand(10,20);
%   tdata = step(Ht,x,10);

%   Copyright 2009-2016 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable, Logical)
        %CombineTrainingData    Output data from all training cells
        %   Set this property to true to output training data from all
        %   training cells. Set this property to false to output training
        %   data separately from leading and trailing training cells. The
        %   default value of this property is true.
        CombineTrainingData = true;
    end

    methods
        function obj = CFARTraining(varargin)
            obj@phased.internal.AbstractTrainingInterface(varargin{:});
        end
    end

    methods (Access = protected)
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function num = getNumOutputsImpl(obj)
            if obj.CombineTrainingData
                num = 1;
            else
                num = 2;
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)  %#ok
            flag = false;
            if index == 1
                flag = false;
            end
            if index == 2
                flag = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~) %#ok
            flag = false;
        end
        
        function validateInputsImpl(~,x,cutidx)
            cond =  ~(isa(x,'float'));
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:invalidInputDataType','X','float');
            end
            cond =  ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:inputMustBeMatrix','X');
            end
            cond =  ~(isa(cutidx,'float'));
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:invalidInputDataType','CUTIDX','float');
            end
            cond =  ~isscalar(cutidx);
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:inputMustBeScalar','CUTIDX');
            end
        end
        
        function varargout = stepImpl(obj,datavec,CUTIdx)
           
            % Assume input data are valid as this is a protected method
            
            NumCells = size(datavec,1);
            
            [FrontTrainingCellIdx, RearTrainingCellIdx] = ...
                get1DTrainingIdx(obj,NumCells,CUTIdx);
                           
            % Retrieve data
            if obj.CombineTrainingData
                varargout{1} = datavec([FrontTrainingCellIdx,RearTrainingCellIdx],:);
            else
                varargout{1} = datavec(FrontTrainingCellIdx,:);
                varargout{2} = datavec(RearTrainingCellIdx,:);
            end
        end
    end

end




