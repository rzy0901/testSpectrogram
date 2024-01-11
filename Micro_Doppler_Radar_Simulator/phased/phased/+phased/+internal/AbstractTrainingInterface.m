classdef (Hidden) AbstractTrainingInterface < matlab.System
%This class is for internal use only. It may be removed in the future.

%AbstractTrainingInterface   Define the AbstractTrainingInterface
%class.

%   Copyright 2009-2011 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable) 

        %NumGuardCells  Number of guard cells
        %   Specify the number of guard cells used in the training as an
        %   even integer. It specifies the total number of cells on both
        %   sides of the cell under test. The default value of this
        %   property is 2, i.e., there is one guard cell at both the front
        %   and back of the cell under test.
        NumGuardCells = 2;
        %NumTrainingCells   Number of training cells
        %   Specify the number of training cells used in the training as an
        %   even integer. Whenever possible, the training cells are equally
        %   divided before and after the cell under test. The default value
        %   of this property is 2, i.e., there is one training cell at both
        %   the front and back of the cell under test.
        NumTrainingCells = 2;
    end
    
    methods

        function set.NumGuardCells(obj,val)
            validateattributes(val,{'double','single'},...
                {'finite','scalar','nonnegative','even'},...
                'phased.internal.TrainingInterface','NumGuardCells');
            obj.NumGuardCells = val;
        end
        
        function set.NumTrainingCells(obj,val)
            validateattributes(val,{'double','single'},...
                {'finite','scalar','positive','even'},...
                'phased.internal.TrainingInterface','NumTrainingCells');
            obj.NumTrainingCells = val;
        end
    end
    
    methods (Access = protected)
        function obj = AbstractTrainingInterface(varargin)
            setProperties(obj, nargin, varargin{:}, 'NumGuardCells', 'NumTrainingCells');
        end
    end

    methods (Access = protected)
        
        function [FrontTrainingCellIdx, RearTrainingCellIdx] = get1DTrainingIdx(obj,NumCells,CUTIdx)
        %get1DTrainingIdx Calculate the training cells indices 
        %   [FrontIndices, RearIndices]=get1DTrainingIdx(htr,Ncells,CUTIdx)
        %   returns the training cell indices in front of the cell under
        %   test (FrontIndices) and after the cell under test
        %   (RearIndices). The calculation is based on number of total
        %   cells NumCells, the cell under test index CUTIdx and the
        %   settings in the training object htr, such as NumGuardCells and
        %   NumTrainingCells.  
        %
        %   Examples: 
        %       % Assume total number of cells is 100, and cell under test
        %       % cell 40.  Let's also assume the NumGuardCells is 4.  Then
        %       % here is what we do:
        %       %   NumTrainingCells     FrontIndex     ReadIndex
        %       %       60                8:37 (30)       43:72 (30)
        %       %       80                1:37 (37)       43:85 (43)
        %       %      100                error out

            OneSideNumGuardCells = obj.NumGuardCells/2;
            NumAvailableCells = NumCells - 2*OneSideNumGuardCells - 1;
            NumRefCells = obj.NumTrainingCells;
            halfNumRefCells = floor(NumRefCells/2);
            guardFrontIdx = CUTIdx - OneSideNumGuardCells;  % front end of guard cells
            guardRearIdx = CUTIdx + OneSideNumGuardCells;   % rear end of guard cells
            
            cond = NumRefCells > NumAvailableCells;
            if cond
                coder.internal.errorIf(cond,...
                   'phased:phased:internal:AbstractTrainingInterface:OutOfBound', NumAvailableCells);
            end
            
            
            NumFrontRemCells = max(guardFrontIdx-1,0);       % number of training cells in front of guard cells
            NumRearRemCells = max(NumCells-guardRearIdx,0);  % number of training cells behind guard cells
            
            if NumFrontRemCells < NumRearRemCells  % front remaining cells are less
                NumFrontTrainingCells = min(halfNumRefCells,...
                    NumFrontRemCells);  % Use half of training cells or all available cells in front
                NumRearTrainingCells = min(NumRefCells-NumFrontTrainingCells,...
                    NumRearRemCells);   % Make up the rest or take all rear cells if not enough
            else   % rear remaining cells are less
                NumRearTrainingCells = min(halfNumRefCells,...
                    NumRearRemCells);  % Use half of training cells or all available rear cells
                NumFrontTrainingCells = min(NumRefCells-NumRearTrainingCells,...
                    NumFrontRemCells); % Make up the rest or take all front cells if not enough
            end
            
            
            assert(NumFrontTrainingCells <  NumRefCells+1); %Help compiler compute upperbounds
            FrontTrainingCellIdx = guardFrontIdx + (-NumFrontTrainingCells:-1);
            assert(NumRearTrainingCells <  NumRefCells+1);
            RearTrainingCellIdx = guardRearIdx + (1:NumRearTrainingCells);
            
        end
        
        function flag = isInputSizeLockedImpl(~,~)
            flag = false;
        end

    end

end




