classdef (Sealed,StrictDefaults) CFARDetector2D < phased.internal.AbstractDetector
%CFARDetector2D   Constant false alarm rate detector for matrix data
%   H = CFARDetector2D creates a constant false alarm rate (CFAR) detector
%   System object, H. This object performs CFAR detection on the matrix
%   input data.
%
%   H = phased.CFARDetector2D(Name,Value) creates a CFAR detector object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,IDX) performs the CFAR detection on the real input data X.
%   X is an MxN matrix or an MxNxP array. Each page of X is an independent
%   2-D signal. Detection is performed for each cell specified in the 2xL
%   matrix IDX. IDX must contain positive integers, and each column of IDX
%   specifies the row and column index of a cell under test (CUT). Each
%   cell specified by IDX must correspond to a training region which lies
%   entirely inside the input X.
%
%   Y contains the detection results. The format of Y depends on how the
%   OutputFormat property is set. By default, the OutputFormat property is
%   set to 'CUT result'.
%
%   When the OutputFormat property is set to 'CUT result', Y is an LxP
%   matrix containing the detection results for each cell under test in
%   IDX. Target detections are indicated with a 1, and a 0 denotes the
%   absence of a target in the cell under test. L is the number of cells
%   specified in IDX and P is the number of independent 2-D signals in X.
% 
%   When the OutputFormat property is set to 'Detection index', Y is a DxQ
%   matrix with each column of Y containing the indices for each detection
%   found in X and the rows of Y defining these indices for each dimension
%   of X. D is the number of dimensions of the input data X. When the
%   NumDetectionsSource property is set to 'Auto', Q is the total number of
%   detections found in the data. When the NumDetectionsSource property is
%   set to 'Property', Q is set to the NumDetections property value. By
%   default, NumDetectionsSource is set to 'Auto'. Columns without
%   detections are set to NaN.
%
%   Y = step(H,X,IDX,K) uses K as the threshold factor used to calculate
%   the detection threshold when you set the ThresholdFactor property to
%   'Input port'. K must be a positive scalar.
%
%   [Y,TH] = step(H,X,IDX) returns additional output TH as the detection
%   threshold for each detection reported in Y when you set the
%   ThresholdOutputPort property to true. When the OutputFormat property is
%   set to 'CUT result', TH has the same dimensionality as Y. When the
%   OutputFormat property is set to 'Detection index', TH is a 1xQ row
%   vector, with the threshold value reported in each column of TH
%   corresponding to the detection reported in the same column of Y.
%   Columns without detections are set to NaN. By default, OutputFormat is
%   set to 'CUT result'.
%
%   [Y,N] = step(H,X,IDX) returns additional output N as the estimated
%   noise power for each detection reported in Y when you set the
%   NoisePowerOutputPort property to true. When the OutputFormat property
%   is set to 'CUT result', N has the same dimensionality as Y. When the
%   OutputFormat property is set to 'Detection index', N is a 1xQ row
%   vector, with the noise power value reported in each column of N
%   corresponding to the detection reported in the same column of Y.
%   Columns without detections are set to NaN. By default, OutputFormat is
%   set to 'CUT result'.
%
%   You can combine optional input and output arguments when their enabling
%   properties are set. Optional inputs and outputs must be listed in the
%   same order as the order of the enabling properties. For example,
%
%   [Y,TH,N] = step(H,X,IDX,K)
%
%   The algorithm used in CFARDetector2D is cell averaging CFAR. Detection
%   is performed in three steps. First, a rectangular band of training
%   cells is identified for the input. The values in these training cells
%   are averaged to form the noise estimate. Second, the noise estimate is
%   multiplied by the threshold factor to form the threshold. Finally, the
%   value in the test cell is compared against the threshold to determine
%   whether a target is present or absent. If the value is greater than the
%   threshold, a target is present, and a logical 1 is returned in the
%   corresponding location in Y.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   CFARDetector2D methods:
%
%   step     - Perform CFAR detection (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create 2D CFAR detector object with same property values
%   isLocked - Locked status (logical)
%
%   CFARDetector2D properties:
%
%   Method                - CFAR algorithm
%   GuardBandSize         - Size in cells of the guard region band
%   TrainingBandSize      - Size in cells of the training region band
%   Rank                  - Rank of order statistic
%   ThresholdFactor       - Threshold factor method
%   ProbabilityFalseAlarm - Probability of false alarm
%   CustomThresholdFactor - Custom threshold factor
%   OutputFormat          - Output reporting format
%   ThresholdOutputPort   - Output detection threshold
%   NoisePowerOutputPort  - Output noise power
%   NumDetectionsSource   - Source of the number of detections
%   NumDetections         - Maximum number of detections
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example 1:
%   %   Perform cell averaging CFAR detection on a Gaussian noise array
%   %   with two pages. In this case there are no targets and the 
%   %   probability of false alarm (PFA) will be estimated for each page. 
%   %   The desired PFA is 0.1. Assume that the data is from a square law 
%   %   detector and no pulse integration is performed. Use a training 
%   %   band of 2 cells in width to estimate the noise level. Use a guard 
%   %   band of 1 cell in width to separate the test cell and the training 
%   %   cells. 
%
%   rs = RandStream.create('mt19937ar','Seed',5);
%   cfar = phased.CFARDetector2D('TrainingBandSize',2,'GuardBandSize',1,...
%                   'ProbabilityFalseAlarm',0.1);
%   N = 10;
%   Ntrial = 100;
%   x = 1/sqrt(2)*(randn(rs,N,N,Ntrial)+1i*randn(rs,N,N,Ntrial));
%   xSq = abs(x).^2;
%   % Perform the detection on two cells at the center of the signal.
%   CUTIdx = [5 6; 5 6];
%   dets = cfar(xSq,CUTIdx);
%   Pfa = sum(dets,2)/Ntrial
%
%   % Example 2:
%   %   Perform CFAR detection on data with a 10-pulse noncoherent 
%   %   integration. Estimate the probability of false alarm (PFA) using a 
%   %   custom threshold factor. The desired PFA is 0.001.
%
%   rs = RandStream.create('mt19937ar','Seed',5);
%   cfar = phased.CFARDetector2D('TrainingBandSize',2,'GuardBandSize',1,...
%             'ThresholdFactor','Custom','CustomThresholdFactor',2.35);
%   N = 10;
%   Ntrial = 10000;
%   xSq = 0;
%   for m = 1:10
%     x = 1/sqrt(2)*(randn(rs,N,N,Ntrial)+1i*randn(rs,N,N,Ntrial));
%     xSq = xSq + abs(x).^2;
%   end
%   % Perform the detection on two cells at the center of the signal.
%   CUTIdx = [5 6; 5 6];
%   dets = cfar(xSq,CUTIdx);
%   Pfa = sum(dets,2)/Ntrial
%
%   See also phased, phased.CFARDetector, phased.MatchedFilter, 
%   npwgnthresh.

%   Copyright 2015-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %GuardBandSize   Size in cells of the guard region band
        %   Specify the number of row and column guard cells on each side
        %   of the cell under test as nonnegative integers. The first
        %   element specifies the guard band size along the row dimension,
        %   and the second along the column dimension. Specifying
        %   GuardBandSize as a scalar is equivalent to specifying a vector
        %   with the same value for both dimensions. The default value of
        %   this property is [1 1], indicating that there is a
        %   one-guard-cell-wide region surrounding each cell under test.
        GuardBandSize = [1 1]
        
        %TrainingBandSize   Size in cells of the training region band
        %   Specify the number of row and column training cells on each
        %   side of the cell under test as nonnegative integers. The first
        %   element specifies the training band size along the row
        %   dimension, and the second along the column dimension.
        %   Specifying TrainingBandSize as a scalar is equivalent to
        %   specifying a vector with the same value for both dimensions.
        %   The default value of this property is [1 1], indicating that
        %   there is a one-training-cell-wide region surrounding the guard
        %   band for each cell under test.
        TrainingBandSize = [1 1]
    end
    
    properties (Access=private)
        pMaximumCellRowIndex;
        pMaximumCellColIndex
    end
    
    properties (Access = private)
        pNumPages = 1
    end
    
    properties (Access = private, Logical)
        pSizeInitialized
    end
    
    methods
        function set.GuardBandSize(obj,val)
            validateattributes( val, { 'double','single' }, ...
              { 'finite', 'nonnan', 'vector','nonnegative','integer'}, '', 'GuardBandSize');
            obj.GuardBandSize = val;
        end        
        
        function set.TrainingBandSize(obj,val)
            validateattributes( val, { 'double','single' }, ...
              { 'finite', 'nonnan', 'vector','nonnegative','integer'}, '', 'TrainingBandSize');
            if ~any(val>0)
                error(message('phased:CFARDetector2D:TrainingCellsCannotZero'));
            end
            obj.TrainingBandSize = val;
        end
        
    end
    
    methods
        function obj = CFARDetector2D(varargin)
            obj = obj@phased.internal.AbstractDetector(varargin{:});
        end
    end
    
    methods (Access = protected)
        
        function processInputSizeChangeImpl(obj,X,~,~)
            sz_x = size(X);
            obj.pMaximumCellRowIndex = sz_x(1);
            obj.pMaximumCellColIndex = sz_x(2);
        end

        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractVarSizeEngine(obj);
        end
        
        function setupImpl(obj,X,varargin)
            setupImpl@phased.internal.AbstractDetector(obj);
            sz_x = size(X);
            obj.pMaximumCellRowIndex = sz_x(1);
            obj.pMaximumCellColIndex = sz_x(2);

            if length(sz_x)== 3
                obj.pNumPages = sz_x(3);
                obj.pNumChannels = sz_x(3); % Number of independent frames
            else
                obj.pNumPages = 1;
                obj.pNumChannels = 1; % Number of indpendent frames
            end
            
            obj.pNumInputChannels = getNumChannels(obj,X);
        end
        
        function validateInputsImpl(obj,x,idx,varargin)
            validateInputsImpl@phased.internal.AbstractDetector(obj,x,idx,varargin{:});
            
            cond =( length(size(x))>3) || (length(size(x))==1) || isempty(x); 
            if cond               
                 coder.internal.errorIf(cond, ...
                     'phased:CFARDetector2D:XMustBe2D3DNonEmpty');
            end
            
            cond = ~ismatrix(idx) || ~isequal(size(idx,1),2) || isempty(idx);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:CFARDetector2D:IDXMustTwoRows');
            end
            
            % Make sure that the training region is not larger
            % than the size of X
            TrainingRegionSize = getTrainingRegionSize(obj);
            cond = any(TrainingRegionSize > size(x(:,:,1)));
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:CFARDetector2D:TrainingRegionLargerX');
            end
            
            validateNumPages(obj,x,obj.pNumPages);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractDetector(obj);
            if isLocked(obj)
                s.pMaximumCellRowIndex = obj.pMaximumCellRowIndex;
                s.pMaximumCellColIndex = obj.pMaximumCellColIndex;
                s.pSizeInitialized = obj.pSizeInitialized;
                s.pNumPages = obj.pNumPages;
            end
        end
        
        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function resetImpl(obj)
            obj.pSizeInitialized = false;
        end
                
        function [yout,varargout] = stepImpl(obj,X,Idx,ThFac)
            
            classtoUse = class(X);
            if ~obj.pSizeInitialized
                processInputSizeChangeImpl(obj,X);
                obj.pSizeInitialized = true;
            end
        
            sigdatatypes.validateIndex(Idx(1,:),'step','Row 1 of Idx',...
                {'double','single'},{'<=',obj.pMaximumCellRowIndex});

            sigdatatypes.validateIndex(Idx(2,:),'step','Row 2 of Idx',...
                {'double','single'},{'<=',obj.pMaximumCellColIndex});
              
            if obj.ThresholdFactor(1) ~= 'I'
                ThFac = obj.pFactor;
            else
                validateattributes(ThFac,{'double','single'},{'positive','scalar'},...
                    'step','K');
            end
            ThFac = cast(ThFac,classtoUse);
            
            % Validate CUTS
            validateCUTS(obj,X,Idx);
            
            numCUT = size(Idx,2);
            th = zeros(numCUT,obj.pNumChannels,classtoUse); 
            noise = zeros(numCUT,obj.pNumChannels,classtoUse); 
            trninds = zeros(numCUT,obj.pNumChannels,class(Idx));
            
            % If X is a row vector, add an extra row. This is needed to
            % produce the intended output (column vector) when using a
            % column vector of linear indices. Since the training cells
            % must be in the same row, the output will not be affected.
            if isrow(X)
                X = [X; X];
            end
            
            for m = 1:numCUT
                if obj.Method(1) == 'C' %'CA'
                    if m == 1
                        trninds = cast(get2DTrainingInds(obj,X,Idx(:,m)),class(Idx));
                        % casting to class of Idx because of codegen error
                    else
                        % Use last training indices to improve performance
                        dIdx = Idx(:,m)-Idx(:,m-1);
                        trninds = trninds + dIdx(1) + size(X,1)*dIdx(2);
                    end
                    noisePowEst = mean(X(trninds),1);
                elseif obj.Method(1) == 'S' %'SOCA'
                    if m == 1  
                        [trnindsl,trnindsr] = get2DTrainingInds(obj,X,Idx(:,m));
                    else
                        dIdx = Idx(:,m)-Idx(:,m-1);
                        trnindsl = trnindsl  + dIdx(1) + size(X,1)*dIdx(2);
                        trnindsr = trnindsr  + dIdx(1) + size(X,1)*dIdx(2);
                    end
                    % Form noise power estimate
                    noisePowEst = min(mean(X(trnindsl),1),mean(X(trnindsr),1));
                elseif obj.Method(1) == 'G' %'GOCA'
                    if m == 1  
                        [trnindsl,trnindsr] = get2DTrainingInds(obj,X,Idx(:,m));
                    else
                        dIdx = Idx(:,m)-Idx(:,m-1);
                        trnindsl = trnindsl  + dIdx(1) + size(X,1)*dIdx(2);
                        trnindsr = trnindsr  + dIdx(1) + size(X,1)*dIdx(2);
                    end
                    % Form noise power estimate
                    noisePowEst = max(mean(X(trnindsl),1),mean(X(trnindsr),1));
                else %'OS'
                    if m == 1
                        trninds = cast(get2DTrainingInds(obj,X,Idx(:,m)),class(Idx));
                        % casting to class of Idx because of codegen error
                    else
                        dIdx = Idx(:,m)-Idx(:,m-1);
                        trninds = trninds + dIdx(1) + size(X,1)*dIdx(2);
                    end
                    % Form noise power estimate
                    temp = sort(X(trninds),1,'ascend');
                    noisePowEst = temp(obj.Rank,:);
                end
                % Form the threshold for this CUT
                th(m,:) = noisePowEst * ThFac;
                noise(m,:) = noisePowEst;
            end
            
            % Compute CUT values to the threshold
            linIdx = sub2ind(size(X),Idx(1,:)',Idx(2,:)');
            xCUTinds = bsxfun(@plus,linIdx,(0:obj.pNumChannels-1)*(size(X,1)*size(X,2)));
            
            % Compare each CUT to the corresponding threshold value
            y = X(xCUTinds) > th;
            
            if strcmp(obj.OutputFormat,'CUT result')
                yout = y;
                thout = th;
                noiseout = noise;
            else % Detection index
                numDet = sum(y(:));
                
                if strcmp(obj.NumDetectionsSource,'Auto')
                    % Auto in Simulink sets the output size to the maximum
                    % possible size, which is the product of the number of
                    % cells under test and the number of independent
                    % signals in the input data (i.e. channels)
                    if isInputDataSizePropagated(obj)
                        numOut = numCUT*obj.pNumChannels;
                    else
                        numOut = numDet;
                    end
                else
                    numOut = obj.NumDetections;
                end
                thout = NaN(1,numOut,classtoUse);
                noiseout = NaN(1,numOut,classtoUse);
                
                numAvail = min(numOut,numDet);
                ind = xCUTinds(y)';
                if obj.pNumChannels == 1
                    yout = NaN(2,numOut,classtoUse);
                    [idx1,idx2] = ind2sub(size(X),ind(1:numAvail));
                    detidx = [idx1;idx2];
                else
                    yout = NaN(3,numOut,classtoUse);
                    [idx1,idx2,idx3] = ind2sub(size(X),ind(1:numAvail));
                    detidx = [idx1;idx2;idx3];
                end
                yout(:,1:numAvail) = detidx;
                tmp = th(y);
                thout(1:numAvail) = tmp(1:numAvail);
                tmp = noise(y);
                noiseout(1:numAvail) = tmp(1:numAvail);
            end
            
            iVarg = 1;
            if obj.ThresholdOutputPort
                varargout{iVarg} = thout;
                iVarg = iVarg+1;
            end
            if obj.NoisePowerOutputPort
                varargout{iVarg} = noiseout;
            end
        end
    end
    
    methods (Access = protected)

          function GuardRegionSize = getGuardRegionSize(obj) 
              if isscalar(obj.GuardBandSize)
                  GuardRegionSize = [2*obj.GuardBandSize + 1 ...
                    2*obj.GuardBandSize + 1];
              else
                  GuardRegionSize = 2*obj.GuardBandSize + 1;
              end
          end
       
          function TrainingRegionSize = getTrainingRegionSize(obj)
              GuardRegionSize = getGuardRegionSize(obj); 
              TrainingRegionSize = 2*obj.TrainingBandSize + GuardRegionSize;
          end
          
          function NumTrainingCells = getNumTrainingCells(obj) 
              GuardRegionSize = getGuardRegionSize(obj);
              TrainingRegionSize = getTrainingRegionSize(obj);
              NumTrainingCells = TrainingRegionSize(1)*TrainingRegionSize(2)...
                -GuardRegionSize(1)*GuardRegionSize(2);
          end
    end
    
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
           props = {...
                'Method',...
                'Rank',...
                'GuardBandSize',...
                'TrainingBandSize',...
                'ThresholdFactor',...
                'ProbabilityFalseAlarm',...
                'CustomThresholdFactor',...
                'OutputFormat',...
                'ThresholdOutputPort',...
                'NoisePowerOutputPort',...
                'NumDetectionsSource',...
                'NumDetections'
                };
            groups = matlab.system.display.Section(...
           'Title','Parameters',...
           'PropertyList',props);
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:CFARDetector2DTitle')),...
                'Text',getString(message('phased:library:block:CFARDetector2DDesc',...
                'CA','GOCA','SOCA','OS')));
        end
    end
    
    methods (Access = protected)
        function str = getIconImpl(obj)
            str = sprintf('%s CFAR 2-D',obj.Method);
        end
        
        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            szCut = propagatedInputSize(obj,2);
            if strcmp(obj.OutputFormat,'CUT result')
                if numel(szX)==3
                    varargout{1} = [szCut(2) szX(3)];
                else
                    varargout{1} = [szCut(2) 1];
                end
                varargout{2} = varargout{1};
                varargout{3} = varargout{1};
            else % Detection index
                numDim = numel(szX);
                if numDim==3
                    numChan = szX(3);
                else
                    numChan = 1;
                end
                if strcmp(obj.NumDetectionsSource,'Auto')
                    numOut = szCut(2)*numChan;
                else
                    numOut = obj.NumDetections;
                end
                varargout{1} = [numDim numOut];
                varargout{2} = [1 numOut];
                varargout{3} = [1 numOut];
            end
        end
    end

    methods (Access = protected)
        function varargout = get2DTrainingInds(obj,X,Idx)
            %get2DTrainingIdx Calculate the 2-D training cells indices
            %   [FrontIndices, RearIndices]=get2DTrainingIdx(htr,Ncells,Idx)
            %   returns the training cell indices in front of the cell under
            %   test (FrontIndices) and after the cell under test
            %   (RearIndices). 

            [ NumRows, NumColumns, ~]  = size(X);

            GuardRegionSize = getGuardRegionSize(obj); 
            TrainingRegionSize = getTrainingRegionSize(obj); 
            
            % Compute Training and Guard cell indices
            indx1_r = linspace(Idx(1)-(TrainingRegionSize(1)-1)/2, ...
              Idx(1)+(TrainingRegionSize(1)-1)/2,TrainingRegionSize(1));     
            indx1_c = linspace(Idx(2)-(TrainingRegionSize(2)-1)/2, ...
              Idx(2)+(TrainingRegionSize(2)-1)/2,TrainingRegionSize(2));     
            indx2_r = linspace(Idx(1)-(GuardRegionSize(1)-1)/2, ...
              Idx(1)+(GuardRegionSize(1)-1)/2,GuardRegionSize(1));
            indx2_c = linspace(Idx(2)-(GuardRegionSize(2)-1)/2, ...
              Idx(2)+(GuardRegionSize(2)-1)/2,GuardRegionSize(2));

            % Generate the subscripts of the Guard and Training Regions
            [indx1_C, indx1_R] = meshgrid(indx1_c,indx1_r);
            [indx2_C, indx2_R] = meshgrid(indx2_c,indx2_r);

            % Convert subscripts to linear indices 
            linIndx1 = sub2ind(size(X),indx1_R(:),indx1_C(:));
            linIndx2 = sub2ind(size(X),indx2_R(:),indx2_C(:));
            linIndxCUT = sub2ind(size(X),Idx(1),Idx(2));

            % Compute indices of training cells, excluding the guard region and
            % CUT. 
            TrainingCellIndx = setdiff(linIndx1,linIndx2); 
            chIdx = 1:obj.pNumChannels;
            if (obj.Method(1) == 'S') || (obj.Method(1) == 'G') %SOCA || GOCA
                %Front training region are cells with indices less than CUT
                FrontTrainingCellIdx = TrainingCellIndx(TrainingCellIndx > linIndxCUT);
                RearTrainingCellIdx  = TrainingCellIndx(TrainingCellIndx < linIndxCUT);
                FrontTrainingCellIdx = repmat(FrontTrainingCellIdx,1,size(chIdx,2)) + ...
                  repmat((chIdx-1)*(NumRows*NumColumns),size(FrontTrainingCellIdx,1),1);
                RearTrainingCellIdx = repmat(RearTrainingCellIdx,1,size(chIdx,2)) + ...
                  repmat((chIdx-1)*(NumRows*NumColumns),size(RearTrainingCellIdx,1),1);
                varargout{1} = FrontTrainingCellIdx;
                varargout{2} = RearTrainingCellIdx;
            else
                TrainingCellIndx = repmat(TrainingCellIndx,1,size(chIdx,2)) + ...
                  repmat((chIdx-1)*(NumRows*NumColumns),size(TrainingCellIndx,1),1);
                varargout{1} = TrainingCellIndx;
            end
        end
        
        function validateCUTS(obj,X,Idx)
            [ NumRows, NumColumns, ~]  = size(X);

            TrainingRegionSize = getTrainingRegionSize(obj);   
            
            % Error if the training region falls outside of x
            outSideRow = any((Idx(1,:)-(TrainingRegionSize(1)-1)/2 < 1)) || ...
                  any((Idx(1,:) + (TrainingRegionSize(1)-1)/2 > NumRows));

            outSideCol = any((Idx(2,:) - (TrainingRegionSize(2)-1)/2 < 1)) || ...
                  any((Idx(2,:) + (TrainingRegionSize(2)-1)/2 > NumColumns));

            cond = outSideRow || outSideCol;
            if cond 
                coder.internal.errorIf(cond, ...
                    'phased:CFARDetector2D:TrainingCellsOutsideX');
            end

        end
    end
end

