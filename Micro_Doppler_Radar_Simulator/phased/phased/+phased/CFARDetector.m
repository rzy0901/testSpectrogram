classdef (Sealed,StrictDefaults) CFARDetector < phased.internal.AbstractDetector
%CFARDetector   Constant false alarm rate (CFAR) detector
%   H = phased.CFARDetector creates a constant false alarm rate (CFAR)
%   detector System object, H. This object performs CFAR detection on the
%   input data.
%
%   H = phased.CFARDetector(Name,Value) creates a CFAR detector object, H,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,IDX) performs the CFAR detection on the real input data X.
%   X can be either a column vector or a matrix. Each row of X is a cell
%   and each column of X is independent data. Detection is performed along
%   each column for the cells specified in IDX. IDX must be a vector of
%   positive integers with each entry specifying the index of a cell under
%   test (CUT).
%
%   Y contains the detection results. The format of Y depends on how the
%   OutputFormat property is set. By default, the OutputFormat property is
%   set to 'CUT result'.
% 
%   When the OutputFormat property is set to 'CUT result', Y is an MxN
%   matrix containing the detection results for each cell under test in
%   IDX. Target detections are indicated with a 1, and a 0 denotes the
%   absence of a target in the cell under test. M is the number of indices
%   specified in IDX and N is the number of independent signals in X.
% 
%   When the OutputFormat property is set to 'Detection index', Y is a DxQ
%   matrix with each column of Y containing the indices for each detection
%   found in X and the rows of Y defining these indices for each dimension
%   of X. D is 1 if X is a vector, 2 if X is a matrix, and 3 if X is a data
%   cube. When the NumDetectionsSource property is set to 'Auto', Q is the
%   total number of detections found in the data. When the
%   NumDetectionsSource property is set to 'Property', Q is set to the
%   NumDetections property value. By default, NumDetectionsSource is set to
%   'Auto'. Columns without detections are set to NaN.
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
%   The algorithm used in CFARDetector is cell averaging CFAR. Detection is
%   performed in three steps. First, the training cells are identified from
%   the input and the values in these training cells are averaged to form
%   the noise estimate. Second, the noise estimate is multiplied by the
%   threshold factor to form the threshold. Finally, the value in the test
%   cell is compared against the threshold to determine whether the target
%   is present or absent. If the value is greater than the threshold, the
%   target is present.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   CFARDetector methods:
%
%   step     - Perform CFAR detection (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create CFAR detector object with same property values
%   isLocked - Locked status (logical)
%
%   CFARDetector properties:
%
%   Method                - CFAR algorithm
%   Rank                  - Rank of order statistic
%   NumGuardCells         - Number of guard cells
%   NumTrainingCells      - Number of training cells
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
%   % Example:
%   %   Perform cell averaging CFAR detection on a given Gaussian noise
%   %   vector with a desired probability of false alarm of 0.1. Assume
%   %   that the data is from a square law detector and no pulse
%   %   integration is performed. Use 50 cells to estimate the noise
%   %   level and 1 cell to separate the test cell and training cells.
%   %   Perform the detection on all cells of input.
%
%   rs = RandStream.create('mt19937ar','Seed',5);
%   cfar = phased.CFARDetector('NumTrainingCells',50,'NumGuardCells',2,...
%                   'ProbabilityFalseAlarm',0.1);
%   N = 1000; x = 1/sqrt(2)*(randn(rs,N,1)+1i*randn(rs,N,1));
%   dresult = cfar(abs(x).^2,1:N);
%   Pfa = sum(dresult)/N
%
%   See also phased, phased.MatchedFilter, phased.TimeVaryingGain,
%   npwgnthresh.

%   Copyright 2009-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %NumGuardCells   Number of guard cells
        %   Specify the number of guard cells used in training as an even
        %   integer. It specifies the total number of cells on both sides
        %   of the cell under test. The default value of this property is 2
        %   indicating that there is one guard cell at both the front and
        %   back of the cell under test.
        NumGuardCells = 2
        %NumTrainingCells   Number of training cells
        %   Specify the number of training cells used in training as an
        %   even integer. Whenever possible, the training cells are equally
        %   divided before and after the cell under test. The default value
        %   of this property is 2 indicating that there is one training
        %   cell at both the front and back of the cell under test.
        NumTrainingCells = 2
    end

    properties (Access=private, Nontunable)
        cTraining;
    end
    
    properties (Access=private)
        pMaximumCellIndex;
    end
    
    properties (Access = private, Logical)
        pSizeInitialized
    end
    
    methods
        
        function set.NumGuardCells(obj,val)
            validateattributes( val, { 'double','single' }, { 'finite', 'nonnan', 'nonnegative', 'even', 'scalar' }, '', 'NumGuardCells');
            obj.NumGuardCells = val;
        end
        
        function set.NumTrainingCells(obj,val)
            validateattributes( val, { 'double','single' }, { 'finite', 'nonnan', 'positive', 'even', 'scalar' }, '', 'NumTrainingCells');
            obj.NumTrainingCells = val;
        end
        
    end
    
    methods
        function obj = CFARDetector(varargin)
            obj = obj@phased.internal.AbstractDetector(varargin{:});
        end
    end
    
    methods (Access = protected)
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractVarSizeEngine(obj);
            release(obj.cTraining);
        end
        
        function resetImpl(obj)
            reset(obj.cTraining);
            obj.pSizeInitialized = false;
        end
        
        function processInputSizeChangeImpl(obj,x,~,~)
            obj.pMaximumCellIndex = size(x,1);
        end

        function setupImpl(obj,X,CUTIdx,ThFac) %#ok<INUSD>
            setupImpl@phased.internal.AbstractDetector(obj);
            obj.pMaximumCellIndex = size(X,1);
            obj.pNumChannels = getNumChannels(obj,X);
            obj.pNumInputChannels = getNumChannels(obj,X);
            obj.pValidatedNumInputChannels = getNumChannels(obj,X);

            if (obj.Method(1) == 'S') || (obj.Method(1) == 'G') %SOCA || GOCA
                obj.cTraining = phased.internal.CFARTraining(...
                    'NumGuardCells',obj.NumGuardCells,...
                    'NumTrainingCells',obj.NumTrainingCells,...
                    'CombineTrainingData',false);
            else
                 obj.cTraining = phased.internal.CFARTraining(...
                    'NumGuardCells',obj.NumGuardCells,...
                    'NumTrainingCells',obj.NumTrainingCells,...
                    'CombineTrainingData',true);
            end
        end
        
        function validateInputsImpl(obj,x,cutidx,varargin)
            validateInputsImpl@phased.internal.AbstractDetector(obj,x,cutidx,varargin{:});
            
            cond = ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeMatrix','X');
            end
            
            validateNumChannels(obj,x);
            
            cond = ~isvector(cutidx) || isempty(cutidx);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeVector','Idx');
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractDetector(obj);
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.cTraining = saveobj(obj.cTraining);
                s.pMaximumCellIndex = obj.pMaximumCellIndex;
                s.pSizeInitialized = obj.pSizeInitialized;
            end
        end
        
        function s = loadSubObjects(obj,s)
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cTraining = phased.internal.CFARTraining.loadobj(s.cTraining);
                    s = rmfield(s,'cTraining');
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
        
        function [yout,varargout] = stepImpl(obj,X,CUTIdx,ThFac)
            
            classtoUse=class(X);
            if ~obj.pSizeInitialized
                processInputSizeChangeImpl(obj,X);
                obj.pSizeInitialized = true;
            end
        
            sigdatatypes.validateIndex(CUTIdx,'step','Idx',...
                {'double','single'},{'<=',obj.pMaximumCellIndex});
            if obj.ThresholdFactor(1) ~= 'I'
                ThFac = obj.pFactor;
            else
                validateattributes(ThFac,{'double','single'},{'positive','scalar'},...
                    'step','K');
            end
            
            NumCUT = numel(CUTIdx);
            th = zeros(NumCUT,obj.pNumChannels,classtoUse);
            noise = zeros(NumCUT,obj.pNumChannels,classtoUse); 
            for m = 1:NumCUT
                if obj.Method(1) == 'C' %'CA'
                    trndata = step(obj.cTraining,X,CUTIdx(m));
                    % Averaging cells
                    noisepowerest = mean(trndata,1);
                elseif obj.Method(1) == 'S' %'SOCA'
                    [trndatal,trndatat] = step(obj.cTraining,X,CUTIdx(m));
                    noisepowerest = min(mean(trndatal,1),mean(trndatat,1));
                elseif obj.Method(1) == 'G' %'GOCA'
                    [trndatal,trndatat] = step(obj.cTraining,X,CUTIdx(m));
                    noisepowerest = max(mean(trndatal,1),mean(trndatat,1));
                else %'OS'
                    trndata = step(obj.cTraining,X,CUTIdx(m));
                    temp = sort(trndata,1,'ascend');
                    noisepowerest = temp(obj.Rank,:);
                end
                % Form threshold
                th(m,:) = noisepowerest * ThFac;
                noise(m,:) = noisepowerest;
            end
            
            xCUTinds = bsxfun(@plus,CUTIdx(:),(0:obj.pNumChannels-1)*size(X,1));
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
                        numOut = NumCUT*obj.pNumChannels;
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
                    yout = NaN(1,numOut,classtoUse);
                    detidx = ind(1:numAvail);
                else
                    yout = NaN(2,numOut,classtoUse);
                    
                    [idx1,idx2] = ind2sub(size(X),ind(1:numAvail));
                    detidx = [idx1;idx2];
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

    methods (Access=protected)
        function NumTrainingCells = getNumTrainingCells(obj) 
            NumTrainingCells = obj.NumTrainingCells;
        end  
    end
    
    methods (Static,Hidden,Access=protected) 
        function groups = getPropertyGroupsImpl
           props = {...
                'Method',...
                'Rank',...
                'NumGuardCells',...
                'NumTrainingCells',...
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
                'Title',getString(message('phased:library:block:CFARDetectorTitle')),...
                'Text',getString(message('phased:library:block:CFARDetectorDesc',...
                'CA','GOCA','SOCA','OS')));
        end
    end

    methods (Access = protected)
        function str = getIconImpl(obj) 
            str = sprintf('%s CFAR',obj.Method);
        end
        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            szCut = propagatedInputSize(obj,2);
            if strcmp(obj.OutputFormat,'CUT result')
                varargout{1} = [max(szCut) szX(2)];
                varargout{2} = varargout{1};
                varargout{3} = varargout{1};
            else % Detection index
                numDim = numel(szX);
                if numDim==3
                    numChan = szX(3);
                else
                    numChan = 1;
                    if szX(2) == 1
                        numDim = 1;
                    end
                end
                if strcmp(obj.NumDetectionsSource,'Auto')
                    numOut = max(szCut)*numChan;
                else
                    numOut = obj.NumDetections;
                end
                varargout{1} = [numDim numOut];
                varargout{2} = [1 numOut];
                varargout{3} = [1 numOut];
            end
        end
    end
end
