classdef (Hidden) AbstractClusterDBSCAN < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.

%AbstractClusterDBSCAN   Define the AbstractClusterDBSCAN class.

%   Copyright 2019 The MathWorks, Inc.

%#codegen

    properties (Nontunable)
        %EpsilonSource          Source of cluster threshold epsilon
        %   Specify the source of the cluster threshold epsilon as one of
        %   'Auto' | 'Property'. When you set the EpsilonSource property to
        %   'Auto', epsilon is set using a k-Nearest Neighbor search. The
        %   k-NN search is calculated with k equal to
        %   (MinNumPoints-1):(MaxNumPoints-1), where the '-1' corresponds
        %   to a cluster consisting of MinNumPoints, including the point
        %   under consideration as a core point. If EpsilonSource is set to
        %   'Property', the cluster threshold is set to the property
        %   Epsilon, and the property MaxNumPoints is inactive. Defaults to
        %   'Property'.
        EpsilonSource = 'Property'
    end

    properties
        %Epsilon                Cluster threshold epsilon
        %   Epsilon is a positive scalar threshold for a neighborhood
        %   search query. It defines a radius around a core point. This
        %   property is active when EpsilonSource is set to 'Property'.
        %   Defaults to 10.
        %
        %   Alternatively, Epsilon can be designated as a 1xP row vector
        %   where P is the number of clustering dimensions in the input
        %   data X. The epsilons for the different clustering dimensions in
        %   the row vector create a hyperellipse search area.
        Epsilon = 10
        %MinNumPoints           Minimum number of points in a cluster
        %   Positive integer used as a threshold to determine whether a
        %   point is a core point. This property defines the minimum number
        %   of points required for a cluster. Defaults to 3 points.
        MinNumPoints = 3 
        %MaxNumPoints           Maximum number of points for 'Auto' epsilon
        %   Positive integer defining the maximum number of points expected
        %   in a cluster. Used to estimate epsilon using a k-Nearest
        %   Neighbor search. This property is only active when
        %   EpsilonSource is set to 'Auto'. Defaults to 10 points.
        MaxNumPoints = 10
    end
    
    properties (Nontunable)
        %EpsilonHistoryLength   Length of cluster threshold epsilon history
        %   Specify the length of the epsilon history to store as a
        %   positive integer. EpsilonHistoryLength set to 1 is memory-less.
        %   Memory-less means that each epsilon estimate is immediately
        %   used and no moving-average smoothing occurs. If
        %   EpsilonHistoryLength is greater than 1, the epsilon value will
        %   be averaged over the history length specified. Defaults to 10.
        EpsilonHistoryLength = 10 
    end

    properties (Nontunable,Logical)
        %EnableDisambiguation   Enable disambiguation of dimensions
        %   A logical that specifies if disambiguation of detections should
        %   occur. If EnableDisambiguation is true, clustering will occur
        %   across boundaries defined by the step input AMBLIMS to ensure
        %   that ambiguous detections are appropriately clustered. The
        %   property AmbiguousIndices must be defined as it specifies the
        %   column indices of X in which the extended clustering should
        %   occur. Up to 2 ambiguous dimensions are permitted. Not
        %   recommended for large data sets. Defaults to false.
        EnableDisambiguation(1,1) = false
    end
   
    properties (Nontunable)
        %AmbiguousDimension     Indices of ambiguous dimensions
        %   If EnableDisambiguation is true, the AmbiguousDimension
        %   property is required as it specifies the column indices of X in
        %   which disambiguation should occur. The AmbiguousDimension
        %   property is either a scalar index for a single ambiguous
        %   dimension in the input data matrix X or a 1x2 length row vector
        %   of indices for 2 ambiguous dimensions in X. The size and order
        %   of AmbiguousDimension must be consistent with the step input
        %   AMBLIMS. This property is inactive when the
        %   EnableDisambiguation property is false. The default assumes
        %   column index 1 of X is the ambiguous dimension.
        %
        %   The AMBLIMS step input may be specified as either a 1x2 vector
        %   for a single ambiguous dimension where the fields are
        %   [MinAmbiguityLimitDimension1, MaxAmbiguityLimitDimension1] or a
        %   2x2 matrix for 2 ambiguous dimensions where the fields are
        %   [MinAmbiguityLimitDimension1, MaxAmbiguityLimitDimension1;
        %   MinAmbiguityLimitDimension2, MaxAmbiguityLimitDimension2].
        AmbiguousDimension = 1
    end

    properties (Access = protected)       
        pEpsilon = 10 
        pEpsilonIsInitialized = false
        pEpsilonHistoryVec
        pEpsIdx = 1
        
        pAmbiguityLimits = [1 1; 1 1]
        pUpdateEpsilon = false 
        
        pAmbiguitiesDetected = [false false]
        
        pFigHandle  
        pFigAxes
        pFigAxesNV
        pFigInitialized = false
        pScatterObjInitialized = false
        pScatterClustersHandle
        pScatterNoiseObjInitialized = false
        pScatterNoiseHandle
        pTextHandle 
        pTitle
    end
    
    properties(Constant, Hidden)
        EpsilonSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        pDefaultTitle = 'Clusters'; 
    end

    methods 
        function obj = AbstractClusterDBSCAN(varargin)
            setProperties(obj, nargin, varargin{:});

            obj.pEpsilonHistoryVec = nan(1,obj.EpsilonHistoryLength);
        end
    end
    
    methods
        function set.Epsilon(obj,value)
            obj.validateEpsilon(value); 
            obj.Epsilon = value;
        end
        function set.MinNumPoints(obj,value)
            validateattributes(value,{'double'},{'scalar','integer','real',...
                'nonempty','finite',...
                'nonnan','positive'},'','MinNumPoints');
            obj.MinNumPoints  = value;
        end
        function set.MaxNumPoints(obj,value)
            validateattributes(value,{'double'},{'scalar','integer','real',...
                'nonempty','finite',...
                'nonnan','positive'},'','MaxNumPoints');
            obj.MaxNumPoints  = value;
        end
        function set.EpsilonHistoryLength(obj,value)
            validateattributes(value,{'double'},{'scalar','integer','real',...
                'nonempty','finite',...
                'nonnan','positive'},'','EpsilonHistoryLength');
            obj.EpsilonHistoryLength = value;
        end
        function set.AmbiguousDimension(obj,value)
            validateattributes(value,{'double'},{'finite','nonnan','nonempty',...
                'positive','row','integer'},'','AmbiguousDimension');
            cond = size(value,2)>2;
            if cond
                coder.internal.errorIf(cond,'phased:phased:tooManyColumns','AmbiguousDimension',2);
            end
            obj.AmbiguousDimension = value;
        end
    end

    methods (Access = protected)
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and system block dialog
            flag = false; 
            if ~obj.EnableDisambiguation
                if strcmp(prop,'AmbiguousDimension')
                    flag = true;
                end
            end
            if strcmp(obj.EpsilonSource,'Auto') && strcmp(prop,'Epsilon')
                flag = true; 
            end
            if strcmp(obj.EpsilonSource,'Property') && ...
                    (strcmp(prop,'MaxNumPoints') || strcmp(prop,'EpsilonHistoryLength'))
                flag = true;
            end
        end
        
        function varargout = clusterdataDBSCAN(obj,x,varargin)     
            if isempty(x) 
                % x is empty; define empty outputs 
                idx = [];
                clusterIDs = []; 
            elseif all(isnan(x))
                % x is all NaNs; define outputs as all noise and assign no
                % cluster IDs
                nDets = size(x,1);
                idx = -1*ones(nDets,1); 
                clusterIDs = nan(1,nDets);
            else
                % x is not empty; cluster 
                validateInputs(obj,x,varargin{:});
                
                % Remove NaNs from consideration
                idxGd  = all(~isnan(x),2);
                
                % Set cluster threshold epsilon
                if strcmp(obj.EpsilonSource,'Auto')
                    if ~obj.pEpsilonIsInitialized || obj.pUpdateEpsilon
                        autoEpsilon(obj,x(idxGd,:));
                    end
                else
                    obj.pEpsilon = obj.Epsilon;
                    obj.pEpsilonIsInitialized = true;
                end
                
                % Roll out grid to handle ambiguous returns
                xCluster = setupAmbiguousClustering(obj,x(idxGd,:));
                
                % Cluster data
                [idxTmp,~] = dbscanHyperellipse(obj,xCluster);
                
                % Update indices based on ambiguous returns
                idxTmp2 = updateAmbigousClustering(obj,xCluster,idxTmp);
                
                % Initialize output values
                nDets = size(x,1);
                idx = -1*ones(nDets,1); 
                
                % Assign idx 
                idx(idxGd) = idxTmp2; 
                
                % Assign cluster IDs for compatibility with other Phased
                % Array System Toolbox System Objects
                clusterIDs = idx.';
                nClusters = max([clusterIDs 0],[],2);
                idxNoise = idx == -1;
                numNoisePoints = sum(double(idxNoise));
                if numNoisePoints > 0 
                    clusterIDsTmp = ones(size(clusterIDs));
                    clusterIDs(idxNoise) = cumsum(clusterIDsTmp(idxNoise),2)+nClusters;
                    clusterIDs(any(isnan(x),2)) = NaN; 
                end
            end
            
            % Set outputs
            if nargout >= 1 
                varargout{1} = idx; 
            end
            if nargout >= 2
                varargout{2} = clusterIDs; 
            end
        end
        
        function validateInputsImpl(obj,x,varargin)
            % Check data 
            cond = ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','double');
            end
            
            % Check conditional inputs
            idx = 1;
            if obj.EnableDisambiguation
                cond = ~isa(varargin{idx},'double');
                idx = idx+1; 
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','AMBLIMS','double');
                end
            end
            if strcmp(obj.EpsilonSource,'Auto')
                cond = ~isa(varargin{idx},'logical');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','UPDATE','double');
                end
            end
        end
             
        function validateInputs(obj,x,varargin)
            % Method for validation of inputs to the step method
            
            [nDets,nDims] = size(x); 
            
            % Verify x
            validateattributes(x,{'double'},{'2d','nonsparse','real',...
                'nonempty'},'','X');
            cond = any(isinf(x(:))); 
            if cond
                coder.internal.errorIf(cond,'phased:clustering:expectedFinite','X');
            end
            
            % Verify epsilon for given x
            if ~strcmp(obj.EpsilonSource,'Auto')
                if ~isscalar(obj.Epsilon)
                    [nRows,nCols] = size(obj.Epsilon);
                    
                    % Check number of columns
                    cond = nCols ~= nDims;
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:NumColumnsMismatch','X','Epsilon');
                    end
                    
                    if ~isrow(obj.Epsilon)
                        % Check number of rows
                        cond = nRows ~= nDets;
                        if cond
                            coder.internal.errorIf(cond,'phased:phased:expectedNumRows','Epsilon',sprintf('%d, %d',1,nDets));
                        end
                    end
                end
            end
            
            % Parse and check property dependent variables
            idx = 1;
            
            % Disambiguation Parameters
            if obj.EnableDisambiguation
                if isscalar(obj.AmbiguousDimension)
                    validateattributes(varargin{idx},{'double'},{'real','finite','nonnan','nonempty',...
                        'ncols',2,'nrows',1},'','AMBLIMS'); 
                else
                    validateattributes(varargin{idx},{'double'},{'real','finite','nonnan','nonempty',...
                        'ncols',2,'nrows',2},'','AMBLIMS'); 
                end
                obj.pAmbiguityLimits = varargin{idx};
                idx = idx+1;
                
                % Check that first column is less than second column
                cond = any(obj.pAmbiguityLimits(:,1) > obj.pAmbiguityLimits(:,2));
                if cond
                    coder.internal.errorIf(cond,'phased:clustering:invalidAmbiguityLimits','AMBLIMS');
                end
               
                % Ambiguity dimensions specified can not be greater than
                % the number of columns in x
                cond = any(obj.AmbiguousDimension>nDims);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:tooFewColumns','X',max(obj.AmbiguousDimension,[],2));
                end
                
                % Check to make sure all ambiguous dimension indices are
                % unique
                if numel(obj.AmbiguousDimension) == 2
                    cond = obj.AmbiguousDimension(1) == obj.AmbiguousDimension(2);
                    if cond
                        coder.internal.errorIf(cond,'phased:clustering:expectedUniqueIndices','AmbiguousDimension');
                    end
                end
                
                % Check that the ambiguity dimensions of x are contained by
                % the ambiguity limits in the corresponding dimension
                for ii = 1:numel(obj.AmbiguousDimension)
                    ambIdx = obj.AmbiguousDimension(ii);
                    greaterThanLim2 = x(:,ambIdx) > obj.pAmbiguityLimits(ii,2);
                    lessThanLim1 = x(:,ambIdx) < obj.pAmbiguityLimits(ii,1);
                    cond = any(lessThanLim1 | greaterThanLim2); 
                    if cond
                        coder.internal.errorIf(cond, 'phased:clustering:expectedWithinLimits','X','AMBLIMS',ambIdx);
                    end
                end
            end
            
            % Auto Epsilon Parameters
            if strcmp(obj.EpsilonSource,'Auto')
                validateattributes(varargin{idx},{'logical'},{'scalar','real',...
                    'nonempty','nonnan'},'','UPDATE');
                obj.pUpdateEpsilon = varargin{idx};
            end
        end

        function icon = getIconImpl(obj) %#ok<MANU>
            % Define icon for System block
            icon = {'DBSCAN' 'Clusterer'}; % Example: text icon
        end
        
        function flag = isInputSizeMutableImpl(obj,index) %#ok<INUSL>
            % Return false if input size cannot change
            % between calls to the System object
            if index == 1
                flag = true;
            else
                flag = false;
            end
        end

        function flag = isInputComplexityMutableImpl(obj,index)  %#ok<INUSD>
            % Return false if input complexity cannot change
            % between calls to the System object
            flag = false;
        end

        function flag = isInputDataTypeMutableImpl(obj,index) %#ok<INUSD>
            % Return false if input data type cannot change
            % between calls to the System object
            flag = false;
        end

        function num = getNumInputsImpl(obj)
            % Define total number of inputs for system with optional inputs
            num = 1; 
            if obj.EnableDisambiguation
                num = num+1; 
            end
            if strcmp(obj.EpsilonSource,'Auto')
                num = num+1;
            end
        end
    
        function loadObjectImpl(obj,s,~)
            % Set properties in object obj to values in structure s
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            
            if isLocked(obj)               
                s.pEpsilon = obj.pEpsilon;
                s.pEpsilonIsInitialized = obj.pEpsilonIsInitialized;
                s.pEpsilonHistoryVec = obj.pEpsilonHistoryVec;
                s.pEpsIdx = obj.pEpsIdx;
                
                s.pAmbiguityLimits = obj.pAmbiguityLimits;
                s.pUpdateEpsilon = obj.pUpdateEpsilon;
                
                s.pAmbiguitiesDetected = obj.pAmbiguitiesDetected;
                
                s.pFigHandle = obj.pFigHandle;
                s.pFigAxes = obj.pFigAxes;
                s.pFigAxesNV = obj.pFigAxesNV;
                s.pFigInitialized = obj.pFigInitialized;
                s.pScatterObjInitialized = obj.pScatterObjInitialized;
                s.pScatterClustersHandle = obj.pScatterClustersHandle;
                s.pScatterNoiseObjInitialized = obj.pScatterNoiseObjInitialized;
                s.pScatterNoiseHandle = obj.pScatterNoiseHandle;
                s.pTextHandle = obj.pTextHandle;
                s.pTitle = obj.pTitle;
            end
        end

        function varargout = getInputNamesImpl(obj)
            % Return input port names for System block
            varargout{1} = 'X';
            idx = 2;
            if obj.EnableDisambiguation
                varargout{idx} = 'AmbLims';
                idx = idx+1;
            end
            if strcmp(obj.EpsilonSource,'Auto')
                varargout{idx} = 'Update';
            end
        end

        function [name,name1] = getOutputNamesImpl(obj) %#ok<MANU>
            % Return output port names for System block
            name = 'Idx';
            name1 = 'Clusters';
        end

        function [out,out1] = getOutputSizeImpl(obj)
            % Return size for each output port
            sz = propagatedInputSize(obj,1);
            out = [sz(1) 1];
            out1 = [1 sz(1)]; 
        end

        function [out,out1] = getOutputDataTypeImpl(obj) %#ok<MANU>
            % Return data type for each output port
            out = "double";
            out1 = out; 
        end

        function [out,out1] = isOutputComplexImpl(obj) %#ok<MANU>
            % Return true for each output port with complex data
            out = false;
            out1 = out; 
        end

        function [out,out1] = isOutputFixedSizeImpl(obj)
            % Return true for each output port with fixed size
            out = propagatedInputFixedSize(obj,1);
            out1 = out; 
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractSampleRateEngine(obj);
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
        end
    end
    
    methods (Static,Hidden,Access=protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:DBSCANClustererTitle')),...
                'Text',getString(message('phased:library:block:DBSCANClustererDesc')));
        end
    end

    methods 
        function varargout = plot(obj,x,idx,varargin)
            %plot    Plot DBSCAN clustering results
            %   HF = plot(H,X,IDX) returns a plot of the DBSCAN clustering
            %   results. X is the NxP data matrix, where N is the number of
            %   points (or detections), and P is the clustering dimensions
            %   such as range, azimuth, and elevation.
            %
            %   IDX is an Nx1 column vector containing cluster IDs that
            %   represents the clustering results of the DBSCAN algorithm.
            %   A value equal to '-1' implies a DBSCAN noise point.
            %   Positive IDX values correspond to clusters that satisfy the
            %   DBSCAN clustering criteria.
            %
            %   Each cluster is assigned a unique color and is tagged with
            %   its corresponding DBSCAN cluster index. Noise points,
            %   points that have an IDX value of '-1', are not labeled and
            %   are assigned a gray color.
            %
            %   Assumes IDX contains positive cluster IDs start at 1 and
            %   continue on to the total number of clusters with no missing
            %   indices.
            %
            %   HF is the figure handle.
            %
            %   HF = plot(...,'Parent',HAX) specifies the target axes HAX
            %   for the cluster results plot.
            %
            %   HF = plot(...'Title',TITLESTR) specifies the title of the
            %   cluster results plot as TITLESTR.
            %
            %   % Example:
            %   %   Cluster data using clusterDBSCAN. Compare results
            %   %   of clustering with 2 Epsilon values using plot method.
            %   
            %   % Create target data
            %   x = [rand(20,2)+12; rand(20,2)+10; rand(20,2)+15];
            %   
            %   % Create DBSCAN system object
            %   clusterRadarData = clusterDBSCAN('MinNumPoints',3);
            %   
            %   % Cluster with Epsilon = 1
            %   clusterRadarData.Epsilon = 1;
            %   idxEpsilon1 = clusterRadarData(x);
            %   
            %   % Cluster with Epsilon = 3
            %   clusterRadarData.Epsilon = 3;
            %   idxEpsilon2 = clusterRadarData(x);
            %   
            %   % Plot clustering results side-by-side. Pass in axes
            %   % handles and titles.
            %   hAx1 = subplot(1,2,1);
            %   plot(clusterRadarData,x,idxEpsilon1,...
            %       'Parent',hAx1,'Title','Epsilon = 1');
            %   hAx2 = subplot(1,2,2);
            %   plot(clusterRadarData,x,idxEpsilon2,...
            %       'Parent',hAx2,'Title','Epsilon = 3');
            %
            %   See also phased, clusterDBSCAN/discoverClusters,
            %   clusterDBSCAN/estimateEpsilon.
            
            narginchk(3,7);
            nargoutchk(0,1);
            
            % Not support for codegen
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                    'phased:Waveform:CodegenNotSupported','plot');
            end
            
            % Validate inputs to method
            parseAndValidatePlotInputs(obj,x,idx,varargin{:}); 
            
            % Check to see if a figure handle has been provided
            if isempty(obj.pFigAxesNV)
                % Check to see if figure was deleted
                if ~ishandle(obj.pFigHandle)
                    obj.pFigInitialized = 0;
                end
                
                % Check to see if scatter object was deleted
                if ~ishandle(obj.pScatterClustersHandle)
                    obj.pScatterObjInitialized = 0;
                end
                
                % Check to see if scatter object was deleted
                if ~ishandle(obj.pScatterNoiseHandle)
                    obj.pScatterNoiseObjInitialized = 0;
                end
            else
                obj.pFigAxes = obj.pFigAxesNV;
                obj.pFigInitialized = 1;
                obj.pFigHandle = obj.pFigAxes.Parent;
                obj.pScatterObjInitialized = 0;
                obj.pScatterNoiseObjInitialized = 0;
                obj.pTextHandle = [];
            end
            
            % Determine input characteristics
            idxNoise = idx == -1; 
            [nDets,nDims] = size(x); 
            idxDets = 1:nDets;
            nClusters = max(idx(~idxNoise),[],1);
            dbscanIDs = unique(idx(~idxNoise));
            
            % Initialize 3D data x and cluster text positions clusterXYZ
            clusterXYZ = zeros(nClusters,3);
            x3Dims = zeros(nDets,3);
            x3Dims(:,1:nDims) = x(:,1:nDims);
            medianVal = abs(median(x(:))); 
            padVal = medianVal*.1;
            padValLabels = medianVal*.05;
            
            % Calculate cluster text positions clusterXYZ
            idxNotChecked = true(nDets,1);
            for iC = 1:nClusters
                thisClusterIdx = idx(idxNotChecked) == dbscanIDs(iC);
                idxTheseDets = idxDets(idxNotChecked);
                thisClusterIdx = idxTheseDets(thisClusterIdx);
                if numel(thisClusterIdx)>1
                    clusterXYZ(iC,:)= max(x3Dims(thisClusterIdx,:),[],1);
                else
                    clusterXYZ(iC,:)= x3Dims(thisClusterIdx,:);
                end
                idxNotChecked(thisClusterIdx) = false;
            end
            
            % Determine cluster colors 
            clusterColors = hsv(nClusters);
            
            % Plot
            if obj.pScatterObjInitialized && obj.pScatterNoiseObjInitialized
                % Noise
                set(obj.pScatterNoiseHandle,...
                    'XData',x3Dims(idxNoise,1),...
                    'YData',x3Dims(idxNoise,2),...
                    'ZData',x3Dims(idxNoise,3)); % Update noise

                % Clusters
                set(obj.pScatterClustersHandle,...
                    'XData',x3Dims(~idxNoise,1),...
                    'YData',x3Dims(~idxNoise,2),...
                    'ZData',x3Dims(~idxNoise,3),...
                    'CData',clusterColors(idx(~idxNoise),:)); % Update clusters
                    
                %Update axes 
                obj.pFigAxes.XLim = [min(x3Dims(:,1))-padVal max(x3Dims(:,1),[],1)+padVal];
                obj.pFigAxes.YLim = [min(x3Dims(:,2))-padVal max(x3Dims(:,2),[],1)+padVal];
                obj.pFigAxes.ZLim = [min(x3Dims(:,3))-padVal max(x3Dims(:,3),[],1)+padVal];
                obj.pFigAxes.XLabel.String = 'Dimension 1';
                obj.pFigAxes.YLabel.String = 'Dimension 2';
                if nDims == 3
                    obj.pFigAxes.ZLabel.String = 'Dimension 3';
                end
                obj.pFigAxes.Title.String = obj.pTitle;
                drawnow limitrate;
                
                % Add text labels
                addClusterLabels(obj,padValLabels,clusterXYZ,dbscanIDs);
            else
                if ~obj.pFigInitialized
                    %Figure is not initialized
                    obj.pFigHandle = figure('Name','Clusters');
                    obj.pFigInitialized = 1;
                    obj.pFigAxes = axes(obj.pFigHandle); 
                else 
                    % Figure is initialized so prep for new axes 
                    if nargin<4
                        clf(obj.pFigHandle);
                        obj.pFigAxes = axes(obj.pFigHandle); 
                    else
                        cla(obj.pFigAxes); 
                    end
                end
                
                % Noise
                obj.pScatterNoiseHandle = line(obj.pFigAxes,...
                    x3Dims(idxNoise,1),x3Dims(idxNoise,2),x3Dims(idxNoise,3),...
                    'Color',[0.8 0.8 0.8],'LineStyle','none',...
                    'Marker','o','MarkerSize',4,....
                    'MarkerFaceColor',[0.8 0.8 0.8],...
                    'MarkerEdgeColor',[0.6 0.6 0.6]);
                hold on;
                
                % Clusters
                obj.pScatterClustersHandle = scatter3(obj.pFigAxes,...
                    x3Dims(~idxNoise,1),...
                    x3Dims(~idxNoise,2),...
                    x3Dims(~idxNoise,3),...
                    20,clusterColors(idx(~idxNoise),:),'filled');
                
                % Designate plots as initialized
                obj.pScatterNoiseObjInitialized = 1;
                obj.pScatterObjInitialized = 1;
                
                % Update axes
                grid(obj.pFigAxes,'on'); 
                obj.pFigAxes.XLim = [min(x3Dims(:,1))-padVal max(x3Dims(:,1),[],1)+padVal];
                obj.pFigAxes.YLim = [min(x3Dims(:,2))-padVal max(x3Dims(:,2),[],1)+padVal];
                obj.pFigAxes.ZLim = [min(x3Dims(:,3))-padVal max(x3Dims(:,3),[],1)+padVal];
                obj.pFigAxes.XLabel.String = 'Dimension 1';
                obj.pFigAxes.YLabel.String = 'Dimension 2';
                obj.pFigAxes.Title.String = obj.pTitle;
                
                % Update view
                if nDims == 3
                    obj.pFigAxes.ZLabel.String = 'Dimension 3';
                    view(obj.pFigAxes,36,15);
                else
                    view(obj.pFigAxes,0,90); % X-Y view
                end
                
                % Add text labels
                addClusterLabels(obj,padValLabels,clusterXYZ,dbscanIDs);
                
                drawnow;
            end
          
            % Define output handle if requested
            if nargout == 1
                varargout{1} = obj.pFigHandle;
            end
        end
    end
    
    methods (Static)
        function varargout = discoverClusters(x,maxEpsilon,minNumPoints)
            %discoverClusters   Discover cluster hierarchy using OPTICS
            %   [ORDER,REACHDIST] = clusterDBSCAN.discoverClusters(...
            %   X,MAXEPSILON,MINNUMPOINTS) assists with the discovery of
            %   the underlying cluster hierarchy in the NxP data matrix X,
            %   where N is the number of points and P is the clustering
            %   dimensions, by implementing the OPTICS (Ordering Points To
            %   Identify the Clustering Structure) algorithm. This method
            %   returns the ORDER, which is a 1xN cluster ordered list of
            %   sample indices and REACHDIST, which is a 1xN vector of
            %   reachability distances per sample. The output of the OPTICS
            %   algorithm is equivalent to density-based clustering with a
            %   broad range of parameter settings and can be used to assist
            %   with the estimation of an appropriate epsilon clustering
            %   threshold for the DBSCAN algorithm.
            %
            %   For example, given 4 Nx1 column-vectors of detected range,
            %   Doppler, azimuth, and elevation values for N radar
            %   detections represented by detRng, detDop, detAz, and detEl,
            %   respectively, the input X can be formed as X = [detRng,
            %   detDop, detAz, detEl].
            %
            %   The OPTICS algorithm is based on an assumption that dense
            %   clusters are entirely contained by less dense clusters, and
            %   it involves processing data in the correct order by
            %   tracking the point density neighborhoods. This is performed
            %   by ordering data points by the shortest reachability
            %   distances, thereby guaranteeing that clusters with higher
            %   density are identified first.
            %
            %   MAXEPSILON is the maximum epsilon threshold to use in the
            %   cluster hierarchy search and is specified as a positive
            %   scalar. The epsilon parameter defines a neighborhood around
            %   a point. Reducing MAXEPSILON results in shorter run times.
            %   Setting MAXEPSILON to inf will identify all cluster
            %   possibilities. MINNUMPOINTS is a positive integer used as a
            %   threshold to determine whether a point is a core point. It
            %   defines the minimum number of points for a cluster. The
            %   OPTICS algorithm is relatively insensitive to parameter
            %   settings, but it is best to err on the side of larger
            %   parameters to ensure good results.
            %
            %   clusterDBSCAN.discoverClusters(...) with no output
            %   arguments creates a bar graph that is a visual
            %   representation of the cluster hierarchy.
            %
            %   The OPTICS reachability plot represents the density-based
            %   clustering structure and is independent of the
            %   dimensionality of the data. The valleys within the OPTICS
            %   plot correspond to primary clusters.
            %
            %   % Example:
            %   %   Create target data. Use the discoverClusters static
            %   %   method to discover the underlying cluster hierarchy.
            %
            %   % Set cluster parameters
            %   maxEpsilon = 10;
            %   minNumPoints = 5;
            %
            %   % Create target data
            %   x = [rand(20,2)+11.5; rand(20,2)+10;...
            %   rand(20,2)+15; 30*rand(5,2)];
            %
            %   % Determine cluster hierarchy
            %   clusterDBSCAN.discoverClusters(x,maxEpsilon,minNumPoints)
            %
            %   See also phased, clusterDBSCAN/estimateEpsilon,
            %   clusterDBSCAN/plot.
            
            %   Copyright 2019 The MathWorks, Inc.
            
            % Reference
            % [1] Ankerst M., Breunig M. M., Kriegel H.-P., and Sander J.
            %    "OPTICS: Ordering Points to Identify the Clustering
            %    Structure." Proc. ACM SIGMOD'99 Int. Conf. on Management
            %    of Data, Philadelphia, PA, ACM Press, 1999, pp. 49-60.
            
            %#codegen
            
            % Check inputs/outputs
            narginchk(3,3);
            if coder.target('MATLAB') % For MATLAB
                nargoutchk(0,2);
            else
                nargoutchk(1,2);
            end
            
            % Validate inputs
            insufficientDataFlag = phased.internal.AbstractClusterDBSCAN.validateInputsDiscoverClusters(x,maxEpsilon,minNumPoints);
            if insufficientDataFlag
                if nargout >= 1
                    varargout{1} = [];
                end
                if nargout >= 2
                    varargout{2} = [];
                end
                return;
            end
            
            % Run OPTICS
            [reachDist,order] = phased.internal.AbstractClusterDBSCAN.optics(x,maxEpsilon,minNumPoints);
            
            % Plot reachability distances
            if nargout == 0
                hFig = figure('Name','Reachability Distances');
                hAxReach = axes(hFig);
                bar(hAxReach,order,reachDist);
                hAxReach.XLabel.String = 'Order Index';
                hAxReach.YLabel.String = '\epsilon';
                hAxReach.Title.String = 'Reachability Distances';
                grid on;
                hAxReach.YLim = [0 maxEpsilon];
            end
            
            % Define outputs
            if nargout >= 1
                varargout{1} = order;
            end
            if nargout == 2
                varargout{2} = reachDist;
            end
        end
        
        function varargout = estimateEpsilon(x,minNumPoints,maxNumPoints)
            %estimateEpsilon    Estimate cluster threshold from data
            %   EPSILON = clusterDBSCAN.estimateEpsilon(...
            %   X,MINNUMPOINTS,MAXNUMPOINTS) returns an estimate of the
            %   DBSCAN (Density-Based Spatial Clustering of Applications
            %   with Noise) clustering threshold epsilon based on the input
            %   NxP data matrix X using a k-Nearest Neighbor (k-NN) search,
            %   where N is the number of points (or detections), and P is
            %   the clustering dimensions.
            %
            %   For example, given 4 Nx1 column-vectors of detected range,
            %   Doppler, azimuth, and elevation values for N radar
            %   detections represented by detRng, detDop, detAz, and detEl,
            %   respectively, the input X can be formed as X = [detRng,
            %   detDop, detAz, detEl].
            %
            %   MINNUMPOINTS defines the minimum number of points expected
            %   in a cluster, and MAXNUMPOINTS defines the maximum number
            %   of points expected in a cluster. DBSCAN defines a core
            %   point as a point P whose epsilon-neighborhood includes the
            %   minimum number of points MINNUMPOINTS, including the core
            %   point. Thus, the k-NN search is calculated with k equal to
            %   (MINNUMPOINTS-1):(MAXNUMPOINTS-1), where the '-1'
            %   corresponds to a cluster consisting of k number of points,
            %   including the point under consideration as a core point.
            %
            %   It is recommended to set MINNUMPOINTS and MAXNUMPOINTS
            %   according to a range of values around the expected number
            %   of detections to ensure the best epsilon estimates. If N <
            %   MINNUMPOINTS, there is insufficient data, and the function
            %   returns EPSILON equal to NaN. If N < MAXNUMPOINTS, the k-NN
            %   search is limited to (MINNUMPOINTS-1):(N-1).
            %
            %   The MINNUMPOINTS property of the DBSCAN algorithm is more
            %   easily set than the clustering threshold Epsilon. The
            %   purpose of MINNUMPOINTS is to smooth the density estimate.
            %   Since a DBSCAN cluster is a maximal set of
            %   density-connected points, it is best to err on the side of
            %   a lower value for MINNUMPOINTS when the expected number of
            %   detections in a cluster is unknown, even though lower
            %   values make the DBSCAN algorithm more susceptible to noise.
            %   A general guideline for MINNUMPOINTS selection is as
            %   follows:
            %      MINNUMPOINTS = 2*P,
            %   where P is the number of clustering dimensions in X. For
            %   challenging data sets with the following qualities:
            %      - highly noisy,
            %      - very large N,
            %      - high dimensionality P,
            %      - or many duplicates
            %   increasing MINNUMPOINTS may improve the DBSCAN clustering
            %   results.
            %
            %   For the input data matrix X, epsilons are calculated as the
            %   knees of the series of k-NN search curves. The estimated
            %   epsilon is given by the average of the estimated knees of
            %   the k-NN curves.
            %
            %   The knees of the curves are estimated as the point farthest
            %   from the line that is created from the first point and last
            %   point in each curve.
            %
            %   Estimating epsilon is computationally intensive. It is not
            %   recommended for large data sets.
            %
            %   clusterDBSCAN.estimateEpsilon(...) with no output arguments
            %   creates a figure with the k-NN search curves and the
            %   epsilon estimate.
            %
            %   % Example:
            %   %   Create target data. Use the estimateEpsilon static 
            %   %   method to calculate an appropriate epsilon threshold 
            %   %   given the input data x.
            %
            %   % Set cluster parameters
            %   minNumPoints = 15;
            %   maxNumPoints = 20;
            %
            %   % Create target data
            %   x = [rand(20,2)+11.5; rand(20,2)+10;...
            %   rand(20,2)+15; 30*rand(5,2)];
            %
            %   % Estimate clustering threshold epsilon
            %   clusterDBSCAN.estimateEpsilon(x,minNumPoints,maxNumPoints)
            %
            % See also phased, clusterDBSCAN/discoverClusters,
            % clusterDBSCAN/plot.

            %   Copyright 2019 The MathWorks, Inc.
            
            % References
            % [1] Ester M., Kriegel H.-P., Sander J., and Xu X. "A
            %     Density-Based Algorithm for Discovering Clusters in Large
            %     Spatial Databases with Noise." Proc. 2nd Int. Conf. on
            %     Knowledge Discovery and Data Mining,  Portland, OR, AAAI
            %     Press, 1996, pp. 226-231.
            % [2] Schubert E., Sander J., Ester M., Kriegel H.-P., and Xu
            %     X. "DBSCAN Revisited, Revisited: Why and How You Should
            %     (Still) Use DBSCAN." ACM Transactions on Database 
            %     Systems, Vol. 42, No. 3, Article 19, 2017.
            
            %#codegen
            
            % Check inputs/outputs
            narginchk(3,3);
            if coder.target('MATLAB') % For MATLAB
                nargoutchk(0,1);
            else
                nargoutchk(1,1);
            end
            
            % Validate inputs
            insufficientDataFlag = phased.internal.AbstractClusterDBSCAN.validateInputsEstimateEpsilon(x,minNumPoints,maxNumPoints);
            if insufficientDataFlag
                if nargout == 1
                    varargout{1} = NaN;
                end
                return;
            end
            
            % Constants
            colorVec = [0.8500 0.3250 0.0980];
            axesPadding = 0.2;
            
            % Initialize figure for convenience plot
            if nargout == 0
                hFig = figure('Name','Estimated Epsilon');
                hAxKnnDist = axes(hFig);
                hold on;
                grid on;
            end
            
            % Calculate expected epsilon values for different k-values
            minNumPoints = min(minNumPoints,size(x,1));
            maxNumPoints = min(maxNumPoints,size(x,1));
            thisEpsEstimate = zeros(1,maxNumPoints-minNumPoints+1);
            D = phased.internal.AbstractClusterDBSCAN.calcPairwiseDist(x,x);
            knn = mink(D,maxNumPoints,2);
            vecLength = numel(knn); 
            nRows = size(knn,1);
            kdist = -1*ones(vecLength,1);
            for ii = minNumPoints:maxNumPoints
                lastIdx = (ii-1)*nRows;
                kdist(1:lastIdx) = reshape(knn(:,2:ii),lastIdx,1); % Remove first column since its the distance to itself
                kDistSort = sort(kdist);
                
                % Find first value that is not -1
                idxPos = vecLength-lastIdx+1;
                
                % Find knee
                idxKnee = phased.internal.AbstractClusterDBSCAN.findKnee(kDistSort(idxPos:end));
                thisEpsEstimate(ii-minNumPoints+1) = kDistSort(idxKnee+idxPos-1);
                
                % Update plot
                if nargout == 0
                    line(hAxKnnDist,(1:lastIdx).',kDistSort(idxPos:end),'Marker','.');
                    if ii == maxNumPoints
                        hEps = line(hAxKnnDist,...
                            idxKnee,thisEpsEstimate(ii-minNumPoints+1),...
                            'Marker','o','MarkerSize',5,...
                            'Color',colorVec,'MarkerFaceColor',colorVec);
                    else
                        line(hAxKnnDist,...
                            idxKnee,thisEpsEstimate(ii-minNumPoints+1),...
                            'Marker','o','MarkerSize',5,...
                            'Color',colorVec,'MarkerFaceColor',colorVec);
                    end
                    plotStr = sprintf('%d-NN',ii-1);
                    hText = text(hAxKnnDist,lastIdx,kDistSort(end),plotStr,...
                        'FontSize',8,'FontWeight','bold',...
                        'EdgeColor',[0 0 0],'BackgroundColor',[1 1 1]);
                    set(hText, 'Clipping', 'on');
                end
            end
            
            % Estimate epsilon
            epsilon = mean(thisEpsEstimate);
            
            % Plot estimated epsilon
            if nargout == 0
                axis(hAxKnnDist,'tight');
                xlims = get(hAxKnnDist,'XLim');
                hLine = line(hAxKnnDist,...
                    [xlims(1) xlims(2)+axesPadding],epsilon*ones(1,2),...
                    'LineStyle','--','Color',colorVec);
                legend([hEps hLine],{'Estimated Epsilon',...
                    'Time-Averaged Epsilon'},'Location','NorthWest');
                plotStr = sprintf('\\epsilon = %.3f',epsilon);
                hText = text(hAxKnnDist,xlims(1)+20,epsilon,plotStr,...
                    'HorizontalAlignment','left','Color',colorVec,...
                    'FontSize',10,'FontWeight','bold',...
                    'EdgeColor',[0 0 0],'BackgroundColor',[1 1 1]);
                set(hText, 'Clipping', 'on');
                
                % Update axes
                xlim(hAxKnnDist,[xlims(1) xlims(2)+axesPadding])
                ylims = get(hAxKnnDist,'YLim');
                ylim(hAxKnnDist,[ylims(1) ylims(2)+axesPadding])
                
                % Label axes
                xlabel('Index')
                ylabel('\epsilon')
                title('Estimated Epsilon')
            end
            
            % Set outputs
            if nargout == 1
                varargout{1} = epsilon;
            end
        end
    end
    
    methods (Access = private)
        function autoEpsilon(obj,x)
            % Estimate epsilon
            idxGd = all(~isinf(x),2); 
            thisEpsEstimate = obj.estimateEpsilon(x(idxGd,:),obj.MinNumPoints,obj.MaxNumPoints);
            
            % If an epsilon estimate was obtained, add to moving-average
            % buffer and take the mean value. Otherwise, use previous (or
            % default) estimate.
            if ~isnan(thisEpsEstimate)
                % Add to history
                if obj.pEpsIdx <= obj.EpsilonHistoryLength
                    obj.pEpsilonHistoryVec(obj.pEpsIdx) = thisEpsEstimate;
                    obj.pEpsIdx = obj.pEpsIdx+1;
                else
                    obj.pEpsilonHistoryVec = [obj.pEpsilonHistoryVec(2:end) thisEpsEstimate];
                end
                epsilon = mean(obj.pEpsilonHistoryVec(~isnan(obj.pEpsilonHistoryVec)));
                
                obj.pEpsilon = epsilon;
                
                if ~obj.pEpsilonIsInitialized
                    obj.pEpsilonIsInitialized = true;
                end
            end
        end
            
        function xUnambig = setupAmbiguousClustering(obj,x)
            % Unwrap the ambiguities. This method unwraps ambiguities if
            % they exist. The existence of ambiguities involves
            % finding detections in the boundary regions defined by the
            % maximum range and maximum Doppler +-
            
            % Codegen initialization
            [nDets,nDims] = size(x);
            
            % Support hyperellipse
            if isscalar(obj.pEpsilon)
                epsilonVec = obj.pEpsilon.*ones(1,nDims); 
            else
                epsilonVec = obj.pEpsilon; 
            end
            
            % Detect ambiguities
            obj.pAmbiguitiesDetected = [false false];
            if obj.EnableDisambiguation
                for ii = 1:numel(obj.AmbiguousDimension)
                    ambIdx = obj.AmbiguousDimension(ii); 
                    boundary1 = (x(:,ambIdx) >= obj.pAmbiguityLimits(ii,1)) & ...
                        (x(:,ambIdx) <= (obj.pAmbiguityLimits(ii,1) + epsilonVec(ambIdx)));
                    boundary2 = (x(:,ambIdx) >= obj.pAmbiguityLimits(ii,2) - epsilonVec(ambIdx)) & ...
                        (x(:,ambIdx) <= obj.pAmbiguityLimits(ii,2));
                    epsLessThanAmbLim = epsilonVec(ambIdx) < sum(abs(obj.pAmbiguityLimits(ii,:)));
                    if any(boundary1 | boundary2) && epsLessThanAmbLim % Detections exist at the edges and epsilon is less than the ambiguity limits
                         % If epsilon is greater than or equal to the
                         % ambiguity limits, there is no reason to unwrap
                         % the data, because all of the data will be
                         % clustered in that dimension.
                        obj.pAmbiguitiesDetected(ii) = true;
                    end
                end
            end
            
            % Set outputs
            if all(obj.pAmbiguitiesDetected) 
                xUnambig = [x; x; x; x]; 
                for ii = 1:numel(obj.AmbiguousDimension)
                    xUnambig(ii*nDets+1:(ii+1)*nDets,obj.AmbiguousDimension(ii)) = x(:,obj.AmbiguousDimension(ii))+sum(abs(obj.pAmbiguityLimits(ii,:)));
                    xUnambig(3*nDets+1:end,obj.AmbiguousDimension(ii)) = x(:,obj.AmbiguousDimension(ii))+sum(abs(obj.pAmbiguityLimits(ii,:)));
                end               
            elseif obj.pAmbiguitiesDetected(1) || obj.pAmbiguitiesDetected(2)
                xUnambig = [x; x]; 
                for ii = 1:numel(obj.AmbiguousDimension)
                    if obj.pAmbiguitiesDetected(ii)
                        xUnambig(nDets+1:end,obj.AmbiguousDimension(ii)) = x(:,obj.AmbiguousDimension(ii))+sum(abs(obj.pAmbiguityLimits(ii,:)));
                    end
                end
            else
                xUnambig = x;
            end
        end

        function [idx,corePts] =  dbscanHyperellipse(obj,x)
            % This is a version of DBSCAN that supports both n-spheres, as
            % well as hyperellipses for the neighborhood search. The
            % n-sphere version is the conventional DBSCAN algorithm, where
            % a single epsilon value is used for all dimensions. The basis
            % of this MATLAB implementation is Algorithm 2 from [2].
            
            % Size of the database
            [nDets,~] = size(x);
            
            if ~isscalar(obj.pEpsilon)
                distances = obj.calcEllipsoidDist(x,x,obj.pEpsilon) <= 1; 
            else
                % Scalar epsilon value for conventional DBSCAN
                distances = obj.calcPairwiseDist(x,x) <= obj.pEpsilon;
            end

            % Create a vector of nans for point labels
            idx = nan(nDets,1);
            
            % 1.) Identify core points
            % Compute neighbors of each point and identify core points.
            corePts = sum(double(distances),1) >= obj.MinNumPoints;
            nCorePts = sum(double(corePts));
            nOtherPts = nDets-nCorePts;
            
            % 2.) Assign core points
            % Join neighboring core points into clusters.
            corePtsLabels = nan(nDets,1);
            if nCorePts >= 1
                clustCount = 1;
                for p = 1:nCorePts
                    if isnan(corePtsLabels(p))
                        corePtsLabels = joinCorePoints(obj, p, nCorePts, corePtsLabels, distances(corePts,corePts), clustCount);
                        clustCount = clustCount+1;
                    end
                end
                idx(corePts) = corePtsLabels(1:nCorePts);
            end
            
            % 3.) Loop over non-core points
            otherPts = ~corePts;
            if any(corePts)
                for p = 1:nOtherPts
                    curP = obj.findFirstIdx(otherPts);
                    % 4.) Assign border points
                    % Add to a neighboring core point if possible. No need to
                    % expand cluster, since it has already been performed.
                    neighbors = distances(corePts,curP);
                    idxFirst = obj.findFirstIdx(neighbors);
                    if neighbors(idxFirst)
                        idx(curP) = corePtsLabels(idxFirst);
                    end
                    otherPts(curP) = false;
                end
            end
            
            % 5.) Assign noise points
            idx(isnan(idx)) = -1;
        end
             
        function corePtsLabels = joinCorePoints(obj, p, nCorePts, corePtsLabels, distances, clustCount)
            % This is a part of the DBSCAN algorithm. This is where the
            % clusters are expanded to form the maximal set of
            % density-connected points. The expansion is performed solely
            % on core points.
            
            seeds = distances(:,p);
            corePtsLabels(isnan(corePtsLabels(1:nCorePts)) & seeds) = clustCount; % Update idx to current cluster count if not already assigned to cluster or noise
            while any(seeds) % Loop over all points in seeds
                curP = obj.findFirstIdx(seeds); % Find first seed under consideration
                neighborsFound = distances(:,curP); % Get neighbors
                addSeed = isnan(corePtsLabels(1:nCorePts)) & neighborsFound; % If not labeled or noise and a neighbor
                seeds(addSeed) = true; % Set neighbors as seeds for consideration
                corePtsLabels(addSeed) = clustCount; % Assign to current cluster count
                seeds(curP) = false; % Remove current point from seed
            end
        end
        
        function idxOut = updateAmbigousClustering(obj,x,idx)
            % Undo the unwrapping of ambiguities.
            
            % Assign appropriate mapping indices based on ambiguities that
            % exist
            nPtsUnambig = size(x,1);
            if all(obj.pAmbiguitiesDetected)
                nPtsOrig = nPtsUnambig/4;
                idxOrig = 1:nPtsOrig;
                idxMapping = [idxOrig idxOrig idxOrig idxOrig];
                idxOut = obj.mapUnambigToAmbig(idx,idxMapping,nPtsOrig);
            elseif any(obj.pAmbiguitiesDetected)
                nPtsOrig = nPtsUnambig/2;
                idxOrig = 1:nPtsOrig;
                idxMapping = [idxOrig idxOrig];
                idxOut = obj.mapUnambigToAmbig(idx,idxMapping,nPtsOrig);
            else
                idxOut = idx;
            end
        end
        
        function parseAndValidatePlotInputs(obj,x,idx,varargin)
            % Validation method for public plot method
            
            % Verify x
            validateattributes(x,{'double'},{'2d','nonsparse','real',...
                'nonempty'},'','X');
            cond = any(isinf(x(:)));
            if cond
                coder.internal.errorIf(cond,'phased:clustering:expectedFinite','X');
            end
            [nDets,nCols] = size(x);
            
            % Number of clustering dimensions must be <= 3
            cond = nCols > 3;
            if cond
                coder.internal.errorIf(cond,'phased:phased:tooManyColumns','X',3);
            end
            
            % Input idx should be a column vector
            validateattributes(idx,{'double'},{'finite','nonsparse',...
                'nonnan','nonempty','vector','integer'},...
                '','idx');
            
            % Number of detections in x (number of rows) should be equal to
            % the length of idx
            nRowsIdx = numel(idx);
            cond = nRowsIdx~=nDets;
            if cond
                coder.internal.errorIf(cond,'phased:phased:NumRowsMismatch','X','IDX');
            end
            
            % Parse
            defaultAxes = [];
            if coder.target('MATLAB')
                % MATLAB parser
                pstruct = inputParser;
                pstruct.FunctionName = '';
                pstruct.addParameter('Parent',defaultAxes);
                pstruct.addParameter('Title',obj.pDefaultTitle);
                pstruct.parse(varargin{:});
                
                % Get values
                obj.pFigAxesNV = pstruct.Results.Parent;
                obj.pTitle = pstruct.Results.Title;
            else
                % Codegen parser
                poptions = struct( ...
                    'CaseSensitivity',false, ...
                    'PartialMatching','none', ...
                    'StructExpand',false, ...
                    'IgnoreNulls',true);
                
                parms = {'Parent','Title'};
                pstruct = coder.internal.parseParameterInputs(parms,poptions,varargin{:});
                
                % Get values
                obj.pFigAxesNV = coder.internal.getParameterValue(pstruct.Parent,defaultAxes,varargin{:});
                obj.pTitle = coder.internal.getParameterValue(pstruct.Title,obj.pDefaultTitle,varargin{:});
            end
            
            % Is this an axes handle?
            if ~isempty(obj.pFigAxesNV)
                if ~ishandle(obj.pFigAxesNV)
                    coder.internal.assert(false, ...
                        'phased:clustering:invalidAxesHandle','Parent');
                end
                
                % Is this a valid axes handle?
                if ~strcmpi(get(obj.pFigAxesNV,'Type'),'axes')
                    coder.internal.assert(false, ...
                        'phased:clustering:invalidAxesHandle','Parent');
                end
            end
            
            % Check title string
            sigdatatypes.checkString('', 'Title', obj.pTitle);
        end
        
        function addClusterLabels(obj,padVal,clusterXYZ,clusterIDs)
            % Add text labels to cluster plots
            nClusters = numel(clusterIDs);
            plotStr = cell(1,nClusters);
            for iC = 1:nClusters
                plotStr{iC} = sprintf('%d',clusterIDs(iC));
            end
            if ishandle(obj.pTextHandle)
                delete(obj.pTextHandle);
            end
            obj.pTextHandle = text(...
                clusterXYZ(:,1)+padVal,clusterXYZ(:,2)+padVal,clusterXYZ(:,3)+padVal,...
                plotStr,'FontSize',6,'HorizontalAlignment','left',...
                'Color',[0 0 0],'FontWeight','bold',...
                'BackgroundColor',[1 1 1],'EdgeColor',[0 0 0]);
            set(obj.pTextHandle,'Clipping','on'); 
        end
    end
    
    methods (Hidden,Static)
        function validateEpsilon(value)
            validateattributes(value,{'double'},{'2d','nonsparse','real',...
                'nonempty','finite','nonnan','positive'},'','Epsilon');
        end
        
        function D = calcPairwiseDist(X,Y)
            % Calculates the pairwise distance using the Euclidean distance
            % metric.
            
            nx = size(X,1);
            ny = size(Y,1);
            D = zeros(nx,ny);
            for i = 1:ny
                del = bsxfun(@minus,X,Y(i,:));
                dsq = sum((del) .^ 2, 2);
                dsq = sqrt(dsq);
                D(:,i) = dsq;
            end
        end
        
        function D = calcEllipsoidDist(X,Y,epsilonMatrix)
            % Calculates ellipsoidal distance. A point lies in an ellipsoid
            % created by a point Y if dsum is <= 1.
            
            nx = size(X,1);
            ny = size(Y,1);
            D = zeros(nx,ny);
            for i = 1:ny
                del = bsxfun(@minus,X,Y(i,:));
                if isrow(epsilonMatrix)
                    del = bsxfun(@rdivide, (del.^2), (epsilonMatrix.^2));
                else
                    del = bsxfun(@rdivide, (del.^2), (epsilonMatrix(i,:).^2));
                end
                dsum = sum(del,2);
                D(:,i) = dsum;
            end
        end
        
        function idx1 = findFirstIdx(x)
            % Find index of first true value. Assumes true value exists.
            % Use in lieu of find for codegen.
            
            idx1 = 1;
            for ii = 1:length(x)
                if x(ii)
                    idx1 = ii;
                    break;
                end
            end
        end

        function idxOut = mapUnambigToAmbig(idx,idxMapping,nPtsOrig)
            % Maps the unambiguous unwrapped detections to the original
            % ambiguous detections.
            
            % Find unique cluster IDs
            clusterIdx = unique(idx);
            
            % Build cluster matrix
            nClusters = length(clusterIdx);
            clusterMatrix = zeros(nPtsOrig,nClusters);
            nPtsInCluster = zeros(1,nClusters);
            for iC = 1:nClusters
                idxThisCluster = idx == iC;
                nPtsInCluster(iC) = sum(double(idxThisCluster));
                clusterMatrix(idxMapping(idxThisCluster),iC) = 1;
            end
            
            % Find largest cluster for each detection
            idxOut = -1*ones(nPtsOrig,1);
            idxOrig = 1:nPtsOrig;
            for iD = 1:nPtsOrig
                assocClustIdx = clusterMatrix(iD,:) == 1;
                if any(assocClustIdx)
                    [~,idxMax] = max(nPtsInCluster(assocClustIdx),[],2);
                    
                    % Find max value within logical array
                    count = 0;
                    for iD2 = 1:nPtsOrig
                        if assocClustIdx(iD2)
                            count = count+1;
                        end
                        if count == idxMax
                            idxMax = iD2;
                            break;
                        end
                    end
                    idxOut(iD) = idxOrig(idxMax);
                end
            end
            
            % Fix cluster numbering
            ambigClustIdx = unique(idxOut);
            clusterCount = 1; 
            for iC = 1:length(ambigClustIdx)
                if ambigClustIdx(iC)~=-1
                    idxOut(idxOut == ambigClustIdx(iC)) = clusterCount;
                    clusterCount = clusterCount+1; 
                end
            end
        end
        
        function insufficientDataFlag = validateInputsDiscoverClusters(x,maxEpsilon,minNumPoints)
            % Validate input x
            validateattributes(x,{'double'},{'2d','nonsparse','real',...
                'nonempty','nonnan','finite'},'','X');
            
            % Validate minimum number of points
            validateattributes(minNumPoints,{'double'},{'scalar','integer','real',...
                'nonempty','finite','nonnan','positive'},'','MINNUMPOINTS');
            
            % Check that x has at least minNumPoints detections. If it doesn't, return
            % to whatever called this function.
            nDets = size(x,1);
            cond = nDets<minNumPoints;
            insufficientDataFlag = false;
            if cond
                insufficientDataFlag = true;
                coder.internal.warning('phased:clustering:insufficientNumObservations','X','MINNUMPOINTS');
            end
            
            % Validate maximum epsilon value
            validateattributes(maxEpsilon,{'double'},{'scalar','real',...
                'nonempty','nonnan','positive'},'','MAXEPSILON');
            
            % minNumPoints >= 2
            cond = minNumPoints<2;
            if cond
                coder.internal.errorIf(cond,'phased:phased:expectedGreaterThanOrEqualTo','MINNUMPOINTS',2);
            end
        end
        
        function [reachDist,order] = optics(x,maxEpsilon,minNumPoints)
            % Performs the OPTICS algorithm. Aids in the discovery of clusters.
            
            % Initialization
            m = size(x,1); % m = number of measurements
            reachDist = inf(1,m); % Reachability distance
            order = zeros(1,m); % Output order
            
            % Calculate core distances
            D = phased.internal.AbstractClusterDBSCAN.calcPairwiseDist(x,x);
            knn = mink(D,minNumPoints,2);
            coreDist = max(knn(:,2:end),[],2); % Remove first column since its the distance to itself
            coreDist(coreDist>maxEpsilon) = inf;
            
            % Check to see if the maximum epsilon value was sufficient
            cond = all(isinf(coreDist));
            if cond
                coder.internal.warning('phased:clustering:reachDistsAreInf','MAXEPSILON');
            end
            
            % Calculate reachability distances
            seeds = true(1,m);
            idxSeed = 1;
            seeds(1) = false;
            order(1) = 1;
            idxOrder = 1;
            mm = zeros(1,m);
            idxSeeds = nan(1,m);
            idxMeas = 1:m; 
            while any(seeds)
                mm(seeds) = max([ones(1,sum(double(seeds)))*coreDist(idxSeed); D(idxSeed,seeds)],[],1);
                ii = false(1,m);
                ii(seeds) = reachDist(seeds) > mm(seeds);
                reachDist(ii) = mm(ii);
                [~,idxSeed] = min(reachDist(seeds),[],2);
                nSeeds = sum(double(seeds),2); 
                idxSeeds(1:nSeeds) = idxMeas(seeds);  
                idxSeed = idxSeeds(idxSeed); 
                idxOrder = idxOrder + 1;
                seeds(idxSeed) = false;
                order(idxOrder) = idxSeed;
            end
            
            % Order seeds
            reachDist(1) = inf;
        end
        
        function insufficientDataFlag = validateInputsEstimateEpsilon(x,minNumPoints,maxNumPoints)
            % Validate input x
            validateattributes(x,{'double'},{'2d','nonsparse','real',...
                'nonempty','nonnan','finite'},'','X');
            
            % Validate minimum number of points
            validateattributes(minNumPoints,{'double'},{'scalar','integer','real',...
                'nonempty','finite','nonnan','positive'},'','MINNUMPOINTS');
            
            % Check that x has at least minNumPoints. If it doesn't, return to
            % whatever called this function.
            nDets = size(x,1);
            cond = nDets<minNumPoints;
            insufficientDataFlag = false;
            if cond
                insufficientDataFlag = true;
                coder.internal.warning('phased:clustering:insufficientNumObservations','X','MINNUMPOINTS');
            end
            
            % Validate maximum number of points
            validateattributes(maxNumPoints,{'double'},{'scalar','integer','real',...
                'nonempty','nonnan','positive'},'','MAXNUMPOINTS');
            
            % minNumPoints must be less than or equal to maxNumPoints
            cond = minNumPoints>maxNumPoints;
            if cond
                coder.internal.errorIf(cond,'phased:phased:expectedLessThanOrEqualTo','MINNUMPOINTS','MAXNUMPOINTS');
            end
            
            % minNumPoints >= 2
            cond = minNumPoints<2;
            if cond
                coder.internal.errorIf(cond,'phased:phased:expectedGreaterThanOrEqualTo','MINNUMPOINTS',2);
            end
        end
        
        function idxK = findKnee(Y)
            % Estimates the knee as the point farthest from the line that
            % connects the first and last points in the curve Y.
            n = length(Y);
            pts = [(1:n).' Y];
            firstPt = pts(1,:);
            lastPt = pts(end,:);
            vecFirstLast = lastPt - firstPt;
            vecFirstLast = vecFirstLast./norm(vecFirstLast);
            vecFromFirst = bsxfun(@minus,pts,firstPt);
            scalarDot = sum(bsxfun(@times,vecFromFirst,vecFirstLast),2);
            vecToLine = vecFromFirst-scalarDot*vecFirstLast;
            distToLine = vecnorm(vecToLine,2,2);
            [~,idxK] = max(distToLine,[],1);
        end
    end
end



