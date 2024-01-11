classdef (Sealed, StrictDefaults) SimulinkClusterDBSCAN < phased.internal.AbstractClusterDBSCAN
%This class is for internal use only. It may be removed in the future.

% clusterDBSCAN    Density-based algorithm for clustering data
%   H = clusterDBSCAN creates a system object, H, that can be used to
%   cluster data using the DBSCAN (Density-Based Spatial Clustering of
%   Applications with Noise) algorithm. The system object clusterDBSCAN is
%   capable of clustering any type of data and has features to assist with
%   epsilon parameter estimation and disambiguation.
%
%   H = clusterDBSCAN(Name,Value) creates a cluster object, H, with the
%   specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   step method syntax:
%
%   [IDX,CLUSTERIDS] = step(H,X) partitions the points in the NxP data
%   matrix X into clusters based on the Epsilon and MinNumPoints
%   properties, where N is the number of points (or detections), and P is
%   the clustering dimensions. The DBSCAN algorithm is capable of
%   clustering any type of data with appropriate MinNumPoints and Epsilon
%   settings.
%
%   Epsilon is a positive scalar that is a threshold for a neighborhood
%   search query. Alternatively, Epsilon can be designated as a 1xP row
%   vector. In the case of the row vector, the epsilons for the different
%   clustering dimensions create a hyperellipse search area.
%
%   MinNumPoints is a positive integer that defines the minimum number of
%   points in a cluster. The purpose of MinNumPoints is to smooth the
%   density estimate. Since a DBSCAN cluster is a maximal set of
%   density-connected points, it is best to err on the side of a lower
%   value for MinNumPoints when the expected number of detections in a
%   cluster is unknown, even though lower values make the DBSCAN algorithm
%   more susceptible to noise. A general guideline for MinNumPoints
%   selection is as follows:
%   MinNumPoints = 2*P,
%   where P is the number of clustering dimensions in X. For challenging
%   data sets with the following qualities:
%      - highly noisy,
%      - very large N,
%      - high dimensionality P,
%      - or many duplicates
%   increasing MinNumPoints may improve the DBSCAN clustering results.
%
%   IDX is an Nx1 column vector containing cluster IDs that represents the
%   clustering results of the DBSCAN algorithm. A value equal to '-1'
%   implies a DBSCAN noise point. Positive IDX values correspond to
%   clusters that satisfy the DBSCAN clustering criteria. DBSCAN is
%   calculated using the Euclidean distance metric.
%
%   CLUSTERIDS is a 1xN row vector of positive integers. Each value in
%   CLUSTERIDS is a unique identifier indicating a hypothetical target
%   cluster. In contrast to the DBSCAN output IDX that assigns noise points
%   an ID of '-1', unique positive cluster IDs are created for all points.
%   The CLUSTERIDS output was created for compatibility with other Phased
%   Array System Toolbox system objects such as phased.RangeEstimator and
%   phased.DopplerEstimator.
%
%   The system object clusterDBSCAN has features to assist with parameter
%   selection and disambiguation. Any combination of the features can be
%   used or not depending on the following property settings and positional
%   inputs to the step method:
%
%   1) Ambiguous detections
%   [IDX,CLUSTERIDS] = step(H,X,AMBLIMS) is a syntax that is enabled if the
%   EnableDisambiguation property is set to true. Clustering will occur
%   across boundaries to ensure that ambiguous detections are appropriately
%   clustered for up to 2 dimensions. The ambiguous columns of X must be
%   defined using the AmbiguousDimension property. AMBLIMS defines the
%   minimum and maximum ambiguity limits in the same units as the
%   AmbiguousDimension columns of the input data matrix X. The AMBLIMS
%   input may be specified as either a 1x2 vector for a single ambiguous
%   dimension where the fields are [MinAmbiguityLimitDimension1,
%   MaxAmbiguityLimitDimension1] or a 2x2 matrix for 2 ambiguous dimensions
%   where the fields are [MinAmbiguityLimitDimension1,
%   MaxAmbiguityLimitDimension1; MinAmbiguityLimitDimension2,
%   MaxAmbiguityLimitDimension2].
%
%   2) Epsilon parameter estimation
%   [IDX,CLUSTERIDS] = step(H,X,UPDATE) is a syntax that is enabled when
%   the EpsilonSource property is set to 'Auto'. This syntax permits the
%   automatic estimation of the clustering threshold epsilon from the input
%   data matrix X using a k-Nearest Neighbor (k-NN) search. This syntax
%   requires the property MaxNumPoints, which defines the maximum number of
%   points expected in a cluster. The k-NN search is calculated with k
%   equal to (MinNumPoints-1):(MaxNumPoints-1), where the '-1' corresponds
%   to a cluster consisting of k number of points, including the point
%   under consideration as a core point. For a single k-NN search curve,
%   set MinNumPoints equal to MaxNumPoints.
%
%   UPDATE is a logical that defines when the step method should update its
%   epsilon estimate. If UPDATE is true, the epsilon threshold is first
%   estimated as the average of the knees of the k-NN search curves. The
%   estimate is then added to a buffer whose length L is defined by the
%   EpsilonHistoryLength property. The final epsilon that is used is
%   calculated as the average of the L-length epsilon history buffer.
%   EpsilonHistoryLength set to 1 is memory-less. Memory-less means that
%   each epsilon estimate is immediately used and no moving-average
%   smoothing occurs. If UPDATE is false a previous epsilon estimate is
%   used. Estimating epsilon is computationally intensive. This syntax is
%   not recommended for large data sets.
%
%   It is recommended to set MinNumPoints and MaxNumPoints according to a
%   range of values around the expected number of detections to ensure the
%   best epsilon estimates. If N < MinNumPoints, there is insufficient
%   data, and a previous epsilon estimate will be used. If N <
%   MaxNumPoints, the k-NN search is limited to (MinNumPoints-1):(N-1).
%
%   [IDX,CLUSTERIDS] = step(H,AMBLIMS,UPDATE) is the full syntax with all
%   enabling properties (EnableDisambiguation is set to true and
%   EpsilonSource is set to 'Auto').
%
%   clusterDBSCAN methods:
%
%   step             - Cluster data in the input NxP data matrix
%   discoverClusters - Discover cluster hierarchy using OPTICS
%   estimateEpsilon  - Estimate cluster threshold from data
%   plot             - Plot DBSCAN clustering results
%   release  - Allow property value and input characteristics changes
%   clone    - Create a clustering object with the same property values
%   isLocked - Locked status (logical)
%   reset    - Reset the clustering object to its initial state
%
%   clusterDBSCAN properties:
%
%   EpsilonSource        - Source of cluster threshold
%   Epsilon              - Cluster threshold
%   MinNumPoints         - Minimum number of points in a cluster
%   MaxNumPoints         - Maximum number of points in a cluster
%   EpsilonHistoryLength - Length of epsilon history
%   EnableDisambiguation - Enable disambiguation of dimensions
%   AmbiguousDimension   - Indices of ambiguous dimensions
%
%   % Example:
%   %   The maximum unambiguous range and Doppler in this example is 20 m
%   %   and 30 Hz, respectively. Load data with the following extended
%   %   targets and false alarms:
%   %      - 1 unambiguous target located at (10, 15)
%   %      - 1 ambiguous target in Doppler located at (10, -30)
%   %      - 1 ambiguous target in range located at (20, 15)
%   %      - 1 ambiguous target in range and Doppler located at (20, 30) 
%   %      - 5 false alarms
%   %   Create a clusterDBSCAN system object and specify that
%   %   clustering should be performed across range and Doppler boundaries
%   %   due to the existence of the ambiguities. Plot the DBSCAN clustering
%   %   results with the plot method. 4 targets are identified.
% 
%   % Load data 
%   load('dataClusterDBSCAN.mat');
%   
%   % Create DBSCAN system object that clusters across range and Doppler
%   % boundaries
%   clusterRadarData = clusterDBSCAN('MinNumPoints',3,'Epsilon',2,...
%       'EnableDisambiguation',true,'AmbiguousDimension',[1 2]); 
%   
%   % Cluster with DBSCAN
%   ambLims = [0 maxRange; ...
%       minDoppler maxDoppler]; 
%   idx = clusterRadarData(x,ambLims);
% 
%   % Plot clustering results
%   plot(clusterRadarData,x,idx); 
%
%   See also phased, phased.CFARDetector, phased.RangeEstimator,
%   phased.DopplerEstimator.
    
%   Copyright 2019 The MathWorks, Inc.

% References
% [1] Ester M., Kriegel H.-P., Sander J., and Xu X. "A Density-Based
%     Algorithm for Discovering Clusters in Large Spatial Databases with
%     Noise." Proc. 2nd Int. Conf. on Knowledge Discovery and Data Mining,
%     Portland, OR, AAAI Press, 1996, pp. 226-231.
% [2] Schubert E., Sander J., Ester M., Kriegel H.-P., and Xu X. "DBSCAN
%     Revisited, Revisited: Why and How You Should (Still) Use DBSCAN." ACM
%     Transactions on Database Systems, Vol. 42, No. 3, Article 19, 2017.
% [3] Wagner T., Feger R., and Stelzer A. "A Fast Grid-Based Clustering
%     Algorithm for Range/Doppler/DoA Measurements." Proc. of the 13th
%     European Radar Confer., London, UK, 2016, pp. 105-108.
% [4] Wagner T., Feger R., and Stelzer A. "Modification of DBSCAN and
%     Application to Range/Doppler/DoA Measurements for Pedestrian
%     Recognition with an Automotive Radar System." Proc. of the 12th
%     European Radar Conf., Paris, France, 2015, pp. 269-272.
    
%#codegen

    properties (Nontunable)
        %BlockOutputs       Define outputs for Simulink block
        %   Specify the output of the Simulink block as 'Index' | 'Cluster
        %   ID' | 'Index and Cluster ID'. 'Index' specifies the output as
        %   an Nx1 column vector containing cluster IDs that represents the
        %   clustering results of the DBSCAN algorithm. A value of '-1'
        %   implies a DBSCAN noise point. Positive indices correspond to
        %   clusters that satisfy the DBSCAN clustering criteria. 'Cluster
        %   ID' specifies the output as a 1xN row vector of positive
        %   integers. Each value is a unique identifier indicating a
        %   hypothetical target cluster. In contrast to the 'Index' output
        %   that assigns noise points a value of '-1', unique positive
        %   cluster IDs are created for all points. Defaults to 'Index and
        %   ID'.
        BlockOutputs = 'Index and Cluster ID'
    end

    properties(Constant, Hidden)
        BlockOutputsSet = matlab.system.StringSet({...
            'Index','Cluster ID','Index and Cluster ID'} );
    end
    
    methods
        function obj = SimulinkClusterDBSCAN(varargin)
            obj@phased.internal.AbstractClusterDBSCAN(varargin{:});
            
            if isscalar(obj.Epsilon)
                validateattributes(obj.Epsilon,{'double'},{'scalar','nonsparse','real',...
                    'nonempty','finite','nonnan','positive'},'','Epsilon');
            else
                validateattributes(obj.Epsilon,{'double'},{'row','nonsparse','real',...
                    'nonempty','finite','nonnan','positive'},'','Epsilon');
            end
        end
    end
    
    methods (Access=protected)
        function num = getNumOutputsImpl(obj)
            if strcmp(obj.BlockOutputs,'Index') || strcmp(obj.BlockOutputs,'Cluster ID')
                num = 1;
            else
                num = 2;
            end
        end
        
        function varargout = getOutputNamesImpl(obj)
            % Return output port names for System block
            if strcmp(obj.BlockOutputs,'Index')
                varargout{1} = 'Idx';
            elseif strcmp(obj.BlockOutputs,'Cluster ID')
                varargout{1} = 'Clusters';
            else
                varargout{1} = 'Idx';
                varargout{2} = 'Clusters';
            end
        end
        
        function varargout = getOutputSizeImpl(obj)
            % Return output port names for System block
            sz = propagatedInputSize(obj,1);
            if strcmp(obj.BlockOutputs,'Index')
                varargout{1} = [sz(1) 1];
            elseif strcmp(obj.BlockOutputs,'Cluster ID')
                varargout{1} = [1 sz(1)];
            else
                varargout{1} = [sz(1) 1];
                varargout{2} = [1 sz(1)];
            end
        end
        
        function varargout = getOutputDataTypeImpl(obj)
            % Return data type for each output port
            if strcmp(obj.BlockOutputs,'Index')
                varargout{1} = "double";
            elseif strcmp(obj.BlockOutputs,'Cluster ID')
                varargout{1} = "double";
            else
                varargout{1} = "double";
                varargout{2} = "double";
            end
        end
        
        function varargout = isOutputComplexImpl(obj)
            % Return true for each output port with complex data
            if strcmp(obj.BlockOutputs,'Index')
                varargout{1} = false;
            elseif strcmp(obj.BlockOutputs,'Cluster ID')
                varargout{1} = false;
            else
                varargout{1} = false;
                varargout{2} = false;
            end
        end
        
        function varargout = isOutputFixedSizeImpl(obj)
            % Return true for each output port with fixed size
            if strcmp(obj.BlockOutputs,'Index')
                varargout{1} = propagatedInputFixedSize(obj,1);
            elseif strcmp(obj.BlockOutputs,'Cluster ID')
                varargout{1} = propagatedInputFixedSize(obj,1);
            else
                varargout{1} = propagatedInputFixedSize(obj,1);
                varargout{2} = propagatedInputFixedSize(obj,1);
            end
        end
        
        function varargout = stepImpl(obj,x,varargin)
            [idx, clusterIDs] = clusterdataDBSCAN(obj,x,varargin{:});
            
            % Set outputs
            if strcmp(obj.BlockOutputs,'Index')
                varargout{1} = idx;
            elseif strcmp(obj.BlockOutputs,'Cluster ID')
                varargout{1} = clusterIDs;
            else
                varargout{1} = idx;
                varargout{2} = clusterIDs;
            end
        end
    end
    
    methods (Hidden,Static)
        function validateEpsilon(value)
            if isscalar(value)
                validateattributes(value,{'double'},{'scalar','nonsparse','real',...
                    'nonempty','finite','nonnan','positive'},'','Epsilon');
            else
                validateattributes(value,{'double'},{'row','nonsparse','real',...
                    'nonempty','finite','nonnan','positive'},'','Epsilon');
            end
        end
    end
end
