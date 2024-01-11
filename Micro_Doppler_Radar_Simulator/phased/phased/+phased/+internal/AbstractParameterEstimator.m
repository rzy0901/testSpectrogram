classdef (Hidden) AbstractParameterEstimator < matlab.System &...
        matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%   Copyright 2016 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %NumEstimatesSource  Source of the number of estimates reported
        %   Specify the source of the number of estimates as one of 'Auto'
        %   | 'Property'. The default is 'Auto'. If you set this property
        %   to 'Auto', the number of estimates is determined from the
        %   number of columns in the input argument DETIDX. If cluster IDs
        %   are provided, the number of estimates is determined from the
        %   number of unique cluster IDs. If you set this property to
        %   'Property', the number of reported estimates is determined by
        %   the value of the NumEstimates property.
        NumEstimatesSource = 'Auto'
    end
    
    properties (Nontunable, PositiveInteger)
        %NumEstimates  Maximum number of estimates
        %   Specify the maximum number of estimates to report as a positive
        %   integer scalar. The default value is 1. This property applies
        %   when you set the NumEstimatesSource property to 'Property'.
        NumEstimates = 1
    end
    
    properties (Nontunable, Logical)
        %ClusterInputPort  Enable cluster ID input
        %   Set this property to true to input cluster IDs as an input
        %   argument. The cluster IDs are used to map detections in the
        %   response data to hypothetical targets, reducing the total
        %   number of reported estimates. The default value of this
        %   property is false.
        ClusterInputPort = false;
    end
    
    properties (Nontunable, Logical)
        %VarianceOutputPort  Output variance for parameter estimates
        %   Set this property to true to output variances for the parameter
        %   estimates. Set this property to false to not output the 
        %   estimated variances. The default value of this property is
        %   false.
        VarianceOutputPort = false;
    end
    
    properties (Nontunable)
        %NoisePowerSource  Source of noise power
        %   Specify how the object determines the noise power values to be
        %   used for variance estimation as one of 'Property' | 'Input
        %   port'. The default value is 'Property'. When set to 'Property',
        %   the NoisePower property must be set to represent the noise
        %   power at the detection locations. When set to 'Input port', the
        %   noise power is specified through an input argument. This
        %   property only applies when you set the VarianceOutputPort
        %   property to true.
        NoisePowerSource = 'Property'
    end

    properties (Nontunable)
        %NoisePower  Noise power
        %   Specify a scalar value representing the noise power at the
        %   detection locations. NOISEPOWER units are linear (not decibels)
        %   and must use the same scale factor as the response data. The
        %   same noise power value is applied to all detections. The
        %   default value of this property is 1. This property only applies
        %   when you set the VarianceOutputPort property to true.
        NoisePower = 1
    end
    
    properties (Access = protected, Hidden=true, PositiveInteger)
        pDimension = 1
    end
    
    properties (Access = protected, Hidden=true)
        pCubeSize = [-1 -1 -1]
    end
    
    properties(Constant, Hidden)
        NumEstimatesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    properties(Constant, Hidden)
        NoisePowerSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    methods (Access = protected, Abstract)
        setParameterDimension(obj)
        var = estimateVariance(obj,snr,paramgrid)
    end
    
    methods 
        function obj = AbstractParameterEstimator(varargin)
            % ParameterEstimator Constructor for phased.ParameterEstimator class
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods
        function set.NoisePower(obj,val)
            validateattributes( val, { 'double','single' }, {'real','scalar','positive'},'','NOISEPOWER');
            cond = any(isinf(val));
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:step:expectedNonInf','NOISEPOWER');
            end
            obj.NoisePower = val;
        end
    end
    
    methods (Access = protected)
                
        function num = getNumOutputsImpl(obj)
            num = 2;
            
            if ~obj.VarianceOutputPort
                num = num-1;
            end
        end
        
        function validateInputsImpl(obj,resp,~,detidx,varargin)
            szResp = size(resp);
            
            cond = ~isa(resp,'float');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','RESP','float');
            end
            
            % resp must have 1 to 3 dimensions
            cond = numel(szResp) > 3 || isempty(resp);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:MustBe1D2D3DNonEmpty','RESP');
            end
            
            % detidx
            if iscolumn(resp)
                validateattributes(detidx,{'double','single'},{'real','row'},'','DETIDX');
            else
                validateattributes(detidx,{'double','single'},{'real','nrows',ndims(resp)},'','DETIDX');
            end
            validateattributes(detidx(~isnan(detidx(:))),{'double','single'},{'integer'},'','DETIDX');
            
            iVarg = 1;
            if obj.VarianceOutputPort && strcmp(obj.NoisePowerSource,'Input port')
                noisepower = varargin{iVarg};
                iVarg = iVarg+1;
                
                % Can be a scalar or a row vector. When a vector, length is
                % dependent on the size of an input argument, so that is
                % checked in step method
                %
                % Can only be positive valued, but during initialization,
                % in Simulink, zeros are propagated to validateInputsImpl,
                % so positive check needs to be done in stepImpl instead.
                validateattributes(noisepower,{'double','single'},{'real','row'},'','NOISEPOWER');
            end
            
            if obj.ClusterInputPort
                clusterIDs = varargin{iVarg};
                
                validateattributes(clusterIDs,{'double','single'},{'real','row','numel',size(detidx,2)},'','CLUSTERIDS');
            end
            
        end
        
        function processInputSizeChangeImpl(obj,resp,~,~,varargin)
            szResp = size(resp);
            if numel(szResp)==3
                obj.pCubeSize(:) = szResp;
            else
                obj.pCubeSize(:) = [szResp -1];
            end
            setParameterDimension(obj);
        end
        
        function setupImpl(obj,resp,~,~,varargin)
            szResp = size(resp);
            if numel(szResp)==3
                obj.pCubeSize(:) = szResp;
            else
                obj.pCubeSize(:) = [szResp -1];
            end
            setParameterDimension(obj);
        end
        
        function flag = isInputComplexityLockedImpl(~,index)
            flag = false;
            
            % Only the first input (the signal response) can be complex
            if index>1
                flag = true;
            end
        end

        function flag = isOutputComplexityLockedImpl(~,~) 
            flag = false;
        end

        function flag = isInputSizeLockedImpl(~,~)
            flag = false;
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            switch prop
                case 'NoisePowerSource'
                    if ~obj.VarianceOutputPort
                        flag = true;
                    end
                case 'NoisePower'
                    if ~obj.VarianceOutputPort || ...
                            strcmp(obj.NoisePowerSource,'Input port')
                        flag = true;
                    end
                case 'NumEstimates'
                    if strcmp(obj.NumEstimatesSource,'Auto')
                        flag = true;
                    end
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.pDimension = obj.pDimension;
                s.pCubeSize = obj.pCubeSize;
            end
        end

        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function [est,var] = stepImpl(obj,resp,paramgrid,detidx,varargin)
            
            classtoUse = class(resp);
            paramgrid = cast(paramgrid,classtoUse);
            szResp = size(resp);
            dim = obj.pDimension;
            dimLen = szResp(dim);
            numDet = size(detidx,2);
            
            % Perform validation on inputs that must be positive. In
            % Simulink, zeros are propagated, so positive validation must
            % be performed inside stepImpl
            
            % Validate detection index to be positive
            cond = any(detidx(:)<=0);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:step:expectedPositive','DETIDX');
            end
                
            iVarg = 1;
            if obj.VarianceOutputPort
                
                % Grab noise power from property or input
                if strcmp(obj.NoisePowerSource,'Property')
                    noisepowerIn = obj.NoisePower;
                else
                    noisepowerIn = varargin{iVarg};
                    iVarg = iVarg+1;
                end
                
                % Perform validation on inputs that must be positive. In
                % Simulink, zeros are propagated, so positive validation
                % must be performed inside stepImpl
                
                % Validate noise power to be positive
                cond = any(noisepowerIn<=0);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:step:expectedPositive','NOISEPOWER');
                end

                cond = any(isinf(noisepowerIn));
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:step:expectedNonInf','NOISEPOWER');
                end
                
                if ~isscalar(noisepowerIn)
                    cond = ~isrow(noisepowerIn) || (numel(noisepowerIn) ~= numDet);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:step:ScalarOrVectorNCols','NOISEPOWER',numDet);
                    end
                end
                
                % Always a vector
                noisepowerVect = NaN(numDet,1,'like',noisepowerIn);
                noisepowerVect(:) = noisepowerIn;
            end
            
            if obj.ClusterInputPort
                clusterIDs = varargin{iVarg};
            else
                clusterIDs = (1:numDet)';
            end
            
            % Remove undefined detection indices
            iGd = ~isnan(detidx(1,:));
            detidx = detidx(:,iGd);
            clusterIDs = clusterIDs(iGd);
            if obj.VarianceOutputPort
                noisepower = noisepowerVect(iGd);
            end
           
            clusters = unique(clusterIDs);
            numClusters = numel(clusters);
            if strcmp(obj.NumEstimatesSource,'Auto')
                numOut = numClusters;
            else
                numOut = obj.NumEstimates;
            end
            
            est = NaN(numOut,1,classtoUse);
            if obj.VarianceOutputPort
                snr = NaN(numOut,1,classtoUse);
                var = NaN(numOut,1,classtoUse);
            end
            
            numEst = min(numClusters,numOut);
            for m = 1:numEst
                thiscluster = find(clusterIDs==clusters(m));
                
                % Find peak in cluster
                switch size(detidx,1)
                    case 1
                        idx = detidx(thiscluster);
                    case 2
                        idx = sub2ind(szResp,detidx(1,thiscluster),detidx(2,thiscluster));
                    case 3
                        idx = sub2ind(szResp,detidx(1,thiscluster),detidx(2,thiscluster),detidx(3,thiscluster));
                end
                
                y = resp(idx);
                [mVal,iMax] = max(abs(y(:)));
                [idx1,idx2,idx3] = ind2sub(szResp,idx(iMax));
                
                % Select cut to estimate parameter along
                switch dim
                    case 2
                        idx2 = idx2+(-1:1)';
                        idx2 = idx2(idx2>0&idx2<=dimLen);
                        idx1 = repmat(idx1,numel(idx2),1);
                        idx3 = repmat(idx3,numel(idx2),1);
                        
                        x = paramgrid(idx2);
                    case 3
                        idx3 = idx3+(-1:1)';
                        idx3 = idx3(idx3>0&idx3<=dimLen);
                        idx1 = repmat(idx1,numel(idx3),1);
                        idx2 = repmat(idx2,numel(idx3),1);
                        
                        x = paramgrid(idx3);
                    otherwise % 1
                        idx1 = idx1+(-1:1)';
                        idx1 = idx1(idx1>0&idx1<=dimLen);
                        idx2 = repmat(idx2,numel(idx1),1);
                        idx3 = repmat(idx3,numel(idx1),1);
                        
                        x = paramgrid(idx1);
                end
                iEst = sub2ind(szResp,idx1,idx2,idx3);
                y = abs(resp(iEst));
                
                % Estimate parameter
                if numel(y)<3 % Centroid interpolation
                    est(m) = sum(x(:).*y(:))/sum(y);
                else % Quadratic interpolation
                    den = 2*y(2)-y(1)-y(3);
                    if den~=0
                        delta = 0.5*(y(3)-y(1))/den;
                        if delta>0.5
                            delta = cast(0.5,classtoUse);%%casting due to codegen datatype mismatch
                        elseif delta<-0.5
                            delta = cast(-0.5,classtoUse);%%casting due to codegen datatype mismatch
                        end
                    else
                        delta = cast(0,classtoUse);%%casting due to codegen datatype mismatch
                    end
                    est(m) = x(2)+delta*(x(2)-x(1));
                end
                
                if obj.VarianceOutputPort
                    snr(m) = mVal^2/noisepower(thiscluster(iMax));
                end
            end
            if obj.VarianceOutputPort
                var(1:numEst) = estimateVariance(obj,snr(1:numEst),paramgrid);
            end            
        end
    end
  
    methods (Access = protected)
        function varargout = getOutputNamesImpl(~)
            varargout = {'Est','Var'};
        end
      
        function varargout = getOutputSizeImpl(obj)
           varargout{1} = [obj.NumEstimates 1];
           varargout{2} = [obj.NumEstimates 1];
        end
        function varargout = isOutputFixedSizeImpl(~)
            varargout{1} = true;
            varargout{2} = true;
        end
        function varargout = getOutputDataTypeImpl(obj)
            varargout{1} = propagatedInputDataType(obj,1);
            varargout{2} = propagatedInputDataType(obj,1);
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = false;
            varargout{2} = false;
        end
    end
end
