classdef (Hidden) AbstractDOA  < phased.internal.AbstractNarrowbandArrayProcessing & ...
        matlab.system.mixin.CustomIcon
    
%This class is for internal use only. It may be removed in the future.

%AbstractDOA Abstract class for narrowband DOA estimation

%   Copyright 2009-2018 The MathWorks, Inc. 
    
%   References:
%   [1] Optimum Array Processing, Part IV of Detection, Estimation, and
%       Modulation Theory, Harry Van Trees 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties(Access = protected, Nontunable)
    %private handle to covariance matrix estimator
    cCovEstimator
end

properties(Access = protected)
    %private handle for number of data snapshots
    pNumSnapshots;
end

methods(Access = protected)

    function obj = AbstractDOA(varargin)
        obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
    end    
    
    function privValidateSensorArray(obj,val) 
        if isSubarraySupported(obj)
            validateattributes( val, { 'phased.internal.AbstractArray',...
                'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
        else
            validateattributes( val, { 'phased.internal.AbstractArray'},...
                { 'scalar' }, '', 'SensorArray');
        end
    end
    
    function flag = isSubarraySupported(obj) %#ok<MANU>
        flag = false;
    end

    function validateInputsImpl(obj,x)
        % check the dimension of the input
        sz_x = size(x); 
        act_chnl = sz_x(2);
        exp_chnl = getDOF(obj.SensorArray);
        cond = act_chnl ~= exp_chnl;
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractDOA:UnmatchedChannel', exp_chnl, act_chnl);
        end
        cond = ~(isa(x,'float')) ;
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','float');
        end
    end
    
    function setupImpl(obj,x)
        processInputSizeChangeImpl(obj,x);
    end
    
    function processInputSizeChangeImpl(obj,x)
        obj.pNumSnapshots = size(x,1);
    end
    
    function flag = isInputComplexityLockedImpl(~,~) 
        flag = false;
    end
                
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
        if isLocked(obj)
            s.cCovEstimator = saveobj(obj.cCovEstimator);
            s.pNumSnapshots = obj.pNumSnapshots;
        end
    end

    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractNarrowbandArrayProcessing(obj,s);
        if isfield(s,'isLocked')
            if s.isLocked
                obj.cCovEstimator = phased.internal.SpatialCovEstimator.loadobj(s.cCovEstimator);
                s = rmfield(s,'cCovEstimator');
            end
        end
    end
end


end

