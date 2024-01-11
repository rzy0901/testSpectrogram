classdef (Hidden) AbstractULADOA < phased.internal.AbstractDOA
%This class is for internal use only. It may be removed in the future.

%AbstractULADOA Abstract class for DOA estimation using ULA

%   Copyright 2009-2013 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Access = protected, Nontunable)
    pSpatialSmoothing = 0;
    pEffChannels;   
end
   
methods (Access = protected)

    function obj = AbstractULADOA(varargin)
       obj = obj@phased.internal.AbstractDOA(varargin{:});
    end
    
    function privValidateSensorArray(~,val)
        %privValidateSensorArray
        %   Certain array operation is limited to certain array geometries.
        %   Each operation can then overload this method to do its own
        %   validation. By default, any array geometry is ok.
        %The array must be a ULA
        cond = ~(isa(val,'phased.ULA') || isa(val,'phased.HeterogeneousULA'));
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractULADOA:InvalidArray');
        end
        
    end
    % in phase mode excitation the effective 
    % number of elements would be the underlying number
    % of elements of the virtual ULA
    function N = getNumEffectiveElements(obj)
        N = getNumElements(obj.SensorArray);
    end
    function validatePropertiesImpl(obj)
        M = getNumEffectiveElements(obj);
        effChannels = M - obj.pSpatialSmoothing;
        cond = effChannels < 2;
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractULADOA:InvalidSpatialSmoothing', M - 2);
        end
    end
    
    function effChan = getEffectiveChannel(obj)
        effChan = getNumEffectiveElements(obj) - obj.pSpatialSmoothing;
    end
    
    function setupImpl(obj,x)
        setupImpl@phased.internal.AbstractDOA(obj,x);
        obj.cCovEstimator = phased.internal.SpatialCovEstimator(...
            'NumSubarrays',obj.pSpatialSmoothing+1);
        obj.pEffChannels = getEffectiveChannel(obj);
    end
    
    function resetImpl(obj)
        resetImpl@phased.internal.AbstractDOA(obj);
        reset(obj.cCovEstimator);
    end
    
    function releaseImpl(obj)
        releaseImpl@phased.internal.AbstractDOA(obj);
        release(obj.cCovEstimator);
    end

    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractDOA(obj);
        s.pSpatialSmoothing = obj.pSpatialSmoothing;
        if isLocked(obj)
            s.pEffChannels = obj.pEffChannels;
        end
    end

    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractDOA(obj,s);
        if isfield(s,'isLocked')
            s = rmfield(s,'isLocked');
        end
    end
    
end
methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl(sensorType)
        if nargin == 1
            groups = getPropertyGroupsImpl@phased.internal.AbstractDOA(sensorType);
        else
            groups = getPropertyGroupsImpl@phased.internal.AbstractDOA('ula');
        end
    end
end
end

