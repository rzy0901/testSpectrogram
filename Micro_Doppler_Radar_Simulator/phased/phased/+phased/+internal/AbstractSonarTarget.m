classdef (Hidden) AbstractSonarTarget < phased.internal.AbstractTarget
%This class is for internal use only. It may be removed in the future.

%ABSTRACTSONARTARGET Define the ABSTRACTSONARTARGET class
% This is an abstract class in support of SONAR target functionality.

%   Copyright 2015-2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties(Access = protected)
        % Private flag for whether pTS needs to be initialized
        pNeedTSInit;
    end

    methods (Access = protected)
        function obj = AbstractSonarTarget(varargin)
            obj = obj@phased.internal.AbstractTarget(varargin{:});
        end
    end

    methods (Access = protected)

        function resetImpl(obj)
            resetImpl@phased.internal.AbstractTarget(obj);
            if obj.pFluctuate
                obj.pNeedTSInit = true;
            else
                obj.pNeedTSInit = false;
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractTarget(obj);
            if isLocked(obj)
                s.pNeedTSInit = obj.pNeedTSInit;
            end
        end
    end
    
end
