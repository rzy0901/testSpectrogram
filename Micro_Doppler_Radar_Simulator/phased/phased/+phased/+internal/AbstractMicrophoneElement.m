classdef (Hidden) AbstractMicrophoneElement < phased.internal.AbstractElement
%This class is for internal use only. It may be removed in the future.

%AbstractMicrophoneElement   Define the AbstractMicrophoneElement class.

%   Copyright 2010-2011 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods (Access = protected)

        function obj = AbstractMicrophoneElement(varargin)
            obj@phased.internal.AbstractElement(varargin{:});
            obj.pIsOutputComplex = false;
        end
    end

end


