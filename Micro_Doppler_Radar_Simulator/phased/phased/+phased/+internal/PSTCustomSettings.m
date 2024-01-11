classdef (Sealed) PSTCustomSettings < handle
    
    %   Copyright 2019 The MathWorks, Inc.

    properties
        % This property indicates whether MATLAB is running on MOTW. This flag
        % is set to true at startup of a MOTW session.
        MOTW = false;
    end

    methods (Access = private)
        function obj = PSTCustomSettings
            % Private constructor, this class can only be instantiated by calling
            % one of its methods.
        end
    end

    methods (Static)
        function singleObj = getInstance
            mlock
            persistent obj
            if isempty(obj) || ~isvalid(obj)
                obj = phased.internal.PSTCustomSettings;
            end
            singleObj = obj;
        end

        function setMOTWFlag(flag)
            obj = phased.internal.PSTCustomSettings.getInstance;
            obj.MOTW = flag;
        end

        function flag = isMOTW
            obj = phased.internal.PSTCustomSettings.getInstance;
            flag = obj.MOTW;
        end
    end
end