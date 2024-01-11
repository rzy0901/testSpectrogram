classdef (Hidden) AbstractRespPattern < sigutils.sorteddisp
%This class is for internal use only. It may be removed in the future.

%AbstractRespPattern Class definition for phased.internal.AbstractRespPattern class
%   This is an abstract class to support basic functionality for response
%   patterns.

%   Copyright 2008-2014 The MathWorks, Inc.

properties (GetAccess = public, SetAccess = protected)
    %Type - Type of the response pattern
    %   Type is a read-only property and must be set at construction time.
    Type
end

properties (Abstract)
    %Pattern - Pattern data of the response pattern
    %   Pattern holds the data of the response pattern.  Depending on the
    %   type of the response pattern, the dimension of Pattern may be
    %   different from type to type.  Therefore, each subclass has to
    %   redefine its own Pattern property.
    Pattern
end

methods (Abstract)
    %PLOT - Plot the pattern
    %   This is an abstract method. Each subclass must override the
    %   implementation.
    plot(obj)
end

methods (Abstract, Access = protected)
    %getPlotOption - get plot option
    %   This is an abstract method. 
    plotoptionobj = getPlotOption(obj,varargin)
end

methods (Static, Hidden)
    function annotatePlotSetup
    %ANNOTATEPLOTSETUP Response pattern plot setting
        [fsize,lblColor,~] = phased.internal.AbstractRespPattern.setAnnotationSizeColor;
        set(gca,'fontsize',fsize)
        set(get(gca,'XLabel'),'Color',lblColor);
        set(get(gca,'YLabel'),'Color',lblColor);
        set(get(gca,'ZLabel'),'Color',lblColor);
        set(get(gca,'Title'),'Color',lblColor);

        % axis tight;
        grid on;
    end
    
    function [fontSize,lblColor,ticklblColor]=setAnnotationSizeColor
    %setAnnotationSizeColor Set the annotation font size and color
        fontSize = 10;
        ticklblColor = [0.4 0.4 0.4];
        lblColor = [0 0 0];
    end
    function dBRespLimited = limitDynamicdBRange(dBResp,dRange)
    % limit dynamic scales in dB responses
        if iscell(dBResp)
            respmax = max(cell2mat(cellfun(@max,dBResp,'UniformOutput',false)));  % maximum response value in dB
            respmin = respmax - dRange;  % minimum response value in dB
            dBRespLimited = dBResp;
            idx = cellfun(@(x)x<respmin, dBResp,'UniformOutput',false);
            for m = 1:numel(dBRespLimited)
                temp = dBRespLimited{m};
                if ~isempty(temp)
                    temp(idx{m}) = respmin;
                end
                dBRespLimited{m} = temp;
            end
        else
            respmax = max(max(dBResp));  % maximum response value in dB
            respmin = respmax - dRange;  % minimum response value in dB
            dBRespLimited = dBResp;
            dBRespLimited(dBResp<respmin) = respmin;
        end
    end
end
 
end
% [EOF]


