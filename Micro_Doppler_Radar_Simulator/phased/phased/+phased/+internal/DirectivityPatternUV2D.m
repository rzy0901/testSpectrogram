classdef (Hidden, Sealed) DirectivityPatternUV2D < phased.internal.AbstractRadiationPatternUV2D
%This class is for internal use only. It may be removed in the future.

%DirectivityPatternUV2D 2D response pattern in U/V space
%   Hresp = phased.internal.DirectivityPatternUV2D('Name',Value,...)
%   returns a 2D magnitude response pattern object Hresp. Hresp can be used
%   to store the 2D spatial response of an array. See properties list below
%   for valid PropertyNames.
%
%   DirectivityPatternUV2D methods:
%       plot          - Plot the 2D response pattern.
%                       pattern
%
%   DirectivityPatternUV2D properties:
%       Type      - '2D Response Pattern in U/V space' (Read only)
%       U         - Response sampling values in U space
%       Pattern   - Magnitude response pattern
%       Freq      - Frequency corresponding to each pattern
%
%   Example:
%       % Construct and plot a 2D isotropic array response pattern.
%       hresp = phased.internal.DirectivityPatternUV2D;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern3D.

%   Copyright 2013 The MathWorks, Inc.
    
    
    methods
        function obj = DirectivityPatternUV2D(varargin)
        %DirectivityPatternUV2D Constructor of phased.internal.DirectivityPatternUV2D class
            obj@phased.internal.AbstractRadiationPatternUV2D(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        function validatePattern(obj,value)
            validateattributes(value,{'numeric'},...
                {'size',[numel(obj.U),numel(obj.Freq)]},...
                sprintf('%s.Pattern',class(obj)),'Pattern');   
        end
        
    end
    
    methods (Access = protected)
        function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
            plotoptionobj = phased.internal.DirectivityPatternUV2DPlotOption(varargin{:});
        end
        
    end


end



% [EOF]
