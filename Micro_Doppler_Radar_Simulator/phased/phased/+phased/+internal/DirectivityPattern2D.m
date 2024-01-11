classdef (Hidden, Sealed) DirectivityPattern2D < phased.internal.AbstractRadiationPattern2D
%This class is for internal use only. It may be removed in the future.

%DirectivityPattern2D 2D response pattern
%   Hresp = phased.internal.DirectivityPattern2D('PropertyName',PropertyValue,...)
%   returns a 2D magnitude response pattern object Hresp. Hresp can be used
%   to store the 2D spatial response of an array. See properties list below
%   for valid PropertyNames.
%
%   DirectivityPattern2D methods:
%       plot          - Plot the 2D response pattern.
%       polar         - Plot the 2D response pattern using polar coordinates.
%       beamwidth     - Calculate beamwidth for the 2D response pattern
%       sidelobeLevel - Calculate sidelobe level for the 2D response
%                       pattern
%
%   DirectivityPattern2D properties:
%       Type      - '2D Response Pattern' (Read only)
%       SliceDir  - Response slicing direction
%       Angle     - Response sampling angles (degrees)
%       Pattern   - Magnitude response pattern
%       Freq      - Frequency corresponding to each pattern
%
%   Example:
%       % Construct and plot a 2D isotropic array response pattern.
%       hresp = phased.internal.DirectivityPattern2D;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern3D.

%   Copyright 2013-2014 The MathWorks, Inc.
    
    methods
        function obj = DirectivityPattern2D(varargin)
            obj@phased.internal.AbstractRadiationPattern2D(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
            plotoptionobj = phased.internal.DirectivityPattern2DPlotOption(varargin{:});
        end
        
        function validatePattern(obj,value)
            validateattributes(value,{'numeric'},...
                {'nrows',numel(obj.Angle)},...
                sprintf('%s.Pattern',class(obj)),'Pattern');   
        end
        
    end
        
end

% [EOF]
