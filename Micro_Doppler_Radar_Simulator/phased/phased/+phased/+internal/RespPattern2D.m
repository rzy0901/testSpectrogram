classdef (Hidden, Sealed) RespPattern2D < phased.internal.AbstractRadiationPattern2D
%This class is for internal use only. It may be removed in the future.

%RespPattern2D 2D response pattern
%   Hresp = phased.internal.RespPattern2D('PropertyName',PropertyValue,...)
%   returns a 2D magnitude response pattern object Hresp. Hresp can be used
%   to store the 2D spatial response of an array. See properties list below
%   for valid PropertyNames.
%
%   RespPattern2D methods:
%       plot          - Plot the 2D response pattern.
%       polar         - Plot the 2D response pattern using polar coordinates.
%       beamwidth     - Calculate beamwidth for the 2D response pattern
%       sidelobeLevel - Calculate sidelobe level for the 2D response
%                       pattern
%
%   RespPattern2D properties:
%       Type      - '2D Response Pattern' (Read only)
%       SliceDir  - Response slicing direction
%       Angle     - Response sampling angles (degrees)
%       Pattern   - Magnitude response pattern
%       Freq      - Frequency corresponding to each pattern
%
%   Example:
%       % Construct and plot a 2D isotropic array response pattern.
%       hresp = phased.internal.RespPattern2D;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern3D.

%   Copyright 2008-2013 The MathWorks, Inc.
    
    methods
        function obj = RespPattern2D(varargin)
            obj@phased.internal.AbstractRadiationPattern2D(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
            plotoptionobj = phased.internal.RespPattern2DPlotOption(varargin{:});
        end
        
    end
        

end

% [EOF]
