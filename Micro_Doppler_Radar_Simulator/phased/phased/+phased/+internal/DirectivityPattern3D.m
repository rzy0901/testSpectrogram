classdef (Hidden,Sealed) DirectivityPattern3D < phased.internal.AbstractRadiationPattern3D
%This class is for internal use only. It may be removed in the future.

%DirectivityPattern3D 3D response pattern
%   Hresp =
%   phased.internal.DirectivityPattern3D('PropertyName',PropertyValue,...)
%   returns a 3D magnitude response pattern object Hresp. Hresp can be used
%   to store the 3D spatial response of an array. See properties list below
%   for valid PropertyNames.
%
%   DirectivityPattern3D methods:
%       plot          - Plot the 3D response pattern.
%       polar         - Plot the 3D response pattern using polar 
%                       coordinates.
%       beamwidth     - Calculate beamwidth for the 3D response pattern
%       sidelobeLevel - Calculate sidelobe level for the 3D response
%                       pattern
%
%   DirectivityPattern3D properties:
%       Type      - '3D Response Pattern'.  This is read-only
%       ElAngle   - Response elevation sampling angles (degrees)
%       AzAngle   - Response azimuth sampling angles (degrees)
%       Pattern   - Magnitude response pattern
%
%   Example:
%       % Construct a 3D isotropic array response pattern and then plot it.
%       hresp = phased.internal.DirectivityPattern3D;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern2D.

%   Copyright 2013 The MathWorks, Inc.

methods
    function obj = DirectivityPattern3D(varargin)
    %DirectivityPattern3D Constructor of phased.internal.DirectivityPattern3D class
        obj@phased.internal.AbstractRadiationPattern3D(varargin{:});
    end
    

end    

methods (Access = protected)
    
    function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
        plotoptionobj = phased.internal.DirectivityPattern3DPlotOption(varargin{:});
        plotoptionobj.PlotType = 'plot';
    end
        
    function validatepattern(obj,value)
        validateattributes(value,{'numeric'},{'nonempty',...
            'size',[numel(obj.ElAngle) numel(obj.AzAngle)]},...
            sprintf('%s.Pattern',class(obj)),'Pattern'); 
    end
end


end
% [EOF]
