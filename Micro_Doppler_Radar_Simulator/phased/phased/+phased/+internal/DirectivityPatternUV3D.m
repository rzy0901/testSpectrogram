classdef (Hidden,Sealed) DirectivityPatternUV3D < phased.internal.AbstractRadiationPatternUV3D
%This class is for internal use only. It may be removed in the future.

%DirectivityPatternUV3D 3D response pattern in UV space
%   Hresp = phased.internal.DirectivityPatternUV3D('Name',Value,...)
%   returns a 3D magnitude response pattern object Hresp. Hresp can be used
%   to store the 3D spatial response of an array. See properties list below
%   for valid PropertyNames.
%
%   DirectivityPatternUV3D methods:
%       plot          - Plot the 3D response pattern.
%
%   DirectivityPatternUV3D properties:
%       Type      - '3D Response UV Pattern'.  This is read-only
%       U         - U grid
%       V         - V grid
%       Pattern   - Magnitude response pattern
%
%   Example:
%       % Construct a 3D isotropic array response pattern and then plot it.
%       hresp = phased.internal.DirectivityPatternUV3D;
%       plot(hresp)
%
%   See also phased.internal, phased.internal.RespPattern3D.

%   Copyright 2013 The MathWorks, Inc.

methods
    function obj = DirectivityPatternUV3D(varargin)

        obj@phased.internal.AbstractRadiationPatternUV3D(varargin{:});
       
        obj.Type = '3D Directivity Pattern in u-v space';

    end
    
end    

methods (Access = protected)
    function plotoptionobj = getPlotOption(obj,varargin)  %#ok<INUSL>
        plotoptionobj = phased.internal.DirectivityPatternUV3DPlotOption(varargin{:});
        plotoptionobj.PlotType = 'plot';
    end
           
    function validatepattern(obj,value)
        validateattributes(value,{'numeric'},{'nonempty',...
            'size',size(obj.U)},...
            sprintf('%s.Pattern',class(obj)),'Pattern');  
    end
    
end


end
% [EOF]