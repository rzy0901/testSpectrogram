classdef ArrayResponseDemo2DPolarScope < matlabshared.scopes.UnifiedSystemScope
% ArrayResponseDemo2DPolarScope  2D polar array response scope

% Copyright 2010-2011 The MathWorks, Inc.
    
    properties
        %   Name Caption to display on scope window
        %   Specify the caption to display on the scope window as any string.
        %   The default value of this property is 'Phased Response Pattern
        %   Scope'. This property is tunable.
        Name = 'Array Response Pattern 2-D Polar Scope';
    end
    
    
    methods
        function this=ArrayResponseDemo2DPolarScope(varargin)
            this@matlabshared.scopes.UnifiedSystemScope(varargin{:});
        end
        function set.Name(this, value)
            setScopeName(this, value);
            this.Name = value;
        end
    end
    
    methods
        function desc = getDescription(~)
            desc = 'Array Response Pattern 2-D Polar Scope';
        end
    end
    methods (Access = protected)
        function hScopeCfg = getScopeCfg (~)
            hScopeCfg = ArrayResponseDemo2DPolarScopeCfg;
        end
    end
end


