classdef ArrayResponseDemo2DPolarScopeCfg < scopeextensions.AbstractSystemObjectScopeCfg
% define the phased scope config class

% Copyright 2010-2011 The MathWorks, Inc.

methods
   function this = ArrayResponseDemo2DPolarScopeCfg(varargin)

      this@scopeextensions.AbstractSystemObjectScopeCfg(varargin{:});
    end
end

methods
    function c = getHiddenTypes(~)
        c = {};
    end
    function appName = getAppName(~)
        %getAppName Returns the simple application name.
        appName = 'Array Response Pattern 2-D Polar Scope';
    end
     
    
    function cfgFile = getConfigurationFile(~)
      cfgFile = 'ArrayResponseDemo2DPolarScope.cfg';
    end          
    
    function helpArgs = getHelpArgs(this,key) %#ok
         helpArgs = [];
    end   

end

end 
  
