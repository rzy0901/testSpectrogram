classdef MatrixViewerSpecification < matlabshared.scopes.SystemObjectScopeSpecification
    %MatrixViewerSpecification   Define the MatrixViewerSpecification class.
    
    %   Copyright 2013 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $  $Date: 2013/05/28 03:41:50 $

    methods

        function this = MatrixViewerSpecification(varargin)
            %MatrixViewerSpecification   Construct the
            %MatrixViewerSpecification class.
            
            mlock;
            
            this@matlabshared.scopes.SystemObjectScopeSpecification(varargin{:});
        end

        function cfgFile = getConfigurationFile(~)
            cfgFile = 'phasedmatrixviewer.cfg';
        end
        
        function args = getHelpArgs(~, ~)
            args = '';
        end
        
        function appName = getAppName(~)
            appName = 'Matrix Viewer';
        end
        
        function hgRoot = getHGRoot(~)
            hgRoot = groot;
        end
        
        function b = isToolbarCompact(~, tag)
            b = strcmp(tag, 'autoscale');
        end
    end
    
    methods (Hidden)
        function b = useMCOSExtMgr(~)
            b = true;
        end
    end
end

% [EOF]
