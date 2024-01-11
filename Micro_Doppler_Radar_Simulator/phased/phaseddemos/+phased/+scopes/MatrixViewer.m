classdef MatrixViewer < matlabshared.scopes.UnifiedSystemScope
    %MatrixViewer   Define the MatrixViewer class.
    
    %   Copyright 2013 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $  $Date: 2013/05/28 03:41:49 $
    
    properties
        Name = 'Matrix Viewer';
    end
    
    properties (Dependent)
        
        XLabel;
        YLabel;
        XStart;
        YStart;
        XScale;
        YScale;
        XInvert;
        YInvert;
        Title;
    end
    
    properties (Hidden, Dependent, SetAccess = private)
        Axes;
    end
    
    methods
        
        function this = MatrixViewer(varargin)
            %MatrixViewer   Construct the MatrixViewer class.
            this@matlabshared.scopes.UnifiedSystemScope();
            
            this.Position = uiscopes.getDefaultPosition([560 420]);

            setProperties(this, nargin, varargin{:});
        end
        
        function set.XLabel(this, newXLabel)
            if ischar(newXLabel)
                setVisualProperty(this, 'XLabel', newXLabel);
            else
                msgObj = MException('Spcuilib:scopes:InvalidLabelValue', ...
                    'XLabel must be a string.');
                throwAsCaller(msgObj);
            end
        end
        function xLabel = get.XLabel(this)
            xLabel = getVisualProperty(this, 'XLabel');
        end
        
        function set.YLabel(this, newYLabel)
            if ischar(newYLabel)
                setVisualProperty(this, 'YLabel', newYLabel);
            else
                msgObj = MException('Spcuilib:scopes:InvalidLabelValue', ...
                    'YLabel must be a string.');
                throwAsCaller(msgObj);
            end
        end
        function yLabel = get.YLabel(this)
            yLabel = getVisualProperty(this, 'YLabel');
        end
        
        function set.XStart(this, newXStart)
            if phased.scopes.MatrixViewerValidator.Start(newXStart)
                setVisualProperty(this, 'XStart', newXStart, true);
            else
                msgObj = MException('Spcuilib:scopes:InvalidStartValue', ...
                    'XStart must be a scalar real number.');
                throwAsCaller(msgObj);
            end
        end
        function XStart = get.XStart(this)
            XStart = getVisualProperty(this, 'XStart', true);
        end
        
        function set.YStart(this, newYStart)
            if phased.scopes.MatrixViewerValidator.Start(newYStart)
                setVisualProperty(this, 'YStart', newYStart, true);
            else
                msgObj = MException('Spcuilib:scopes:InvalidStartValue', ...
                    'YStart must be a scalar real number.');
                throwAsCaller(msgObj);
            end
        end
        function YStart = get.YStart(this)
            YStart = getVisualProperty(this, 'YStart', true);
        end
        
        function set.XScale(this, newXScale)
            if phased.scopes.MatrixViewerValidator.Scale(newXScale);
                setVisualProperty(this, 'XScale', newXScale, true);
            else
                msgObj = MException('Spcuilib:scopes:InvalidScaleValue', ...
                    'XScale must be a scalar real number greater than zero.');
                throwAsCaller(msgObj);
            end
        end
        function XScale = get.XScale(this)
            XScale = getVisualProperty(this, 'XScale', true);
        end
        
        function set.YScale(this, newYScale)
            if phased.scopes.MatrixViewerValidator.Scale(newYScale);
                setVisualProperty(this, 'YScale', newYScale, true);
            else
                msgObj = MException('Spcuilib:scopes:InvalidScaleValue', ...
                    'YScale must be a scalar real number greater than zero.');
                throwAsCaller(msgObj);
            end
        end
        function YScale = get.YScale(this)
            YScale = getVisualProperty(this, 'YScale', true);
        end
        
        function set.XInvert(this, newXInvert)
            if phased.scopes.MatrixViewerValidator.Scale(newXInvert);
                setVisualProperty(this, 'XInvert', newXInvert);
            else
                msgObj = MException('Spcuilib:scopes:InvalidInvertValue', ...
                    'XInvert must be a scalar logical.');
                throwAsCaller(msgObj);
            end
        end
        function xInvert = get.XInvert(this)
            xInvert = getVisualProperty(this, 'XInvert');
        end
        
        function set.YInvert(this, newYInvert)
            if phased.scopes.MatrixViewerValidator.Scale(newYInvert);
                setVisualProperty(this, 'YInvert', newYInvert);
            else
                msgObj = MException('Spcuilib:scopes:InvalidInvertValue', ...
                    'YInvert must be a scalar logical.');
                throwAsCaller(msgObj);
            end
        end
        function yInvert = get.YInvert(this)
            yInvert = getVisualProperty(this, 'YInvert');
        end
        
        function set.Title(this, newTitle)
            setVisualProperty(this, 'Title', newTitle);
        end
        function title = get.Title(this)
            title = getVisualProperty(this, 'Title');
        end
        
        function hAxes = get.Axes(this)
            hAxes = this.pFramework.Visual.Axes;
        end
    end
    
    methods (Access = protected)
        function cfg = getScopeCfg(~)
            cfg = phased.scopes.MatrixViewerSpecification;
        end
    end
end

% [EOF]
