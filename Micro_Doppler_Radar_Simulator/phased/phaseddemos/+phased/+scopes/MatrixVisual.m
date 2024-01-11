classdef MatrixVisual < matlabshared.scopes.visual.AxesVisual
    %MatrixVisual   Define the MatrixVisual class.
    
    %   Copyright 2013 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $  $Date: 2013/05/28 03:41:52 $
    
    properties (SetAccess = protected)
        Image;
    end
    
    properties (Access = private)
        XStart = 1;
        YStart = 1;
        XScale = 1;
        YScale = 1;
    end
    
    properties (Access = protected)
        DataSourceChangedListener;
    end
    
    methods
        function this = MatrixVisual(varargin)
            %MatrixVisual   Construct the MatrixVisual class.
            
            mlock;
            
            this@matlabshared.scopes.visual.AxesVisual(varargin{:});
            
            this.DataSourceChangedListener = addlistener(this.Application, ...
                'DataSourceChanged', @(h, ev) onDataSourceChanged(this));
            
            this.XStart = evalPropertyValue(this, 'XStart');
            this.YStart = evalPropertyValue(this, 'YStart');
            this.XScale = evalPropertyValue(this, 'XScale');
            this.YScale = evalPropertyValue(this, 'YScale');
        end
        
        function [primary, secondary] = getPlotNavigationAxes(~)
            primary = 'C';
            secondary = 'XY';
        end
        
        function extents = getDataExtents(this, ~)
            maxDims = this.Application.DataSource.getMaxDimensions;
            
            extents.X = [-0.5 maxDims(2)-0.5] * this.XScale + this.XStart;
            extents.Y = [-0.5 maxDims(1)-0.5] * this.YScale + this.YStart;
            
            c = get(this.Image, 'CData');
            extents.C = [min(c(:)) max(c(:))];
        end
        
        function update(this)
            d = getRawData(this.Application.DataSource, 1);
            set(this.Image, 'CData', double(d));
        end
        
        function varargout = validateSource(this, hSource)
            if nargin < 2
                hSource = this.Application.DataSource;
            end
            
            b = true;
            exception = MException.empty;
            
            if getNumInputs(hSource) ~= 1
                b = false;
                exception = MException('signal:MatrixViewer:InvalidNumInputs', ...
                    'The Matrix Viewer can only support a single input.');
            else
                maxDims = getMaxDimensions(hSource);
                if numel(maxDims) > 2 && prod(maxDims(3:end)) > 1
                    b = false;
                    exception = MException('signal:MatrixViewer:InvalidDimensions', ...
                        'The Matrix Viewer cannot support more than 2 dimensions.');
                end
            end
            
            if nargout
                varargout = {b, exception};
            elseif ~b
                throw(exception);
            end
        end
    end
    
    methods (Static)
        function h = getPropertySet()
            h = getPropertySet@matlabshared.scopes.visual.AxesVisual( ...
                'XStart',       'string', '1', ...
                'YStart',       'string', '1', ...
                'XScale',       'string', '1', ...
                'YScale',       'string', '1', ...
                'XInvert',      'bool',   false, ...
                'YInvert',      'bool',   true, ...
                'XLabel',       'string', 'X', ...
                'YLabel',       'string', 'Y', ...
                'Title',        'string', '');
        end
    end
    
    methods (Hidden)
        function setup(this, hVisParent)
            
            % Ignore HG2 warning about DrawMode.
            w = warning('off', 'MATLAB:hg:WillBeRemovedReplaceWith');
            c = onCleanup(@() warning(w));
            
            setup@matlabshared.scopes.visual.AxesVisual(this, hVisParent);
            
            fgColor = [175 175 175]/255;
            set(this.Axes, 'Box', 'on', ...
                'Layer', 'top', ...
                'XColor', fgColor, ...
                'YColor', fgColor, ...
                'Units', 'Pixel');
            
            this.Image = image('Parent', this.Axes, ...
                'CDataMapping', 'Scaled');
            
            xlabel(this.Axes, getPropertyValue(this, 'XLabel'), 'Color', fgColor);
            ylabel(this.Axes, getPropertyValue(this, 'YLabel'), 'Color', fgColor);
            title(this.Axes,  getPropertyValue(this, 'Title'),  'Color', fgColor);
            
            set(hVisParent, 'ResizeFcn', @(~,~) resizeFunction(this, hVisParent));
            
            resizeFunction(this, hVisParent);
            
            propertyChanged(this, 'XInvert');
            propertyChanged(this, 'YInvert');
        end
        
        function propertyChanged(this, eventData)
            
            if ~ischar(eventData)
                eventData = eventData.AffectedObject.Name;
            end
            
            switch eventData
                case 'XScale'
                    this.XScale = evalPropertyValue(this, 'XScale');
                    updateXProperties(this);
                case 'YScale'
                    this.YScale = evalPropertyValue(this, 'YScale');
                    updateYProperties(this);
                case 'XStart'
                    this.XStart = evalPropertyValue(this, 'XStart');
                    updateXProperties(this);
                case 'YStart'
                    this.YStart = evalPropertyValue(this, 'YStart');
                    updateYProperties(this);
                case 'XInvert'
                    this.Axes.XDir = convertInvertToDirection(this, 'XInvert');
                case 'YInvert'
                    this.Axes.YDir = convertInvertToDirection(this, 'YInvert');
                case 'XLabel'
                    xlabel(this.Axes, getPropertyValue(this, 'XLabel'));
                case 'YLabel'
                    ylabel(this.Axes, getPropertyValue(this, 'YLabel'));
                case 'Title'
                    title(this.Axes, getPropertyValue(this, 'Title'));
            end
        end
    end
    
    methods (Access = protected)
        
        function hInstall = createGUI(~)
            
            mVideoTools = uimgr.uimenugroup('VideoTools', -inf); %, mColormap);
            hInstall = uimgr.Installer({ ...
                mVideoTools, 'Base/Menus/Tools'});
        end
        
        function onDataSourceChanged(this)
            maxDims = getMaxDimensions(this.Application.DataSource);
            
            % Set a default blank image when changing the source.
            set(this.Image, 'CData', zeros(maxDims(1:2)));
            
            % Update the X & Y Properties based on the new size of the data
            % inputs.
            updateXProperties(this);
            updateYProperties(this);
        end
        
        function updateXProperties(this)
            
            maxDims = getMaxDimensions(this.Application.DataSource);
            
            if all(maxDims == 0)
                return;
            end
            xStart  = this.XStart;
            xScale  = this.XScale;
            
            % Calculate the xlimits and xdata based on the size of the data
            % from the source and the scaling and starting factors.
            this.Axes.XLim   = [-0.5 maxDims(2)-0.5] * xScale + xStart;
            this.Image.XData = [0    maxDims(2)-1  ] * xScale + xStart;
            
            % Reset the zoom state since we're changing the limits.
            resetplotview(this.Axes, 'SaveCurrentView');
        end
        
        function updateYProperties(this)
            
            maxDims = getMaxDimensions(this.Application.DataSource);
            if all(maxDims == 0)
                return;
            end
            yStart  = this.YStart;
            yScale  = this.YScale;
            
            % Calculate the ylimits and ydata based on the size of the data
            % from the source and the scaling and starting factors.
            this.Axes.YLim   = [-0.5 maxDims(1)-0.5] * yScale + yStart;
            this.Image.YData = [0    maxDims(1)-1  ] * yScale + yStart;
            
            % Reset the zoom state since we're changing the limits.
            resetplotview(this.Axes, 'SaveCurrentView');
        end
    end
end

function resizeFunction(this, hVisParent)

parentPosition = getpixelposition(hVisParent);

% Have the axes take up all of the available space.  Do not use normalized
% here because the insets need to be set in pixels.
this.Axes.OuterPosition = [0 0 parentPosition(3:4)];

% Keep the axes inset small to minimize wasted space.
this.Axes.LooseInset    = [1 1 15 15];

end

% -------------------------------------------------------------------------
function direction = convertInvertToDirection(this, propName)

if getPropertyValue(this, propName)
    direction = 'reverse';
else
    direction = 'normal';
end
end

% [EOF]
