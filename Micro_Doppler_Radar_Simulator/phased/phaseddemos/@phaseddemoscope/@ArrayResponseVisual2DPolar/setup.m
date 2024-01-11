function setup(this, hVisualParent)
%SETUP  Create basic ui elements for the visualization

% Copyright 2010-2011 The MathWorks, Inc.

setupAxes(this, hVisualParent);

%% for 2D polar
if this.polarplot
    setpolaraxes(this);
    
    % Set the YTicks.
    set(this.Axes,'YLim',[-1 1]);
    set(this.Axes,'XLim',[-1 1]);
    set(this.Axes,'Visible','off');

    ylbl = 'Normalized Response Pattern ';
    hylbl = get(this.Axes,'YLabel');
    set(hylbl,'String',ylbl);
    set(hylbl,'Visible','on');
else
    % Set the X and Y axes labels.
    ylbl = 'Normalized Response Pattern';
    set(get(this.Axes,'YLabel'),'String',ylbl);
    xlabel(this.Axes, 'Angle (degrees)');
    
    % Set the YTicks.
    set(this.Axes,'YLim',[-180 180]);
    set(this.Axes,'XLim',[0 1]);
end
end
% [EOF]
