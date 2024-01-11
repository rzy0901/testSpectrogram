classdef (Sealed,StrictDefaults) IntensityScope < matlabshared.scopes.UnifiedSystemScope
%IntensityScope Display scrolling intensity vectors or arrays
%   IS = phased.IntensityScope returns a System object, IS, that can show a
%   scrolling plot of intensity vectors versus time.
%
%   IS = phased.IntensityScope('Name', Value, ...) returns a IntensityScope
%   System object, IS, with each specified property name set to the
%   specified value. You can specify name-value pair arguments in any order
%   as (Name 1, Value 1, ..., Name N, Value N).
%
%   Step method syntax:
%
%   step(IS, X) displays the real signal, X, in the IntensityScope figure.
%   X is an N by M matrix where N is the number of intensity bins in an
%   intensity vector and M is the number of intensity vectors. N should be
%   greater than 1 and M is equal or greater than 1. Each column of the
%   matrix represents an intensity vector from successive times. The time
%   between intensity vectors is specified in the TimeResolution property.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   IntensityScope methods:
%
%   step      - Display signal in the IntensityScope figure (see above)
%   release   - Allow property value and input characteristics changes, and
%               release IntensityScope resources
%   clone     - Create IntensityScope object with same property values
%   isLocked  - Display locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>     - Clear IntensityScope figure
%   show      - Turn on visibility of IntensityScope figure
%   hide      - Turn off visibility of IntensityScope figure
%   isVisible - Return visibility of IntensityScope figure (logical)
%
%   IntensityScope properties:
%
%   Name                 - Caption to display on the IntensityScope window
%   XResolution          - X-axis sample spacing
%   XOffset              - X-axis display offset
%   XLabel               - X-axis label
%   Title                - Display title
%   TimeResolution       - Time resolution    
%   TimeSpan             - Time span
%   IntensityUnits       - Intensity units
%   Position             - Scope window position in pixels
%   ReducePlotRate       - Reduce plot rate to improve performance
%
%   % Example: 
%   %   Plot a range-time intensity for three simulated targets.
%    
%   % Function to simulate received range bins
%   txpow = 200; gain = 1e8; std = 5;
%   rangePow = @(bins,range) ...
%       gain.*exp(-0.5*((bins-range)./std).^2).* ...
%       txpow./(range.^4)./(sqrt(2*pi).*std);
%
%   % Setup the scope
%   scope = phased.IntensityScope( ...
%               'Name','Range-Time Intensity Scope',...
%               'XLabel','Range (m)', ...
%               'XResolution',1, ... %Assume range bin is 1 meter
%               'TimeResolution',0.05,'TimeSpan',5, ...
%               'IntensityUnits','Watts');
%   
%   % Create range bins for three targets, two moving in opposite 
%   % directions
%   x = 0:1000;
%   ranges(:,1) = 250:10:850;
%   ranges(:,2) = 850:-10:250;
%   ranges(:,3)  = 400;
%   for k = 1:size(ranges,1)
%       y = ((rangePow(x,ranges(k,1))+ ...
%             rangePow(x,ranges(k,2)) + ...
%             rangePow(x,ranges(k,3))));
%       scope(y.');
%       pause(.1);
%   end
%
%   See also phased.ScenarioViewer

%   Copyright 2015-2018 The MathWorks, Inc.

properties(Nontunable)
    %Name Caption to display on the IntensityScope window
    %   Specify the caption to display on the IntensityScope window as any
    %   string. The default value of this property is 'Intensity Scope'.
    Name = 'Intensity Scope';
    %XResolution X-axis sample spacing
    %   Specify the spacing between samples along the X-axis as a finite
    %   numeric scalar. The default value of this property is 1.
    XResolution = 1;
    %XOffset
    %   Specify the offset to apply to the x-axis. This property is a
    %   numeric scalar with default value 0.
    XOffset = 0;    
end
properties (Nontunable, Dependent)
    %XLabel X-axis label
    %   Specify the x-axis label as a string. The default value of this
    %   property is ''.
    XLabel
    %Title Display title
    %   Specify the display title as a string. The default value of this
    %   property is ''.
    Title
    %TimeResolution Time resolution
    %   Specify the time resolution of each intensity line as a positive
    %   scalar in seconds. The default is 1 ms.
    TimeResolution
    %TimeSpan Time span
    %   Specify the time span of the intensity display as a positive
    %   scalar in seconds. The default is 100 ms.
    TimeSpan
    %IntensityUnits Intensity units
    %   Specify the label in the colorbar as a string. This label is
    %   usually used to indicate units of the displayed intensity. The
    %   default is 'dB'.
    IntensityUnits
end
properties (Dependent)
    %ReducePlotRate Reduce plot rate to improve performance
    %   Set this property to false to update the viewer every time the step
    %   method is called. This will negatively impact performance. The
    %   default is true. This property is tunable.
    ReducePlotRate;
end

%properties (Access = private)
%    pNumOfSamples = 0;
%end

methods
    
    function obj = IntensityScope(varargin)
        
        obj@matlabshared.scopes.UnifiedSystemScope();
        obj.Position = uiscopes.getDefaultPosition([800 450]);
        setProperties(obj, nargin, varargin{:});
        
    end
    function set.Name(obj, value)
        setScopeName(obj, value);
        obj.Name = value;
    end
    function set.TimeResolution(obj, value)
        validateattributes(value,{'double'}, {'positive','real','scalar','finite','nonnan'},'','TimeResolution');
        setVisualProperty(obj, 'TimeResolution', value, true);
    end
    function value = get.TimeResolution(obj)
        value = getVisualProperty(obj, 'TimeResolution', true);
    end
    function set.TimeSpan(obj, value)
        validateattributes(value,{'double'}, {'positive','real','scalar','finite','nonnan'},'','TimeSpan');
        setVisualProperty(obj, 'HistoryDepth', value, true);
    end
    function value = get.TimeSpan(obj)
        value = getVisualProperty(obj, 'HistoryDepth', true);
    end
    function set.XLabel(obj, value)
        val = convertStringsToChars(value);
        if  ~ischar(val)
            error(message('phased:scopes:InvalidString','XLabel'));
        end
        setVisualProperty(obj, 'XLabel', val, false);
    end
    function value = get.XLabel(obj)
        value = getVisualProperty(obj, 'XLabel', false);
    end
    function set.Title(obj, value)
        val = convertStringsToChars(value);
        if  ~ischar(val)
            error(message('phased:scopes:InvalidString','Title'));
        end
        setVisualProperty(obj, 'Title', val, false);
    end
    function value = get.Title(obj)
        value = getVisualProperty(obj, 'Title', false);
    end
    function set.IntensityUnits(obj, value)
        val = convertStringsToChars(value);
        if  ~ischar(val)
            error(message('phased:scopes:InvalidString','IntensityUnits'));
        end
        setVisualProperty(obj, 'IntensityUnit', val, false);
    end
    function value = get.IntensityUnits(obj)
        value = getVisualProperty(obj, 'IntensityUnit', false);
    end
    function set.ReducePlotRate(obj, value)
        validateattributes(value,{'logical'}, {'scalar'},'','ReducePlotRate');
        setPropertyValue(obj.getSource, 'LockSynchronous', ~value);
    end
    function value = get.ReducePlotRate(obj)
        value = ~getPropertyValue(obj.getSource, 'LockSynchronous');
    end    
end
methods (Access = protected)
    function h = getScopeCfg(~)
        h = phased.scopes.IntensitySystemScopeSpecification('AppName','Intensity Scope');
    end
    function num = getNumInputsImpl(~)
        num = 1;
    end
    function flag = isInputSizeLockedImpl(~,~)
        %flag = false;
        flag = true;
    end
    function setupImpl(this, x)
        %this.pNumOfSamples = size(x,1);
        xData  = ['(0:' num2str(size(x,1)-1) ')*' ...
            num2str(this.XResolution) '+' num2str(this.XOffset)];
        setVisualProperty(this, 'XData',xData, false);
        setupImpl@matlabshared.scopes.UnifiedSystemScope(this, x);
        %don't buffer data see g1195827
        %set the property values in the config file instead of here
        %setPropertyValue(this.getFramework.DataSource,'LockSynchronous',false);        
    end
    function validateInputsImpl(~,x)
        validateattributes(x,{'double'},{'nonempty','2d','finite','real'},'step','X');
        sz_x = size(x);
        cond = sz_x(1) < 2;
        if cond
            coder.internal.errorIf(cond,'phased:scopes:InvalidInputDimension');
        end
        %         cond = this.pNumOfSamples ~= 0 && this.pNumOfSamples ~= sz_x(1);
        %         % Cannot vary the number of samples, only number of frames is
        %         % allowed to vary
        %         if cond
        %             coder.internal.errorIf(cond,'phased:scopes:InvalidInputDimension');
        %         end
    end
    
end
methods (Hidden)
    function st = getInputSampleTime(obj)
        % Return frame time
        dm = getMaxDimensions(obj.pSource,1);
        st = obj.TimeResolution*dm(2);
    end
end
methods (Static, Hidden)
    function flag = isAllowedInSystemBlock(~)
        flag = false;
    end
end
end
