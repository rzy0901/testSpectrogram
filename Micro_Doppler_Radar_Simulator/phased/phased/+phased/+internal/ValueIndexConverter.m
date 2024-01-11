classdef (Hidden, StrictDefaults) ValueIndexConverter < matlab.System & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%ValueIndexConverter   Convert value to index
%   H = phased.internal.ValueIndexConverter creates a System object, H,
%   which converts the value to its corresponding index on a uniform grid.
%
%   H = phased.internal.ValueIndexConverter(Name,Value) creates a value to
%   index converter System object, H, with the specified property Name set
%   to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   IDX = step(H,X) returns in IDX the indices correspond to the values in
%   X. If a value falls between two grids, the next grid is returned. X is
%   a vector and dX is a positive scalar. IDX has the same dimension as X.
%
%   IDX = step(H,X,dX) specifies the step size of the grid in a positive
%   scalar, dX. This syntax is applicable when you set the
%   GridStepSizeSource property to 'Input port'.
%
%   IDX = step(H,X,X0) specifies the starting value of the grid, X0, as a 
%   scalar. This syntax is applicable when you set the
%   GridInitialValueSource property to 'Input port'.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   IDX = step(H,X,dX,X0)
%   
%   % Example:
%   %   Find the indices of the time stamps, 1.25 and 1.5 seconds in a 
%   %   sampling grid with a sample rate of 10 Hz.
%
%   fs = 10; x = [1.25 1.5]; dt = 1/fs;
%   myConverter = phased.internal.ValueIndexConverter('GridStepSize',dt);
%   idx = step(myConverter,x)
%
%   See also phased, val2ind.

%   Copyright 2014 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %GridStepSizeSource   Source of grid step size
    %   Specify the source of grid's step size using one of 'Property' |
    %   'Input port', where the default is 'property'.
    GridStepSizeSource = 'Property'
    %GridInitialValueSource     Source of initial value of the grid
    %   Specify the source of grid's initial value using one of
    %   'Property' | 'Input port', where the default is 'Property'.
    GridInitialValueSource = 'Property'
end

properties (Nontunable)
    %GridStepSize   Grid step 
    %   Specify the grid size as a positive scalar. This property is only
    %   applicable when you set the GridStepSizeSource property to
    %   'Property'. The default is 1.
    GridStepSize = 1
    % GridInitialValue      Initial value of the grid
    %   The initial value of the grid as a real scalar. This property is
    %   only applicable when you set the GridInitialValueSource property to
    %   'Property'. The default is 0.
    GridInitialValue = 0
end

properties (Constant, Hidden)
    GridStepSizeSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    GridInitialValueSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
end

properties (Access = private, Nontunable)
    % pIsGridStepSizeInput  flag for source of step size of grid
    pIsGridStepSizeInput
    % pIsGridInitialValueInput  flag for source of initial value of grid
    pIsGridInitialValueInput
end

methods
    function set.GridStepSize(obj,val)
        validateattributes(val, {'numeric'}, {'finite','nonnan','positive','scalar'},...
            '','GridStepSize');
        obj.GridStepSize = val;
    end
    function set.GridInitialValue(obj,val)
        validateattributes(val, {'numeric'}, {'finite','nonnan','real','scalar'},...
            '','GridInitialValue');
        obj.GridInitialValue = val;
    end
end

methods
    function obj = ValueIndexConverter(varargin)
        setProperties(obj, nargin, varargin{:});
    end
end

methods (Access=protected)

    function setupImpl(obj,~,~,~)
        obj.pIsGridStepSizeInput = ...
            ~strcmp(obj.GridStepSizeSource,'Property');
        obj.pIsGridInitialValueInput = ...
            ~strcmp(obj.GridInitialValueSource,'Property');
    end
    
    function idx = stepImpl(obj, value_in,delta,GridStartVal)
        if obj.pIsGridStepSizeInput
            deltaval = delta;
        else
            deltaval = obj.GridStepSize;
        end
        if obj.pIsGridInitialValueInput
            if obj.pIsGridStepSizeInput
                startval = GridStartVal;
            else
                startval = delta;
            end
        else
            startval = obj.GridInitialValue;
        end
        cond = any(value_in < startval);
        if cond
            coder.internal.errorIf(cond,'phased:val2ind:OutOfRange', ...
                sprintf( '%5.2f', startval ));
        end
        idx = phased.internal.val2ind(value_in,deltaval,startval);
            
    end


    function num = getNumInputsImpl(obj)
        num = 1;
        if ~strcmp(obj.GridStepSizeSource,'Property')
            num = num+1;
        end
        if ~strcmp(obj.GridInitialValueSource,'Property')
            num = num+1;
        end
    end

    function validateInputsImpl(obj,value_in,delta,GridStartVal)
        validateattributes(value_in, {'numeric'}, {'real','vector'},'','X');
        if ~strcmp(obj.GridStepSizeSource,'Property')
            validateattributes(delta, {'numeric'}, {'finite','nonnan','scalar'},...
                '','dX');
        end
        if ~strcmp(obj.GridInitialValueSource,'Property')
            if ~strcmp(obj.GridStepSizeSource,'Property')
                validateattributes(GridStartVal, {'numeric'}, {'real','scalar'},...
                    '','X0');
            else
                validateattributes(delta, {'numeric'}, {'real','scalar'},...
                    '','X0');
            end
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        s.pIsGridStepSizeInput = obj.pIsGridStepSizeInput;
        s.pIsGridInitialValueInput = obj.pIsGridInitialValueInput;
    end

    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        if strcmp(prop,'GridStepSize') && ...
                strcmp(obj.GridStepSizeSource,'Input port')
            flag = true;
        end
        if strcmp(prop,'GridInitialValue') && ...
                strcmp(obj.GridInitialValueSource,'Input port')
            flag = true;
        end
    end
end

methods (Access = protected) %for Simulink
    
    function varargout = getInputNamesImpl(obj) 
        if strcmp(obj.GridStepSizeSource,'Property')
            if strcmp(obj.GridInitialValueSource,'Property')
                varargout = {'X'};
            else
                varargout = {'X','X0'};
            end
        else
            if strcmp(obj.GridInitialValueSource,'Property')
                varargout = {'X','dX'};
            else
                varargout = {'X','dX','X0'};
            end
        end

    end

    function varargout = getOutputNamesImpl(obj) %#ok<MANU>
        varargout = {'Idx'};
    end
    
    function varargout = getOutputSizeImpl(obj)
        varargout{1} = propagatedInputSize(obj,1);
    end
    function varargout = isOutputFixedSizeImpl(obj)
        varargout{1} = propagatedInputFixedSize(obj, 1);
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj, 1);
    end
    function varargout = isOutputComplexImpl(obj) %#ok<MANU>
        varargout{1} = false;
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Value -> Index');
    end        
end
methods (Static,Hidden,Access=protected)
  function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:ValueIndexConverterTitle')),...
          'Text',getString(message('phased:library:block:ValueIndexConverterDesc')));
  end
end
 
end