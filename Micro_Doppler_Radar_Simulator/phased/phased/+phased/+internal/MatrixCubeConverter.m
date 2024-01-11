classdef (Hidden,StrictDefaults) MatrixCubeConverter < matlab.System & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%MatrixCubeConverter   Convert matrix to cube
%   H = phased.internal.MatrixCubeConverter creates a System object, H,
%   which converts a data matrix to a data cube.
%
%   H = phased.internal.MatrixCubeConverter(Name,Value) creates a value to
%   index converter System object, H, with the specified property Name set
%   to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   XC = step(H,X) returns in XC the data cube obtained by reshaping the
%   data matrix, X. The size of X is PxQ and the size of XC is MxQxN, where
%   N is specified in the NumPages property of H. P must be the same as
%   M*N.
%
%   % Example:
%   %   Reshape a 6x2 matrix to a 3x2x2 cube.
%
%   x = [ones(3,2);2*ones(3,2)];
%   myConverter = phased.internal.MatrixCubeConverter('NumPages',2);
%   xc = step(myConverter,x)
%
%   See also phased, val2ind.

%   Copyright 2014 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %NumPages Number of pages
    %   Specify the number of pages in the cube as a positive integer.
    NumPages = 1
end

methods
    function set.NumPages(obj,val)
        validateattributes(val, {'numeric'}, {'finite','integer','scalar'},...
            '','NumPages');
        obj.NumPages = val;
    end
end

methods
    function obj = MatrixCubeConverter(varargin)
        setProperties(obj, nargin, varargin{:});
    end
end

methods (Access=protected)

    function xc = stepImpl(obj,x) 
        sz_x = size(x);
        N = obj.NumPages;
        xc = permute(reshape(x,[sz_x(1)/N N sz_x(2)]),[1 3 2]);
    end

    function validateInputsImpl(obj,x) 
        validateattributes(x, {'numeric'}, {'2d'},'','X');
        N = obj.NumPages;
        cond = (rem(size(x,1),N)~=0);
        if cond
            coder.internal.errorIf(cond,'phased:phased:numRowsNotIntegerMultiple','X',N);
        end
    end
    
end

methods (Access = protected) %for Simulink
    
    function varargout = getInputNamesImpl(obj)  %#ok<MANU>
        varargout = {'X'};

    end

    function varargout = getOutputNamesImpl(obj) %#ok<MANU>
        varargout = {'XC'};
    end
    
    function varargout = getOutputSizeImpl(obj)
        sz_x = propagatedInputSize(obj,1);
        N = obj.NumPages;
        varargout{1} = [sz_x(1)/N sz_x(2) N];
    end
    function varargout = isOutputFixedSizeImpl(obj)
        varargout{1} = propagatedInputFixedSize(obj, 1);
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj, 1);
    end
    function varargout = isOutputComplexImpl(obj)
        varargout{1} = propagatedInputComplexity(obj,1);
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Matrix -> Cube');
    end        
end

methods (Static,Hidden,Access=protected)
  function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:MatrixCubeConverterTitle')),...
          'Text',getString(message('phased:library:block:MatrixCubeConverterDesc')));
  end
end
 
end