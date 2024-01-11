classdef MatrixViewerValidator < matlabshared.scopes.Validator
    %MatrixViewerValidator   Define the MatrixViewerValidator class.

    %   Copyright 2013 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $  $Date: 2013/05/28 03:41:51 $

    properties (Constant)
        Scale  = @(val) ~(val <= 0 || ~isscalar(val) || isnan(val) || isinf(val) || ~isreal(val));
        Start  = @(val) ~(~isscalar(val) || isnan(val) || isinf(val) || ~isreal(val));
        Invert = @(val) ~(~islogical(val) || ~isscalar(val));
    end
end

% [EOF]
