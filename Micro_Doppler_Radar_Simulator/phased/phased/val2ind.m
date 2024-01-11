function idx_out = val2ind(value_in,delta,GridStartVal)
%val2ind   Convert the grid value to grid index
%   IDX = val2ind(VAL,GRIDDELTA) returns the index IDX corresponding to the
%   value VAL in a uniformly sampled grid.  The starting value of the grid
%   is 0 and the increment between the grid samples is GRIDDELTA. If VAL is
%   a vector, IDX is also a vector of the same size containing the indices
%   of the elements in VAL.
%
%   IDX = val2ind(VAL,GRIDDELTA,GRIDSTART) specifies the grid starting
%   value GRIDSTART.
%
%   If the distance between VAL and GRIDSTART is not an integer multiple of
%   GRIDDELTA, the index of next larger grid value is returned.
%
%   % Example:
%   %   A signal is sampled at 25 Hz. Find out which sample is received 
%   %   at 0.4 seconds.
%
%   Fs = 25;
%   ind = val2ind(0.4,1/Fs)
%
%   See also phased, linspace, unigrid.

%   Copyright 2008-2014 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,3,nargin);

if nargin < 3
    GridStartVal = 0;
end

validateattributes(delta, {'numeric'}, {'finite','nonnan','scalar','positive'},...
    'val2ind','DELTA');
validateattributes(GridStartVal, {'numeric'}, {'finite','nonnan','real','scalar'},...
    'val2ind','GRIDSTART');
validateattributes(value_in, {'numeric'}, {'finite','nonnan','real','vector'},'val2ind','VAL');

cond = any(value_in(:).' < GridStartVal);
if cond
    coder.internal.errorIf(cond,'phased:val2ind:OutOfRange', ...
        sprintf( '%5.2f', GridStartVal ));
end

idx_out = phased.internal.val2ind(value_in,delta,GridStartVal);


% [EOF]
