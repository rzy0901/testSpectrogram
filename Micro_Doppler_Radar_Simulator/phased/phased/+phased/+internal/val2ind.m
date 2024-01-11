function idx_out = val2ind(value_in,delta,GridStartVal)
%This function is for internal use only. It may be removed in the future.

%val2ind   Convert the grid value to grid index
%   IDX = phased.internal.val2ind(VAL,GRIDDELTA) returns the index IDX
%   corresponding to the value VAL in a uniformly sampled grid.  The
%   starting value of the grid is 0 and the increment between the grid
%   samples is GRIDDELTA. If VAL is a vector, IDX is also a vector of the
%   same size containing the indices of the elements in VAL.
%
%   IDX = phased.internal.val2ind(VAL,GRIDDELTA,GRIDSTART) specifies the
%   grid starting value GRIDSTART.
%
%   If the distance between VAL and GRIDSTART is not an integer multiple of
%   GRIDDELTA, the index of next larger grid value is returned.
%
%   % Example:
%   %   A signal is sampled at 25 Hz. Find out which sample is received 
%   %   at 0.4 seconds.
%
%   Fs = 25;
%   ind = phased.internal.val2ind(0.4,1/Fs)
%
%   See also phased, linspace, unigrid.

%   Copyright 2014 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>
inputIsColumn = iscolumn(value_in);
value = value_in(:).';

idx = (value-GridStartVal)/delta;
for m = 1:numel(idx)
    if abs(idx(m)-round(idx(m))) <= 10*eps(idx(m))
        idx(m) = round(idx(m));
    end
end
% roundidx = find(abs(idx-round(idx)) <= 10*eps(idx));
% if ~isempty(roundidx)
%     idx(roundidx) = round(idx(roundidx));
% else
%     idx = idx; %#ok codegen requires this to constant fold
% end
idx =  ceil(idx) + 1;

if inputIsColumn
    idx_out = idx.';
else
    idx_out = idx;
end


% [EOF]
