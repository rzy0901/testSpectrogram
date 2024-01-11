function L = fspl(R,lambda)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2015 The MathWorks, Inc.

%fspl     Free space path loss
%   L = fspl(R,LAMBDA) returns the free space path loss L (in dB) suffered
%   by a signal with wavelength LAMBDA (in meters) when it is propagated in
%   free space for a distance of R (in meters). R can be a length-M vector
%   and LAMBDA can be a length-N vector. L has the same dimensionality as
%   MxN. Each element in L is the free space path loss for the
%   corresponding propagation distance specified in R.
%
%   Note that the best case is lossless so the loss is always greater than
%   or equal to 0 dB.
%
%   % Example:
%   %   Calculate the free space loss for a signal whose wavelength is 30
%   %   cm. The signal is propagated for 1 km.
%
%   L = phased.internal.fspl(1000,0.3)
%
%   See also phased, phased.FreeSpace.


%   Reference
%   [1] John Proakis, Digital Communications, 4th Ed., McGraw-Hill, 2001

%#codegen 
%#ok<*EMCA

L = 4*pi*R(:)*(1./lambda(:).');
% L(L<1) = 1;
L = validateLoss(L);
L = mag2db(L);

end

function L = validateLoss(L)
    for m = 1:numel(L)
        if L(m)<1
            L(m)=1;  % Loss cannot be less than 1
        end
    end
end



% [EOF]
