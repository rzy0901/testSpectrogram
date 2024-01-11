function currentclutter = privDopplerShift(time, timeoffset,  sampRate, patchDoppler, clutterpatch)
 %PRIVDOPPLERSHIFT
%   Computes the Doppler Shift. For use with gpuArray/arrayfun

%   Copyright 2012 The MathWorks, Inc.

    t  = time/sampRate + timeoffset;
    currentclutter = clutterpatch * exp(1i*2*pi*patchDoppler*t);
    
    

end