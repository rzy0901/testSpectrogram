function [azOut, elOut] = angleAxesConversion(azin, elin, axisAz, axisEl)
%This function is for internal use only. It may be removed in the future.

%ANGLEAXESCONVERSION converts angles to new rotated coordinate system
%  [AZOUT, ELOUT] = ANGLEAXESCONVERSION(AZIN, ELIN, AXISAZ, AXISEL) takes an
%  azimuth, AXISAZ, and elevation, AXISEL, rotation to an original
%  coordinate systems axes and applies that rotation to the angles pairs defined
%  by all combinations of azimuth angles AZIN and elevation ELIN which are
%  defined in the original coordinate system.
%  The new azimuth and elevation angles in the new (rotated) coordinate
%  system are output as the matrices AZOUT and ELOUT.
%
%  All angles are in radians.
%
%  See also phased, global2localcoord.

%  Copyright 2012 The MathWorks, Inc.
%
    azin = reshape( azin, 1, 1, []);
    elin = reshape( elin, [], 1, 1);
            

    cosAxisAz = cos(axisAz);
    sinAxisAz = sin(axisAz);
    cosAxisEl = cos(axisEl);
    sinAxisEl = sin(axisEl);
    
    [azOut, elOut] = arrayfun(@angBasisChange, azin, elin, cosAxisAz, sinAxisAz, cosAxisEl, sinAxisEl);
    
end

function [azo, elo] = angBasisChange( azin, elin, cA, sA, cE, sE)

%convert azin elin to rectangular
z = sin(elin);
coselev = cos(elin);
x = coselev * cos(azin);
y = coselev * sin(azin);

%Apply rotation
xnew = cA*cE*x + sA*cE*y + sE*z;
ynew = -sA*x + cA*y;
znew = -cA*sE*x -sA*sE*y + cE*z;


%Convert back to spherical
hypotxy = sqrt(xnew*xnew + ynew*ynew);

elo = atan2(znew,hypotxy);
azo = atan2(ynew,xnew);


end

