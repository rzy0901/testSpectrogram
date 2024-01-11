function sv = privSteeringVec(posx, posy, posz, freq,c,azang, elang)
%PRIVSTEERINGVEC
%   Computes the steering vector. For use with gpuArray/arrayfun

%   Copyright 2012 The MathWorks, Inc.
  
  % angles defined in local coordinate system
  % angles expected in radians.



cdel = -cos(elang);
sdel = sin(elang);
cdaz = cos(azang);
sdaz = sin(azang);


incidentdirx = cdel*cdaz;
incidentdiry = cdel*sdaz;
incidentdirz = -sdel;

tau = posx*incidentdirx + ...
      posy*incidentdiry + ...
      posz*incidentdirz;

tau = tau/c;

sv = exp(-1i*2*pi*freq*tau);


