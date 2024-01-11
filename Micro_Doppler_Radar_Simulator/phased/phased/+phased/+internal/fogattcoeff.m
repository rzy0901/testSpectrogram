function Kl = fogattcoeff(f,T)

% f: frequency (Hz) row vector
% T: temperature (K) scalar
% Kl: attenuation coefficient (db/km)/(g/m^3)

fGHz = f/1e9;

theta = 300/T;
epsilon0 = 77.66+103.3*(theta-1);
epsilon1 = 0.0671*epsilon0;
epsilon2 = 3.52;

fpGHz = 20.20-146*(theta-1)+316*(theta-1).^2;
fsGHz = 39.8*fpGHz;

epsilon2p = fGHz.*(epsilon0-epsilon1)./(fpGHz.*(1+(fGHz./fpGHz).^2)) + ...
    fGHz.*(epsilon1-epsilon2)./(fsGHz.*(1+(fGHz./fsGHz).^2));
epsilon1p = (epsilon0-epsilon1)./(1+(fGHz./fpGHz).^2) + ...
    (epsilon1-epsilon2)./(1+(fGHz./fsGHz).^2) + epsilon2;

eta = (2+epsilon1p)./epsilon2p;
Kl = 0.819*fGHz./(epsilon2p.*(1+eta.^2));
