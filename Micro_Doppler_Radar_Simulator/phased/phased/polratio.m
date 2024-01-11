function P = polratio(fv)
%polratio Polarization ratio
%   P = polratio(FV) returns the polarization ratio, P, of the field FV. FV
%   must be a 2-row matrix. Each column of FV represents the field in the
%   form of [Eh;Ev], where Eh is the horizontal component of the field and
%   Ev is the vertical component of the field. P is a row vector whose
%   number of columns is the same as the number of columns in FV. Each
%   entry in P contains the polarization ratio of the corresponding field
%   in FV.
%
%   The polarization ratio is defined as Ev/Eh.
%
%   % Example:
%   %   Calculate the polarization ratio of a 45-degree linear 
%   %   polarization.
%
%   fv = [1;1];
%   p = polratio(fv)
%
%   See also phased, polellip, pol2circpol, circpol2pol, stokes.

%   Copyright 2012 The MathWorks, Inc.

%   References:
%   [1] Harold Mott, Polarization in Antennas and Radar, John Wiley & Sons,
%   1986

%#codegen
%#ok<*EMCA>

validateattributes(fv,{'double'},{'finite','2d','nrows',2},...
    'polratio','FV');

Eh = fv(1,:);
Ev = fv(2,:);

EhZeroIdx = (Eh==0);
if any(EhZeroIdx&(Ev==0))
    coder.internal.errorIf(any(EhZeroIdx&(Ev==0)),'phased:phased:zeroColumns','FV');
end

if isreal(fv)
    P = inf(1,size(fv,2));
else
    P = complex(inf(1,size(fv,2)));
end

idx = ~EhZeroIdx;
P(idx) = Ev(idx)./Eh(idx);

% [EOF]
