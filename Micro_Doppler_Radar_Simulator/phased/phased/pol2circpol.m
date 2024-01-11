function cfv = pol2circpol(fv)
%pol2circpol   Linear to circular polarization representation conversion
%   CFV = pol2circpol(FV) converts the linear polarization representation,
%   FV, of a field to the corresponding circular polarization
%   representation, CFV.
%
%   FV is either a row vector or a 2-row matrix. If FV is a matrix, each
%   column of FV represents the field in the form of [Eh;Ev], where Eh is
%   the horizontal component of the field and Ev is the vertical component
%   of the field. If FV is a vector, each entry in FV represents the
%   polarization ratio defined as Ev/Eh.
%
%   CFV has the same dimension as FV. If FV is a matrix, CFV is also a
%   2-row matrix. Each column of CFV contains the circular polarization
%   representation, [El;Er], of the corresponding field in FV. El is the
%   left circular component of the field and Er is the right circular
%   component of the field. If FV is a vector, then CFV is also a vector.
%   Each entry in CFV represents the circular polarization ratio, defined
%   as Er/El, of the corresponding polarization ratio in FV.
%
%   % Examples:
%   
%   % Example 1:
%   %   Calculate the circular polarization representation of a left 
%   %   circular polarization field.
%
%   fv = [1;1i];
%   cfv = pol2circpol(fv)
%
%   % Example 2:
%   %   Calculate the circular polarization ratio corresponding to a right 
%   %   circular polarization field. The field is described by the
%   %   polarization ratio.
%
%   fv = [1;-1i];
%   P = fv(2)/fv(1);
%   q = pol2circpol(P)
%
%   See also phased, polellip, polratio, circpol2pol, stokes.

%   Copyright 2012-2013 The MathWorks, Inc.

%   References:
%   [1] Harold Mott, Polarization in Antennas and Radar, John Wiley & Sons,
%   1986
%   [2] IEEE Std 149-1979, IEEE Standard Test Procedures for Antennas, 1979

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(1,1,nargin);

Nrow = size(fv,1);

cond = (Nrow ~= 1) && (Nrow ~= 2);

if cond
    coder.internal.errorIf(cond, ...
           'phased:polarization:invalidPolarizationInput','FV');
end

if Nrow == 1
    validateattributes(fv,{'double'},{'nonnan','row'},'pol2circpol','P');
    cfv = complex(-1*ones(size(fv)));
    idx = ~isinf(fv);
    cfv(idx) = (1-1i*fv(idx))./(1+1i*fv(idx));  
    % According to IEEE std 149-1979, circular polarization is Er/El,
    % reciprocal of q in Mott88, (2.85). Also see Mott92 pg 133.
    cfv = 1./cfv;
else %Nrow == 2
    validateattributes(fv,{'double'},{'finite','2d','nrows',2},...
        'pol2circpol','FV');
    if any((fv(1,:)==0)&(fv(2,:)==0))
        coder.internal.errorIf(any((fv(1,:)==0)&(fv(2,:)==0)),...
            'phased:phased:zeroColumns','FV');
    end
    M = 1/sqrt(2)*[1 -1i;1 1i];
    cfv = M*fv;
end



% [EOF]
