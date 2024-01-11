function fv = circpol2pol(cfv)
%circpol2pol   Circular to linear polarization representation conversion
%   FV = circpol2pol(CFV) converts the circular polarization
%   representation, CFV, of a field to the corresponding linear
%   polarization representation, FV.
%
%   CFV is either a row vector or a 2-row matrix. If CFV is a matrix, each
%   column of CFV represents the field in the circular component form of
%   [El;Er], where El is the left circular component of the field and Er is
%   the right circular component of the field. If CFV is a vector, each
%   entry in CFV represents the circular polarization ratio defined as
%   Er/El.
%
%   FV has the same dimension as CFV. If CFV is a matrix, FV is also a
%   2-row matrix. Each column of FV contains the linear polarization
%   representation, [Eh;Ev], of the corresponding field in CFV. Eh is the
%   horizontal component of the field and Er is the vertical component of
%   the field. If CFV is a vector, then FV is also a vector. Each entry in
%   FV represents the polarization ratio, defined as Ev/Eh, of the
%   corresponding circular polarization ratio in FV.
%
%   % Examples:
%   
%   % Example 1:
%   %   Calculate the linear polarization representation of a horizontal
%   %   polarization field.
%
%   cfv = [1;1];
%   fv = circpol2pol(cfv)
%
%   % Example 2:
%   %   Calculate the polarization ratio corresponding to a right circular
%   %   polarization field. The field is described by the circular
%   %   polarization ratio.
%
%   cfv = [0;1];
%   q = cfv(2)/cfv(1);
%   P = circpol2pol(q)
%
%   See also phased, polellip, polratio, pol2circpol, stokes.

%   Copyright 2012 The MathWorks, Inc.

%   References:
%   [1] Harold Mott, Polarization in Antennas and Radar, John Wiley & Sons,
%   1986
%   [2] IEEE Std 149-1979, IEEE Standard Test Procedures for Antennas, 1979

%   Copyright 2012-2013 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(1,1,nargin);

[Nrow, Ncol] = size(cfv);

cond = (Nrow ~= 1) && (Nrow ~= 2);

if cond
    coder.internal.errorIf(cond,...
        'phased:polarization:invalidPolarizationInput','FV');
end

if Nrow == 1
    validateattributes(cfv,{'double'},{'nonnan','row'},'circpol2pol','q');
    
    % According to IEEE std 149-1979, circular polarization is Er/El,
    % reciprocal of q in Mott88, (2.85). Also see Mott92 pg 133.
    cfv = 1./cfv;
    idx1 = (cfv==-1);
    idx2 = isinf(cfv);
    
    fv = complex(inf(Nrow,Ncol));
    fv(idx2) = 1i;
    idx = ~(idx1|idx2);
    fv(idx) = -1i*(1-cfv(idx))./(1+cfv(idx));
else %Nrow == 2
    validateattributes(cfv,{'double'},{'finite','2d','nrows',2},...
        'circpol2pol','CFV');
    if any((cfv(1,:)==0)&(cfv(2,:)==0))
        coder.internal.errorIf(any((cfv(1,:)==0)&(cfv(2,:)==0)),...
            'phased:phased:zeroColumns','CFV');
    end
    M = 1/sqrt(2)*[1 1;1i -1i];
    fv = M*cfv;
end

% [EOF]
