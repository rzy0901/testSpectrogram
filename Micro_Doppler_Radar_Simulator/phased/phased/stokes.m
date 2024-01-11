function varargout = stokes(fv_in)
%stokes   Stokes parameters
%   S = stokes(FV) returns the Stokes parameters of the polarization field
%   FV. FV can be either a row vector or a 2-row matrix containing the
%   linear polarization representation of the field.
%
%   If FV is a matrix, each column in FV represents the field in the form
%   of [Eh;Ev], where Eh is the horizontal component of the field and Ev is
%   the vertical component of the field. If FV is a vector, each entry in
%   FV represents the polarization ratio, Ev/Eh, of a unit-power field. 
%
%   S is a 4-row matrix whose number of columns is the same as the number
%   of columns in FV. Each column of S contains Stokes parameters of the
%   corresponding field in FV.
%
%   stokes(FV) normalizes the fields in FV to unit power and then plots
%   them on the Poincare sphere. Each field generates to a point on the
%   sphere.
%
%   % Examples:
%   
%   % Example 1:
%   %   Calculate Stokes parameter for a left circularly polarized field 
%   %   and a horizontally polarized field.
%
%   fv = [1 1i;1 0].';
%   s = stokes(fv)
%
%   % Example 2:
%   %   Plot a left circularly polarized field and a horizontally polarized
%   %   field on the Poincare sphere.
%
%   fv = [1 1i;1 0].';
%   stokes(fv)
%
%   See also phased, polellip, polratio, pol2circpol, circpol2pol.

%   Copyright 2012-2013 The MathWorks, Inc.

%   References:
%   [1] Harold Mott, Polarization in Antennas and Radar, John Wiley & Sons,
%   1986
%   [2] Warren L. Stutzman, Polarization in Electromagnetic Systems, Artech
%   House, 1993

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(1,1,nargin);

[Nrow, Ncol] = size(fv_in);

cond = (Nrow ~= 1) && (Nrow ~= 2);

if cond
    coder.internal.errorIf(cond, ...
          'phased:polarization:invalidPolarizationInput','FV');
end

if Nrow == 1
    validateattributes(fv_in,{'double'},{'nonnan','row'},'polellip','P');
    fv = [ones(1,Ncol);fv_in];
    fv(1,isinf(fv_in)) = 0;
    fv(2,isinf(fv_in)) = 1;
    fvmag = hypot(fv(1,:),fv(2,:));
    fv = bsxfun(@rdivide,fv,fvmag);
else %Nrow == 2
    validateattributes(fv_in,{'double'},{'finite','2d','nrows',2},...
        'polellip','FV');
    fv = fv_in;
    if any((fv(1,:)==0)&(fv(2,:)==0))
        coder.internal.errorIf(any((fv(1,:)==0)&(fv(2,:)==0)),...
            'phased:phased:zeroColumns','FV');
    end
end

normtemp = real(fv.*conj(fv));
crosstemp = 2*fv(2,:).*conj(fv(1,:));

G = [sum(normtemp);-diff(normtemp);real(crosstemp);imag(crosstemp)];

if nargout
    varargout{1} = G;
else
    if ~isempty(coder.target)
        coder.internal.assert(false,'phased:phased:invalidCodegenOutput','stokes');
    end
    [px,py,pz] = sphere(20);
    hsurf = surf(px,py,pz,'Tag','Sphere');
    set(hsurf,'CData',ones(size(get(hsurf,'ZData'))));
    set(hsurf,'LineStyle','--','LineWidth',0.5,...
        'EdgeColor',[0.4 0.4 0.4],'EdgeAlpha',0.5,'FaceAlpha',0.5);
    daspect([1 1 1]);
    %normalize to unit sphere
    Px = G(2,:)./G(1,:);
    Py = G(3,:)./G(1,:);
    Pz = G(4,:)./G(1,:);
    hold on;
    plot3(Px,Py,Pz,'.','MarkerSize',20,'MarkerFaceColor','r',...
        'Tag','StokesPoints');
    plot3([1 1.2],[0 0],[0 0],'Tag','HAxis');
    text(1.3,0,0,'H','Tag','HPolLabel');
    plot3([-1 -1.2],[0 0],[0 0],'Tag','VAxis');
    text(-1.3,0,0,'V','Tag','VPolLabel');
    plot3([0 0],[0 0],[-1.3 -1],'Tag','RCAxis');
    text(0,0,-1.3,'RC','Tag','RCPolLabel');
    plot3([0 0],[0 0],[1 1.3],'Tag','LCAxis');
    text(0,0,1.3,'LC','Tag','LCPolLabel');
    hold off;
    title(getString(message('phased:polarization:PoincareSphere')),...
        'Tag','PoincareTitle');
    axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
    axis off;
    view(150,30);
end


% [EOF]
