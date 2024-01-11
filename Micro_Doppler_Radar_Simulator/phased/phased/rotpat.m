function rpat = rotpat(pat,az,el,rotax,expval_in)
%rotpat     Rotate pattern
%   RPAT = rotpat(PAT,AZ,EL,ROTAX) rotates the pattern PAT so its boresight
%   is aligned with the x axis of a new local coordinate system specified
%   in ROTAX.
%
%   AZ is a length-P vector specifying the azimuth angles (in degrees) at
%   which the pattern, PAT, is sampled. EL is a length-Q vector specifying
%   the elevation angles (in degrees) at which the pattern, PAT, is
%   sampled. PAT is a QxP matrix or a QxPxL array. When PAT is a QxP
%   matrix, it represents a 3D pattern sampled at specified azimuth and
%   elevation angles. When PAT is a QxPxL array, each page of PAT is a 3D
%   pattern sampled at specified azimuth and elevation angles. 
%
%   PAT is defined in an orthonormal coordinate system whose axes are noted
%   as x, y, and z. The azimuth angle is measured from x axis toward y axis
%   and the elevation angle is measured from x-y plane toward z axis. The
%   boresight of the PAT is considered to be aligned with x axis.
%
%   ROTAX represents the rotated coordinate system. ROTAX can be either a
%   3x3 matrix or a 3x3xL array. If ROTAX is a matrix, its columns
%   represent the orthonormal axes, noted as u, v, and w, of the rotated
%   coordinate system. Each axis is represented in the form of [x; y; z]
%   using PAT's original coordinate system. All specified patterns are then
%   rotated so their boresight directions are now along the u axis. If
%   ROTAX is a 3x3xL array, each page represents the rotated coordinate
%   system and patterns are rotated to be aligned with the corresponding u
%   axes.
%
%   RPAT = rotpat(...,EXPVAL) specifies a scalar as the extrapolated value
%   when AZ and EL do not cover the entire 3D space. In general, consider
%   setting EXPVAL to 0 if the pattern is specified in the linear scale or
%   -inf if the pattern is specified in dB scale. The default value of
%   EXPVAL is 0.
%
%   % Example:
%   %   Rotate the polarized pattern of a short dipole antenna.
%
%   % compute the 3D pattern of a short dipole antenna
%   ant = phased.ShortDipoleAntennaElement;
%   el = -90:90;
%   az = -180:180;
%   pat_h = zeros(numel(el),numel(az),'like',1+1i);
%   pat_v = pat_h;
%   fc = 3e8;
%   for m = 1:numel(el)
%       temp = step(ant,fc,[az;el(m)*ones(1,numel(az))]);
%       pat_h(m,:) = temp.H;
%       pat_v(m,:) = temp.V;
%   end
% 
%   % specify the rotated axes and rotate the pattern
%   newax = rotx(45)*roty(45);
%   pat2_h = rotpat(pat_h,az,el,newax);
%   pat2_v = rotpat(pat_v,az,el,newax);
% 
%   ant2 = phased.CustomAntennaElement(...
%       'SpecifyPolarizationPattern',true,...
%       'HorizontalMagnitudePattern',mag2db(abs(pat2_h)),...
%       'HorizontalPhasePattern',rad2deg(angle(pat2_h)),...
%       'VerticalMagnitudePattern',mag2db(abs(pat2_v)),...
%       'VerticalPhasePattern',rad2deg(angle(pat2_v)));
%   pattern(ant2,fc,'Type','Power');
%
%   See also phased, rotx, roty, rotz.

%   Copyright 2018 The MathWorks, Inc.

%#codegen

narginchk(4,5);

if nargin < 5
    expval = 0;
else
    validateattributes(expval_in,{'double','single'},{'scalar','nonnan','nonempty'},'rotpat','EXPVAL');
    expval = expval_in;
end

sigdatatypes.validateAngle(az,'rotpat','AZ',{'double','single'},{'vector','increasing','>=',-180,'<=',180});
sigdatatypes.validateAngle(el,'rotpat','EL',{'double','single'},{'vector','increasing','>=',-90,'<=',90});

p = numel(az);
q = numel(el);
validateattributes(pat,{'double','single'},{'3d','nrows',q,'ncols',p,'nonempty','nonnan'},'rotpat','PAT');

validateattributes(rotax,{'double','single'},{'3d','nrows',3,'ncols',3,'nonempty','nonnan',...
    'real','finite'});

Lp = size(pat,3);
Lax = size(rotax,3);
cond = Lp>1 && Lax>1 && Lax~=Lp;
if cond
    coder.internal.errorIf(cond,'phased:phased:invalidPageNumbers','ROTAX',Lp);
end   

for m = 1:Lax
    tempax = rotax(:,:,m);
    cond = norm(tempax'*tempax-eye(3))>sqrt(eps);
    if cond
        coder.internal.errorIf(cond,'phased:phased:expectedOrthonormalAxes','ROTAX');
    end
end

if Lp == 1 && Lax > 1
    rpat = zeros([size(pat) Lax],'like',pat);
else
    rpat = zeros(size(pat),'like',pat);
end

for m = 1:numel(el)
    gcoord = [az(:).';el(m)*ones(1,numel(az));ones(1,numel(az))];
    if Lax == 1
        lccord = phased.internal.global2localcoord(gcoord,'ss',zeros(3,1),rotax);
        for n = 1:size(pat,3)
            rpat(m,:,n) = interp2(az,el,pat(:,:,n),lccord(1,:),lccord(2,:),'linear',expval);
        end
    else 
        for n = 1:Lax
            lccord = phased.internal.global2localcoord(gcoord,'ss',zeros(3,1),rotax(:,:,axidx));
            if L == 1
                temppat = pat;
            else
                temppat = pat(:,:,n);
            end
            rpat(m,:,n) = interp2(az,el,temppat,lccord(1,:),lccord(2,:),'linear',expval);
        end
    end
end



