function patazel = azelcut2pat(azcut,elcut)
%azelcut2pat    Interpolate azimuth and elevation cuts to form 3D pattern
%   PAT = azelcut2pat(ACUT,ECUT) returns a 3D response pattern, PAT, from
%   its azimuth cut, ACUT (in dB), and elevation cut, ECUT (in dB).
%
%   ACUT specifies the azimuth cut of the 3D pattern as a 1xQ vector or a
%   LxQ matrix, where Q is the number of azimuth angles and L is the number
%   of frequencies. The default value of this property is a 1x361 row
%   vector with all elements equal to 0.
%
%   ECUT specifies the elevation cut of the 3D pattern as a 1xP vector or a
%   LxP matrix, where P is the number of elevation angles and L is the
%   number of the frequencies. The default value of this property is a
%   1x181 row vector with all elements equal to 0.
%
%   PAT is a a PxQxL array that contains the interpolated 3D response
%   pattern (in dB). The pattern interpolation is given as
%
%       PAT(az,el) = ACUT(az)+ECUT(el)
%
%   where az and el represent azimuth and elevation angles, respectively.
%
%   % Example:
%   %   Form a custom antenna from the azimuth and elevation cuts of its
%   %   pattern and assign it to a custom antenna element.
%
%   az = -180:180;
%   acut = mag2db(ones(size(az)));
%   el = -90:90;
%   ecut = mag2db(cosd(el));
%   pat = azelcut2pat(acut,ecut);
%   ant = phased.CustomAntennaElement('AzimuthAngles',az,...
%       'ElevationAngles',el,'MagnitudePattern',pat,...
%       'PhasePattern',zeros(size(pat)));
%   pattern(ant,3e8);
%
%   See also phased, rotpat.

%#codegen

validateattributes(azcut,{'double'},{'2d','real','nonempty','nonnan'},'azelcut2pat','ACUT');
validateattributes(elcut,{'double'},{'2d','real','nonempty','nonnan'},'azelcut2pat','ECUT');
cond = size(azcut,1)~=size(elcut,1);
if cond
    coder.internal.errorIf(cond,'phased:phased:NumRowsMismatch','ACUT','ECUT');
end

patazel = bsxfun(@plus,permute(elcut,[2 3 1]),permute(azcut,[3 2 1]));
