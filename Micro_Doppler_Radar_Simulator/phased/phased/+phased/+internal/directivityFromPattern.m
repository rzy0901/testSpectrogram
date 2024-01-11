function d = directivityFromPattern(pat,az,el)
%This function is for internal use only. It may be removed in the future.

%directivityFromPattern     Compute directivity from pattern
% D = phased.interna.directivityFromPattern(PAT,AZ,EL) computes the 3D
% directivity pattern, d (in dBi), from the 3D field pattern PAT. PAT is
% defined across angles defined in azimuth and elevation angles (in
% degrees) given by AZ and EL. Note that PAT is assumed to cover the entire
% 3D space. 
%
% PAT can be either a matrix or a struct. If PAT is a matrix, it represents
% the field pattern. PAT has a dimension of PxQ where P is the number of
% angles in EL and Q is the number of angles in AZ. D has the same
% dimensions as PAT. Both AZ and EL is assumed to be uniformly sampled. If
% PAT is a struct, then it has two fields, H and V, and each field contains
% a PXQ matrix representing the field pattern in the corresponding
% polarization. The D will have two fields too where each field contains
% the directivity pattern in the corresponding polarization. 

%   Copyright 2016 The MathWorks, Inc.

%#codegen
daz = deg2rad(mean(diff(az)));
del = deg2rad(mean(diff(el)));
elr = deg2rad(el(:));
pat2 = abs(pat).^2;
intpat2 = phased.internal.integratePattern(pat2,elr,daz,del);
d = phased.internal.normalizeIntegratedPower(pat2,intpat2,true);
d = pow2db(d);

