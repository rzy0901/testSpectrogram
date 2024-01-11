function polchan = polchanmat(fv_t,fv_r,pos_t,ax_t,pos_r,ax_r,pos_s,ax_s,scat_s)
%This function is for internal use only. It may be removed in the future.

%polchan Polarized scattering channel
%   POLCHAN =
%   phased.internal.polchanmat(FVT,FVR,POST,AXT,POSR,AXR,POSS,AXS,SCATS)
%   returns an NrxNt channel matrix POLCHAN where Nr and Nt are number of
%   elements in receive array and number of elements in transmit array,
%   respectively. FVT (Ntx1 in H and V fields) is polarized steering vector
%   at transmit array. FVR (Ntx1 in H and V fields) is polarized steering
%   vector at receive array. POSR (3x1), POST (3x1), and POSS (3x1) are
%   positions for receive array, transmit array, and scatterers,
%   respectively. AXR (3x3), AXT (3x3), and AXS (3x3) are local coordinates
%   axes for receive array, transmit array, and scatterers, respectively.
%   SCATS (2x2) is the scattering matrix.

%   Copyright 2018 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% convert field in tx to global coordinates

scoord_t = phased.internal.global2localcoord(pos_s,'rs',pos_t,ax_t);
fv_t_local = sph2cartvec([fv_t;zeros(1,size(fv_t,2))],...    % no radial components
    scoord_t(1),scoord_t(2)); 
fv_t_global = phased.internal.local2globalvec(fv_t_local,ax_t);

% convert to target incident field
fv_sin_local = phased.internal.global2localvec(fv_t_global,ax_s);
tcoord_s = phased.internal.global2localcoord(pos_t,'rs',pos_s,ax_s);
fv_sin_sph = cart2sphvec(fv_sin_local,tcoord_s(1),tcoord_s(2));
fv_sin_hv = fv_sin_sph(1:2,:);

% apply scattering
fv_sout_hv = scat_s*fv_sin_hv;

% convert scattered field to global coordinates
rcoord_s = phased.internal.global2localcoord(pos_r,'rs',pos_s,ax_s);
fv_sout_sph = [fv_sout_hv;zeros(1,size(fv_sout_hv,2))];
fv_sout_local = sph2cartvec(fv_sout_sph,rcoord_s(1),rcoord_s(2));
fv_sout_global = phased.internal.local2globalvec(fv_sout_local,ax_s);

% % convert rx to global
% scoord_r = phased.internal.global2localcoord(pos_s,'rs',pos_r,ax_r);
% fv_r_local = sph2cartvec([fv_r;zeros(1,size(fv_r,2))],...
%     scoord_r(1),scoord_r(2));
% fv_r_global = phased.internal.local2globalvec(fv_r_local,ax_r);
% 
% % collapse to form channel matrix
% polchan = fv_r_global.'*fv_sout_global;

% convert to rx local coordinates
fv_rin_local = phased.internal.global2localvec(fv_sout_global,ax_r);
scoord_r = phased.internal.global2localcoord(pos_s,'rs',pos_r,ax_r);
fv_rin_sph = cart2sphvec(fv_rin_local,scoord_r(1),scoord_r(2));
fv_rin_hv = fv_rin_sph(1:2,:);

% collapse with receiver
polchan = fv_r.'*fv_rin_hv;


