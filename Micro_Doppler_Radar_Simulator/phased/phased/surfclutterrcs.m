function rcs = surfclutterrcs(nrcs,r,thetaa,thetae,thetag,tau,c)
%surfclutterrcs Surface clutter radar cross section
%   RCS = surfclutterrcs(NRCS,R,BEAMAZ,BEAMEL,GRAZ,TAU) returns the radar
%   cross section (RCS) of a clutter patch that is of range R (in meters)
%   away from the radar system. The clutter patch's normalized radar cross
%   section (in m^2/m^2) is specified in NRCS.
%
%   BEAMAZ and BEAMEL are the radar's azimuth and elevation beamwidths (in
%   degrees) corresponding to the clutter patch, respectively. GRAZ is
%   the grazing angle of the clutter patch that is relative to the radar.
%   TAU is the pulse width (in seconds) of the transmitted signal.
%
%   The calculation automatically determines whether the surface clutter
%   area is beam limited or pulse limited based on the given parameters.
%
%   RCS = surfclutterrcs(NRCS,R,BEAMAZ,BEAMEL,GRAZ,TAU,C) specifies the
%   propagation speed in C (in m/s). The default value of C is the speed
%   of light.
%
%   % Example: 
%   %   Calculate the RCS of a clutter patch and estimate the 
%   %   clutter-to-noise  ratio at the receiver. Assume that the patch has 
%   %   an NRCS of 1 m^2/m^2 and is 1000 meters away from the radar system.
%   %   The azimuth and elevation beamwidths are 1 degree and 3 degrees, 
%   %   respectively. The grazing angle is 10 degrees. The pulse width is 
%   %   10 microseconds. The radar is operated at a wavelength of 1 cm with
%   %   a peak power of 5 kw.
%
%   nrcs = 1; rng = 1000; bw_az = 1; bw_el = 3; graz = 10;
%   tau = 10e-6; lambda = 0.01; ppow = 5000;
%   rcs = surfclutterrcs(nrcs,rng,bw_az,bw_el,graz,tau);
%   cnr = radareqsnr(lambda,rng,ppow,tau,'rcs',rcs)
%
%   See also grazingang, radareqsnr, surfacegamma.

%   Copyright 2011 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(6,7,nargin);

if nargin < 7
    c = physconst('lightspeed');
end
eml_assert_no_varsize(1:nargin,nrcs,r,thetaa,thetae,thetag,tau,c);
sigdatatypes.validateArea(nrcs,'surfclutterrcs','NRCS',{'scalar'});
sigdatatypes.validateDistance(r,'surfclutterrcs','R',{'scalar'});
sigdatatypes.validateAngle(thetaa,'surfclutterrcs','BEAMAZ',...
    {'scalar','nonnegative'});
sigdatatypes.validateAngle(thetae,'surfclutterrcs','BEAMEL',...
    {'scalar','nonnegative'});
sigdatatypes.validateAngle(thetag,'surfclutterrcs','GRAZ',...
    {'scalar','nonnegative'});
sigdatatypes.validateDuration(tau,'surfclutterrcs','TAU',{'scalar'});
sigdatatypes.validateSpeed(c,'surfclutterrcs','C',{'scalar','positive'});

limitLength = min(r*phased.internal.deg2rad(thetae)*...
                  cotd(thetag),...  % beam limited
                  c*tau/2);    % pulse limited
A = r*phased.internal.deg2rad(thetaa)*limitLength*secd(thetag);
rcs = nrcs*A;




% [EOF]
