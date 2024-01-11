function [estpos, esttaper] = pilotcalib(nompos_in, x, pilotang_in, nomtaper, uncerts)
% pilotcalib Array calibration using pilot sources
%   ESTPOS = pilotcalib(NOMPOS, X, PILOTANG) returns the estimated element
%   positions of the sensor array in ESTPOS. NOMPOS represents the nominal
%   positions of the sensor array before calibration. X represents the
%   signals received from the array in the presence of the pilot sources.
%   PILOTANG represents the directions of each of the pilot sources. A
%   minimum of 3 pilot sources are required in this case.
%
%   NOMPOS specifies the locations of elements in the sensor array in the
%   unit of signal wavelength. All elements in the sensor array are assumed
%   to be isotropic. NOMPOS must be a 3xN matrix where N is the number of
%   elements in the sensor array. Each column of NOMPOS represents the
%   [x;y;z] coordinates of the corresponding element. The nominal position
%   of all sensors must be within half a wavelength of their true positions
%   for the calibration to be successful.
%
%   X represents the signals received from the array when pilot sources are
%   transmitting. It is assumed that pilot sources transmit independently
%   and a narrowband plane wave signal is received in the absence of
%   multipath. X must be an LxNxM matrix where L is the number of snapshots
%   received by the array for each pilot source and M is the number of
%   pilot sources.
%
%   PILOTANG represents the directions of the pilot sources and must be a
%   2xM matrix where M is the number of pilot sources. Each column
%   specifies the direction in [azimuth; elevation] in degrees.  The
%   azimuth angle must be between -180 and 180 degrees and the elevation
%   angle must be between -90 and 90 degrees. The azimuth angle is defined
%   in the xy plane; it is the angle measured from the x axis, which is
%   also the array normal direction, toward the y axis. The elevation angle
%   is defined as the angle from the xy plane toward the z axis.
%   Calibration sources must span over different azimuth and elevation
%   angles.
%
%   ESTPOS represents the estimated locations of the elements in the sensor
%   array, specified in the unit of signal wavelength. ESTPOS is a 3xN
%   matrix where each column of ESTPOS represents the [x;y;z] coordinates
%   of the corresponding element.
%
%   [ESTPOS, ESTTAPER] = pilotcalib(NOMPOS, X, PILOTANG) returns ESTTAPER
%   as the estimated taper acting on the array. Hence, the taper associated
%   with the array as well as the array location is estimated. ESTTAPER is
%   an Nx1 vector with each element representing the taper of the
%   corresponding element. In this case, it is assumed that the nominal
%   taper is a vector of ones. A minimum of 4 pilot sources are required.
%
%   [ESTPOS, ESTTAPER] = pilotcalib(NOMPOS, X, PILOTANG, NOMTAPER)
%   specifies NOMTAPER as the nominal taper acting on the array. The taper
%   associated with the array as well as the array location is estimated.
%   NOMTAPER must be an Nx1 vector with each element representing the
%   nominal taper of the corresponding element. A minimum of 4 pilot
%   sources are required.
%
%   [ESTPOS, ESTTAPER] = pilotcalib(NOMPOS, X, PILOTANG, NOMTAPER, UNCERTS)
%   specifies UNCERTS as the configuration settings to use for calibrating
%   the array. UNCERTS must be a 4x1 vector with elements having a value of
%   0 or 1 to denote which uncertainties should be estimated. The vector
%   takes the form of [x;y;z;taper] where setting x, y or z to 1 will
%   estimate uncertainties associated with the x, y or z dimension of the
%   array. Setting taper to 1 will estimate the taper. For example, setting
%   UNCERTS to [0;1;1;1] will estimate uncertainties in the y and z axes of
%   the array and the taper. The number of pilot sources must be at least
%   equal to the number of 1's in the vector UNCERTS.
%
%   % Example:
%   %   Consider a 6-element isotropic ULA with half-wavelength sensor
%   %   spacing which has been geometrically perturbed in the x-y plane and
%   %   contains an unknown taper. Perform pilot calibration on the array 
%   %   using 3 pilot sources at (azimuth, elevation) angles in degrees of 
%   %   (-60,0), (10,0) and (40,0). For the calibration process, an SNR of 
%   %   40dB is assumed and 10000 snapshots of data are recorded for each 
%   %   pilot.
%
%   NominalElementPositions = [zeros(1,6);-1.25:0.5:1.25;zeros(1,6)];
%   NominalTaper = ones(6,1);
%   L=10000;
%   N=length(NominalTaper);
%   ncov = 0.001;
%   PilotAng = [-60,10,40; 0, 0, 0];
%   M = size(PilotAng,2);
%   TrueElementPositions = NominalElementPositions + ...
%       0.01*[zeros(3,1),[randn(2,N-1);zeros(1,N-1)]];
%   TrueTaper = NominalTaper + [0; randn(N-1,1).*exp(1j*randn(N-1,1))];
%   for m = 1:M
%       X(:,:,m) = sensorsig(TrueElementPositions,L,PilotAng(:,m),ncov,...
%           'Taper',TrueTaper);
%   end
%   [estpos,esttaper] = pilotcalib(NominalElementPositions,X,PilotAng,...
%       NominalTaper, [1;1;0;1]);
%
%   % Check error after calibration
%   posRMSE = sqrt(sum(sum((TrueElementPositions - estpos).^2)/N))
%   taperMagRMSE = sqrt(sum((abs(TrueTaper) - abs(esttaper)).^2)/N)
%   taperPhaseRMSE = sqrt(sum((angle(TrueTaper) - angle(esttaper)).^2)/N)

%   Copyright 2014 The MathWorks, Inc.

%   References
%   [1] N. Fistas and A. Manikas, "A New General Global Array Calibration
%   Method", IEEE Proceedings of ICASSP, Vol. IV, pp. 73-76, April 1994.

%#ok<*EMCA>
%#codegen

if nargin <4
    nomtaper = ones(size(x,2),1);
end
if nargin<5
    if nargout<2
        uncerts = [1;1;1;0];
    else
        uncerts = ones(4,1);
    end
end

esttaper = nomtaper;

L = size(x,1); % Number of snapshots
N = size(x,2); % Number of sensors
M = size(x,3); % Number of pilots

validateattributes(uncerts,{'double'},{'binary','size',[4 1]},...
    'pilotcalib','UNCERTS');

if size(nompos_in,1) == 1
    nompos = [zeros(1,size(nompos_in,2));nompos_in;zeros(1,size(nompos_in,2))];
elseif size(nompos_in,1) == 2;
    nompos = [zeros(1,size(nompos_in,2));nompos_in];
else
    nompos = nompos_in;
end
sigdatatypes.validate3DCartCoord(nompos,'pilotcalib','NOMPOS',{'size',[3 N]});
estpos = nompos;


if size(pilotang_in,1) == 1
    pilotang = [pilotang_in;zeros(1,size(pilotang_in,2))];
else
    pilotang = pilotang_in;
end
sigdatatypes.validateAzElAngle(pilotang,'pilotcalib','PILOTANG',{'size',[2 M]});

validateattributes(nomtaper,{'double'},{'size',[N 1]},'pilotcalib','NOMTAPER');

NumUncerts = sum(uncerts);
cond = (M<NumUncerts);
if cond
    coder.internal.errorIf(cond,'phased:system:array:InsufficientPilots',NumUncerts);
end

cond = (uncerts(3) == 1) && (length(unique(pilotang(2,:)))<2);
if cond
    coder.internal.errorIf(cond,'phased:system:array:CalibInZ')
end

cond = (sum(uncerts(1:2)) > 0) && (length(unique(pilotang(1,:)))<2);
if cond
    coder.internal.errorIf(cond,'phased:system:array:CalibInXY')
end

Es = complex(zeros(N,M));
for m=1:M
    Rxx = x(:,:,m).'*conj(x(:,:,m))/L;
    [E,D]=eig(Rxx);
    d = diag(D);
    [~,di] = sort(d);
    E = E(:,di);
    Es(:,m)=(E(:,end)/E(1,end));
end

PilotKI = [cosd(pilotang(1,:)).*cosd(pilotang(2,:));...
        sind(pilotang(1,:)).*cosd(pilotang(2,:));sind(pilotang(2,:))];

if any(uncerts(1:3)) % If sensor location uncertainties

    PilotKI_t = PilotKI;

    % Shift the reference point to sensor 1 (the known sensor) for the
    % purposes of the estimation of uncertainties
    nompos_t = nompos - nompos(:,1)*ones(1,size(nompos,2));
    
    if uncerts(4) % If taper uncertainties
        Es_divided=Es(:,2:end)./(Es(:,1)*ones(1,M-1));
        angle_Es_t = angle(Es_divided);
        PilotKI_t = PilotKI(:,2:end)-PilotKI(:,1)*ones(1,M-1);
    else
        angle_Es_t = angle(Es./(nomtaper*ones(1,M)));
    end
        
    K = round(-((1/(2*pi))*angle_Es_t)+(nompos_t.'*PilotKI_t));
    
    angle_Es_t = angle_Es_t + 2*pi*K - 2*pi*nompos_t.'*PilotKI_t;

    PilotKI_t = PilotKI_t(logical(uncerts(1:3)),:);
    
    uncert_est = ((1/(2*pi))*angle_Es_t*pinv(PilotKI_t)).';
    
    estpos(logical(uncerts(1:3)),:) = nompos(logical(uncerts(1:3)),:) + ...
        uncert_est;
    
end
    
if uncerts(4) 
    
    Taper_est_vec=Es./exp(1j*2*pi*estpos.'*PilotKI);
    Taper_est_vec=Taper_est_vec./(ones(N,1)*Taper_est_vec(1,:));
	esttaper=((sum(Taper_est_vec,2)/M)*nomtaper(1));

end