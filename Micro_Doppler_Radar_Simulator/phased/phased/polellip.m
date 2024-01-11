function [tau,epsilon,ar,rs] = polellip(fv)
%polellip   Polarization ellipse
%   TAU = polellip(FV) returns the tilt angle (in degrees) of the
%   polarization ellipse of the field specified in FV. FV can be either a
%   row vector or a 2-row matrix containing the linear polarization
%   representation of the field.
%
%   If FV is a matrix, each column in FV represents the field in the form
%   of [Eh;Ev], where Eh is the horizontal component of the field and Ev is
%   the vertical component of the field. If FV is a vector, each entry in
%   FV represents the polarization ratio, Ev/Eh, of a unit-power field. 
%
%   TAU is a row vector whose number of columns matches the number of
%   columns in FV. Each entry in TAU represents the tilt angle of the
%   polarization ellipse associated with the corresponding field specified
%   in FV. TAU is defined as the angle between the horizontal axis and the
%   major axis of the polarization ellipse, and is within the range of [-90
%   90].
%
%   [TAU,EPSILON] = polellip(FV) also returns the ellipticity angle (in
%   degrees) of the polarization ellipse in EPSILON. The ellipticity angle
%   is measured from major axis, and is within the range of [-45 45].
%
%   [TAU,EPSILON,AR] = polellip(FV) also returns the axial ratio of the
%   polarization ellipse in AR. The AR is defined as the ratio of the major
%   axis to the minor axis.
%
%   [TAU,EPSILON,AR,RS] = polellip(FV) also returns the rotation sense of
%   the polarization ellipse in RS. The RS is a string that describes the
%   polarization. It can be one of 'Linear', 'Left Circular', 'Right
%   Circular', 'Left Elliptical', 'Right Elliptical'.
%
%   polellip(FV) plots the polarization ellipse of the field specified in
%   FV. This syntax requires FV to have only 1 column. 
%
%   % Examples:
%
%   % Example 1:
%   %   Calculate the tilt angle and ellipticity angle of a 45 degree 
%   %   linear polarization field.
%
%   fv = [1;1];
%   [tau,epsilon] = polellip(fv)
%
%   % Example 2:
%   %   Calculate the tilt angles, ellipticity angles, axial ratios, and
%   %   rotation senses of both left and right circular polarization 
%   %   fields.
%
%   fv = [1 1i;1 -1i].';
%   [tau,epsilon,ar,rs] = polellip(fv)
%
%   % Example 3:
%   %   Plot the polarization ellipse of an elliptically polarized field. 
%   %   The field is specified by a polarization ratio.
%   
%   fv = 1+1i;
%   polellip(fv)
%
%   See also phased, polratio, pol2circpol, circpol2pol, stokes.

%   Copyright 2012-2013 The MathWorks, Inc.

%   References:
%   [1] Harold Mott, Polarization in Antennas and Radar, John Wiley & Sons,
%   1986
%   [2] Warren L. Stutzman, Polarization in Electromagnetic Systems, Artech
%   House, 1993
%   [3] IEEE Std 145-1983, IEEE Standard Definitions of Terms for Antennas,
%   1983
%   [4] IEEE Std 149-1979, IEEE Standard Test Procedures for Antennas, 1979

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(1,1,nargin);

[Nrow, Ncol] = size(fv);

cond = (Nrow ~= 1) && (Nrow ~= 2);

if cond
    coder.internal.errorIf(cond, ...
        'phased:polarization:invalidPolarizationInput','FV');
end 

if Nrow == 1
    validateattributes(fv,{'double'},{'nonnan','row'},'polellip','P');
    Eh = ones(1,Ncol);
    Ev = fv;
    Eh(isinf(fv)) = 0;
    Ev(isinf(fv)) = 1;
    Emag = hypot(Eh,Ev);
    Eh = Eh./Emag;
    Ev = Ev./Emag;
else %Nrow == 2
    validateattributes(fv,{'double'},{'finite','2d','nrows',2},...
        'polellip','FV');
    Eh = fv(1,:);
    Ev = fv(2,:);
    if any((Eh==0)&(Ev==0))
        coder.internal.errorIf(any((Eh==0)&(Ev==0)), ...
            'phased:phased:zeroColumns','FV');
    end
end

if ~nargout && Ncol>1
    if Nrow == 1
        in_name = 'P';
    else
        in_name = 'FV';
    end
    coder.internal.errorIf(~nargout && Ncol>1, ...
        'phased:phased:invalidColumnNumbers',in_name,1);
end

ExAbs = abs(Eh);
EyAbs = abs(Ev);
ExPhi = angle(Eh);
EyPhi = angle(Ev);
phi = EyPhi-ExPhi;

tol = sqrt(eps);

if nargout
    
    % Initialize
    tau = zeros(1,Ncol);
    epsilon = zeros(1,Ncol);
    ar = inf(1,Ncol);          % IEEE Standard 145-1983, ar = major/minor
    
    idx1 = (ExAbs==0);
    idx2 = (EyAbs==0);
    idx3 = ~(idx1|idx2);
    
    % --------------
    % special linear
    % --------------
    
    % vertical
    tau(idx1) = 90;
    % epsilon(idx1) = 0; % true by default
    % ar(idx1) = inf;    % true by default
    
    % horizontal
    % tau(idx2) = 0;     % true by default    
    % epsilon(idx2) = 0; % true by default
    % ar(idx2) = 0;      % true by default
    
    % ----------
    % general
    % ----------
    
    EyAbs3 = EyAbs(idx3);
    ExAbs3 = ExAbs(idx3);
    phi3 = phi(idx3);
    
    % tilt angle
    
    tanalpha = EyAbs3./ExAbs3;
    
    % tilt angle should be between -pi/2 to pi/2 while the above
    % equation returns value between -pi/4 to pi/4. The tilt angle
    % (orientation) is defined as the angle between x axis and major
    % axis. However, the above tau is between x axis and the closest
    % axis. Therefore, depending on Eh and Ev, the angle needs to be
    % adjusted.
    
    tan2alpha = 2*tanalpha./(1-tanalpha.^2);
    tau3 = atan(tan2alpha.*cos(phi3))/2;
    convert_idx = (EyAbs3>ExAbs3);
    tau3(convert_idx) = tau3(convert_idx)-sign(tau3(convert_idx))*pi/2;
    tau(idx3) = phased.internal.rad2deg(tau3);
    
    % ellipticity
    % ellipticity is measured from major axis
    sin2alpha = 2*tanalpha./(1+tanalpha.^2);
    temp = sin2alpha.*sin(phi3);
    temp(temp>1) = 1;
    temp(temp<-1) = -1;
    epsilon3 = asin(temp)/2;
    % The result is between -pi/4 to pi/4
    epsilon3 = phased.internal.rad2deg(epsilon3);
    epsilon(idx3) = epsilon3;
    
    % axis ratio
    % epsilon negative -> right-handed. However, according to IEEE standard
    % 145-1983, axis ratio should be positive when right-handed, so we need
    % to negate the result.
    ar3 = -tand(epsilon3);
    ar3(abs(epsilon3)<=tol)=0;
    ar(idx3) = 1./ar3;
    
    % rotation sense
    if nargout == 4
        cond = isempty(coder.target);
        if ~cond
            coder.internal.assert(cond,'phased:phased:invalidCodegenOutputWithOutputs','polellip',4);
        end
        rs = cell(1,Ncol);
        
        % linear
        lin_idx = (abs(epsilon)<=tol);
        rs(lin_idx) = cellstr(getString(message(...
            'phased:polarization:linear')));
        
        % circular
        left_idx = (epsilon>tol);
        right_idx = (epsilon<-tol);
        circ_idx = (abs(abs(ar)-1)<tol);
        lc_idx = (left_idx&circ_idx);
        rc_idx = (right_idx&circ_idx);
        rs(lc_idx) = cellstr(getString(message(...
            'phased:polarization:leftCircular')));
        rs(rc_idx) = cellstr(getString(message(...
            'phased:polarization:rightCircular')));
        
        % elliptical
        ellip_idx = ~(lin_idx|circ_idx);
        le_idx = (left_idx&ellip_idx);
        re_idx = (right_idx&ellip_idx);
        rs(le_idx) = cellstr(getString(message(...
            'phased:polarization:leftElliptical')));
        rs(re_idx) = cellstr(getString(message(...
            'phased:polarization:rightElliptical')));
        
    end  
    
else  
    cond = isempty(coder.target);
    if ~cond
        coder.internal.assert(cond,'phased:phased:invalidCodegenOutput','polellip');
    end
    if ishold
        holdstatus = true;
    else
        holdstatus = false;
        cla;
        hold on;
    end
   
    if Ev == 0
        plot([-ExAbs ExAbs],[0 0],'LineWidth',2,'Tag','HPolLine');
        arrow_x1 = [-ExAbs/3 -ExAbs/4 -ExAbs/3];
        arrow_y1 = [0.05 0 -0.05];
        arrow_x2 = -arrow_x1;
        arrow_y2 = arrow_y1;
        plot([arrow_x1 nan arrow_x2], [arrow_y1 nan arrow_y2],...
            'LineWidth',2,'Tag','ArrowHead');
    elseif Eh == 0
        plot([0 0],[-EyAbs EyAbs],'LineWidth',2,'Tag','VPolLine');
        arrow_x1 = [-0.05 0 0.05];
        arrow_y1 = [-EyAbs/3 -EyAbs/4 -EyAbs/3];
        arrow_x2 = arrow_x1;
        arrow_y2 = -arrow_y1;
        plot([arrow_x1 nan arrow_x2], [arrow_y1 nan arrow_y2],...
            'LineWidth',2,'Tag','ArrowHead');
    else
        EhField = linspace(-ExAbs,ExAbs,1000);
        a = 1/(EyAbs^2);
        b = -2*cos(phi)/(ExAbs*EyAbs).*EhField;
        c = EhField.^2/(ExAbs^2)-sin(phi)^2;
        temp1 = -b;
        temp2 = sqrt(abs(b.^2-4.*a.*c));
        % use abs because numbers could be close and numeric error could lead to negative numbers
        EvFieldP = (temp1 + temp2)./(2*a);
        EvFieldN = (temp1 - temp2)./(2*a);
        plot([EhField fliplr(EhField)],...
            [EvFieldP fliplr(EvFieldN)],'LineWidth',2,'Tag','PolEllipse','Color','b');
        
        % plot propagation direction
        plot(0,0,'Marker','o','MarkerSize',5,'MarkerFaceColor','k',...
            'Color','k','Tag','ZAxisDot');
        plot(0,0,'Marker','o','MarkerSize',10,'MarkerFaceColor','none',...
            'Color','k','Tag','ZAxisCircle');
        text(ExAbs/20,0,'z','Tag','ZAxisLabel');

        % plot axes
        % depending on the tiltness, the box limit is either in
        % horizontal or vertical
        [taup,epsilonp] = polellip([Eh;Ev]);
        if EyAbs > ExAbs
            plot(EyAbs/tand(taup)*[-1 1],EyAbs*[-1 1],...
                'Tag','MajorAxis','Color','b');
            plot(ExAbs*[-1 1],ExAbs/tand(taup)*[1 -1],...
                'Tag','MinorAxis','Color','b');
        else
            plot(ExAbs*[-1 1],ExAbs*tand(taup)*[-1 1],...
                'Tag','MajorAxis','Color','b');
            plot(EyAbs*tand(taup)*[1 -1],EyAbs*[-1 1],...
                'Tag','MinorAxis','Color','b');
        end
        
        % plot rotation arrow
        baseidx = 250;
        if epsilonp > tol
            % left
            [arrow_x,arrow_y] = arrowpoints(...
                EhField(baseidx:-100:baseidx-100),...
                EvFieldP(baseidx:-100:baseidx-100));
            plot(arrow_x,arrow_y,...
                'LineWidth',2,'Tag','ArrowHead','Color','b');
        elseif epsilonp < -tol
            % right
            [arrow_x,arrow_y] = arrowpoints(...
                EhField(baseidx:100:baseidx+100),...
                EvFieldP(baseidx:100:baseidx+100));
            plot(arrow_x,arrow_y,...
                'LineWidth',2,'Tag','ArrowHead','Color','b');
        else
            % linear
            [arrow_x1,arrow_y1] = arrowpoints(...
                EhField(baseidx:-100:baseidx-100),...
                EvFieldP(baseidx:-100:baseidx-100));
            baseidx = 750;
            [arrow_x2,arrow_y2] = arrowpoints(...
                EhField(baseidx:100:baseidx+100),...
                EvFieldP(baseidx:100:baseidx+100));
            plot([arrow_x1 nan arrow_x2], [arrow_y1 nan arrow_y2],...
                'LineWidth',2,'Tag','ArrowHead');

        end
    
    end
       
    if ~holdstatus
        hold off;
    end

    title(getString(message('phased:polarization:polarizationEllipse')),...
        'Tag','EllipseTitle');
    xlabel(getString(message('phased:polarization:Eh')),...
        'Tag','EllipseXLabel');
    ylabel(getString(message('phased:polarization:Ev')),...
        'Tag','EllipseYLabel');
    axis equal;
    grid on;

end
    
function [x,y] = arrowpoints(arrow_u,arrow_v)
%arrowpoints Compute the points used to draw arrow

diffu = diff(arrow_u);
line_slope = diff(arrow_v)/diffu;
pline_x0 = mean(arrow_u);
pline_y0 = mean(arrow_v);
arrow_span = norm(pline_x0-arrow_u(1),pline_y0-arrow_v(1))*tand(30);
cos_alpha = 1/sqrt(1+line_slope^2);
sin_alpha = line_slope/sqrt(1+line_slope^2);
xdist = arrow_span*sin_alpha;
ydist = arrow_span*cos_alpha;
x = [pline_x0-xdist arrow_u(1) pline_x0+xdist];
y = [pline_y0+ydist arrow_v(1) pline_y0-ydist];

% [EOF]
