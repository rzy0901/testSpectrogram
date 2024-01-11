
function plotGratingLobeDiagramPlanar(RowElementSpacing,ColumnElementSpacing,Lattice,f,ang,c)
%This function is for internal use only. It may be removed in the future.

%plotGratingLobeDiagram  Plot grating lobe diagram in sine space (u v coordinates)
%   RowElementSpacing,ColumnElementSpacing are defined in m (see doc
%   phased.URA)
%   Lattice can be 'Rectangular' or 'Triangular'
%   Frequency, f, defined in Hz.
%   Steering angle, ang, can be
%     either a 2x1 vector or a scalar. If angle is a vector, the
%     form is [azimuth; elevation] (in degrees) where the azimuth
%     angle must be between [-180 180] and the elevation angle
%     must be between [-90 90]. If ANGLE is a scalar, angle
%     specifies the azimuth angle and the corresponding elevation
%     angle is assumed to be 0 degrees.
%   Propagation speed , C,is defined in m/s.

%   Copyright 2013-2016 The MathWorks, Inc.

% Validate inputs
sigdatatypes.validateFrequency(f,'plotGratingLobeDiagram',...
    'FREQ',{'scalar'});
sigdatatypes.validateSpeed(c,'plotGratingLobeDiagram','C',{'scalar'});
if isscalar(ang)
    ang = [ang;0];
end
sigdatatypes.validateAzElAngle(ang,'plotGratingLobeDiagram',...
    'ANGLE',{'vector'});
az = ang(1);
el = ang(2);

lambda = c/f;
f0 = f;
u0 = cosd(el)*sind(az);
u_center = f0/f*u0;
v0 = sind(el);
v_center = f0/f*v0;

gl_u_spacing = lambda/ColumnElementSpacing;
if strcmp(Lattice,'Rectangular')
    gl_v_spacing = lambda/RowElementSpacing;
else
    gl_v_spacing = lambda/(2*RowElementSpacing);
end

nu = 10;
nl = -10;
u_grat_grid = (nl:nu)*gl_u_spacing;%
v_grat_grid = (nl:nu)*gl_v_spacing;

cla reset;
hold on;

theta = 0:5:360;
c_x = cosd(theta);
c_y = sind(theta);

% plot visible region patch
hgreen = patch(c_x,c_y,'g','Tag','green_patch');

% plot grating lobe region
alphaVal = 0.8;

for m = 1:(nu-nl+1)
    xref = u_grat_grid(m);
    for n = 1:(nu-nl+1)
        yref = v_grat_grid(n);
        if strcmp(Lattice,'Rectangular')
            if xref ~= 0 || yref ~= 0
                tag_patch = sprintf('red_patch_%d_%d',m,n);
                hred = patch(xref+c_x,yref+c_y,[1 0.54 0.54],...
                    'EdgeColor',[1 0.64 0.64],'FaceAlpha',alphaVal,'Tag',tag_patch);
            end
        else
            if (xref ~= 0 || yref ~= 0) && ~rem(2*nl+m+n-2,2)
                tag_patch = sprintf('red_patch_%d_%d',m,n);
                hred = patch(xref+c_x,yref+c_y,[1 0.54 0.54],...
                    'EdgeColor',[1 0.64 0.64],'FaceAlpha',alphaVal,'Tag',tag_patch);
            end
        end
    end
end

% plot grating lobe point
[u_plot,v_plot] = meshgrid(u_grat_grid,v_grat_grid);
if ~strcmp(Lattice,'Rectangular')
    u_plot(2:2:end) = nan;
    v_plot(2:2:end) = nan;
end
h = plot(u_center+u_plot(:),v_center+v_plot(:),'ko','MarkerSize',10,'Tag','GratingLobe_plot');
haxes = get(h,'Parent');
% plot main lobe point
hmainlobe = plot(u_center,v_center,'ko', 'MarkerEdgeColor','k',...
    'MarkerFaceColor','k','MarkerSize',10,'Tag','MainLobe_plot');
% plot a unit cirle edge
plot(c_x,c_y,'k','Tag','Edge_plot');
% text information of scan area
if strcmp(Lattice,'Rectangular')
    if (gl_u_spacing < 1)||(gl_v_spacing < 1)
        arrayspan_str = sprintf('%s\n%s',...
            getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
            getString(message('phased:system:array:AlwaysGratingLobe')));
    elseif (gl_u_spacing >= 2)&&(gl_v_spacing >= 2)
        arrayspan_str = sprintf('%s\n%s',...
            getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
            getString(message('phased:system:array:NoGratingLobe')));
    else
        apu = min((gl_u_spacing-1),1);
        apv = min((gl_v_spacing-1),1);
        arrayspan_str = sprintf('%s\nU: [%.2f %.2f] (Az: [%.1f %.1f] deg)\nV: [%.2f %.2f] (El: [%.1f %.1f] deg)',...
            getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
            -apu,apu, -asind(apu),asind(apu),-apv,apv, -asind(apv),asind(apv));
    end
else
    % triangular lattice
    if ((gl_u_spacing >= 2)&&(gl_v_spacing >= 1))||...
            ((gl_v_spacing >= 2)&&(gl_u_spacing >= 1))
        arrayspan_str = sprintf('%s\n%s',...
            getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
            getString(message('phased:system:array:NoGratingLobe')));
    elseif (gl_v_spacing < 0.5)||...
            (gl_u_spacing < 0.5)||...
            ((gl_u_spacing < 1)&&(gl_v_spacing < 1)&&((gl_u_spacing^2+gl_v_spacing^2)<1))
        % (gl_u_spacing^2+gl_v_spacing^2) is the distance from the center
        % to the grating lobe. when it is <=1 the four diagonal corner grating 
        % lobes enter visible area (unit circle)
        % The grating lobe located at the top and bottom enter are at a distance
        % of 2*gl_v_spacing. When gl_v_spacing < 0.5 those lobes enter the visible 
        % area. 
        % Similarly when gl_u_spacing <0.5 the  left and right grating lobes enter 
        % visible area.
        arrayspan_str = sprintf('%s\n%s',...
            getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
            getString(message('phased:system:array:AlwaysGratingLobe')));
    else
        if (gl_u_spacing < 1)&&(gl_v_spacing < 1)&&~(((gl_u_spacing^2+gl_v_spacing^2)<1))
            apu = min((2*gl_u_spacing-1),(gl_u_spacing-sqrt(1-gl_v_spacing^2)));
            apv = min((2*gl_v_spacing-1),(gl_v_spacing-sqrt(1-gl_u_spacing^2)));
            % Here we are trying to determine which edge is closer to the center of
            % the visible area. (2*gl_u_spacing-1) calculates the distance from the
            % center to the edge of the left (or right lobes).
            % (gl_u_spacing-sqrt(1-gl_v_spacing^2)) calculates the distance from the
            % center to the intersection of the upper left and lower right lobes.
            % min((2*gl_u_spacing-1),(gl_u_spacing-sqrt(1-gl_v_spacing^2))) will
            % give us the closer distance
            % try [0.52 1.6] lambda spacing and then try [0.52 1.2] lambda spacing
            % in sensorArrayAnalyzer to visualize those lobes.
        elseif ((gl_u_spacing >= 2)&&(gl_v_spacing < 1)&&(gl_v_spacing >= 0.5))||...
                ((gl_v_spacing >= 2)&&(gl_u_spacing < 1)&&(gl_u_spacing >= 0.5))||...
                ((gl_u_spacing < 2)&&(gl_u_spacing >= 1)&&(gl_v_spacing < 2)&&(gl_v_spacing >= 1))
            apu = min((2*gl_u_spacing-1),1);
            apv = min((2*gl_v_spacing-1),1);
        elseif (gl_v_spacing < 1)&&(gl_u_spacing < 2)&&(gl_u_spacing >= 1)
            apu = min((gl_u_spacing-sqrt(1-gl_v_spacing^2)),1);
            apv = 2*gl_v_spacing-1;
            % For now the following code will never be hit since the ratio of 
            % gl_u_spacing and gl_v_spacing is fexed to either 1 or 2/sqrt(3)
            % 
            %elseif (gl_u_spacing < 1)&&(gl_v_spacing < 2)&&(gl_v_spacing >= 1)
            %  apu = 2*gl_u_spacing-1;
            %  apv = min((gl_v_spacing-sqrt(1-gl_u_spacing^2)),1);
        end
        arrayspan_str = sprintf('%s\nU: [%.2f %.2f] (Az: [%.1f %.1f] deg)\nV: [%.2f %.2f] (El: [%.1f %.1f] deg)',...
            getString(message('phased:system:array:GratingLobeScanRegionHeader')),...
            -apu,apu, -asind(apu),asind(apu),-apv,apv, -asind(apv),asind(apv));
    end
end
xlabel('U','Tag','ULabel');
ylabel('V','Tag','VLabel')
title(getString(message('phased:system:array:GratingLobeDiagram','U-V')),'Tag','GratingLobeTitle');
set(haxes,'Color',get(gcf,'Color'));
hleg = legend([hmainlobe,h,hgreen,hred],...
    {getString(message('phased:system:array:MainLobeLegend')),...
    getString(message('phased:system:array:GratingLobeLegend')),...
    getString(message('phased:system:array:GratingLobeFreeAreaLegend')),...
    getString(message('phased:system:array:GratingLobeAreaLegend'))},...
    'Location','NorthEastOutside','Tag','GratingLobeLegend','AutoUpdate','off');
set(hleg,'FontSize',9);
insetdist = get(haxes,'TightInset');
text('Unit','Normalized','Position',[0.52 -3*insetdist(2)],'String',arrayspan_str,...
    'Color',[.501 .501 .501],'FontSize',9,'Tag','GratingLobeAnnotation');
hold off;
axis equal
axis([-3 3 -3 3]);
set(haxes,'OuterPosition',[0.1 0.25 0.8 0.7]);
zoom reset;
set(hleg,'Position',[0.2 0.1 0.28 0.15]);
end

