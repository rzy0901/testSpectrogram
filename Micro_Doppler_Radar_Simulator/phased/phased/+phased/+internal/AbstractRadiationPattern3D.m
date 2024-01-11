classdef (Hidden) AbstractRadiationPattern3D < phased.internal.AbstractRespPattern3D
%This class is for internal use only. It may be removed in the future.

%   Copyright 2013-2017 The MathWorks, Inc.

properties
    %AzAngle - Azimuth angle
    %   AzAngle is a column vector containing the azimuth angles (in
    %   degrees) at which the response pattern are sampled.
    AzAngle
    %ElAngle - Elevation angle
    %   ElAngle is a column vector containing the elevation angles (in
    %   degrees) at which the response pattern are sampled.
    ElAngle
    %Pattern - Response pattern
    %   Pattern is a matrix containing the samples of the magnitude
    %   response pattern at each corresponding (ElAngle,AzAngle) pair. The
    %   number of rows of Pattern must match the number of elements in
    %   ElAngle and the number of columns of Pattern must match the number
    %   of elements in AzAngle. The default pattern is an isotropic
    %   pattern.
    Pattern
end

methods (Abstract, Access = protected)
    validatepattern(obj,value)
end

methods
    function obj = AbstractRadiationPattern3D(varargin)
    %AbstractRadiationPattern3D Constructor of phased.internal.AbstractRadiationPattern3D class
        ElAngle = -90:90; %#ok<*PROP>
        AzAngle = -180:180;
        Pattern = ones(numel(ElAngle),numel(AzAngle));
        sigutils.pvparse(varargin{:});
        obj.ElAngle = ElAngle;
        obj.AzAngle = AzAngle;
        obj.Pattern = Pattern;
        obj.Type = '3D Radiation Pattern';
    end
    
    function varargout = polar(obj,varargin)
        %POLAR Plot response pattern in spherical-polar coordinates
        %   polar(Hresp) plots the response pattern Hresp in
        %   spherical-polar coordinates.
        %
        %   polar(Hresp, 'Units', UNIT) plots the response pattern using
        %   the unit specified in UNIT. UNIT can be any of the following:
        %   ['mag' | 'power' | {'db'}].
        %
        %   polar(..., 'NormalizeResp', NFLAG) plots the normalized
        %   response pattern if NFLAG is true. NFLAG can be any of the
        %   following: [true | {false}].
        %
        %   polar(..., 'Title', TITLE) uses TITLE as the title of the
        %   resulting plot.
        %
        %   polar(..., 'PlotPlaneCircles',PCFLAG) uses PCFLAG to enable
        %   plotting circles on three major planes. PCFLAG must be a
        %   logical scalar. The default value of PCFLAG is false.
        %
        %   h = polar(...) returns a handle to the surface object.
        %
        %   Example: 
        %       % Construct and plot a 3D isotropic array response pattern. 
        %       hresp = phased.internal.RespPattern3D; 
        %       polar(hresp);
        %
        %   See also phased.internal.RespPattern3D
            
            plotoption = getPlotOption(obj,varargin{:});
            plotoption.PlotType = 'polar';
            Response = privRespPattern(plotoption,obj.Pattern);
            dBRange = 50;
            if strcmpi(plotoption.Units,'db')
                Response = phased.internal.AbstractRadiationPattern3D.limitDynamicdBRange(Response,dBRange);
            end
            ZLbl = getRespLabel(plotoption);
            
            [RespRow,RespCol] = size(Response);

            % Setting up Az surface
            Theta = repmat(obj.AzAngle(:)',RespRow,1);
            
            % Setting up El surface
            Rho = repmat(obj.ElAngle(:),1,RespCol);
            
            % Correcting for dB conversions
            if ~isempty(strfind(ZLbl,'dB'))
                MaxResp = max(max(Response));
                MinResp = min(min(Response));
                RespRange = MaxResp - MinResp;
                if RespRange == 0
                    RespOffset = 1;
                else
                    % if the dynamic range is less than 50 dB, plot 20%
                    % more than the dynamic range, otherwise, move the
                    % minimum response to the center of polar plot
                    
                    if RespRange < dBRange
                        MinResp = max(MinResp-0.2*RespRange,...
                            MaxResp-dBRange);
                    end
                    RespOffset = -MinResp;
                end
            else
                RespOffset = 0;
            end
            
            % Creating Cartesian equivalent in radians
            [X,Y,Z] = sph2cart(deg2rad(Theta),deg2rad(Rho),Response+RespOffset);
            
            % Plotting surface in a new figure
            surfHdl = surf(X,Y,Z,'FaceColor','Interp');
            set(surfHdl,'cdata',Response,'LineStyle','none','EdgeColor','none','FaceAlpha',1,'Tag','3D polar plot');
            axHdl = get(surfHdl,'parent');
            set(axHdl,'DataAspectRatio',[1 1 1]);
            view_az = 135; view_el = 20;
            view([view_az view_el]);
            axis vis3d
            axis off;
            hold on;
            zHdl = zoom;
            zHdl.setAxes3DPanAndZoomStyle(axHdl,'camera');
            
            % Get current positions of axes
            xmax = max(max(abs(X)));
            ymax = max(max(abs(Y)));
            zmax = max(max(abs(Z)));
            maxdata = max([xmax ymax zmax]);
            xmax = max(maxdata/4,xmax);  % ensure the axis span is at least maxdata/4 
            ymax = max(maxdata/4,ymax);
            zmax = max(maxdata/4,zmax);
            if isfinite(xmax) && xmax~=0
                set(get(surfHdl,'parent'),'XLim',1.3*[-xmax xmax]);
            end
            if isfinite(ymax) && ymax~=0
                set(get(surfHdl,'parent'),'YLim',1.3*[-ymax ymax]);
            end
            if isfinite(zmax) && zmax~=0
                set(get(surfHdl,'parent'),'ZLim',1.3*[-zmax zmax]);
            end
            
            axisfactor = 1.2;
            XPos=axisfactor*xmax;
            YPos=axisfactor*ymax;
            ZPos=axisfactor*zmax;
            
            % Create pseudo axes and mark ticks
            plot3( [0,XPos],[0,0],[0,0],'k','LineWidth',1.5,'Tag','XAxis' );
            text(1.3*XPos,0,0,sprintf('x\nAz 0\nEl 0'),'Tag','Az0El0Label');
            plot3( [0,0],[0,YPos],[0,0],'k','LineWidth',1.5,'Tag','YAxis' );
            text(0,1.15*YPos,0,sprintf('y\nAz 90\nEl 0'),'Tag','Az90El0Label');
            plot3( [0,0],[0,0],[0,ZPos],'k','LineWidth',1.5,'Tag','ZAxis' );
            text(0,0,1.15*ZPos,sprintf('z\nAz 0\nEl 90'),'Tag','Az0El90Label');
            
            % Create circules
            if plotoption.PlotPlaneCircles
                circfactor = 1.1;
                XCircR=circfactor*maxdata;
                YCircR=circfactor*maxdata;
                ZCircR=circfactor*maxdata;
                circ_ang = 0:360;
                plot3(ZCircR*cosd(circ_ang),ZCircR*sind(circ_ang),...
                    zeros(size(circ_ang)),'b','LineWidth',2,'Tag','XYCircle');
                plot3(YCircR*cosd(circ_ang),zeros(size(circ_ang)),...
                    YCircR*sind(circ_ang),'g','LineWidth',2,'Tag','XZCircle');
                plot3(zeros(size(circ_ang)),XCircR*cosd(circ_ang),...
                    XCircR*sind(circ_ang),'r','LineWidth',2,'Tag','YZCircle');
            end
            
            % Create az/el arrows
            arrowang = 22.8;
            arrowlen = 0.25*maxdata;
            arrowheadlen = 0.25*arrowlen;
            XPos = 0.95*XPos;
            plot3(XPos+[0 0],[0 arrowlen],[0 0],'k','LineWidth',2,'Tag','AzArrowBody');
            plot3(XPos+[0 0],[0 0],[0 arrowlen],'k','LineWidth',2,'Tag','ElArrowBody');
            fill3(XPos+[0 -arrowheadlen*tand(arrowang) arrowheadlen*tand(arrowang)],...
                arrowlen+[0 -arrowheadlen -arrowheadlen],zeros(1,3),...
                'k','Tag','AzArrowHead');
            text(XPos,1.1*arrowlen,0,'az','Tag','AzArrowLabel');
            fill3(XPos+zeros(1,3),[0 -arrowheadlen*tand(arrowang) arrowheadlen*tand(arrowang)],...
                arrowlen+[0 -arrowheadlen -arrowheadlen],...
                'k','Tag','ElArrowHead');
            text(XPos,0,1.2*arrowlen,'el','Tag','ElArrowLabel');
            
            hold off;
            
            % Hide data cursor in polar mode
            hfig = get(get(surfHdl,'Parent'),'Parent');
            if ishghandle(hfig) && strcmp(get(hfig,'Type'),'figure')
                hdcm = datacursormode(hfig);
                set(hdcm,'Enable','off');
                set(hdcm,'UpdateFcn',@(dummy,event_obj)getString(...
                    message('phased:apps:arrayapp:DataCursorNotSupported')));
                hdcb = findall(hfig,'Tag','Exploration.DataCursor');
                set(hdcb,'Visible','off');
            end
            
            surfHdl.FaceLighting = 'gouraud';
            surfHdl.AmbientStrength = 0.5;
            surfHdl.DiffuseStrength = 0.8;
            surfHdl.SpecularStrength = 0.9;
            surfHdl.SpecularExponent = 25;
            surfHdl.BackFaceLighting = 'unlit';
            camlight infinite;
            material dull;
            title(plotoption.Title);
            
            % Setting label of colorbar based on response's (z-axis) label
            clrbarHdl = colorbar('peer',get(surfHdl,'parent'));
            set(clrbarHdl,'Tag','Colorbar for 3D polar plot')
            [fsize,lblColor] = phased.internal.AbstractRespPattern.setAnnotationSizeColor;
            ylabel(clrbarHdl,ZLbl,'fontsize',fsize,'Color',lblColor);
            
            % Adjust camera's view angle to ensure that the title is not
            % clipped
            vecTgt = campos(axHdl)-camtarget(axHdl);
            vecTitle = campos(axHdl)-get(get(axHdl,'Title'),'Position');
            angTitle = 2*acosd(vecTgt*vecTitle'/(norm(vecTgt)*norm(vecTitle)));
            camva(axHdl,max(angTitle, camva(axHdl)));
            
            if nargout == 1,
                varargout{1} = surfHdl;
            end
        end
end    

methods (Access = protected)
    function xlbl = getXLabel(obj,~)  %#ok<INUSD>
    %getXLabel Return x-axis label of the plot
        xlbl = 'Azimuth Angle (degrees)';
    end
    
    function ylbl = getYLabel(obj,~)  %#ok<INUSD>
    %getYLabel Return y-axis label of the plot
        ylbl = 'Elevation Angle (degrees)';
    end
    
    function x = getXData(obj,~)
    %getXData Get x-axis data
        x = obj.AzAngle;
    end
    
    function y = getYData(obj,~)
    %getYData Get y-axis data
        y = obj.ElAngle;
    end
    
    function annotatePlot(obj,plotoptionobj)
        annotatePlot@phased.internal.AbstractRespPattern3D(obj,plotoptionobj);
        resplbl = getRespLabel(plotoptionobj);
        zlabel(resplbl);
        hcbar = colorbar;
        [fsize,lblColor] = phased.internal.AbstractRespPattern.setAnnotationSizeColor;
        ylabel(hcbar,resplbl,'fontsize',fsize,'Color',lblColor);
    end
    
end


methods (Access = protected)
    function sortedList = getSortedPropDispList(this)  %#ok<MANU>
        % Get the sorted list of the properties to be displayed. 
        sortedList = {'Type','ElAngle','AzAngle','Pattern'};
    end
    
end

methods 
    function set.Pattern(obj,value)
        validatepattern(obj,value)
        obj.Pattern = value;
    end
    
    function set.AzAngle(obj,value)
        sigdatatypes.validateAngle(value,...
            sprintf('%s.AzAngle',class(obj)),'AzAngle',...
            {'vector','>=',-180,'<=',180});
       if length(value) == 1
            error(message('phased:phased:internal:RespPattern3D:AzScalar'));
        end
        obj.AzAngle = value;
    end

    function set.ElAngle(obj,value)
        sigdatatypes.validateAngle(value,...
            sprintf('%s.ElAngle',class(obj)),'ElAngle',...
            {'vector','>=',-90,'<=',90});
        if length(value) == 1
            error(message('phased:phased:internal:RespPattern3D:ElScalar'));
        end
        obj.ElAngle = value;
    end
end

end

% [EOF]
