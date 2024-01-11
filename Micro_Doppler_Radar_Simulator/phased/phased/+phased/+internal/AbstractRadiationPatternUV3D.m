classdef (Hidden) AbstractRadiationPatternUV3D < phased.internal.AbstractRespPattern3D
%This class is for internal use only. It may be removed in the future.

%   Copyright 2013 The MathWorks, Inc.

properties
    %Pattern - Response pattern
    %   Pattern is a matrix containing the samples of the magnitude
    %   response pattern at each corresponding (u, v) pair. 
    %   The default pattern is all zero.
    Pattern
    % U - U grid
    %   U is a matrix containing the grid of U coordinates in the U/V
    %   space.
    U
    % V - V grid
    %   V is a matrix containing the grid of V coordinates in the U/V
    %   space.
    V
end

methods (Abstract, Access = protected)
    validatepattern(obj,value)
end

methods
    function obj = AbstractRadiationPatternUV3D(varargin)

        u_def = -1:0.01:1;
        v_def = -1:0.01:1;
        [u_plot,v_plot]=meshgrid(u_def,v_def);
        U=u_plot; %#ok<*PROP>
        V=v_plot;
        Pattern = zeros(size(u_plot));
        sigutils.pvparse(varargin{:});
        obj.U = U;
        obj.V = V;
        obj.Pattern = Pattern;
       
    end
    
end    

methods (Access = protected)
    function xlbl = getXLabel(obj,~)  %#ok<INUSD>
    %getXLabel Return x-axis label of the plot
        xlbl = 'U';
    end
    
    function ylbl = getYLabel(obj,~)  %#ok<INUSD>
    %getYLabel Return y-axis label of the plot
        ylbl = 'V';
    end
    
    function x = getXData(obj,~)
    %getXData Get x-axis data
          x = obj.U;
    end
    
    function y = getYData(obj,~)
    %getYData Get y-axis data
          y = obj.V;
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
        sortedList = {'Type','U','V','Pattern'};
    end
    
end

methods 
    function set.Pattern(obj,value)
        validatepattern(obj,value);
        obj.Pattern = value;
    end
    
     function set.U(obj,value)                                             
     validateattributes(value,{'numeric'},...
       {'2d','>=',-1,'<=',1},...
        sprintf('%s.U',class(obj)),'U');                               
        obj.U = value; 
     end  
     
     function set.V(obj,value)                                             
      validateattributes(value,{'numeric'},...
         {'2d','>=',-1,'<=',1},...
        sprintf('%s.V',class(obj)),'V');                                
        obj.V = value; 
     end  
        
                                                                        
    
end

end
% [EOF]
