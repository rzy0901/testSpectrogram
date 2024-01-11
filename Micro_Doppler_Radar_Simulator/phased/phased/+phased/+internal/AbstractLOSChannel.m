classdef (Hidden) AbstractLOSChannel < phased.internal.AbstractFreeSpace
%This class is for internal use only. It may be removed in the future.

%   Copyright 2015-2016 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable, Logical)
        %SpecifyAtmosphere  Specify atmosphere parameters
        %   Set this property to true to specify atmosphere parameters and
        %   consider the loss due to atmosphere gases. Set this property to
        %   false to ignore atmosphere effect in propagation. The default
        %   value of this property is false.
        SpecifyAtmosphere = false
    end
    
    properties (Nontunable)
        %Temperature    Temperature (degrees Celsius)
        %   Specify the temperature of the propagation channel in degrees
        %   Celsius as a scalar. The default value of this property is 15.
        %   This property only applies when you set the SpecifyAtmosphere
        %   property to true.
        Temperature = 15
        %DryAirPressure     Dry air pressure (Pa)
        %   Specify the dry air pressure of the propagation channel in
        %   Pascal (Pa) as a scalar. The default value of this property is
        %   101325 (1 atmosphere). This property only applies when you set
        %   the SpecifyAtmosphere property to true.
        %
        %   The ITU atmosphere gas model is valid between 1 to 1000 GHz.
        DryAirPressure = 101325
        %WaterVapourDensity Water vapour density (g/m^3)
        %   Specify the water vapour density of the propagation channel in
        %   g/m^3 as a scalar. The default value of this property is 7.5.
        %   This property only applies when you set the SpecifyAtmosphere
        %   property to true.
        %
        %   The ITU atmosphere gas model is valid between 1 to 1000 GHz.
        WaterVapourDensity = 7.5
        %LiquidWaterDensity   Liquid water density (g/m^3)
        %   Specify the liquid water density of the propagation channel in
        %   g/m^3 as a scalar. The default value of this property is 0. The
        %   liquid water density is used to describe the characteristics of
        %   fog or cloud. Typical values for water density are 0.05 for
        %   medium fog and 0.5 for thick fog. This property only applies
        %   when you set the SpecifyAtmosphere property to true.
        %
        %   The ITU fog model is valid between 10 to 1000 GHz.
        LiquidWaterDensity = 0
        %RainRate   Rain rate (mm/h)
        %   Specify the rain rate of the propagation channel in mm/h as a
        %   scalar. The default value of this property is 0, indicating no
        %   rain. This property only applies when you set the
        %   SpecifyAtmosphere property to true.
        %
        %   The ITU rain model is valid between 1 to 1000 GHz.
        RainRate = 0
    end
    
    methods
        function set.Temperature(obj,value)
            obj.Temperature = sigdatatypes.validateCFTemperature(value,...
                '','Temperature',{'scalar','>=',-273.15});
        end
        
        function set.DryAirPressure(obj,value)
            obj.DryAirPressure = sigdatatypes.validatePressure(value,...
                '','DryAirPressure',{'scalar','positive'});
        end
        
        function set.WaterVapourDensity(obj,value)
            obj.WaterVapourDensity = sigdatatypes.validateDensity(value,...
                '','WaterVapourDensity',{'scalar'});
        end
        
        function set.LiquidWaterDensity(obj,value)
            obj.LiquidWaterDensity = sigdatatypes.validateDensity(value,...
                '','LiquidWaterDensity',{'scalar'});
        end
        
        function set.RainRate(obj,value)
            obj.RainRate = sigdatatypes.validateSpeed(value,...
                '','RainRate',{'scalar'});
        end
    end
    
    methods (Access = protected)
        function obj = AbstractLOSChannel(varargin)
            obj@phased.internal.AbstractFreeSpace(varargin{:});
        end
    end

    methods (Access = protected)

        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractFreeSpace(obj, prop);
            if ~obj.SpecifyAtmosphere && ...
                    (strcmp(prop,'DryAirPressure') || ...
                    strcmp(prop,'WaterVapourDensity') || ...
                    strcmp(prop,'LiquidWaterDensity') || ...
                    strcmp(prop,'RainRate') || ...
                    strcmp(prop,'Temperature') )
                flag = true;
            end
        end
                
        function tau = computePolarizationTiltAngle(obj,x,startLoc,endLoc)
            switch obj.pValidFields
                case 'XYZ'
                    % the propagation direction is determined by start and ending
                    % locations so the polarization of H and V is defined along
                    % that axis
                    propaxis = bsxfun(@minus,endLoc,startLoc);
                    % normalize
                    propaxis = bsxfun(@rdivide,propaxis,sqrt(sum(abs(propaxis).^2,1)));
                    azel = phased.internal.dirvec2azel(propaxis);
                    polvec = complex(zeros(3,size(azel,2)));
                    % assume x is a 3-row matrix, each column is a vector defining
                    % the field in Cartesian
                    for m = 1:size(azel,2)
                        polvec(:,m) = cart2sphvec(x(:,m),azel(1,m),azel(2,m));
                    end
                    polvec_nz = any(polvec(1:2,:));
                    tau = zeros(1,size(azel,2));
                    tau(polvec_nz) = polellip(polvec(1:2,polvec_nz));
                case 'HV'
                    % assume x is a 2-row vector
                    tau = polellip(x);
                case 'XYZHV'
                    % assume x is a 5-row vector
                    % use HV to compute the polarization. Note that the
                    % customer is responsible for matching the polarization
                    % in XYZ with HV
                    tau = polellip(x(4:5,:));
            end
                    
            % the tilt angle for the return trip is the negative of the
            % tilt angle for the forward trip, because to propagate back,
            % the horizontal component needs to be reversed. For the case
            % where the polarization gets changed by the target, then the
            % two way propagation should not be used since the assumption
            % is no longer true
            
        end
        
    end
    
    methods (Static,Hidden)
        function el = computeElevationAngle(startLoc,endLoc,propdistance)
            % Assume xy plane is ground
            el = asind(bsxfun(@minus,endLoc(3,:),startLoc(3,:))./propdistance);  
            % the return trip el is just the negative of the el of forward
            % trip
        end
        
    end
    
end