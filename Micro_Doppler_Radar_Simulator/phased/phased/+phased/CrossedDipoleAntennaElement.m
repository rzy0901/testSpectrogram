classdef (Sealed,StrictDefaults) CrossedDipoleAntennaElement < phased.internal.AbstractPolarizedAntennaElement
%CrossedDipoleAntennaElement  Crossed dipole antenna
%   H = phased.CrossedDipoleAntennaElement creates a crossed dipole antenna
%   System object, H. This object models a crossed dipole antenna.
%
%   H = phased.CrossedDipoleAntennaElement(Name,Value) creates a crossed
%   dipole antenna object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   A crossed dipole antenna is formed by two orthogonal short dipole
%   antennas, and the main beam is along the 0 degree azimuth and 0 degree
%   elevation in the antenna's local coordinate system.
%
%   The 0 degree azimuth and 0 degree elevation is considered to be the
%   main response axis of the antenna. When placed in a linear or a
%   rectangular array, the main response axis is aligned with the array
%   normal.
%
%   Step method syntax:
%
%   RESP = step(H,FREQ,ANGLE) returns the antenna voltage response, RESP,
%   given the antenna's operating frequency FREQ (in Hz) and the directions
%   specified in ANGLE (in degrees). FREQ is a row vector of length L and
%   ANGLE can be either a row vector of length M or a 2xM matrix.
%
%   RESP is a structure containing two fields, H and V. H represents the
%   antenna's response in horizontal polarization and V represents the
%   antenna's response in vertical polarization. Each field contains an MxL
%   matrix whose columns contain the responses of the antenna element in
%   the indicated polarization, at angles specified in ANGLE, and at
%   corresponding frequencies specified in FREQ.
%
%   When ANGLE is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in [azimuth; elevation] form. The azimuth angle
%   should be between [-180 180] degrees and the elevation angle should be
%   between [-90 90] degrees. If ANGLE is a length M row vector, each
%   element specifies a direction's azimuth angle and the corresponding
%   elevation angle is assumed to be 0.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   CrossedDipoleAntennaElement methods:
%
%   step                  - Output the response of the antenna element
%   release               - Allow property name and input characteristics
%                           changes
%   clone                 - Create a crossed dipole antenna object with 
%                           same property values
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the element is capable of 
%                           simulating polarization
%   directivity           - Compute element directivity
%   pattern               - Plot element response pattern
%   patternAzimuth        - Plot azimuth pattern
%   patternElevation      - Plot elevation pattern
%
%   CrossedDipoleAntennaElement properties:
%   
%   FrequencyRange        - Operating frequency range
%
%   % Examples:
%
%   % Example 1:
%   %   Plot the radiation pattern of a crossed dipole antenna.
%
%   element = phased.CrossedDipoleAntennaElement;
%   fc = 3e8;
%   pattern(element,fc,0,-90:90,'CoordinateSystem','polar');
%
%   % Example 2:
%   %   Compute the polarization response of a crossed dipole antenna
%   %   at the direction of 30 degrees azimuth and 0 degrees elevation.
%
%   element = phased.CrossedDipoleAntennaElement;
%   fc = 3e8; ang = [30;0];
%   resp = element(fc,ang)
%
%   See also phased, phased.IsotropicAntennaElement,
%   phased.ShortDipoleAntennaElement, phased.CustomAntennaElement.

%   Copyright 2012-2016 The MathWorks, Inc.

%   References:
%   [1] Harold Mott, Antennas for Radar and Communications, A Polarimetric
%   Approach, John Wiley & Sons, 1992

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties(Nontunable)

        %FrequencyRange     Operating frequency range (Hz)
        %   Specify the operating frequency range (in Hz) of the antenna
        %   element as a 1x2 row vector in the form of [LowerBound
        %   HigherBound]. The default value of this property is [0 1e20].
        %   The antenna element has no response outside the specified
        %   frequency range.
        FrequencyRange = [0 1e20]
    end
    
    methods
        function set.FrequencyRange(obj,val)
            validateattributes( val, { 'double' },...
                {'nonempty','finite','size',[ 1, 2 ],'nonnegative'},'','FrequencyRange');
            cond = val(1) > val(2);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:element:InvalidFrequencyRange');
            end
            obj.FrequencyRange = val;
        end
    end
    
    methods
        function obj = CrossedDipoleAntennaElement(varargin)
            obj@phased.internal.AbstractPolarizedAntennaElement(varargin{:});
        end
    end
    methods(Hidden)
        function cl = clonecg(obj)
            cl = phased.CrossedDipoleAntennaElement(...
                'FrequencyRange',obj.FrequencyRange);      
            if ~isPolarizationEnabled(obj)
                disablePolarization(cl);
            end

        end
    end

    methods (Access = protected)
        
        function g = getHResponse(obj,freq,ang)  %#ok<INUSL>
            g = -1i*cosd(ang(1,:).');
            g = repmat(g,1,numel(freq))*sqrt(3/2);
        end
        
        function g = getVResponse(obj,freq,ang)  %#ok<INUSL>
            g = -cosd(ang(2,:).')+...
                1i*sind(ang(2,:).').*sind(ang(1,:).');
            g = repmat(g,1,numel(freq))*sqrt(3/2);
        end
        
        function H = getFrequencyResponse(obj,freq)  %#ok<INUSL>
            H = ones(size(freq,2),1);
        end
        
        function frange = getFrequencyRange(obj)
            frange = obj.FrequencyRange;
        end
        
    end
    
    methods (Access = {?phased.internal.AbstractElement,...
            ?phased.internal.AbstractSensorOperation})
        function [ppat,pfreq,angmax] = getPowerPattern(obj,azang,elang)
            Naz = numel(azang);
            Nel = numel(elang);
            ppat = struct('H',zeros(Nel,Naz),'V',zeros(Nel,Naz));
            
            tempant = clonecg(obj);
            release(tempant);
            
            % compute power pattern
            pfreq = mean(obj.FrequencyRange);
            for m = 1:Nel;
                temp = step(tempant,pfreq,[azang;elang(m)*ones(1,Naz)]);
                ppat.H(m,:) = abs(temp.H);
                ppat.V(m,:) = abs(temp.V);
            end
            ppat.H = (ppat.H).^2;
            ppat.V = (ppat.V).^2;
            
            angmax = zeros(2,1);
        end
    end
    
    methods (Access = protected)
        
        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
    end
    methods (Hidden)
        function epat = getgpuElemResponse(obj, az, el, freq) %
            %This method is used by the array to compute the response of
            %the antenna element at all combinations of Azimuth, az, and
            %Elevation, el, angles. This is only used in a GPU simulation.
            if ~isempty(coder.target)
                coder.internal.assert(false, 'phased:element:NoGPUCodegen');
            end 
            %output size and input shape
            if isvector(az)
                az = reshape(az, 1, 1,[]);
                el = reshape(el, [], 1, 1);
                epat = gpuArray.zeros(numel(el), 1, numel(az));
            else
                epat = gpuArray.zeros(size(el));
            end
            
            
            %Note for clutter freq is always scalar so there is no repmat.
            
            %Check that the frequency is valid
            isVal = searchValidFrequency(obj,freq);
            if ~all(isVal),
                return;
            end
            
            epat = privgpuCrossedDipoleAntElemResp(az,el, epat);
        end
    end
end


% [EOF]
