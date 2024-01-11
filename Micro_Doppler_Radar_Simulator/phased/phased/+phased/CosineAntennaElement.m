classdef (Sealed,StrictDefaults) CosineAntennaElement < phased.internal.AbstractAntennaElement
%CosineAntennaElement Cosine antenna
%   H = phased.CosineAntennaElement creates a cosine antenna System object,
%   H. This object models an antenna element whose response is cosine
%   raised to a certain power in both azimuth and elevation direction.
%
%   H = phased.CosineAntennaElement(Name,Value) creates a cosine antenna
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   The cosine pattern is given by
%       P = cos(theta_az)^pow_az*cos(theta_el)^pow_el
%   where theta_az and theta_el are azimuth and elevation angles,
%   respectively. The exponent pow_az and pow_el are real non-negative
%   numbers. Note that a cosine antenna has no response at the back.
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
%   ANGLE can be either a row vector of length M or a 2xM matrix. RESP is
%   an MxL matrix whose columns contain the responses of the antenna
%   element at angles specified in ANGLE at corresponding frequencies
%   specified in FREQ.
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
%   CosineAntennaElement methods:
%
%   step                  - Output the response of the antenna element
%   release               - Allow property name and input characteristics
%                           changes
%   clone                 - Create an isotropic antenna object with same 
%                           property values
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the element is capable of 
%                           simulating polarization
%   directivity           - Compute element directivity
%   pattern               - Plot element response pattern
%   patternAzimuth        - Plot azimuth pattern
%   patternElevation      - Plot elevation pattern
%
%   CosineAntennaElement properties:
%
%   FrequencyRange - Operating frequency range
%   CosinePower    - Exponent of cosine pattern
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a cosine pattern antenna and plot its elevation response.
%   %   Assume the antenna can work between 800 MHz and 1.2 GHz and the
%   %   operating frequency is 1 GHz.
%
%   element = phased.CosineAntennaElement('FrequencyRange',[800e6 1.2e9]);
%   fc = 1e9;
%   pattern(element,fc,-180:180,0,'CoordinateSystem','polar');
%
%   % Example 2:
%   %   Find the response of the above antenna at the boresight.
%
%   element = phased.CosineAntennaElement('FrequencyRange',[800e6 1.2e9]);
%   fc = 1e9; ang = [0;0];
%   resp = element(fc,ang)
%
%   See also phased, phased.IsotropicAntennaElement,
%   phased.CustomAntennaElement, phased.ULA, phased.URA,
%   phased.ConformalArray.

%   Copyright 2008-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %FrequencyRange     Operating frequency range (Hz)
        %   Specify the operating frequency range (in Hz) of the antenna
        %   element as a 1x2 row vector in the form of [LowerBound
        %   HigherBound]. The default value of this property is [0 1e20].
        %   The antenna element has no response outside the specified
        %   frequency range.
        FrequencyRange = [0 1e20];
    end
    
    properties(Dependent, Nontunable)
        %CosinePower    Exponent of cosine pattern
        %   Specify the exponent of cosine pattern as a scalar or a 1x2
        %   vector. All specified values must be real non-negative numbers.
        %   When you set CosinePower to a scalar, both azimuth direction
        %   cosine pattern and elevation direction cosine pattern are
        %   raised to the specified value. When you set CosinePower to a
        %   1x2 vector, the first element is the exponent for the azimuth
        %   direction cosine pattern and the second element is the exponent
        %   for the elevation direction cosine pattern. The default value
        %   of this property is [1.5 1.5].
        CosinePower 
    end
    
    properties(Access=private, Nontunable)
        pAzCosinePower = 1.5;
        pElCosinePower = 1.5;
    end
        
    
    methods
        function set.FrequencyRange(obj,val)
            validateattributes( val, { 'double' }, { 'nonempty', 'finite',...
                'size', [ 1, 2 ], 'nonnegative' }, '', 'FrequencyRange');
            cond = val(1) > val(2);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:phased:element:InvalidFrequencyRange');
            end
            obj.FrequencyRange = val;
        end
        function set.CosinePower(obj,valArg)
            if isscalar(valArg)
                val = [valArg valArg];
            else
                val = valArg;
            end
            validateattributes( val, { 'double' }, { 'nonempty', 'finite',...
                'size', [ 1, 2 ], 'nonnegative' }, '', 'CosinePower');
            obj.pAzCosinePower = val(1);
            obj.pElCosinePower = val(2);
        end
        function val = get.CosinePower(obj)
            val = [obj.pAzCosinePower obj.pElCosinePower];
        end
    end
    
    methods
        function obj = CosineAntennaElement(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
    methods(Hidden)
        function clCos = clonecg(obj)
            clCos = phased.CosineAntennaElement(...
                'FrequencyRange',obj.FrequencyRange, ...
                'CosinePower',obj.CosinePower);      
        end
    end

    methods (Access = protected)
        function g = getSpatialResponse(obj,freq,ang) 
            g = zeros(size(ang,2),1);
            az_ang = ang(1,:);
            el_ang = ang(2,:);
            validAngleIdx = ((az_ang>=-90 & (az_ang<=90)));
            g(validAngleIdx) = cosd(az_ang(validAngleIdx)).^obj.pAzCosinePower.*...
                cosd(el_ang(validAngleIdx)).^obj.pElCosinePower;
            g = g*ones(1,numel(freq));
        end
        
        function H = getFrequencyResponse(obj,freq) %#ok
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
            ppat = zeros(Nel,Naz);
            
            tempant = clonecg(obj);
            release(tempant);
            
            % compute power pattern
            pfreq = mean(obj.FrequencyRange);
            for m = 1:Nel
                ppat(m,:) = step(tempant,pfreq,[azang;elang(m)*ones(1,Naz)]);
            end
            ppat = abs(ppat).^2;
            
            angmax = zeros(2,1);
        end
    end
    
    methods (Access = protected)
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractAntennaElement(obj);
            s.pAzCosinePower = obj.pAzCosinePower;
            s.pElCosinePower = obj.pElCosinePower;
        end
        
        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
    end
    
    methods (Hidden)
        function epat= getgpuElemResponse(obj, az, el, freq)
            %This method is used by the array to compute the response of
            %the antenna element at all combinations of Azimuth, az, and
            %Elevation, el, angles. This is only used in a GPU simulation.
            %Note that az and el are in radians.
            if(isempty(coder.target))
                if searchValidFrequency(obj,freq)
                    epat = privgpuCosAntElemResp(az, el, ...
                        obj.pAzCosinePower, obj.pElCosinePower);
                else
                    if isvector(az)
                        epat = gpuArray.zeros(numel(el), 1, numel(az));
                    else
                        epat = gpuArray.zeros(size(el));
                    end
                end
            else
                coder.internal.assert(false, 'phased:element:NoGPUCodegen');
            end
        end
    end

    methods (Static, Hidden, Access = protected)        
        function groups = getPropertyGroupsImpl
            groups = matlab.system.display.Section('phased.CosineAntennaElement', ...
                            'DependOnPrivatePropertyList',...
                            {'CosinePower'});            
        end
    end
end

