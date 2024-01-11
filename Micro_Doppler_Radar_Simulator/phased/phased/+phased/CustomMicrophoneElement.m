classdef (Sealed,StrictDefaults) CustomMicrophoneElement < phased.internal.AbstractMicrophoneElement
%CustomMicrophoneElement Custom microphone
%   H = phased.CustomMicrophoneElement creates a custom microphone system
%   object, H. This object models a custom microphone element.
%
%   H = phased.CustomMicrophoneElement(Name,Value) creates a custom
%   microphone object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The 0 degree azimuth and 0 degree elevation is considered to be the
%   main response axis of the microphone. When placed in a linear or a
%   rectangular array, the main response axis is aligned with the array
%   normal.
%
%   Step method syntax:
%
%   RESP = step(H,FREQ,ANGLE) returns the microphone's magnitude response,
%   RESP, at frequencies specified in FREQ (in Hz) and directions specified
%   in ANGLE (in degrees). FREQ is a row vector of length L and ANGLE can
%   be either a row vector of length M or a 2xM matrix. RESP is an MxL
%   matrix whose columns contain the responses of the microphone element at
%   the M angles specified in ANGLE and the L frequencies specified in
%   FREQ.
%
%   When ANGLE is a 2xM matrix, each column of the matrix specifies the
%   direction in [azimuth; elevation] form. The azimuth angle should be
%   between [-180 180] degrees and the elevation angle should be between
%   [-90 90] degrees. If ANGLE is a length M row vector, each element
%   specifies a direction's azimuth angle and the corresponding elevation
%   angle is assumed to be 0.
%
%   The total response of a custom microphone element is a combination of
%   its frequency response and spatial response. Both responses are
%   calculated using nearest neighbor interpolation and then multiplied to
%   form the total response. When the PolarPatternFrequencies property
%   value is nonscalar, the object specifies multiple polar patterns. In
%   this case, the interpolation uses the polar pattern closest to the
%   specified frequency.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   CustomMicrophoneElement methods:
%
%   step                  - Output the response of the microphone
%   release               - Allow property name and input characteristics 
%                           changes
%   clone                 - Create an omnidirectional microphone object 
%                           with same property values
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the element is capable of 
%                           simulating polarization
%   directivity           - Compute element directivity
%   pattern               - Plot element response pattern
%   patternAzimuth        - Plot azimuth pattern
%   patternElevation      - Plot elevation pattern
%
%   CustomMicrophoneElement properties:
%
%   FrequencyVector         - Operating frequency vector
%   FrequencyResponse       - Frequency responses
%   PolarPatternFrequencies - Polar pattern frequencies
%   PolarPatternAngles      - Polar pattern angles
%   PolarPattern            - Polar pattern
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a custom Cardioid microphone and plot its azimuth   
%   %   response. Assume the operating frequency is 500 Hz.
%
%   element = phased.CustomMicrophoneElement;
%   element.PolarPatternFrequencies = [500 1000];
%   element.PolarPattern = mag2db([...
%           0.5+0.5*cosd(element.PolarPatternAngles);...
%           0.6+0.4*cosd(element.PolarPatternAngles)]);
%   fc = 500;
%   pattern(element,fc,-180:180,0,'CoordinateSystem','polar');
%
%   % Example 2:
%   %   Find the response of the above microphone in the directions of 
%   %   [0;0] and [40;50] for the operating frequencies of 500, 1500, and 
%   %   2000 Hz.
%
%   element = phased.CustomMicrophoneElement;
%   element.PolarPatternFrequencies = [500 1000];
%   element.PolarPattern = mag2db([...
%           0.5+0.5*cosd(element.PolarPatternAngles);...
%           0.6+0.4*cosd(element.PolarPatternAngles)]);
%   fc = [500 1500 2000]; ang = [0 0;40 50]';
%   resp = element(fc,ang)
%
%   See also phased, phased.OmnidirectionalMicrophoneElement, phased.ULA,
%   phased.URA, phased.ConformalArray.

%   Copyright 2010-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %FrequencyVector    Operating frequency vector (Hz)
        %   Specify the frequencies (in Hz) where the frequency responses
        %   of element are measured as a vector. The elements of the vector
        %   must be increasing. The default of this property is [0 1e20].
        %   The microphone element has no response outside the specified
        %   frequency range.
        FrequencyVector = [0 1e20];
        %FrequencyResponse  Frequency responses (dB)
        %   Specify the frequency responses (in dB) measured at the
        %   frequencies defined in the FrequencyVector property as a row
        %   vector. The length of the vector must equal the length of the
        %   frequency vector specified in the FrequencyVector property. The
        %   default value of this property is [0 0].
        FrequencyResponse = [0 0];
        %PolarPatternFrequencies    Polar pattern frequencies (Hz)
        %   Specify the measuring frequencies (in Hz) of the polar patterns
        %   as a length-M row vector. The measuring frequencies must be
        %   within the frequency range specified in the FrequencyVector
        %   property. The default value of this property is 1e3.
        PolarPatternFrequencies = 1e3;
        %PolarPatternAngles Polar pattern angles (deg)
        %   Specify the measuring angles (in degrees) of the polar patterns
        %   as a length-N row vector. The angles are measured from the
        %   central pickup axis of the microphone, and must be within [-180
        %   180]. The default value of this property is -180:180.
        PolarPatternAngles = -180:180;  
    end
    
    properties (Nontunable, Dependent)
        %PolarPattern   Polar pattern (dB)
        %   Specify the polar patterns of the microphone element as an MxN
        %   matrix. M is the number of measuring frequencies specified in
        %   the PolarPatternFrequencies property. N is the number of
        %   measuring angles specified in the PolarPatternAngles property.
        %   Each row of the matrix represents the magnitude of the polar
        %   pattern (in dB) measured at the corresponding frequency
        %   specified in the PolarPatternFrequencies property and
        %   corresponding angles specified in the PolarPatternAngles
        %   property.
        %
        %   The pattern is assumed to be measured in the azimuth plane
        %   where the elevation angle is 0 and where the central pick up
        %   axis is assume to be 0 degrees azimuth and 0 degrees elevation.
        %   The polar pattern is assumed to be symmetric around the central
        %   axis and therefore the microphone's response pattern in 3D
        %   space can be constructed from the polar pattern. The default
        %   value of this property is an omnidirectional pattern with 0 dB
        %   response everywhere.
        PolarPattern 
    end
    
    properties (Access = private, Nontunable)
        pPolarPattern = ones(1,361)
    end

    methods
        function set.FrequencyVector(obj,value)
            validateattributes( value, { 'double' }, { 'nonempty', ...
                'finite', 'row', 'nonnegative' }, '', 'FrequencyVector');
            cond =  length(value) < 2;
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:NotEnoughSamples', 'FrequencyVector');
            end
            cond =  any(diff(value)<0);
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:NotIncreasing');
            end
            obj.FrequencyVector = value;
        end
        
        function set.FrequencyResponse(obj,value)
            validateattributes(value,{'double'},{'real','row','nonnan'},...
                'phased.CustomMicrophoneElement','FrequencyResponse');
            cond =  length(value) < 2;
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:NotEnoughSamples', 'FrequencyResponse');
            end
            obj.FrequencyResponse = value;
        end
            
        function set.PolarPattern(obj,value)
            validateattributes( value, { 'double' }, { '2d', 'real', 'nonempty', 'nonnan' }, '', 'PolarPattern');
            obj.pPolarPattern = phased.internal.dbtomag(value);
        end
        
        function value = get.PolarPattern(obj)
            value = phased.internal.magtodb(obj.pPolarPattern);
        end

        function set.PolarPatternAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.CustomMicrophoneElement','PolarPatternAngles',...
                {'row','>=',-180,'<=',180});
            cond =  length(value) < 2;
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:NotEnoughSamples', 'PolarPatternAngles');
            end
            cond =  any(diff(value) <= 0);
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:phased:CustomMicrophoneElement:InvalidPolarPatternAngles');
            end
            obj.PolarPatternAngles = value;
        end

        function set.PolarPatternFrequencies(obj,value)
            validateattributes( value, { 'double' }, { 'nonempty', 'row', 'nonnegative', 'finite' }, '', 'PolarPatternFrequencies');
            obj.PolarPatternFrequencies = value;
        end
    end
    
    methods

        function obj = CustomMicrophoneElement(varargin)
            obj@phased.internal.AbstractMicrophoneElement(varargin{:});
        end

    end
    methods(Hidden)
        function cl = clonecg(obj)
            cl = phased.CustomMicrophoneElement(...
                'FrequencyVector',obj.FrequencyVector, ...
                'FrequencyResponse',obj.FrequencyResponse, ...
                'PolarPatternFrequencies',obj.PolarPatternFrequencies, ...
                'PolarPatternAngles',obj.PolarPatternAngles, ...
                'PolarPattern',obj.PolarPattern);
        end
    end

    methods(Access = protected)
        function validatePropertiesImpl(obj)
            cond =  numel(obj.FrequencyVector) ~= numel(obj.FrequencyResponse);
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:FrequencyResponseMismatch');
            end
            same_freq = find(diff(obj.FrequencyVector)==0);
            cond =  ...
                any(obj.FrequencyResponse(same_freq) ~= obj.FrequencyResponse(same_freq+1));
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:InvalidFrequencyResponse');
            end
            frange = getFrequencyRange(obj);
            pfreq = obj.PolarPatternFrequencies;
            cond =  any((pfreq<frange(1))|(pfreq>frange(2)));
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:FrequencyOutOfBound');
            end
            M = numel(pfreq);
            N = numel(obj.PolarPatternAngles);
            [pM, pN] = size(obj.PolarPattern);
            cond =  ((pM~=M)||(pN~=N));
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:element:PatternMismatch', 'PolarPattern', M, N);
            end
            
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractMicrophoneElement(obj);
            % starting from R2014a, we save linear scale pattern in
            % pPolarPattern to improve the performance, avoiding doing
            % db conversion all the time. 
            s.pPolarPattern = obj.pPolarPattern;
        end
        
        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
    end

    methods (Access = protected)
        function g = getSpatialResponse(obj,freq,angle)
            pfreq = obj.PolarPatternFrequencies;
            ppattern = obj.pPolarPattern;
            freq_idx = getClosestElementIndex(pfreq.',freq);
            axisangle = acosd(cosd(angle(1,:)).*cosd(angle(2,:))).';
            [unique_freq_idx, ~, f2p_idx] = unique(freq_idx);
            interp_ppattern = interp1(obj.PolarPatternAngles,...
                ppattern(unique_freq_idx,:).',...
                axisangle,'nearest',0);
            g = interp_ppattern(:,f2p_idx);
                
        end
        
        function H = getFrequencyResponse(obj,freq)
            
            H = interp1(obj.FrequencyVector,...
                phased.internal.dbtomag(obj.FrequencyResponse),...
                freq.','nearest',0);
        end
        
        function frange = getFrequencyRange(obj)
            frange = obj.FrequencyVector(1,[1 end]);
        end
    end

    methods (Access = {?phased.internal.AbstractElement,...
            ?phased.internal.AbstractSensorOperation})
        function [ppat,pfreq,angmax] = getPowerPattern(obj,azang,elang)
            Naz = numel(azang);
            Nel = numel(elang);
            
            % compute power pattern
            ppattern = obj.pPolarPattern;
            ppatternang = obj.PolarPatternAngles;
            
            freq = union(obj.FrequencyVector,obj.PolarPatternFrequencies);
            
            % compute distinct spatial patterns
            Nfreq = size(ppattern,1);
            if Nfreq == 1
                ppats = zeros(Nel,Naz);
                for m = 1:Nel
                    axisangle = acosd(cosd(azang).*cosd(elang(m))).';
                    ppats(m,:) = interp1(ppatternang,...
                        ppattern.',axisangle,'nearest',0);
                end
            else
                ppats = zeros(Nel,Naz,Nfreq);
                for m = 1:Nel
                    axisangle = acosd(cosd(azang).*cosd(elang(m))).';
                    for n = 1:Nfreq
                        ppats(m,:,n) = interp1(ppatternang,...
                            ppattern(n,:).',axisangle,'nearest',0);
                    end
                end
            end
            ppats = abs(ppats).^2;
            
            patc = ppats;
            [rmax,idxc] = max(patc,[],2);
            rmax = squeeze(rmax);
            idxc = squeeze(idxc);
            [~,idxr] = max(rmax,[],1);
            angmaxs = zeros(2,Nfreq);
            for m = 1:Nfreq
                angmaxs(2,m) = elang(idxr(m));
                angmaxs(1,m) = azang(idxc(idxr(m),m));
            end
            
            % compute all frequencies
            
            fr = db2pow(interp1(obj.FrequencyVector,obj.FrequencyResponse,...
                freq,'nearest'));
            
            if Nfreq == 1
                ppat = bsxfun(@times,ppats,permute(fr,[1 3 2]));
                angmax = repmat(angmaxs,1,numel(freq));
            else
                sridx = interp1(obj.PolarPatternFrequencies,1:Nfreq,...
                    freq,'nearest','extrap');
                ppat = bsxfun(@times,ppats(:,:,sridx),permute(fr,[1 3 2]));
                angmax = angmaxs(:,sridx);
            end
            
            pfreq = freq;
                
        end
    end  
    
    methods (Static, Hidden, Access = protected)        
        function groups = getPropertyGroupsImpl
            p1 = matlab.system.display.internal.Property(...
                'PolarPattern','Description', 'Polar pattern (dB)', ...
                'UseClassDefault',false,'Default', 'zeros(1,361)');
            groups = matlab.system.display.Section(...
                            'PropertyList',...
                            {'FrequencyVector',...
                            'FrequencyResponse',...
                            'PolarPatternFrequencies',...
                            'PolarPatternAngles',...
                            p1},...
                            'DependOnPrivatePropertyList',...
                            {'PolarPattern'});   
                       
        end
    end
end

function idx = getClosestElementIndex(a,b)
%getNearestIndex    Find the index of closest element
%   IDX = getClosestElementIndex(A,B) returns the indices IDX for the
%   elements of A who is closest to the elements specified in B. A is a
%   column vector and B is a row vector. IDX is a row vector whose length
%   matches the length of B.
%
%   % Example:
%   %   Find indices for elements in [3 5 12 6 9] that are closest to
%   %   1 and 7.
%   idx = getClosestElementIndex([3 5 12 6 9],[1 7]);

[~,idx] = min(abs(ones(numel(a),1)*b-a*ones(1,numel(b))),[],1);
end

