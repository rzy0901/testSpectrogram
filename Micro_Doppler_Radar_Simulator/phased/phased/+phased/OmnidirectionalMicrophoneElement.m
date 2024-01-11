classdef (Sealed,StrictDefaults) OmnidirectionalMicrophoneElement < phased.internal.AbstractMicrophoneElement
%OmnidirectionalMicrophoneElement Omnidirectional microphone
%   H = phased.OmnidirectionalMicrophoneElement creates an omnidirectional
%   microphone System object, H. This object models an omnidirectional
%   microphone element whose response is one in all directions.
%
%   H = phased.OmnidirectionalMicrophoneElement(Name,Value) creates an
%   omnidirectional microphone object, H, with the specified property Name
%   set to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
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
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   OmnidirectionalMicrophoneElement methods:
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
%   OmnidirectionalMicrophoneElement properties:
%
%   FrequencyRange - Operating frequency range
%   BackBaffled    - Baffle the back of the element
%
%   % Examples:
%
%   % Example 1:
%   %   Construct an omnidirectional microphone and plot its azimuth   
%   %   response. Assume the microphone can work between 20 Hz and 2 kHz
%   %   and the operating frequency is 200 Hz.
%
%   microphone = phased.OmnidirectionalMicrophoneElement(...
%           'FrequencyRange',[20 2e3]);
%   fc = 200;
%   pattern(microphone,fc,-180:180,0,'CoordinateSystem','polar');
%
%   % Example 2:
%   %   Find the response of the above microphone at the boresight for the
%   %   operating frequencies of 200, 300, and 400 Hz.
%
%   microphone = phased.OmnidirectionalMicrophoneElement(...
%           'FrequencyRange',[20 2e3]);
%   fc = [200 300 400]; ang = [0;0];
%   resp = step(microphone,fc,ang)
%
%   See also phased, phased.CustomMicrophoneElement, phased.ULA, 
%   phased.URA, phased.ConformalArray.

%   Copyright 2010-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %FrequencyRange     Operating frequency range (Hz)
        %   Specify the operating frequency range (in Hz) of the microphone
        %   element as a 1x2 row vector in the form of [LowerBound
        %   HigherBound]. The default value of this property is [0 1e20].
        %   The microphone element has no response outside the specified
        %   frequency range.
        FrequencyRange = [0 1e20];
    end

    properties (Nontunable, Logical) 
        %BackBaffled    Baffle the back of the element
        %   Set this property to true to baffle the back of the microphone
        %   element. Set this property to false to not baffle the back of
        %   the microphone element. When the back of the microphone element
        %   is baffled, the microphone responses to all azimuth angles
        %   beyond +/- 90 degrees from the broadside (0 degree azimuth and
        %   elevation) are 0. The default value of this property is false.
        BackBaffled = false;
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
    end
    
    methods
        function obj = OmnidirectionalMicrophoneElement(varargin)
            obj@phased.internal.AbstractMicrophoneElement(varargin{:});
        end
    end
    methods
        function d = directivity(obj,freq,ang)
        %directivity  Compute element directivity
        %   D = directivity(H,FREQ,ANGLE) computes the directivity (in dBi)
        %   of the element for the directions specified in ANGLE (in
        %   degrees) and frequencies specified in FREQ (in Hz). FREQ is a
        %   row vector of length L and ANGLE can be either a row vector of
        %   length M or a 2xM matrix. D is an MxL matrix whose columns
        %   contain the directivity of the element at angles specified in
        %   ANGLE at corresponding frequencies specified in FREQ.
        %  
        %   When ANGLE is a 2xM matrix, each column of the matrix specifies
        %   the direction in space in [azimuth; elevation] form. The
        %   azimuth angle should be between [-180 180] degrees and the
        %   elevation angle should be between [-90 90] degrees. If ANGLE is
        %   a length M row vector, each element specifies a direction's
        %   azimuth angle and the corresponding elevation angle is assumed
        %   to be 0.
        %
        %   % Examples:
        %
        %   % Example 1:
        %   %   Compute the directivity of an omnidirectional microphone 
        %   %   at 300 Hz toward the boresight.
        %
        %   myMic = phased.OmnidirectionalMicrophoneElement;
        %   d = directivity(myMic,300,0)
        %
        %   % Example 2:
        %   %   Compute the directivity of a back baffled omnidirectional 
        %   %   microphone at 10 kHz toward 30 degrees azimuth and 10 
        %   %   degrees elevation.
        %
        %   myMic = phased.OmnidirectionalMicrophoneElement(...
        %           'BackBaffled',true);
        %   d = directivity(myMic,10e3,[30;10])
        %
        %   See also phased, phased.ArrayResponse.
        
            phased.internal.narginchk(3,3,nargin);
            
            sigdatatypes.validateFrequency(freq,'directivity','FREQ',...
                {'row'});
            
            sigdatatypes.validateAngle(ang,'directivity','ANGLE');
            if isrow(ang)
                angIn = [ang;zeros(size(ang))];
            else
                angIn = ang;
            end
            sigdatatypes.validateAzElAngle(angIn,'directivity','ANGLE');
            
            M = size(angIn,2);
            L = numel(freq);
            
            if obj.BackBaffled
                d = 2*ones(M,L);
                idx = phased.internal.isAngleAtBack(angIn);
                d(idx,:) = 0;
            else
                d = ones(M,L);
            end

            d = pow2db(d);
        end
    end

    methods (Hidden)
        function cl = clonecg(obj)
            cl = phased.OmnidirectionalMicrophoneElement(...
                'FrequencyRange',obj.FrequencyRange, ...
                'BackBaffled',obj.BackBaffled);      
        end
    end
    methods (Access = protected)
        function g = getSpatialResponse(obj,freq,ang) 
            g = ones(size(ang,2),1);
            if obj.BackBaffled
                g = phased.internal.backbaffle(g,ang);
            end
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
            ppat = ones(Nel,Naz);
            if obj.BackBaffled
                ppat(bsxfun(@and,abs(elang(:)-90)>sqrt(eps),...
                    (azang<(-90-sqrt(eps))) | (azang>(90+sqrt(eps))))) = 0;
            end
            angmax = zeros(2,1);
            pfreq = mean(obj.FrequencyRange);
        end
        
        function flag = isDirectivityKnown(obj) %#ok<MANU>
            flag = true;
        end
    end
end


