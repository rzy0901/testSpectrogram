classdef (Abstract) AbstractIsotropicSonarElement < phased.internal.AbstractElement
%This class is for internal use only. It may be removed in the future.

%AbstractIsotropicSonarElement   Define the AbstractIsotropicSonarElement class.

%   Copyright 2016 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %FrequencyRange     Operating frequency range (Hz)
        %   Specify the operating frequency range (in Hz) of the
        %   element as a 1x2 row vector in the form of [LowerBound
        %   HigherBound]. The default value of this property is [0 100e6].
        %   The sonar element has no response outside the specified
        %   frequency range.
        FrequencyRange = [0 100e6];
    end
    
    properties (Nontunable, Logical) 
        %BackBaffled    Baffle the back of the element
        %   Set this property to true to baffle the back of the element.
        %   Set this property to false to not baffle the back of the
        %   hydrophone element. When the back of the hydrophone element is
        %   baffled, the hydrophone responses to all azimuth angles beyond
        %   +/- 90 degrees from the broadside (0 degree azimuth and
        %   elevation) are 0. The default value of this property is false.
        BackBaffled = false;
    end
   
    properties (Access='protected', Nontunable)
        pFrequencyVector
        pResponse
    end
   
    methods
        function obj=AbstractIsotropicSonarElement(varargin)
           obj=obj@phased.internal.AbstractElement(varargin{:});
        end
    end
    methods (Access = protected)
        function y = stepImpl(obj,FREQ,ANG)
            y = stepImpl@phased.internal.AbstractElement(obj,FREQ,ANG);
            y = y.*repmat(getResponse(obj,FREQ),size(ANG,2),1); 
        end
    end
    
    methods
        function set.FrequencyRange(obj,val)
            validateattributes( val, { 'double' }, ...
                {'finite','nonnegative','nondecreasing','row','numel',2},'','FrequencyRange');
            cond = val(1) > val(2);
            if cond
                coder.internal.errorIf(cond,...
                        'phased:phased:element:InvalidFrequencyRange');
            end
            obj.FrequencyRange = val;
        end
    end

    methods
        function d = directivity(obj,freq,ang)
        %directivity  Compute element directivity index
        %   D = directivity(H,FREQ,ANGLE) computes the directivity index
        %   (in dBi) of the element for the directions specified in ANGLE
        %   (in degrees) and frequencies specified in FREQ (in Hz). FREQ is
        %   a row vector of length L and ANGLE can be either a row vector
        %   of length M or a 2xM matrix. D is an MxL matrix whose columns
        %   contain the directivity index of the element at angles
        %   specified in ANGLE at corresponding frequencies specified in
        %   FREQ.
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
        %   %   Compute the directivity index of an isotropic hydrophone 
        %   %   at 30 kHz toward the boresight.
        %
        %   myHydro = phased.IsotropicHydrophone;
        %   d = directivity(myHydro,30e3,0)
        %
        %   % Example 2:
        %   %   Compute the directivity index of a back baffled isotropic
        %   %   projector at 1 kHz toward 30 degrees azimuth and 10 degrees 
        %   %   elevation.
        %
        %   myProj = phased.IsotropicProjector('BackBaffled',true);
        %   d = directivity(myProj,1e3,[30;10])
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
    
    methods (Access = protected)
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractElement(obj);
            if isLocked(obj)
                s.pFrequencyVector = obj.pFrequencyVector;
                s.pResponse = obj.pResponse;
            end
        end
        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
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

        function H = getResponse(obj,FREQ) 
            H = interp1(obj.pFrequencyVector,obj.pResponse,FREQ,'linear','extrap');
            H(FREQ<obj.FrequencyRange(1) | FREQ>obj.FrequencyRange(2)) = 0;
        end

        function H = getFrequencyResponse(obj,freq)   %#ok<INUSL>
            H = ones(size(freq,2),1);
        end

        function frange = getFrequencyRange(obj)
            frange = obj.FrequencyRange;
        end
    end
    
    methods ( Access = {?phased.internal.AbstractElement,...
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

            
