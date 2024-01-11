classdef (Hidden, Sealed, StrictDefaults) AntennaAdapter < phased.internal.AbstractPolarizedAntennaElement
%This class is for internal use only. It may be removed in the future.

%AntennaAdapter   Define the antenna adapter for Antenna Toolbox elements
%   H = phased.internal.AntennaAdapter creates an antenna adapter System
%   object, H. This object serves as an adapter to use elements created in
%   Antenna Toolbox in Phased Array System Toolbox.
%
%   H = phased.TwoRayChannel(Name,Value) returns an antenna adapter object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
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
%   IsotropicAntennaElement methods:
%
%   step                  - Output the response of the antenna element
%   release               - Allow property name and input characteristics
%                           changes
%   clone                 - Create an antenna adapter object with same 
%                           property values
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the element is capable of 
%                           simulating polarization
%   directivity           - Compute element directivity
%
%   IsotropicAntennaElement properties:
%
%   Antenna - Handle to the antenna from Antenna Toolbox
%
%   % Example:
%   %   Construct an adapter for dipole and find its response at boresight.
%   %   Assume the operating frequency is 1 GHz.
%
%   ha = phased.internal.AntennaAdapter('Antenna',dipole);
%   fc = 1e9; ang = [0;0];
%   resp = step(ha,fc,ang)
%
%   See also phased, phased.CosineAntennaElement,
%   phased.CustomAntennaElement, phased.ULA, phased.URA,
%   phased.ConformalArray.

%   Copyright 2014-2016 The MathWorks, Inc.

    properties (Nontunable)
        %Antenna    Handle to antenna
        %   A scalar handle to antennas from Antenna Toolbox.
        Antenna
    end
    
    properties (Access=private)
        pComputedFreq
        pComputedAng
        pComputedField
        pMaxResponse
        pAzimuthAngle = -180:5:180;
        pElevationAngle = -90:5:90;
    end
    
    methods
        function obj = AntennaAdapter(varargin)
            obj@phased.internal.AbstractPolarizedAntennaElement(varargin{:});
        end
    end
    
    methods (Access=protected)
        
        function resp = getHResponse(obj,freq,angle)
            updateComputedField(obj,freq);
            resp = interpolateResponse(obj,obj.pComputedField.H,angle);
        end
        
        function resp = getVResponse(obj,freq,angle)
            updateComputedField(obj,freq);
            resp = interpolateResponse(obj,obj.pComputedField.V,angle);           
        end
        
        function frange = getFrequencyRange(obj) %#ok<MANU>
            frange = [0 1e20];
        end
        
        function resp = getFrequencyResponse(obj,freq)  %#ok<INUSL>
            resp = ones(size(freq,2),1);
        end
        
        function d = getDirectivity(obj,freq,ang)
            d = zeros(size(ang,2),numel(freq));
            for m = 1:numel(freq)
                d_temp = pattern(obj.Antenna,freq(m),ang(1,:),ang(2,:));
                d(:,m) = diag(d_temp);
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractPolarizedAntennaElement(obj);
            s.Antenna = saveobj(obj.Antenna);
            if isLocked(obj)
                s.pComputedFreq = obj.pComputedFreq;
                s.pComputedAng = obj.pComputedAng;
                s.pComputedField = obj.pComputedField;
                s.pMaxResponse = obj.pMaxResponse;
                s.pAzimuthAngle = obj.pAzimuthAngle;
                s.pElevationAngle = obj.pElevationAngle;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) 
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked) %#ok<INUSD>
            obj.Antenna = em.EmStructures.loadobj(s.Antenna);
            s = rmfield(s,'Antenna');
        end
        
        function resp = stepImpl(obj,freq,ang)
            resp = stepImpl@phased.internal.AbstractPolarizedAntennaElement(obj,freq,ang);
            maxresp = getMaxResponse(obj,freq); % row
            if isstruct(resp)
                resp.H = bsxfun(@rdivide,resp.H,maxresp);
                resp.V = bsxfun(@rdivide,resp.V,maxresp);
            else
                resp = bsxfun(@rdivide,resp,maxresp);
            end
                
        end
        
   end
    
    methods (Access = {?phased.internal.AbstractElement,...
            ?phased.internal.AbstractArray})
        %flag to make sure whether the array normal and the element normal
        %are aligned or not. In all element cases the normals are aligned,
        %For the short dipole case the element normal and element normal
        %are not aligned.
        function isAligned = isElementNormalArrayNormalAligned(obj)  %#ok<MANU>
            isAligned = false;
        end
    end
    
    methods (Access = {?phased.internal.AbstractElement,...
            ?phased.internal.AbstractSensorOperation})
        
        function [ppat,pfreq,angmax] = getPowerPattern(obj,azang,elang) %#ok<INUSD,STOUT>

            % no-op, the directivity computation will resort to Antenna
            % Toolbox if needed
        end
    end
    
    methods (Hidden)
            
        function epat = getgpuElemResponse(obj, az, el, freq)
            %This method is used by the phased.gpu.ConstantGammaClutter to
            %compute the response of the antenna element at all
            %combinations of Azimuth, az, and Elevation, el, angles.
            %Note that az and el are in radians.
            
            if ~isempty(coder.target)
               coder.internal.assert(false, 'phased:element:NoGPUCodegen');
            end  
            updateComputedField(obj,freq);
            %output size and input shape
            if isvector(az)
                epat_H = zeros(numel(el), 1, numel(az));
                epat_V = zeros(numel(el), 1, numel(az));
                el = el(:);
                az = az(:).';
            else
                epat_H = zeros(size(el));
                epat_V = zeros(size(el));
            end
            
            %Note for clutter freq is always scalar so there is no repmat.
            
            %For elements from Antenna Toolbox, frequency is always valid
            %and frequency response is always 1.
            
            az = phased.internal.rad2deg(az);
            el = phased.internal.rad2deg(el);
            
            for m = 1:size(el,1)
                if isvector(az)
                    ang = [az;el(m)*ones(size(az))];
                else
                    ang = [az(m,:);el(m,:)];
                end
                epat_H(m,:) = interpolateResponse(obj,obj.pComputedField.H,gather(ang));
                epat_V(m,:) = interpolateResponse(obj,obj.pComputedField.V,gather(ang));
            end
            epat_max = getMaxResponse(obj,freq); % row
            epat_H = bsxfun(@rdivide,epat_H,epat_max);
            epat_V = bsxfun(@rdivide,epat_V,epat_max);
            
            epat = gpuArray(hypot(epat_H,epat_V));
        end
        
    end
    
    methods (Access=private)
        function pat_max = getMaxResponse(obj,freq)
            if isempty(obj.pMaxResponse) || ~isequal(obj.pComputedFreq,freq)
                updateComputedField(obj,freq);
                pat_max = obj.pMaxResponse;
            else
                pat_max = obj.pMaxResponse;
            end
        end
        
        function updateComputedField(obj,freq)
            if isempty(obj.pComputedFreq) || ~isequal(obj.pComputedFreq,freq) 
                az = obj.pAzimuthAngle;
                el = obj.pElevationAngle;
                hpat = complex(zeros(numel(el),numel(az)),numel(freq));
                vpat = complex(zeros(numel(el),numel(az)),numel(freq));
                elvecbase = ones(size(az));
                for f = 1:numel(freq)
                    for m = 1:numel(el)
                        ang = [az;el(m)*elvecbase];
                        [pat_H,pat_V] = computeHVField(obj,freq(f),ang);
                        hpat(m,:,f) = pat_H;
                        vpat(m,:,f) = pat_V;
                    end
                end
                pat_abs = hypot(hpat,vpat);
                pat_max_temp = max(max(pat_abs,[],1),[],2);
                pat_max = pat_max_temp(:).';  % change to row
                obj.pMaxResponse = pat_max;
                obj.pComputedFreq = freq;
                obj.pComputedField.H = hpat;
                obj.pComputedField.V = vpat;
            end
        end
        
        function [pat_H,pat_V] = computeHVField(obj,freq,angle)
            pat = step(obj.Antenna,freq,angle);
            [hvec,vvec] = phased.internal.azel2vec(angle);
            num_ang = size(angle,2);
            pat_H = reshape(sum(bsxfun(@times,pat,hvec)),num_ang,[]);
            pat_V = reshape(sum(bsxfun(@times,pat,vvec)),num_ang,[]);
        end
        
        function resp = interpolateResponse(obj,pat,angle)
            numfreq = size(pat,3);
            resp = zeros(size(angle,2),numfreq,'like',1+1i);
            az = obj.pAzimuthAngle;
            el = obj.pElevationAngle;
            for m = 1:numfreq
                resp(:,m) = interpolatePattern(az,el,...
                    pat(:,:,m),angle(1,:),angle(2,:));
            end
        end
    end
        
end

function pat = interpolatePattern(az, el, pattern, az_q, el_q)
    pat = interp2(az,el,pattern,az_q, el_q,'nearest',0);
end

