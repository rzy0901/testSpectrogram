classdef (Sealed,StrictDefaults) CustomAntennaElement < phased.internal.AbstractPolarizedAntennaElement
%CustomAntennaElement Custom antenna
%   H = phased.CustomAntennaElement creates a custom antenna System object,
%   H. This object models an antenna element with a custom response
%   pattern. The default custom antenna element has an isotropic response
%   in space.
%
%   H = phased.CustomAntennaElement(Name,Value) creates a custom antenna
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   The custom antenna object uses nearest neighbor interpolation to
%   estimate the response of the antenna at a given direction. Therefore,
%   to avoid interpolation errors, the custom response pattern should cover
%   azimuth angles for [-180 180] degrees and elevation angles for [-90 90]
%   degrees.
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
%   If you set the SpecifyPolarizationPattern property to false, RESP is an
%   MxL matrix whose columns contain the responses of the antenna element
%   at angles specified in ANGLE, and at corresponding frequencies
%   specified in FREQ. If you set the SpecifyPolarizationPattern property
%   to true, RESP is a structure containing two fields, H and V. H
%   represents the antenna's response in horizontal polarization and V
%   represents the antenna's response in vertical polarization. Each field
%   contains an MxL matrix whose columns contain the responses of the
%   antenna element in the indicated polarization, at angles specified in
%   ANGLE, and at corresponding frequencies specified in FREQ.
%
%   When ANGLE is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in [azimuth; elevation] form. The azimuth angle
%   should be between [-180 180] degrees and the elevation angle should be
%   between [-90 90] degrees. If ANGLE is a length M row vector, each
%   element specifies a direction's azimuth angle and the corresponding
%   elevation angle is assumed to be 0.
%
%   The total response of a custom antenna element is a combination of its
%   frequency response and spatial response. Both responses are calculated
%   using nearest neighbor interpolation and then multiplied together to
%   form the total response.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   CustomAntennaElement methods:
%
%   step                  - Output the response of the antenna element
%   release               - Allow property name and input characteristics
%                           changes
%   clone                 - Create a custom antenna object with same 
%                           property values
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the element is capable of 
%                           simulating polarization
%   directivity           - Compute element directivity
%   pattern               - Plot element response pattern
%   patternAzimuth        - Plot azimuth pattern
%   patternElevation      - Plot elevation pattern
%
%   CustomAntennaElement properties:
%
%   FrequencyVector            - Operating frequency vector
%   FrequencyResponse          - Frequency responses
%   AzimuthAngles              - Azimuth angles 
%   ElevationAngles            - Elevation angles 
%   SpecifyPolarizationPattern - Specify pattern in polarization
%   MagnitudePattern           - Magnitude pattern 
%   PhasePattern               - Phase pattern 
%   HorizontalMagnitudePattern - Horizontal magnitude pattern
%   HorizontalPhasePattern     - Horizontal phase pattern
%   VerticalMagnitudePattern   - Vertical magnitude pattern
%   VerticalPhasePattern       - Vertical phase pattern
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a user defined antenna with cosine pattern and plot
%   %   its elevation response. Assume the antenna works at 1 GHz.
%   %   The user defined pattern is omnidirectional in the azimuth 
%   %   direction but has a cosine pattern in the elevation direction.
%
%   element = phased.CustomAntennaElement;
%   element.AzimuthAngles = -180:180;
%   element.ElevationAngles = -90:90;
%   element.MagnitudePattern = ...
%       mag2db(repmat(cosd(element.ElevationAngles)',1,...
%       numel(element.AzimuthAngles)));
%   fc = 1e9;
%   pattern(element,fc,-180:180,0,'CoordinateSystem','polar');
%
%   % Example 2:
%   %   Find the response of the above antenna at the boresight.
%
%   element = phased.CustomAntennaElement;
%   element.AzimuthAngles = -180:180;
%   element.ElevationAngles = -90:90;
%   element.MagnitudePattern = ...
%       mag2db(repmat(cosd(element.ElevationAngles)',1,...
%       numel(element.AzimuthAngles)));
%   fc = 1e9; ang = [0;0];
%   resp = element(fc,ang)
%
%   % Example 3:
%   %   Find the polarization response of a custom antenna
%
%   element = phased.CustomAntennaElement(...
%       'SpecifyPolarizationPattern',true);
%   element.VerticalMagnitudePattern = repmat(mag2db(...
%       cosd(element.ElevationAngles(:))),1,numel(element.AzimuthAngles));
%   fc = 1e9; ang = [0;0];
%   resp = element(fc,ang)
%
%   % Example 4:
%   %   Plot the radiation pattern of a frequency dependent antenna at 400
%   %   and 800 MHz. The antenna has a cosine pattern at 300 MHz and a
%   %   raised cosine pattern at 1 GHz.
%
%   element = phased.CustomAntennaElement('FrequencyVector',[3e8 1e9]);
%   az = -180:180; el = -90:90;
%   pat1 = cosd(el')*ones(size(az));
%   pat2 = cosd(el').^1.5*ones(size(az));
%   pat = cat(3,pat1,pat2);
%   element.MagnitudePattern = mag2db(pat);
%   element.PhasePattern = zeros(size(pat));
%   pattern(element,[4e8 8e8],-180:180,0,'CoordinateSystem','polar');
%
%   See also phased, phased.IsotropicAntennaElement,
%   phased.CosineAntennaElement, phased.ShortDipoleAntennaElement, 
%   phased.ULA, phased.URA, phased.ConformalArray.

%   Copyright 2008-2017 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %FrequencyVector    Operating frequency vector (Hz)
        %   Specify the frequencies (in Hz) where the frequency responses
        %   of element are measured as a row vector. The elements of the
        %   vector must be increasing. The default of this property is [0
        %   1e20]. The antenna element has no response outside the
        %   specified frequency range.
        FrequencyVector = [0 1e20];
        %FrequencyResponse  Frequency responses (dB)
        %   Specify the frequency responses (in dB) measured at the
        %   frequencies defined in FrequencyVector property as a row
        %   vector. The length of the vector must equal to the length of
        %   the frequency vector specified in the FrequencyVector property.
        %   The default value of this property is [0 0].
        FrequencyResponse = [0 0];
        %AzimuthAngles   Azimuth angles (deg)
        %   Specify the azimuth angles (in degrees) as a length P vector.
        %   These are the azimuth angles where the custom pattern is
        %   evaluated. P must be greater than 2. The default value of this
        %   property is -180:180.
        AzimuthAngles = -180:180;
        %ElevationAngles   Elevation angles (deg)
        %   Specify the elevation angles (in degrees) as a length Q vector.
        %   These are the elevation angles where the custom pattern is
        %   evaluated. Q must be greater than 2. The default value of this
        %   property is -90:90.
        ElevationAngles = -90:90;
    end

    properties (Logical, Nontunable)
        %SpecifyPolarizationPattern Specify pattern in polarization
        %   Set this property to true to specify pattern in H and V
        %   polarizations. Set this property to false to specify a combined
        %   radiation pattern. The default value of this property is false.
        SpecifyPolarizationPattern = false
    end
    
    properties (Access = private, Nontunable)
        %pRadiationPattern - private property which holds the pattern info
        pMagnitudePattern = ones(181,361)
        pPhasePattern = zeros(181,361)
        pRadiationPattern 
        pHorizontalPattern  
        pVerticalPattern 
        pHorizontalMagnitudePattern = ones(181,361)
        pVerticalMagnitudePattern = ones(181,361)
        pHorizontalPhasePattern = zeros(181,361)
        pVerticalPhasePattern = zeros(181,361)
    end
    
    properties (Access = private)
        pPatternsInitialized = false;
    end
       
    
    properties (Access = private, Logical, Nontunable)
        pFrequencyDependentPattern
    end
    
    properties (Dependent, Hidden, Nontunable)
        %RadiationPattern   Radiation pattern (dB)
        %   Specify the 3D custom field magnitude pattern (in dB) as a QxP
        %   matrix or a QxPxL array, where Q is the number of elements
        %   presented in the ElevationAngles property, P is the number of
        %   elements presented in the AzimuthAngles property, and L is the
        %   number of elements in the FrequencyVector property. The default
        %   value of this property is a 181x361 matrix with all elements
        %   equal to 0.
        %
        %   If RadiationPattern is a matrix, the specified pattern is
        %   applied to all frequencies specified in the FrequencyVector
        %   property. If RadiationPattern is a 3-dimensional array, each
        %   page of the array specifies the pattern for the corresponding
        %   frequency specified in the FrequencyVector property.
        %
        %   If the pattern contains NaN at certain azimuth and elevation
        %   directions, these NaNs are converted to -inf, indicating that
        %   there is no response at those directions.
        %
        %   This property applies when the SpecifyPolarizationPattern
        %   property is false.
        RadiationPattern 
    end
    
    properties (Dependent, Nontunable)
        %MagnitudePattern   Magnitude pattern (dB)
        %   Specify the 3D custom field magnitude pattern (in dB) as a QxP
        %   matrix or a QxPxL array, where Q is the number of elements
        %   presented in the ElevationAngles property, P is the number of
        %   elements presented in the AzimuthAngles property, and L is the
        %   number of elements in the FrequencyVector property. The default
        %   value of this property is a 181x361 matrix with all elements
        %   equal to 0.
        %
        %   If MagnitudePattern is a matrix, the specified pattern is
        %   applied to all frequencies specified in the FrequencyVector
        %   property. If MagnitudePattern is a 3-dimensional array, each
        %   page of the array specifies the pattern for the corresponding
        %   frequency specified in the FrequencyVector property.
        %
        %   If the pattern contains NaN at certain azimuth and elevation
        %   directions, these NaNs are converted to -inf, indicating that
        %   there is no response at those directions.
        %
        %   This property applies when the SpecifyPolarizationPattern
        %   property is false.
        MagnitudePattern 
        %PhasePattern   Phase pattern (dB)
        %   Specify the 3D custom field phase pattern (in dB) as a QxP
        %   matrix or a QxPxL array, where Q is the number of elements
        %   presented in the ElevationAngles property, P is the number of
        %   elements presented in the AzimuthAngles property, and L is the
        %   number of elements in the FrequencyVector property. The default
        %   value of this property is a 181x361 matrix with all elements
        %   equal to 0.
        %
        %   If PhasePattern is a matrix, the specified pattern is
        %   applied to all frequencies specified in the FrequencyVector
        %   property. If PhasePattern is a 3-dimensional array, each
        %   page of the array specifies the pattern for the corresponding
        %   frequency specified in the FrequencyVector property.
        %
        %   This property applies when the SpecifyPolarizationPattern
        %   property is false.
        PhasePattern 
        %HorizontalMagnitudePattern   Horizontal magnitude pattern (dB)
        %   Specify the 3D custom field magnitude pattern (in dB) in
        %   horizontal polarization as a QxP matrix or a QxPxL array, where
        %   Q is the number of elements presented in the ElevationAngles
        %   property, P is the number of elements presented in the
        %   AzimuthAngles property, and L is the number of elements in the
        %   FrequencyVector property. The default value of this property is
        %   a 181x361 matrix with all elements equal to 0.
        %
        %   If HorizontalMagnitudePattern is a matrix, the specified
        %   pattern is applied to all frequencies specified in the
        %   FrequencyVector property. If HorizontalMagnitudePattern is a
        %   3-dimensional array, each page of the array specifies the
        %   pattern for the corresponding frequency specified in the
        %   FrequencyVector property.
        %
        %   If the pattern contains NaN at certain azimuth and elevation
        %   directions, these NaNs are converted to -inf, indicating that
        %   there is no response at those directions.
        %
        %   This property applies when the SpecifyPolarizationPattern
        %   property is true.
        HorizontalMagnitudePattern
        %HorizontalPhasePattern   Horizontal phase pattern (deg)
        %   Specify the 3D custom field phase pattern (in degrees) in
        %   horizontal polarization as a QxP matrix or a QxPxL array, where
        %   Q is the number of elements presented in the ElevationAngles
        %   property, P is the number of elements presented in the
        %   AzimuthAngles property, and L is the number of elements in the
        %   FrequencyVector property. The default value of this property is
        %   a 181x361 matrix with all elements equal to 0.
        %
        %   If HorizontalPhasePattern is a matrix, the specified pattern is
        %   applied to all frequencies specified in the FrequencyVector
        %   property. If HorizontalPhasePattern is a 3-dimensional array,
        %   each page of the array specifies the pattern for the
        %   corresponding frequency specified in the FrequencyVector
        %   property.
        %
        %   This property applies when the SpecifyPolarizationPattern
        %   property is true.
        HorizontalPhasePattern
        %VerticalMagnitudePattern   Vertical magnitude pattern (dB)
        %   Specify the 3D custom field magnitude pattern (in dB) in
        %   vertical polarization as a QxP matrix or a QxPxL array, where Q
        %   is the number of elements presented in the ElevationAngles
        %   property, P is the number of elements presented in the
        %   AzimuthAngles property, and L is the number of elements in the
        %   FrequencyVector property. The default value of this property is
        %   a 181x361 matrix with all elements equal to 0.
        %
        %   If VerticalMagnitudePattern is a matrix, the specified
        %   pattern is applied to all frequencies specified in the
        %   FrequencyVector property. If VerticalMagnitudePattern is a
        %   3-dimensional array, each page of the array specifies the
        %   pattern for the corresponding frequency specified in the
        %   FrequencyVector property.
        %
        %   If the pattern contains NaN at certain azimuth and elevation
        %   directions, these NaNs are converted to -inf, indicating that
        %   there is no response at those directions.
        %
        %   This property applies when the SpecifyPolarizationPattern
        %   property is true.
        VerticalMagnitudePattern
        %VerticalPhasePattern   Vertical phase pattern (deg)
        %   Specify the 3D custom field phase pattern (in degrees) in
        %   vertical polarization as a QxP matrix or a QxPxL array, where Q
        %   is the number of elements presented in the ElevationAngles
        %   property, P is the number of elements presented in the
        %   AzimuthAngles property, and L is the number of elements in the
        %   FrequencyVector property. The default value of this property is
        %   a 181x361 matrix with all elements equal to 0.
        %
        %   If VerticalPhasePattern is a matrix, the specified pattern is
        %   applied to all frequencies specified in the FrequencyVector
        %   property. If VerticalPhasePattern is a 3-dimensional array,
        %   each page of the array specifies the pattern for the
        %   corresponding frequency specified in the FrequencyVector
        %   property.
        %
        %   This property applies when the SpecifyPolarizationPattern
        %   property is true.
        VerticalPhasePattern
    end
    
    methods
        function obj = CustomAntennaElement(varargin)
            obj@phased.internal.AbstractPolarizedAntennaElement(varargin{:});
        end
    end
    methods(Hidden)
        function cl = clonecg(obj)
            if obj.SpecifyPolarizationPattern
                cl = phased.CustomAntennaElement(...
                    'FrequencyVector',obj.FrequencyVector, ...
                    'FrequencyResponse',obj.FrequencyResponse, ...
                    'AzimuthAngles',obj.AzimuthAngles, ...
                    'ElevationAngles',obj.ElevationAngles, ...
                    'SpecifyPolarizationPattern',obj.SpecifyPolarizationPattern, ...
                    'HorizontalMagnitudePattern',obj.HorizontalMagnitudePattern, ...
                    'HorizontalPhasePattern',obj.HorizontalPhasePattern, ...
                    'VerticalMagnitudePattern',obj.VerticalMagnitudePattern, ...
                    'VerticalPhasePattern',obj.VerticalPhasePattern);
                if ~isPolarizationEnabled(obj)
                    disablePolarization(cl);
                end
            else
                 cl = phased.CustomAntennaElement(...
                    'FrequencyVector',obj.FrequencyVector, ...
                    'FrequencyResponse',obj.FrequencyResponse, ...
                    'AzimuthAngles',obj.AzimuthAngles, ...
                    'ElevationAngles',obj.ElevationAngles, ...
                    'SpecifyPolarizationPattern',obj.SpecifyPolarizationPattern, ...
                    'MagnitudePattern',obj.MagnitudePattern, ...
                    'PhasePattern',obj.PhasePattern);
            end
        end
    end

    methods (Access = protected)
        function g = getSpatialResponse(obj,freq,ang)
            if obj.SpecifyPolarizationPattern
                g = getSpatialResponse@...
                    phased.internal.AbstractPolarizedAntennaElement(...
                    obj,freq,ang);
            else
                g = getSinglePattern(obj,freq,ang,'Magnitude');
            end
        end
        
        function g = getHResponse(obj,freq,ang)
            g = getSinglePattern(obj,freq,ang,'Horizontal');
        end
        
        function g = getVResponse(obj,freq,ang)
            g = getSinglePattern(obj,freq,ang,'Vertical');
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
            
            [azq,elq] = meshgrid(azang,elang);
            
            % compute power pattern
            if obj.SpecifyPolarizationPattern
                Nfreq = size(obj.pHorizontalMagnitudePattern,3);
                if Nfreq == 1
                    ppat = struct('H',zeros(Nel,Naz),'V',zeros(Nel,Naz));
                    ppat.H = abs(interp2(obj.AzimuthAngles(:).',...
                        obj.ElevationAngles(:),...
                        obj.pHorizontalMagnitudePattern,...
                        azq,elq,'nearest',0)).^2;
                    ppat.V = abs(interp2(obj.AzimuthAngles(:).',...
                        obj.ElevationAngles(:),...
                        obj.pVerticalMagnitudePattern,...
                        azq,elq,'nearest',0)).^2;
                    
                    % max direction irrelevant to frequency response
                    patc = hypot(ppat.H,ppat.V); 
                    
                    freqr = permute(db2pow(obj.FrequencyResponse),[1 3 2]);
                    ppat.H = bsxfun(@times,ppat.H,freqr);
                    ppat.V = bsxfun(@times,ppat.V,freqr);
                    
                else
                    ppat = struct('H',zeros(Nel,Naz,Nfreq),'V',zeros(Nel,Naz,Nfreq));
                    for m = 1:Nfreq
                        ppat.H(:,:,m) = abs(interp2(obj.AzimuthAngles(:).',...
                            obj.ElevationAngles(:),...
                            obj.pHorizontalMagnitudePattern(:,:,m),...
                            azq,elq,'nearest',0)).^2;
                        ppat.V(:,:,m) = abs(interp2(obj.AzimuthAngles(:).',...
                            obj.ElevationAngles(:),...
                            obj.pVerticalMagnitudePattern(:,:,m),...
                            azq,elq,'nearest',0)).^2;
                    end
                    
                    % max direction irrelevant to frequency response
                    patc = hypot(ppat.H,ppat.V);
                    
                    freqr = permute(db2pow(obj.FrequencyResponse),[1 3 2]);
                    ppat.H = bsxfun(@times,ppat.H,freqr);
                    ppat.V = bsxfun(@times,ppat.V,freqr);
                    
               end
                
                
            else
                Nfreq = size(obj.pMagnitudePattern,3);
                if Nfreq == 1
                    ppat = abs(interp2(obj.AzimuthAngles(:).',...
                        obj.ElevationAngles(:),obj.pMagnitudePattern,...
                        azq,elq,'nearest',0)).^2;
                    
                    % max direction irrelevant to frequency response
                    patc = ppat;
                    
                    freqr = permute(db2pow(obj.FrequencyResponse),[1 3 2]);
                    ppat = bsxfun(@times,ppat,freqr);
                    
                else
                    ppat = zeros(Nel,Naz,Nfreq);
                    for m = 1:Nfreq
                        ppat(:,:,m) = abs(interp2(obj.AzimuthAngles(:).',...
                            obj.ElevationAngles(:),...
                            obj.pMagnitudePattern(:,:,m),...
                            azq,elq,'nearest',0)).^2;
                    end
                    
                    % max direction irrelevant to frequency response
                    patc = ppat;
                    
                    freqr = permute(db2pow(obj.FrequencyResponse),[1 3 2]);
                    ppat = bsxfun(@times,ppat,freqr);
                    
                end
                
            end
            
            [rmax,idxc] = max(patc,[],2);
            rmax = squeeze(rmax);
            idxc = squeeze(idxc);
            [~,idxr] = max(rmax,[],1);
            angmax = zeros(2,Nfreq);
            for m = 1:Nfreq
                angmax(2,m) = elang(idxr(m));
                angmax(1,m) = azang(idxc(idxr(m),m));
            end
            
            pfreq = obj.FrequencyVector;
                
        end
    end

    methods
        function set.AzimuthAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.CustomAntennaElement','AzimuthAngles',...
                {'vector','>=',-180,'<=',180});
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'AzimuthAngles');
            end
            obj.AzimuthAngles = value;
        end

        function set.ElevationAngles(obj,value)
            sigdatatypes.validateAngle(value,...
                'phased.CustomAntennaElement','ElevationAngles',...
                {'vector','>=',-90,'<=',90});
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'ElevationAngles');
            end
            obj.ElevationAngles = value;
        end
        
        function value = get.RadiationPattern(obj)
            if coder.target('MATLAB')
                warning(message('phased:element:ObsoletePropertyByTwoProperties',...
                    'RadiationPattern','MagnitudePattern','PhasePattern'));
            end
            value = phased.internal.magtodb(obj.pMagnitudePattern);
        end

        function value = get.MagnitudePattern(obj)
            value = phased.internal.magtodb(obj.pMagnitudePattern);
        end

        function value = get.PhasePattern(obj)
            value = phased.internal.rad2deg(obj.pPhasePattern);
        end

        function value = get.HorizontalMagnitudePattern(obj)
            value = phased.internal.magtodb(obj.pHorizontalMagnitudePattern);
        end

        function value = get.VerticalMagnitudePattern(obj)
            value = phased.internal.magtodb(obj.pVerticalMagnitudePattern);
        end

        function value = get.HorizontalPhasePattern(obj)
            value = phased.internal.rad2deg(obj.pHorizontalPhasePattern);
        end

        function value = get.VerticalPhasePattern(obj)
            value = phased.internal.rad2deg(obj.pVerticalPhasePattern);
        end

        function set.RadiationPattern(obj,value)
            coder.extrinsic('phased.CustomAntennaElement.removeNAN');
            if coder.target('MATLAB')
                warning(message('phased:element:ObsoletePropertyByTwoProperties',...
                    'RadiationPattern','MagnitudePattern','PhasePattern'));
            end
            if numel(size(value)) == 3
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles) numel(obj.FrequencyVector)];
            else
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles)];
            end
            validateattributes(value,{'numeric'},...
                {'real','nonempty','size',pattern_size},...
                'phased.CustomAntennaElement','RadiationPattern');

            obj.pMagnitudePattern = ...
                coder.internal.const(obj.removeNAN(value));
            obj.pPhasePattern = ...
                coder.internal.const(zeros(size(value)));
        end
        
        function set.MagnitudePattern(obj,value)
            coder.extrinsic('phased.CustomAntennaElement.removeNAN');
            if numel(size(value)) == 3
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles) numel(obj.FrequencyVector)];
            else
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles)];
            end
            validateattributes(value,{'numeric'},...
                {'nonempty','real','size',pattern_size},...
                'phased.CustomAntennaElement','MagnitudePattern');
            obj.pMagnitudePattern =...
                coder.internal.const(obj.removeNAN(value));
        end
        
        function set.HorizontalMagnitudePattern(obj,value)
            coder.extrinsic('phased.CustomAntennaElement.removeNAN');
            if numel(size(value)) == 3
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles) numel(obj.FrequencyVector)];
            else
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles)];
            end
            validateattributes(value,{'numeric'},...
                {'nonempty','real','size',pattern_size},...
                'phased.CustomAntennaElement','HorizontalMagnitudePattern');
            obj.pHorizontalMagnitudePattern =...
                coder.internal.const(obj.removeNAN(value));
        end
        
        function set.VerticalMagnitudePattern(obj,value)
            coder.extrinsic('phased.CustomAntennaElement.removeNAN');
            if numel(size(value)) == 3
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles) numel(obj.FrequencyVector)];
            else
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles)];
            end
            validateattributes(value,{'numeric'},...
                {'nonempty','real','size',pattern_size},...
                'phased.CustomAntennaElement','VerticalMagnitudePattern');

            obj.pVerticalMagnitudePattern = ...
                coder.internal.const(obj.removeNAN(value));
        end
        
        function set.PhasePattern(obj,value)
            if numel(size(value)) == 3
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles) numel(obj.FrequencyVector)];
            else
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles)];
            end
            validateattributes(value,{'numeric'},...
                {'nonempty','real','size',pattern_size},...
                'phased.CustomAntennaElement','PhasePattern');
            obj.pPhasePattern = phased.internal.deg2rad(value);
        end
        
        function set.HorizontalPhasePattern(obj,value)
            if numel(size(value)) == 3
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles) numel(obj.FrequencyVector)];
            else
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles)];
            end
            validateattributes(value,{'numeric'},...
                {'nonempty','real','size',pattern_size},...
                'phased.CustomAntennaElement','HorizontalPhasePattern');
            obj.pHorizontalPhasePattern = phased.internal.deg2rad(value);
        end
        
        function set.VerticalPhasePattern(obj,value)
            if numel(size(value)) == 3
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles) numel(obj.FrequencyVector)];
            else
                pattern_size = [numel(obj.ElevationAngles)...
                    numel(obj.AzimuthAngles)];
            end
            validateattributes(value,{'numeric'},...
                {'nonempty','real','size',pattern_size},...
                'phased.CustomAntennaElement','VerticalPhasePattern');
            obj.pVerticalPhasePattern = phased.internal.deg2rad(value);
        end
        
        function set.FrequencyVector(obj,value)
            validateattributes( value, { 'double' }, ...
                {'nonempty','finite','row','nonnegative' }, '', 'FrequencyVector');
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'FrequencyVector');
            end
            cond = any(diff(value)<0);
            if cond
                coder.internal.errorIf(cond,'phased:element:NotIncreasing');
            end
            obj.FrequencyVector = value;
        end
        
        function set.FrequencyResponse(obj,value)
            validateattributes(value,{'double'},{'real','row','nonnan'},...
                'phased.CustomAntennaElement','FrequencyResponse');
            cond = length(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:element:NotEnoughSamples', 'FrequencyResponse');
            end
            obj.FrequencyResponse = value;
        end
            

    end
    
    methods(Access = protected)
        function validatePropertiesImpl(obj)
            cond = length(obj.FrequencyVector) ~= length(obj.FrequencyResponse);
            if cond
                coder.internal.errorIf(cond,'phased:element:FrequencyResponseMismatch');
            end
            same_freq = find(diff(obj.FrequencyVector)==0);
            cond = any(obj.FrequencyResponse(same_freq) ~= obj.FrequencyResponse(same_freq+1));
            if cond
                coder.internal.errorIf(cond,'phased:element:InvalidFrequencyResponse');
            end
            if obj.SpecifyPolarizationPattern 
                size_hm = size(obj.HorizontalMagnitudePattern);
                size_vm = size(obj.VerticalMagnitudePattern);
                size_hp = size(obj.HorizontalPhasePattern);
                size_vp = size(obj.VerticalPhasePattern);
                cond = ~isequal(size_hm,size_vm) ;
                if cond
                    coder.internal.errorIf(cond,'phased:phased:sizeMismatch',...
                        'HorizontalMagnitudePattern','VerticalMagnitudePattern');
                end
                cond = ~isequal(size_vm,size_hp) ;
                if cond
                    coder.internal.errorIf(cond,'phased:phased:sizeMismatch',...
                        'HorizontalPhasePattern','VerticalMagnitudePattern');
                end
                cond = ~isequal(size_hp,size_vp);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:sizeMismatch',...
                        'HorizontalPhasePattern','VerticalPhasePattern');
                end
            else
                size_m = size(obj.MagnitudePattern);
                size_p = size(obj.PhasePattern);
                cond = ~isequal(size_p,size_m) ;
                if cond
                    coder.internal.errorIf(cond,'phased:phased:sizeMismatch',...
                        'PhasePattern','MagnitudePattern');
                end               
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractPolarizedAntennaElement(obj);
            % starting from R2013a, we save linear scale pattern in
            % pRadiationPattern to improve the performance, avoiding doing
            % db conversion all the time. However, to avoid backwards
            % compatibility issue, we will save pRadiationPattern always as
            % dB and then convert it to linear scale when it gets loaded.
            s.pRadiationPattern = mag2db(obj.pRadiationPattern);
            s.pMagnitudePattern = obj.pMagnitudePattern;
            s.pPhasePattern = obj.pPhasePattern;
            s.pHorizontalMagnitudePattern = obj.pHorizontalMagnitudePattern;
            s.pHorizontalPhasePattern = obj.pHorizontalPhasePattern;
            s.pHorizontalPattern = obj.pHorizontalPattern;
            s.pVerticalMagnitudePattern = obj.pVerticalMagnitudePattern;
            s.pVerticalPhasePattern = obj.pVerticalPhasePattern;
            s.pVerticalPattern = obj.pVerticalPattern;
            s.pPatternsInitialized = obj.pPatternsInitialized;
            s.SupportComplexPattern = true;
        end
        
        function loadObjectImpl(obj,s,~)
            if isfield(s,'pRadiationPattern')
                pat = db2mag(s.pRadiationPattern);
                obj.pRadiationPattern = pat;
                s = rmfield(s,'pRadiationPattern');
                if ~isfield(s,'SupportComplexPattern')
                    s.pMagnitudePattern = abs(pat);
                    s.pPhasePattern = zeros(size(pat));
                else
                    s = rmfield(s,'SupportComplexPattern');
                end
            end
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if (obj.SpecifyPolarizationPattern) 
                if strcmp(prop, 'MagnitudePattern') || ...
                        strcmp(prop, 'PhasePattern')
                    flag = true;
                end
            else
                if strcmp(prop, 'HorizontalMagnitudePattern') || ...
                        strcmp(prop, 'VerticalMagnitudePattern') || ...
                        strcmp(prop, 'HorizontalPhasePattern') || ...
                        strcmp(prop, 'VerticalPhasePattern')
                    flag = true;
                end
            end
        end
        
        function setupImpl(obj,freq,ang,varargin)
            setupImpl@phased.internal.AbstractPolarizedAntennaElement(...
                obj,freq,ang);
            initializePatterns(obj);
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractPolarizedAntennaElement(obj);
            obj.pPatternsInitialized = false;
        end
        
    end
    
    methods (Static, Hidden, Access = protected)        
        function groups = getPropertyGroupsImpl
            p1m = matlab.system.display.internal.Property(...
                'MagnitudePattern','Description', 'Magnitude pattern (dB)', ...
                'UseClassDefault',false,'Default', 'zeros(181,361)');
            p1p = matlab.system.display.internal.Property(...
                'PhasePattern','Description', 'Phase pattern (deg)', ...
                'UseClassDefault',false,'Default', 'zeros(181,361)');
            p2 = matlab.system.display.internal.Property(...
                'HorizontalMagnitudePattern','Description', 'Horizontal magnitude pattern (dB)', ...
                'UseClassDefault',false,'Default', 'zeros(181,361)');
            p3 = matlab.system.display.internal.Property(...
                'HorizontalPhasePattern','Description', 'Horizontal phase pattern (deg)', ...
                'UseClassDefault',false,'Default', 'zeros(181,361)');
            p4 = matlab.system.display.internal.Property(...
                'VerticalMagnitudePattern','Description', 'Vertical magnitude pattern (dB)', ...
                'UseClassDefault',false,'Default', 'zeros(181,361)');
            p5 = matlab.system.display.internal.Property(...
                'VerticalPhasePattern','Description', 'Vertical phase pattern (deg)', ...
                'UseClassDefault',false,'Default', 'zeros(181,361)');
            dSpecifyPolarizationPattern = ...
                matlab.system.display.internal.Property('SpecifyPolarizationPattern', ...
                                                        'IsGraphical', false);

            groups = matlab.system.display.Section(...
                            'PropertyList',...
                            {'FrequencyVector',...
                            'FrequencyResponse',...
                            'AzimuthAngles',...
                            'ElevationAngles',...
                            dSpecifyPolarizationPattern,...
                            p1m,p1p,p2,p3,p4,p5},...
                            'DependOnPrivatePropertyList',...
                            {'MagnitudePattern',...
                            'PhasePattern',...
                            'HorizontalMagnitudePattern',...
                            'HorizontalPhasePattern',...
                            'VerticalMagnitudePattern',...
                            'VerticalPhasePattern'});   
                       
        end
    end
    
    methods (Access=private)
        function initializePatterns(obj) 
            coder.extrinsic('phased.CustomAntennaElement.combinePattern');
            %in codegen initializePatterns can only be called once
            if  ~isempty(coder.target) || ~(obj.pPatternsInitialized)
                if obj.SpecifyPolarizationPattern
                    obj.pFrequencyDependentPattern = ...
                        (numel(size(obj.HorizontalMagnitudePattern)) == 3);
                    
                    obj.pHorizontalPattern = coder.internal.const(...
                        obj.combinePattern( ...
                        obj.pHorizontalMagnitudePattern, ...
                        obj.pHorizontalPhasePattern));
                    
                    obj.pVerticalPattern = coder.internal.const(...
                        obj.combinePattern( ...
                        obj.pVerticalMagnitudePattern, ...
                        obj.pVerticalPhasePattern));
                    
                else
                    obj.pFrequencyDependentPattern = ...
                        (numel(size(obj.MagnitudePattern)) == 3);
                    
                    obj.pRadiationPattern = coder.internal.const(...
                        obj.combinePattern( ...
                        obj.pMagnitudePattern,...
                        obj.pPhasePattern));
                end
                obj.pPatternsInitialized = true;
            end
        end
    end   
    methods (Hidden)    
        function epat = getgpuElemResponse(obj, az, el, freq)
            %This method is used by the phased.gpu.ConstantGammaClutter to
            %compute the response of the antenna element at all
            %combinations of Azimuth, az, and Elevation, el, angles.
            %Note that az and el are in radians.
            initializePatterns(obj);
            
            if(isempty(coder.target))
                
                %Compute the az,el grid for interp2
                [azg, elg] = meshgrid(reshape(deg2rad(obj.AzimuthAngles), [],1), ...
                        reshape(deg2rad(obj.ElevationAngles), 1,[]) );
               
                %If az and el are vectors, reshape them for singleton expansion    
                if isvector(az)
                    azCPU = reshape(gather(az), [], 1);
                    elCPU = reshape(gather(el), 1, []);
                
                else
                    azCPU = gather(az);
                    elCPU = gather(el);
                end
                
                %Assume freq is a scalar. Find nearest frequency if needed.
                if obj.pFrequencyDependentPattern
                    pfreq = obj.FrequencyVector;
                    [~,freq_idx] = min(abs(bsxfun(@minus,...
                        pfreq(:),freq(:).')),[],1);
                else
                    freq_idx = ':';
                end
                
                if obj.SpecifyPolarizationPattern
                    h = interpolatePatternRadians(azg, elg, ...
                        obj.pHorizontalPattern(:,:,freq_idx), ...
                        azCPU, elCPU);
                    v = interpolatePatternRadians(azg, elg, ...
                        obj.pVerticalPattern(:,:,freq_idx), ...
                        azCPU, elCPU);
                    g = hypot(h,v);
                else
                    g = interpolatePatternRadians(azg, elg, ...
                        obj.pRadiationPattern, azCPU, elCPU);
                end
                
                
                H = getFrequencyResponse(obj, freq);
                
                %Multiply the frequency and spatial responses
                epat = bsxfun(@times, H, g);
                
                %If the data came in as a vector, permute to the proper
                %convention of el x position x az
                if isvector(azCPU)
                    epat = permute(epat, [1,3,2]);
                end
                
            else
                coder.internal.assert(false, 'phased:element:NoGPUCodegen');
            end
        end
    end

    methods (Access = private)
        function g = getSinglePattern(obj,freq,ang,option)
            if strcmp(option,'Horizontal')
                ppattern = obj.pHorizontalPattern;
            elseif strcmp(option,'Vertical')
                ppattern = obj.pVerticalPattern;
            else  % magnitude pattern
                ppattern = obj.pRadiationPattern;
            end
            if obj.pFrequencyDependentPattern
                pfreq = obj.FrequencyVector;
                [~,freq_idx] = min(abs(bsxfun(@minus,...
                    pfreq(:),freq(:).')),[],1);
                [unique_freq_idx, ~, f2p_idx] = unique(freq_idx);
                numUniqFreq = numel(unique_freq_idx);
                interp_ppattern = complex(zeros(size(ang,2),numUniqFreq));
                for m = numUniqFreq :-1:1
                    interp_ppattern(:,m) = interpolatePattern(...
                        obj.AzimuthAngles,obj.ElevationAngles,...
                        ppattern(:,:,unique_freq_idx(m)),...
                        ang(1,:), ang(2,:));
                end
                g = interp_ppattern(:,f2p_idx);

            else
                g = interpolatePattern(...
                    obj.AzimuthAngles,obj.ElevationAngles,...
                    ppattern,...
                    ang(1,:), ang(2,:));
                g = g*ones(1,numel(freq));
            end
        end
    end

    methods 
        function flag = isPolarizationCapable(obj)  
            flag = obj.SpecifyPolarizationPattern;
        end
    end
    methods (Static, Hidden)
        function value = removeNAN(value)
        %replace NaN's with -Inf and convert to dB
            value(isnan(value)) = -inf; 
            value = db2mag(value);
        end
        function pattern = combinePattern(magPattern,phasePattern)
            pattern = magPattern.*exp(1i*phasePattern);
        end
    end
end

function pat = interpolatePattern(az,el,pattern,az_q,el_q)
    pat = interpolatePatternRadians(phased.internal.deg2rad(az), phased.internal.deg2rad(el), ...
        pattern, phased.internal.deg2rad(az_q(:)), phased.internal.deg2rad(el_q(:)));
end

function pat = interpolatePatternRadians(azr, elr, pattern, az_qr, el_qr)
 pat = interp2(azr,elr,pattern,...
        az_qr, el_qr,'nearest',0);
end


