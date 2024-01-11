classdef (Hidden) AbstractElement < matlab.System
%This class is for internal use only. It may be removed in the future.

%AbstractElement   Define the AbstractElement class.

%   Copyright 2009-2016 The MathWorks, Inc.
%     

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    properties (Access = protected,Logical,Nontunable)
        % whether to output the combine pattern. For polarization capable
        % antennas, turning this to true to use combined pattern at each
        % antenna, i.e., discard polarization information.
        pOutputCombinedPattern = false
        % whether to initialize output as complex or real
        pIsOutputComplex = true;        
    end
    
    methods (Access = protected)

        function obj = AbstractElement(varargin)
            %AbstractElement   Construct the AbstractElement class.
            setProperties(obj, nargin, varargin{:});
        end
        
        function integratePattern(obj,angstep)
            if nargin < 2
                angstep = [1 1];
            end
            [obj.pTotalPowerIntensity,obj.pMaximumDirectivityAngle] = ...
                calcIntegratedPower(obj,angstep);
        end

    end
    
    methods (Abstract, Access = protected)
        resp = getSpatialResponse(obj,freq,angle);
        %getSpatialResponse Returns the spatial response of the
        %element at a given direction for a given freq that is inside the
        %element's frequency range.
        resp = getFrequencyResponse(obj,freq);
        %getFrequencyResponse Returns the frequency response of the element
        %at a given frequency.
        frange = getFrequencyRange(obj);
        %getFrequencyRange Returns the frequency range of the element.
    end
    
    methods (Abstract, Access = {?phased.internal.AbstractElement,...
            ?phased.internal.AbstractSensorOperation})
        [ppat,pfreq,maxdir] = getPowerPattern(obj,azang,elang)
        %getPowerPattern Returns the power pattern and maximum direction.
        %For multiple frequency patterns, it returns a 3D matrix of ppat
        %and a matrix of maxdir whose columns represents direction. Pattern
        %frequency is also returned in pfreq as a vector. azang and elang
        %are row vectors containing sampling angles along azimuth and
        %elevation.
    end
    
    methods (Access = {?phased.internal.AbstractElement,...
            ?phased.internal.AbstractArray})
        %flag to make sure whether the array normal and the element normal
        %are aligned or not. In all element cases the normals are aligned,
        %For the short dipole case the element normal and element normal
        %are not aligned.
        function isAligned = isElementNormalArrayNormalAligned(obj)  %#ok<MANU>
            isAligned = true;
        end
    end
    
    methods (Access = {?phased.internal.AbstractElement,...
            ?phased.internal.AbstractSensorOperation});
        function flag = isDirectivityKnown(obj) %#ok<MANU>
            flag = false;
        end
    end
    
    methods (Access = protected)
        function num = getNumInputsImpl(obj)  %#ok<MANU>
            num = 2;
        end
        
        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end

        function validateInputsImpl(~,freq,angle)
            
            cond =  ~isrow(freq) || isempty(freq);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeRowVector','FREQ');
            end

            cond =  ~isa(freq,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','FREQ','double');
            end

            cond =  ~isreal(freq);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:system:element:InvalidFrequency');
            end

            sz_angle = size(angle);
            cond =  ~ismatrix(angle) || isempty(angle);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeMatrix','ANGLE');
            end

            cond =  sz_angle(1) > 2;
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:system:element:NeedTwoRows');
            end

            cond =  ~isreal(angle);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:system:element:InvalidAngle');
            end

            cond =  ~isa(angle,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','ANGLE','double');
            end
                  
        end
        
        function flag = isInputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
    end
    
    methods (Access = protected)
        
        function resp = stepImpl(obj,freq,angArg, varargin)
            
            [ang, num_ang] = validateInputValues(obj,freq,angArg);
            
            resp = getOutputResponse(obj,freq,ang,'combined',num_ang);

        end
        
        function d = getDirectivity(obj,freq,ang)
            myDirectivity = phased.internal.Directivity(...
                'Sensor',clone(obj));
            d = step(myDirectivity,freq,ang);
        end
        
        function resp = getOutputResponse(obj,freq,ang,option,num_ang)
            isValidFreq = searchValidFrequency(obj,freq);
            validfreq = freq(isValidFreq);           
            if obj.pIsOutputComplex
                resp = complex(zeros(num_ang,numel(freq)));
            else
                resp = zeros(num_ang,numel(freq));                
            end            
            if ~isempty(validfreq)
                validResp = (ones(num_ang,1)*...
                    getFrequencyResponse(obj,validfreq).').*...
                    getOutputSpatialResponse(obj,validfreq,ang,option);

                resp(:,isValidFreq) = validResp;
            end
        end
        
        function resp = getOutputSpatialResponse(obj,freq,ang,~)
            resp = getSpatialResponse(obj,freq,ang);
        end
        
        function [ang, num_ang] = validateInputValues(obj,freq,angArg) %#ok<INUSL>
            cond = any(freq < 0);
            if cond
                coder.internal.errorIf(cond,'phased:step:expectedNonnegative', 'FREQ');
            end
            
            ang_size = size(angArg);
            num_ang = ang_size(2);
            
            if ang_size(1) == 1
                ang = [angArg; zeros(1,num_ang)];
            else
                ang = angArg;
            end
            cond = any(ang(1,:)>180) || any(ang(2,:)>90);
            if cond
                coder.internal.errorIf(cond, ...
                                       'phased:step:AzElNotLessEqual');
            end
        
            cond = any(ang(1,:)<-180) || any(ang(2,:)<-90);
            if cond
                coder.internal.errorIf(cond,...
                                       'phased:step:AzElNotGreaterEqual');
            end
            
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.pOutputCombinedPattern = obj.pOutputCombinedPattern;
            s.pIsOutputComplex = obj.pIsOutputComplex;
        end
        
    end
    
    methods(Hidden)
        function num = getDOF(obj) %#ok<MANU>
        %getDOF  Degree of freedom of the array
        %   N = getDOF(H) returns the degree of freedom (DOF), N, of the
        %   element H.
        %
        %   % Example:
        %   %   Construct an isotropic antenna and then obtain its degree 
        %   %   of freedom.
        %
        %   ha = phased.IsotropicAntennaElement;
        %   dof = getDOF(ha)
            
            num = 1;
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
        %   %   Compute the directivity of a short dipole at 300 MHz toward
        %   %   the boresight.
        %
        %   myAnt = phased.ShortDipoleAntennaElement;
        %   d = directivity(myAnt,3e8,0)
        %
        %   % Example 2:
        %   %   Compute the directivity of a cosine antenna at 1 GHz toward
        %   %   30 degrees azimuth and 10 degrees elevation.
        %
        %   myAnt = phased.CosineAntennaElement;
        %   d = directivity(myAnt,1e9,[30;10])
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

            d = getDirectivity(obj,freq,angIn);
        end
    end
    
    methods 
        
        function varargout = pattern(obj,freq,varargin)
        %pattern Plot element response pattern
        %   pattern(H,FREQ) plots the 3D element directivity pattern (in
        %   dBi). The operating frequency is specified in FREQ (in Hz) as a
        %   positive scalar.
        %
        %   pattern(H,FREQ,AZ) specifies the azimuth angles (in degrees) as
        %   a vector in AZ and plots the 3D element directivity pattern (in
        %   dBi) within the specified azimuth angles. 
        %
        %   pattern(H,FREQ,AZ,EL) specifies the elevation angles (in
        %   degrees) as a vector in EL and plots the element directivity
        %   pattern (in dBi) based on the value of FREQ, AZ, and EL. At
        %   least one of FREQ, AZ, and EL must be a scalar. If FREQ is a
        %   scalar, and both AZ and EL are vectors, the 3D pattern at the
        %   specified frequency in FREQ is plotted in the region of azimuth
        %   and elevation angles specified by AZ and EL. If FREQ and AZ are
        %   vectors and EL is a scalar, the patterns at specified
        %   frequencies in FREQ are plotted along azimuth angles specified
        %   in AZ with fixed elevation angle specified in EL. If FREQ and
        %   EL are vectors and AZ is a scalar, the patterns at specified
        %   frequencies in FREQ are plotted along elevation angles
        %   specified in EL with fixed azimuth angle specified in AZ.
        %
        %   pattern(H,FREQ,Name,Value), pattern(H,FREQ,AZ,Name,Value),
        %   pattern(H,FREQ,AZ,EL,Name,Value) plots the element response
        %   pattern with the specified parameter Name set to the specified
        %   Value. You can specify additional name-value pair arguments in
        %   any order as (Name1,Value1,...,NameN,ValueN). The parameter
        %   Names are
        %       
        %        CoordinateSystem: Specify the coordinate system of the 
        %                          pattern using one of | 'polar' | 
        %                          'rectangular' | 'uv' |. The default 
        %                          value is 'polar'.
        %                          If CoordinateSystem is 'polar' or
        %                          'rectangular', AZ and EL specifies the
        %                          azimuth and elevation angles (in 
        %                          degrees), respectively. AZ must be
        %                          between -180 and 180 and EL must be
        %                          between -90 and 90.
        %                          If CoordinateSystem is 'uv', AZ and EL
        %                          specifies the u and v coordinates,
        %                          respectively. Both AZ and EL must be
        %                          between -1 and 1.
        %                    Type: Specify the type of the pattern as one
        %                          of 'directivity' | 'efield' | 'power' |
        %                          'powerdb', where the default is
        %                          'directivity'. 
        %                          If Type is 'directivity', the pattern is
        %                          the directivity measured in dBi. If Type
        %                          is 'efield', the pattern is the
        %                          field pattern of the sensor. If Type is
        %                          'power', the pattern is the power
        %                          pattern of the sensor, i.e., the
        %                          square of the field pattern. If Type is
        %                          'powerdb', the pattern is the power
        %                          pattern measured in dB.
        %               Normalize: Set this to true to normalize the
        %                          response pattern. Set this to false to
        %                          not normalize the response pattern. The
        %                          default value is true. This parameter is
        %                          not applicable when you set Type to
        %                          'directivity'.
        %               PlotStyle: Specify the plotting style when multiple
        %                          frequencies are specified as one of
        %                          'overlay' | 'waterfall', where 'Overlay'
        %                          is the default value.
        %            Polarization: Specify the pattern polarization using 
        %                          one of 'combined' | 'H' | 'V', where the
        %                          default is 'combined'. This parameter is
        %                          only applicable when the sensor is
        %                          polarization capable and when the Type
        %                          is not 'directivity'. 
        %                          If Polarization is 'combined', the H and
        %                          V polarization pattern is combined. If
        %                          Polarization is 'H', only the H
        %                          polarization is considered. If
        %                          Polarization is 'V', only the V
        %                          polarization is considered.
        %
        %   [PAT,AZ_ANG,EL_ANG] = pattern(...) returns the pattern in PAT.
        %   AZ_ANG is a length-Q row vector. If Coordinate parameter is
        %   'uv', it contains sampling values in U. Otherwise, it contains
        %   sampling angles (in degrees) in azimuth. EL_ANG is a length-P
        %   row vector. If Coordinate parameter is 'uv', it contains
        %   sampling values in V. Otherwise, it contains sampling angles
        %   (in degrees) in elevation. PAT is a PxQ matrix whose elements
        %   are the patterns at corresponding sampling points specified by
        %   AZ_ANG and EL_ANG.
        %
        %   % Example:
        %   %   Plot the azimuth cut response of an isotropic antenna along 
        %   %   0 elevation using a polar plot. Assume the operating 
        %   %   frequency is 1 GHz.
        %
        %   ha = phased.IsotropicAntennaElement;
        %   pattern(ha,1e9,-180:180,0)
        %
        %   See also phased, phased.ArrayResponse.
            
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                    'phased:Waveform:CodegenNotSupported','pattern');
            end
            
            narginchk(2,inf);
            nargoutchk(0,3);
            validateattributes(obj,...
                {'phased.internal.AbstractElement'},...
                {'scalar'},'pattern','H');
            
            [fc,~,plotArgs] = phased.internal.parsePatternInputs(...
                false,obj,freq,varargin{:});
            
            if nargout
                plotdata = phased.internal.plotRadiationPattern(...
                    false,obj,fc,plotArgs{:});
                varargout{1} = plotdata.resp;
                varargout{2} = plotdata.az;
                varargout{3} = plotdata.el;
            else
                phased.internal.plotRadiationPattern(...
                    false,obj,fc,plotArgs{:});
            end
           
        end
        
        function varargout = patternAzimuth(obj,freq,varargin)
        %patternAzimuth Plot azimuth pattern
        %   patternAzimuth(H,FREQ) plots directivity pattern (in dBi) along
        %   azimuth at 0 degree elevation angle. The operating frequency is
        %   specified in FREQ (in Hz) as a positive scalar.
        %
        %   patternAzimuth(H,FREQ,EL) specifies elevation angles (in
        %   degrees) in EL as a length-Q row vector.
        %
        %   patternAzimuth(H,FREQ,EL,Name,Value), plots the array response
        %   pattern with the specified parameter Name set to the specified
        %   Value. You can specify additional name-value pair arguments in
        %   any order as (Name1,Value1,...,NameN,ValueN). The parameter
        %   Names are
        %       
        %                    Type: Specify the type of the pattern as one
        %                          of 'directivity' | 'efield' | 'power' |
        %                          'powerdb', where the default is
        %                          'directivity'. 
        %                          If Type is 'directivity', the pattern is
        %                          the directivity measured in dBi. If Type
        %                          is 'efield', the pattern is the
        %                          field pattern of the sensor. If Type is
        %                          'power', the pattern is the power
        %                          pattern of the sensor, i.e., the
        %                          square of the field pattern. If Type is
        %                          'powerdb', the pattern is the power
        %                          pattern measured in dB.
        %                 Azimuth: A length-P row vector specifying the 
        %                          azimuth angles where the pattern is
        %                          calculated. The default value is
        %                          -180:180.
        %
        %   PAT = patternAzimuth(...) returns the pattern in PAT. PAT is a
        %   PxQ matrix whose elements are the patterns at corresponding
        %   sampling points specified by Azimuth parameter and EL input
        %   argument.
        %
        %   % Example:
        %   %   Plot the azimuth cut response of an cosine antenna 
        %   %   along 0 and 10 degrees elevation. Assume the operating 
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.CosineAntennaElement;
        %   patternAzimuth(ha,3e8,[0 10])
        %
        %   See also phased, phased.ArrayResponse.
            
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                    'phased:Waveform:CodegenNotSupported','pattern');
            end
            
            narginchk(2,inf);
            nargoutchk(0,1);
            validateattributes(obj,...
                {'phased.internal.AbstractElement'},...
                {'scalar'},'pattern','H');
            
            [fc,~,plotArgs] = phased.internal.parsePatternAzimuthInputs(...
                false,obj,freq,varargin{:});
            
            if nargout
                plotdata = phased.internal.plotRadiationPattern(...
                    false,obj,fc,plotArgs{:});
                varargout{1} = plotdata.resp;
            else
                phased.internal.plotRadiationPattern(...
                    false,obj,fc,plotArgs{:});
            end
           
        end
    
        function varargout = patternElevation(obj,freq,varargin)
        %patternElevation Plot elevation pattern
        %   patternElevation(H,FREQ) plots directivity pattern (in dBi)
        %   along elevation at 0 dgree azimuth angle (in degrees). The
        %   operating frequency is specified in FREQ (in Hz) as a positive
        %   scalar.
        %
        %   patternElevation(H,FREQ,AZ) specifies azimuth angles (in
        %   degrees) in AZ as a length-Q row vector.
        %
        %   patternElevation(H,FREQ,AZ,Name,Value), plots the array
        %   response pattern with the specified parameter Name set to the
        %   specified Value. You can specify additional name-value pair
        %   arguments in any order as (Name1,Value1,...,NameN,ValueN). The
        %   parameter Names are
        %       
        %                    Type: Specify the type of the pattern as one
        %                          of 'directivity' | 'efield' | 'power' |
        %                          'powerdb', where the default is
        %                          'directivity'. 
        %                          If Type is 'directivity', the pattern is
        %                          the directivity measured in dBi. If Type
        %                          is 'efield', the pattern is the
        %                          field pattern of the sensor. If Type is
        %                          'power', the pattern is the power
        %                          pattern of the sensor, i.e., the
        %                          square of the field pattern. If Type is
        %                          'powerdb', the pattern is the power
        %                          pattern measured in dB.
        %               Elevation: A length-P row vector specifying the 
        %                          elevation angles where the pattern is
        %                          calculated. The default value is
        %                          -90:90.
        %
        %   PAT = patternElevation(...) returns the pattern in PAT. PAT is
        %   a PxQ matrix whose elements are the patterns at corresponding
        %   sampling points specified by Elevation parameter and AZ input
        %   argument.
        %
        %   % Example:
        %   %   Plot the elevation cut response of a cosine antenna
        %   %   along 0 and 10 degrees azimuth. Assume the operating 
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.CosineAntennaElement;
        %   patternElevation(ha,3e8,[0 10])
        %
        %   See also phased, phased.ArrayResponse.
            
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                    'phased:Waveform:CodegenNotSupported','pattern');
            end
            
            narginchk(2,inf);
            nargoutchk(0,1);
            validateattributes(obj,...
                {'phased.internal.AbstractElement'},...
                {'scalar'},'pattern','H');
            
            [fc,~,plotArgs] = phased.internal.parsePatternElevationInputs(...
                false,obj,freq,varargin{:});
            
            if nargout
                plotdata = phased.internal.plotRadiationPattern(...
                    false,obj,fc,plotArgs{:});
                varargout{1} = plotdata.resp;
            else
                phased.internal.plotRadiationPattern(...
                    false,obj,fc,plotArgs{:});
            end
           
        end
        
    end
    
    methods (Hidden)

        function varargout = plotResponse(obj,freq,varargin)
        %plotResponse Plot element response
        %   plotResponse(H,FREQ) plots the element response pattern along
        %   the azimuth cut where the elevation angle is 0. The operating
        %   frequency is specified in FREQ (in Hz).
        %
        %   plotResponse(H,FREQ,Name,Value) plots the element response with
        %   the specified parameter Name set to the specified Value. You
        %   can specify additional name-value pair arguments in any order
        %   as (Name1,Value1,...,NameN,ValueN). The parameter Names are
        %       
        %                  Format: Format of the plot, using one of 
        %                          | 'Line' | 'Polar' | 'UV' |. The default
        %                          value is 'Line'.
        %                 RespCut: Cut of the response. If Format is 'Line'
        %                          or 'Polar', RespCut uses one of | 'Az' |
        %                          'El' | '3D' |, where the default
        %                          value is 'Az'.  If Format is 'UV',  
        %                          RespCut uses one of | 'U' | '3D' |,
        %                          where the default is 'U'.
        %                          If you set RespCut to 'Az' or 'El' or 
        %                          'U', FREQ can be either a scalar or a 
        %                          row vector. If you set RespCut to '3D', 
        %                          FREQ must be a scalar.
        %                CutAngle: A scalar to specify the cut angle. This
        %                          parameter is only applicable when
        %                          RespCut is set to 'Az' or 'El'. If
        %                          RespCut is 'Az', the CutAngle must be
        %                          between -90 and 90. If RespCut is 'El',
        %                          the CutAngle must be between -180 and
        %                          180.
        %            Polarization: Polarization response, using one of 
        %                          |'None' | 'Combined' | 'H' | 'V' |,
        %                          indicating the response without
        %                          polarization, the combined response, the
        %                          horizontal polarization response and the
        %                          vertical polarization response,
        %                          respectively. The default value is
        %                          'None'. Setting this parameter to other
        %                          than 'None' requires the element to be
        %                          capable of simulating polarization. This
        %                          option is not applicable when you set
        %                          Unit to 'dbi'.
        %                    Unit: The unit of the plot, using one of
        %                          | 'db' | 'mag' | 'pow' | 'dbi' |. The
        %                          default value is 'db'. If unit is 'mag',
        %                          the plot is the field pattern; if unit
        %                          is 'pow', the plot is the power pattern;
        %                          if unit is 'db', the plot is the power
        %                          pattern in dB scale; and if unit is
        %                          'dbi', the plot is the directivity
        %                          pattern.
        %       NormalizeResponse: Set this to true to normalize the
        %                          response pattern. Set this to false to
        %                          not normalize the response pattern. The
        %                          default value is true. This option is
        %                          not applicable when you set Unit to
        %                          'dbi'.
        %             OverlayFreq: A logical flag to specify whether 
        %                          pattern cuts will be overlaid in a 2-D
        %                          line plot or plotted against frequency
        %                          in a waterfall plot. The default value
        %                          is true. This parameter only applies
        %                          when Format is not 'Polar' and RespCut
        %                          is not '3D'. When OverlayFreq is false,
        %                          FREQ must be a vector with at least two
        %                          entries.
        %           AzimuthAngles: A row vector to specify the Azimuthal
        %                          angles and the spacing between these 
        %                          angles over which you can visualize
        %                          the radiation pattern. This parameter is
        %                          only applicable when the RespCut is 'Az'
        %                          or '3D' and the Format is 'Line' or 
        %                          'Polar'. The range of the Azimuthal angles
        %                          specified should be between -180 and 180. 
        %                          The default value is -180:180.
        %         ElevationAngles: A row vector to specify the Elevation
        %                          angles and the spacing between these 
        %                          angles over which you can visualize the
        %                          radiation pattern. This parameter is only
        %                          applicable when the RespCut is 'El' or
        %                          '3D' and the Format is 'Line' or 'Polar'.
        %                          The range of the Elevation angles
        %                          specified should be between -90 and 90. 
        %                          The default value is -90:90.
        %                   UGrid: A row vector to specify the U coordinates
        %                          in the U/V Space over which you want to 
        %                          visualize the radiation pattern. This 
        %                          parameter is applicable only when the
        %                          Format is 'UV'. The range of the U
        %                          coordinates specified should be between
        %                          -1 and 1. The default value is -1:0.01:1.
        %                   VGrid: A row vector specifying the V coordinates
        %                          in the U/V Space over which you want to 
        %                          visualize the radiation pattern. This 
        %                          parameter is applicable only when the 
        %                          Format is 'UV' and the RespCut is '3D'.
        %                          The range of the V coordinates specified
        %                          should be between -1 and 1. The default
        %                          value is -1:0.01:1.
        %
        %   % Example:
        %   %   Plot the azimuth cut response of an isotropic antenna along 
        %   %   0 elevation using a line plot. Assume the operating 
        %   %   frequency is 1 GHz.
        %
        %   ha = phased.IsotropicAntennaElement;
        %   plotResponse(ha,1e9)
        %
        %   See also phased, phased.ArrayResponse.

            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','plotResponse');
            end
            
            narginchk(2,inf);
            validateattributes(obj,...
                {'phased.internal.AbstractElement'},...
                {'scalar'},'plotResponse','H');
            
            [freq,~,plotArgs] = phased.internal.parsePlotResponseInputs(...
                false,obj,freq,varargin{:});
            
            if nargout
                plotdata = phased.internal.plotRadiationPattern(...
                    false,obj,freq,plotArgs{:});
                varargout{1} = plotdata.fighandle;
            else
                phased.internal.plotRadiationPattern(...
                    false,obj,freq,plotArgs{:});
            end
           
        end
    end
    
    methods (Hidden, Sealed)
        function [isValidFreq, freqRange] = searchValidFrequency(obj,freq)
            freqRange = getFrequencyRange(obj);
            isValidFreq = ~((freq < freqRange(1)) | (freq > freqRange(2)));
        end
    end
    
    methods 
        function flag = isPolarizationCapable(obj) 
        %isPolarizationCapable Indicate if the element is capable of 
        %simulating polarization
        %   F = isPolarizationCapable(H) returns the flag, F, which
        %   indicates whether the element H is capable of simulating
        %   polarization.
        %
        %   % Example:
        %   %   Determine whether the IsotropicAntennaElement is capable of
        %   %   simulating polarization
        %   
        %   h = phased.IsotropicAntennaElement;
        %   f = isPolarizationCapable(h)
        %
        %   See also phased.
            flag = false;
        end
        
    end    
    
    methods (Hidden)
        function flag = isPolarizationEnabled(obj) 
        %isPolarizationEnabled Indicate if the element is enabled to 
        %simulate polarization
        %   F = isPolarizationEnabled(H) returns the flag, F, which
        %   indicates whether the element H is enabled to simulate
        %   polarization.
        %
        %   % Example:
        %   %   Determine whether the IsotropicAntennaElement is capable of
        %   %   simulating polarization
        %   
        %   h = phased.IsotropicAntennaElement;
        %   f = isPolarizationEnabled(h)
        %
        %   See also phased.
            flag = isPolarizationCapable(obj) && ~obj.pOutputCombinedPattern;
        end
        
    end    
    
    methods (Hidden) %, Access = {?phased.SteeringVector, ?phased.internal.AbstractArray, ?phased.SteeringVector})
        function disablePolarization(obj)
            if isa(obj,'phased.internal.AbstractPolarizedAntennaElement')
                obj.pOutputCombinedPattern = true;
            end
        end
        
        function newObj = cloneSensor(obj)
            newObj = clone(obj);
            release(newObj);
        end
        
        function flag = isElementFromAntenna(obj) %#ok<MANU>
            flag = false;
        end
    end
    methods (Static,Hidden)
        function a = getAlternateBlock
            a = 'phasedtxrxlib/Narrowband Receive Array';
        end
    end    
end


