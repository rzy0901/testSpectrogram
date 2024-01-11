classdef (Hidden) AbstractArray < matlab.System
%This class is for internal use only. It may be removed in the future.

%AbstractArray   Define the AbstractArray class.

%   Copyright 2009-2017 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Abstract)
        %Taper Array taper
        %   Need to document and process according to the array.
        Taper
    end
    
    properties (Access = protected, Nontunable)
        cElement
        pNumAngles
        pNumFreqs
    end
    
    properties (Access = protected)
        pTaper
    end
    
    properties (Access = protected, Nontunable, Logical)
        pAzimuthOnly;
    end
    
    methods (Access = protected)
        
        function obj = AbstractArray(varargin)
            %AbstractArray   Construct the AbstractArray class.
            setProperties(obj, varargin{1}, varargin{2:end});
        end
        
    end
    
    methods 
        
        function varargout = pattern(obj,freq,varargin)
        %pattern Plot array pattern
        %   pattern(H,FREQ) plots the 3D array directivity pattern (in
        %   dBi). The operating frequency is specified in FREQ (in Hz) as a
        %   positive scalar.
        %
        %   pattern(H,FREQ,AZ) specifies the azimuth angles (in degrees) as
        %   a vector in AZ and plots the 3D array directivity pattern (in
        %   dBi) within the specified azimuth angles. 
        %
        %   pattern(H,FREQ,AZ,EL) specifies the elevation angles (in
        %   degrees) as a vector in EL and plots the array directivity
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
        %   pattern(H,FREQ,AZ,EL,Name,Value) plots the array response
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
        %        PropagationSpeed: A positive scalar specifying the
        %                          signal propagation speed (in m/s) in the
        %                          medium. The default value is
        %                          physconst('lightspeed').
        %                 Weights: A length N column vector or an NxM
        %                          matrix containing the weights applied to
        %                          the array. N is the number of elements
        %                          in the array and M is the number of
        %                          frequencies supplied in FREQ if FREQ is 
        %                          a vector. If Weights is a column vector,
        %                          the same weights will be applied to each
        %                          frequency. If Weights is a matrix, each
        %                          column of weight values will be applied
        %                          to the corresponding frequency in FREQ.
        %                          If Weights is a matrix and FREQ is a
        %                          scalar, each column of weight values
        %                          will be applied to the same frequency in
        %                          FREQ.                          
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
        %   %   Plot the azimuth cut response of a uniform linear array
        %   %   along 0 elevation using a line plot. Assume the operating
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.ULA;
        %   pattern(ha,3e8,-180:180,0,'CoordinateSystem','rectangular')
        %
        %   See also phased, phased.ArrayResponse.
            
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                    'phased:Waveform:CodegenNotSupported','pattern');
            end
            
            narginchk(2,inf);
            nargoutchk(0,3);
            validateattributes(obj,...
                {'phased.internal.AbstractArray'},...
                {'scalar'},'pattern','H');
            
            [fc,c,plotArgs] = phased.internal.parsePatternInputs(...
                true,obj,freq,varargin{:});
            
            if nargout
                plotdata = phased.internal.plotRadiationPattern(...
                    true,obj,fc,c,plotArgs{:});
                varargout{1} = plotdata.resp;
                varargout{2} = plotdata.az;
                varargout{3} = plotdata.el;
            else
                phased.internal.plotRadiationPattern(...
                    true,obj,fc,c,plotArgs{:});
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
        %        PropagationSpeed: A positive scalar specifying the
        %                          signal propagation speed (in m/s) in the
        %                          medium. The default value is
        %                          physconst('lightspeed').
        %                 Weights: A length N column vector or an NxM
        %                          matrix containing the weights applied to
        %                          the array. N is the number of elements
        %                          in the array and M is the number of
        %                          frequencies supplied in FREQ if FREQ is 
        %                          a vector. If Weights is a column vector,
        %                          the same weights will be applied to each
        %                          frequency. If Weights is a matrix, each
        %                          column of weight values will be applied
        %                          to the corresponding frequency in FREQ.
        %                          If Weights is a matrix and FREQ is a
        %                          scalar, each column of weight values
        %                          will be applied to the same frequency in
        %                          FREQ.         
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
        %   %   Plot the azimuth cut response of a uniform linear array
        %   %   along 0 and 10 degrees elevation. Assume the operating 
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.ULA;
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
                {'phased.internal.AbstractArray'},...
                {'scalar'},'pattern','H');
            
            [fc,c,plotArgs] = phased.internal.parsePatternAzimuthInputs(...
                true,obj,freq,varargin{:});
            
            if nargout
                plotdata = phased.internal.plotRadiationPattern(...
                    true,obj,fc,c,plotArgs{:});
                varargout{1} = plotdata.resp;
            else
                phased.internal.plotRadiationPattern(...
                    true,obj,fc,c,plotArgs{:});
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
        %        PropagationSpeed: A positive scalar specifying the
        %                          signal propagation speed (in m/s) in the
        %                          medium. The default value is
        %                          physconst('lightspeed').
        %                 Weights: A length N column vector or an NxM
        %                          matrix containing the weights applied to
        %                          the array. N is the number of elements
        %                          in the array and M is the number of
        %                          frequencies supplied in FREQ if FREQ is 
        %                          a vector. If Weights is a column vector,
        %                          the same weights will be applied to each
        %                          frequency. If Weights is a matrix, each
        %                          column of weight values will be applied
        %                          to the corresponding frequency in FREQ.
        %                          If Weights is a matrix and FREQ is a
        %                          scalar, each column of weight values
        %                          will be applied to the same frequency in
        %                          FREQ.         
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
        %   %   Plot the elevation cut response of a uniform linear array
        %   %   along 0 and 10 degrees azimuth. Assume the operating 
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.ULA;
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
                {'phased.internal.AbstractArray'},...
                {'scalar'},'pattern','H');
            
            [fc,c,plotArgs] = phased.internal.parsePatternElevationInputs(...
                true,obj,freq,varargin{:});
            
            if nargout
                plotdata = phased.internal.plotRadiationPattern(...
                    true,obj,fc,c,plotArgs{:});
                varargout{1} = plotdata.resp;
            else
                phased.internal.plotRadiationPattern(...
                    true,obj,fc,c,plotArgs{:});
            end
           
        end
        
    end
    
    methods
        
        function y = collectPlaneWave(obj,x,ang,freq,c)
        %collectPlaneWave Simulate received plane waves
        %   Y = collectPlaneWave(H,X,ANGLE) returns the received signals at
        %   the sensor array H when the input signals X arrive at the array
        %   from the directions specified in ANGLE.
        %
        %   X must be an M column matrix where M is the number of signals
        %   arriving at the array. Each column of X represents an
        %   individual incoming plane wave signal.
        %
        %   ANGLE can be either a 2xM matrix or a 1xM vector. If ANGLE is a
        %   2xM matrix, the columns are the directions of arrival of the
        %   corresponding signal in X. Each column of ANGLE is in the form
        %   [azimuth; elevation] (in degrees) where the azimuth angle must
        %   be between [-180 180] and the elevation angle must be between
        %   [-90 90]. If ANGLE is a 1xM vector, each entry in ANGLE
        %   specifies the azimuth angle and the corresponding elevation
        %   angle is assumed to be 0.
        %
        %   Y is an N column matrix where N is the number of elements in
        %   the array. Each column of Y is the received signal at the
        %   corresponding element with all incoming signals combined.
        %
        %   Y = collectPlaneWave(H,X,ANGLE,FREQ) uses FREQ (in Hz) as the
        %   incoming signal's carrier frequency. FREQ must be a scalar. The
        %   default value of FREQ is 3e8.
        %
        %   Y = collectPlaneWave(H,X,ANGLE,FREQ,C) uses C (in m/s) as the
        %   signal's propagation speed. C must be a scalar. The default
        %   value of C is the speed of light.
        %
        %   collectPlaneWave modulates the input signal with a phase
        %   corresponding to the delay caused by the direction of arrival.
        %   It does not account for the response of individual elements in
        %   the array.
        %
        %   % Example:
        %   %   Simulate the received signal at a 4-element ULA when two
        %   %   signals are arriving from 10 degrees and 30 degrees
        %   %   azimuth. Both signals have an elevation angle of 0 degrees.
        %   %   Assume the propagation speed is the speed of light and the
        %   %   carrier frequency of the signal is 100 MHz.
        %
        %   ha = phased.ULA(4);
        %   y = collectPlaneWave(ha,randn(4,2),[10 30],1e8,...
        %       physconst('LightSpeed'));
        %
        %   See also phased, phased.Collector.
            
            
            phased.internal.narginchk(3,5,nargin);
            if nargin < 5
                c = physconst('LightSpeed');
            elseif ~isempty(coder.target)
                %error out if c is not constant during codegen
                coder.internal.prefer_const(c);
                coder.internal.assert(eml_is_const(c), ...
                    'phased:collectPlaneWave:propSpeedNotConstCG');
            end
            if nargin < 4
                freq = 3e8;
            end
            validateattributes(x,{'double'},{'2d','nonempty','finite'},...
                'collectPlaneWave','X');
            Nchannel = size(x,2);
            sz_ang = size(ang);
            cond = size(ang,2) ~= Nchannel || sz_ang(1) > 2 || numel(sz_ang) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:collectPlaneWave:incorrectSize', ...
                                       Nchannel, Nchannel, sz_ang( 1 ), sz_ang( 2 ));
            end
            if isa(ang,'double') && (sz_ang(1) == 1)
                angv = [ang; zeros(1,Nchannel)];
            else
                angv = ang;
            end
            sigdatatypes.validateAzElAngle(angv,'collectPlaneWave',...
                'ANGLE');
            sigdatatypes.validateFrequency(freq,'collectPlaneWave',...
                'FREQ',{'scalar'});
            sigdatatypes.validateSpeed(c,'collectPlaneWave','C',{'scalar'});
            validateProperties(obj);

            if isempty(coder.target)
                privArray = clone(obj);
                release(privArray);
            else
                privArray = clonecg(obj);
            end


            %Include Taper but not response
            privHstv = phased.SteeringVector('SensorArray',privArray,...
                'PropagationSpeed',c,'IncludeElementResponse',false);
            sv = step(privHstv,freq,angv);
            lTaper = getTaper(obj);
            y = x*(bsxfun(@times,sv,lTaper)).';
        end
    end
    
    methods (Hidden)
        
        function varargout = plotResponse(obj,freq,v,varargin)
        %plotResponse Plot array response
        %   plotResponse(H,FREQ,V) plots the array response pattern along
        %   the azimuth cut where the elevation angle is 0. The operating
        %   frequency is specified in FREQ (in Hz). The propagation speed
        %   is specified in V (in m/s).
        %
        %   plotResponse(H,FREQ,V,Name,Value,...) plots the array response
        %   with the specified parameter Name set to the specified Value.
        %   You can specify additional name-value pair arguments in any
        %   order as (Name1,Value1,...,NameN,ValueN). The parameter Names
        %   are
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
        %                          180. The default value is 0.
        %            Polarization: Polarization response, using one of 
        %                          |'None' | 'Combined' | 'H' | 'V' |,
        %                          indicating the response without
        %                          polarization, the combined response, the
        %                          horizontal polarization response and the
        %                          vertical polarization response,
        %                          respectively. The default value is
        %                          'None'. Setting this parameter to other
        %                          than 'None' requires the array to be
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
        %                 Weights: A length N column vector or an NxM
        %                          matrix containing the weights applied to
        %                          the array. N is the number of elements
        %                          in the array and M is the number of
        %                          frequencies supplied in FREQ if FREQ is 
        %                          a vector. If Weights is a column vector,
        %                          the same weights will be applied to each
        %                          frequency. If Weights is a matrix, each
        %                          column of weight values will be applied
        %                          to the corresponding frequency in FREQ.
        %                          If Weights is a matrix and FREQ is a
        %                          scalar, each column of weight values
        %                          will be applied to the same frequency in
        %                          FREQ.                          
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
        %                           
        %
        %   % Example:
        %   %   Plot the azimuth cut response of a uniform linear array
        %   %   along 0 elevation using a line plot. Assume the operating
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.ULA;
        %   plotResponse(ha,3e8,physconst('LightSpeed'))
        %
        %   See also phased, phased.ArrayResponse.
            
            if isempty(coder.target)
                narginchk(3,inf);
                validateattributes(obj,{'phased.internal.AbstractArray'},...
                                   {'scalar'},'plotResponse','H');
            
              
                [freq,v,plotArgs] = phased.internal.parsePlotResponseInputs(...
                    true,obj,freq,v,varargin{:});

                if nargout
                    plotdata = phased.internal.plotRadiationPattern(...
                        true,obj,freq,v,plotArgs{:});
                    varargout{1} = plotdata.fighandle;
                else
                    phased.internal.plotRadiationPattern(...
                        true,obj,freq,v,plotArgs{:});
                end
            else
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','plotResponse');
            end
            
        end
    end
    
    methods
        
        function varargout = viewArray(obj,varargin)
        % viewArray     View array geometry 
        %   viewArray(H) plots the geometry of array specified in H. 
        %
        %   viewArray(H,Name,Value,...) plots the geometry of the array
        %   with the specified parameter Name set to the specified Value.
        %   You can specify additional name-value pair arguments in any
        %   order as (Name1,Value1,...,NameN,ValueN). The parameter Names
        %   are
        %
        %   ShowNormals: A logical flag to specify whether to show the 
        %                normal directions of elements in the array. The 
        %                default value is false.
        %
        %     ShowTaper: A logical flag to specify whether to change the  
        %                element color brightness in proportion to the    
        %                element Taper magnitude. The default value is 
        %                false.
        %
        %     ShowIndex: A vector specifying the element indices to show in
        %                the figure. You can also specify 'All' to show
        %                indices of all elements of the array or 'None' to
        %                suppress indices. The default value is 'None'.
        %
        %         Title: A string specifying the title of the plot. The
        %                default value is 'Array Geometry'.
        %
        %   % Examples:
        %   
        %   % Example 1:
        %   %   Display the geometry of a 4x4 URA.
        %
        %   ha = phased.URA(4);
        %   viewArray(ha);
        %
        %   % Example 2:
        %   %   Display the geometry of 6-element ULA and show the indices
        %   %   for the first and third elements.
        %
        %   ha = phased.ULA(6);
        %   viewArray(ha,'ShowIndex',[1 3]);
        %
        %   % Example 3:
        %   %   Display the geometry of a 4x4 URA and display the normal
        %   %   directions and indices for all elements.
        %
        %   ha = phased.URA(4);
        %   viewArray(ha,'ShowNormals',true,'ShowIndex','All');
        %
        %   See also phased, phased.ArrayResponse
            
          if ~isempty(coder.target)
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','viewArray');
          end
          
            %Validating the attributes of input object.
            validateattributes(obj,{'phased.internal.AbstractArray'},...
                {'scalar'},'','H');
            
            ShowTaper     = false;
            ShowPattern   = false;
            ShowNormals   = false;
            ShowIndex     = 'None';
            Title         = 'Array Geometry';
            AxesHandle    = [];
            elem_pos      = getElementPosition(obj) ;
            N             = getNumElements(obj);
            
            sigutils.pvparse(varargin{:});
            
            defaultPosRGB = [0.04 0.51 0.78];
            defaultNormalRGB = [0.84 0.16 0];
            
            if ShowNormals
                validateattributes(ShowNormals,{'logical'},{'scalar' },...
                    'viewArray','ShowNormals'); %#ok<UNRCH>
            end
            
            if ShowTaper
                validateattributes(ShowTaper,{'logical'},{'scalar' },...
                    'viewArray','ShowTaper'); %#ok<UNRCH>
                validateProperties(obj);
            end
            
            if ShowPattern
                validateattributes(ShowPattern,{'logical'},{'scalar' },...
                    'viewArray','ShowPattern'); %#ok<UNRCH>
                validateProperties(obj);
                if ~isa(obj,'phased.internal.AbstractHeterogeneousArray')
                   %ignore Showpattern                   
                    ShowPattern = 0;
                end
            end
            
            if ~isempty(ShowIndex)
                if ischar(ShowIndex)
                    ShowIndex = validatestring(ShowIndex,{'All','None'},...
                        'viewArray','ShowIndex');
                    if strcmp(ShowIndex,'All')
                        ShowIndex = 1:N;
                    else
                        ShowIndex = [];
                    end
                else
                    sigdatatypes.validateIndex(ShowIndex,'viewArray',...
                        'ShowIndex',{'vector','<=',N});
                end
            end
            
            validateattributes(Title,{'char'},{},'viewArray','Title');
            
            hold_status = ishold;
            % ClnupObj      = onCleanup(@()hold('off'));
            if ~hold_status
                if ~isempty(AxesHandle)
                  tag = get(AxesHandle,'Tag');
                  cla(AxesHandle,'reset');
                  set(AxesHandle,'Tag',tag,'Visible','off');
                  hold(AxesHandle,'on');
                else
                  clf;
                  hold on;
                end
            end
            
            % compute normal direction
            normal_pos_sph = deg2rad(getElementNormal(obj));
            [normalX, normalY, normalZ] = sph2cart(normal_pos_sph(1,:),...
                   normal_pos_sph(2,:),1); %#ok<ASGLU>
           
            % compute axis span
            elem_pos_max   = max(elem_pos,[],2);
            elem_pos_min   = min(elem_pos,[],2);
            axis_length = max(elem_pos_max - elem_pos_min);
                      
            XPos  = axis_length/2;
            YPos  = axis_length/2;
            ZPos  = axis_length/2;
            
            % Determine the marker and font size
            if size(elem_pos,2) > 1
                thresh = sqrt(sum((elem_pos(:,1)-elem_pos(:,2)).^2))/(2*YPos);
            else  % single element
                thresh = 1;
            end
            
            if thresh<0.012
                sz_marker = 36;
                sz_font   = 6;
            elseif thresh>0.21
                sz_marker = 108;
                sz_font = 8;
            else
                sz_marker = 72;
                sz_font   = 8;
            end
            
            if ShowPattern
                ColorOrder = [...
                    perms([0 0.25 0.5]);...
                    perms([0 0.25 0.75]);...
                    perms([0 0.5 0.75]);...
                    perms([0.25 0.5 0.75]);...
                    perms([0.25 0.5 1]);...
                    perms([0.25 0.75 1]);...
                    perms([0.5 0.75 1])];
                rs = rng(6);
                ColorOrder = ColorOrder(randperm(size(ColorOrder,1)),:);
                rng(rs);
                ColorOrder = [get(0,'DefaultAxesColorOrder');ColorOrder];
                ColorOrderIdx = mod(obj.ElementIndices(:)-1,length(ColorOrder))+1;
                default_color =  ColorOrder(ColorOrderIdx,:); 
            else
                default_color     = repmat(defaultPosRGB,...
                    numel(elem_pos(1,:)),1);
            end
            
            if ShowTaper
                lTaper = abs(getTaper(obj)); %#ok<UNRCH>
                normTaper = lTaper/max(lTaper);
                default_color = default_color.*repmat(normTaper,1,3);
            end
            
            %If AxesHandle passes set gca to AxesHandle
            if ~isempty(AxesHandle) 
                axes(AxesHandle); 
            end
        
            if ShowPattern
                % Create dummy elements for legend      
                legendLimit = min(length(obj.ElementSet),10);
                legendColor = ColorOrder(1:legendLimit,:);
                hlegends = gobjects(legendLimit,1);
                for m = 1:legendLimit
                    legend_pos = elem_pos(:,1)';
                    hlegends(m) = scatter3(...
                        legend_pos(:,1),legend_pos(:,2),legend_pos(:,3),...
                        sz_marker,legendColor(m,:),'filled',...
                        'Tag','Elements',...
                        'MarkerEdgeColor','black');
                end

                legendstr = cell(legendLimit,1);
                for i = 1:legendLimit
                    legendstr{i} = sprintf('Pattern %d',i);
                end
                % in painter's mode we get a warning that rgb is not supported
                set(gcf, 'Renderer', 'zbuffer');
                legend(hlegends,legendstr{:},'AutoUpdate','off');
            end
            hElements = scatter3(...
                elem_pos(1,:)',elem_pos(2,:)',elem_pos(3,:)',...
                sz_marker,default_color,'filled',...
                'Tag','Elements',...
                'MarkerEdgeColor','black');
            set(hElements,'DeleteFcn',{@deleteAssociatedText});
            if nargout
                varargout{1} = hElements;
            end
                      
            % show indices
            if ~isempty(ShowIndex)
                if size(elem_pos,2)>1
                    min_xdis = min(sqrt(sum(...
                        (elem_pos-circshift(elem_pos,[0 -1])).^2)));
                else  % single element
                    min_xdis = 0.15;
                end
                for m = 1:length(ShowIndex)
                    distfactor = 9;
                    text(elem_pos(1,ShowIndex(m))+min_xdis/distfactor,...
                        elem_pos(2,ShowIndex(m))+min_xdis/distfactor,...
                        elem_pos(3,ShowIndex(m))+min_xdis/distfactor,...
                        num2str(ShowIndex(m)),'FontSize',sz_font);
                end
            end
            
            % show normals
            if ShowNormals
                hQuiver = quiver3(elem_pos(1,:)',elem_pos(2,:)',elem_pos(3,:)',...
                    normalX',normalY',normalZ',0.5,...
                    'Color',defaultNormalRGB,...
                    'LineWidth',1,...
                    'Tag','Normals',...
                    'MaxHeadSize',0.1); %#ok<UNRCH>
            end
            
            origFormat = get(0, 'format');
            format('shortG');
            plot_align = 'right'; 
            
            if ~(isa(obj,'phased.ConformalArray') || ...
                    isa(obj,'phased.HeterogeneousConformalArray'))
                if isa(obj,'phased.UCA')
                    [rpv, rpu] = convert2engstrs(obj.Radius,1);
                    [spv, spu] = convert2engstrs(getElementSpacing(obj,'arc'),1);
                    arPl = getArrayPlaneOrAxis(obj);
                    arrayspan_str = {...
                           ['Radius = ' rpv ' ' rpu 'm'],... 
                           'Element Spacing: ',...
                           [' Arc = ' spv ' ' spu 'm'],...
                           ['Array Plane: ' arPl ' plane']};
                    span = abs(elem_pos_max - elem_pos_min);
                elseif isa(obj,'phased.ULA') || isa(obj,'phased.HeterogeneousULA')
                    spac = obj.ElementSpacing;
                    aper = N*spac(1);
                    [apv, apu] = convert2engstrs(aper);
                    [spv, spu] = convert2engstrs(spac(1));
                    arAx = getArrayPlaneOrAxis(obj);
                    arrayspan_str = {...
                           'Aperture Size: ',...
                           ['    ' upper(arAx) ' axis = ' apv ' ' apu 'm'],... 
                           'Element Spacing: ',...
                           ['    \Delta ' arAx ' = ' spv ' ' spu 'm'],...
                           ['Array Axis: ' upper(arAx) ' axis']};
                else
                    spac = obj.ElementSpacing;
                    arPl = getArrayPlaneOrAxis(obj);
                    if isa(obj,'phased.URA')
                        siz = obj.Size;
                    else % phased.HeterogeneousURA
                        siz = size(obj.ElementIndices);
                    end
                    aper = siz.*spac;
                    [apv, apu] = convert2engstrs(aper);
                    [spv, spu] = convert2engstrs(spac);
                    arrayspan_str = {...
                           'Aperture Size: ',...
                           ['    ' arPl(1) ' axis = ' apv{2} ' ' apu 'm'],...
                           ['    ' arPl(2) ' axis = ' apv{1} ' ' apu 'm'],...
                           'Element Spacing: ',...
                           ['    \Delta ' lower(arPl(1)) ' = ' spv{2} ' ' spu 'm'],...
                           ['    \Delta ' lower(arPl(2)) ' = ' spv{1} ' ' spu 'm']}; 
              
                end
            else
                span = abs(elem_pos_max - elem_pos_min);
                [spanv, spanu] = convert2engstrs(span,3);
                arrayspan_str = {...
                       'Array Span:    ',...
                       ['    X axis = ' spanv(1,:) ' ' spanu 'm'],...
                       ['    Y axis = ' spanv(2,:) ' ' spanu 'm'],...
                       ['    Z axis = ' spanv(3,:) ' ' spanu 'm']};
            end
            
            % scaling
            axis_scale = 2;
            axis_pos_min = elem_pos_min-axis_length/10*axis_scale;
            
            % plot axes           
            line([axis_pos_min(1),(XPos/4)+axis_pos_min(1)],...
                [axis_pos_min(2),axis_pos_min(2)],...
                [axis_pos_min(3),axis_pos_min(3)],...
                'Color','k','LineWidth',1,'Tag','X Axis','Clipping','off');
            line([axis_pos_min(1),axis_pos_min(1)],...
                [axis_pos_min(2),(YPos/4)+axis_pos_min(2)],...
                [axis_pos_min(3),axis_pos_min(3)],...
                'Color','k','LineWidth',1,'Tag','Y Axis','Clipping','off');
            line([axis_pos_min(1),axis_pos_min(1)],...
                [axis_pos_min(2),axis_pos_min(2)],...
                [axis_pos_min(3),(ZPos/4)+axis_pos_min(3)],...
                'Color','k','LineWidth',1,'Tag','Z Axis','Clipping','off');
            
            text(1.11*(XPos/4)+axis_pos_min(1),...
                axis_pos_min(2),axis_pos_min(3),...
                'x','Tag','X Label','FontSize',8);
            text(axis_pos_min(1),...
                1.11*(YPos/4)+axis_pos_min(2),axis_pos_min(3),...
                'y','Tag','Y Label','FontSize',8);
            text(axis_pos_min(1),...
                axis_pos_min(2),1.12*(ZPos/4)+axis_pos_min(3),...
                'z','Tag','Z Label','FontSize',8);
            
            figureRGB = [0.94 0.94 0.94];
            axisRGB = [0.49 0.49 0.49];
            % Don't change colors or add annotations when AxesHandle
            % passed
            if isempty(AxesHandle)
                set(gcf,'Color',figureRGB);
                set(gca,'XColor',axisRGB,...
                    'YColor',axisRGB,...
                    'ZColor',axisRGB,...
                    'Color',figureRGB,...
                    'Position',[0 0 1 1]);

                xlabel('x axis (Az 0 El 0) -->','Color',axisRGB);
                ylabel('y axis -->','Color',axisRGB);
                zlabel('z axis (Az 0 El 90) --->','Color',axisRGB);

                annotatebox = ...
                    annotation('textbox',[0.7 0 0.2 0.25],'String',arrayspan_str,...
                    'FitBoxToText','on','LineStyle','none','FontSize',8,...
                    'HorizontalAlignment',plot_align,'Tag','Annotation',...
                    'Color',[0.501 0.501 0.501]);
                titlebox = annotation('textbox',[0 0.75 1 0.2],...
                    'String',Title,'FitBoxToText','off',...
                    'FontWeight','demi',...
                    'LineStyle','none','HitTest','off',...
                    'HorizontalAlignment','center','Tag','Title');
                setappdata(hElements,'annotatebox',annotatebox);
                setappdata(hElements,'titlebox',titlebox);
            end
            
            view(45, 45);
            axis vis3d;
            box off;
            grid off;
            axis off;
            if(~isempty(AxesHandle))
                hold(AxesHandle,'off');
            else
                hold off;
            end
            daspect([1 1 1])
            camlight('infinite');
            
            if isa(obj,'phased.URA') || isa(obj,'phased.HeterogeneousURA')
                if ~ShowNormals
                    % Show planar arrays in their plane
                    view(90,0)
                else
                    view(60,45)
                end
            elseif isa(obj,'phased.ConformalArray') || ...
                     isa(obj,'phased.UCA') || ...
                    isa(obj,'phased.HeterogeneousConformalArray')
                if ~ShowNormals
                    % Show planar arrays in their plane
                    if span(1)==0 && span(2)~=0 && span(3)~=0
                        view(90,0),
                    elseif span(2)==0 && span(1)~=0 && span(3)~=0
                        view(0,0)
                    elseif span(3)==0 && span(1)~=0 && span(2)~=0
                        view(0,90)
                    end
                end
                if isa(obj,'phased.UCA')
                    arNr = obj.ArrayNormal;
                    hold on;
                    r = obj.Radius;
                    t = -180:180;
                    if strcmpi(arNr,'x')
                        plot3(zeros(1,numel(t)),r*cosd(t),r*sind(t),'Color',[.3 .3 .3],'LineStyle','--');
                        hold off;
                    elseif strcmpi(arNr,'y')
                        plot3(r*cosd(t),zeros(1,numel(t)),r*sind(t),'Color',[.3 .3 .3],'LineStyle','--');
                        hold off;
                    else
                        plot3(r*cosd(t),r*sind(t),zeros(1,numel(t)),'Color',[.3 .3 .3],'LineStyle','--');
                        hold off;
                    end
                end            
            end

            % clean up
            format(origFormat);
            if(~isempty(AxesHandle))
                if hold_status 
                    hold(AxesHandle,'on');
                else
                    hold(AxesHandle,'off');
                end
            else
                if hold_status
                    hold on;
                else
                    hold off;
                end
            end
        end
    end
    
    methods
        function d = directivity(obj,freq,ang,varargin)
        %directivity Compute array directivity
        %   D = directivity(H,FREQ,ANGLE) computes the directivity (in dBi)
        %   of the array for the directions specified in ANGLE (in degrees)
        %   and frequencies specified in FREQ (in Hz). FREQ is a row vector
        %   of length L and ANGLE can be either a row vector of length M or
        %   a 2xM matrix. D is an MxL matrix whose columns contain the
        %   directivity of the array at angles specified in ANGLE at
        %   corresponding frequencies specified in FREQ.
        %  
        %   When ANGLE is a 2xM matrix, each column of the matrix specifies
        %   the direction in space in [azimuth; elevation] form. The
        %   azimuth angle should be between [-180 180] degrees and the
        %   elevation angle should be between [-90 90] degrees. If ANGLE is
        %   a length M row vector, each element specifies a direction's
        %   azimuth angle and the corresponding elevation angle is assumed
        %   to be 0.
        %
        %   D = directivity(H,FREQ,ANGLE,Name,Value,...) computes the
        %   directivity with specified parameter Name set to the specified
        %   Value. You can specify additional name-value pair arguments in
        %   any order as (Name1,Value1,...,NameN,ValueN). The parameter
        %   Names are
        %
        %       PropagationSpeed: A positive scalar to specify the
        %                         propagation speed (in m/s). The default
        %                         value is the speed of light.
        %                Weights: A length-N column vector or an NxL
        %                         matrix containing the weights applied to
        %                         the array. N is the number of elements
        %                         in the array and L is the number of
        %                         frequencies supplied in FREQ if FREQ is 
        %                         a vector. If Weights is a column vector,
        %                         the same weights will be applied to each
        %                         frequency. If Weights is a matrix, each
        %                         column of weight values will be applied
        %                         to the corresponding frequency in FREQ.
        %
        %   % Examples:
        %
        %   % Example 1:
        %   %   Compute the directivity of a 10-element ULA with quarter
        %   %   wavelength spacing. The element is assumed to be isotropic.
        %
        %   c = 3e8; fc = 3e8; lambda = c/fc;
        %   myArray = phased.ULA(10,lambda/4);
        %   d = directivity(myArray,fc,0,'PropagationSpeed',c)
        %
        %   % Example 2:
        %   %   Compute the directivity of a ULA of 10 isotropic element
        %   %   with quarter wavelength spacing. The array is steered 
        %   %   toward end-fire direction.
        %
        %   c = 3e8; fc = 3e8; lambda = c/fc; ang = [90;0];
        %   myArray = phased.ULA(10,lambda/4);
        %   w = steervec(getElementPosition(myArray)/lambda,ang);
        %   d = directivity(myArray,fc,ang,'PropagationSpeed',c,...
        %           'Weights',w)
        %
        %   See also phased, phased.ArrayResponse.
        
            phased.internal.narginchk(3,7,nargin);
            defaultPropagationSpeed = physconst('lightspeed');
            defaultWeights = ones(getDOF(obj),1);
            if coder.target('MATLAB')
                p = inputParser;
                addParameter(p,'PropagationSpeed',defaultPropagationSpeed);
                addParameter(p,'Weights',defaultWeights);

                parse(p,varargin{:});
                c = p.Results.PropagationSpeed;
                w = p.Results.Weights;

            else
                parms = struct('PropagationSpeed',uint32(0), ...
                            'Weights',uint32(0));

                pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
                c = eml_get_parameter_value(pstruct.PropagationSpeed,...
                    defaultPropagationSpeed,varargin{:});
                w = eml_get_parameter_value(pstruct.Weights,...
                    defaultWeights,varargin{:});
            end
            
            sigdatatypes.validateFrequency(freq,'directivity','FREQ',...
                {'row'});
            
            sigdatatypes.validateSpeed(c,'directivity',...
                'PropagationSpeed',{'scalar'});
            if ~coder.target('MATLAB')
                %error out if c is not constant during codegen
                coder.internal.prefer_const(c);
                coder.internal.assert(eml_is_const(c), ...
                    'phased:collectPlaneWave:propSpeedNotConstCG');
            end
                        
            sigdatatypes.validateAngle(ang,'directivity','ANGLE');
            if isrow(ang)
                angIn = [ang;zeros(size(ang))];
            else
                angIn = ang;
            end
            sigdatatypes.validateAzElAngle(angIn,'directivity','ANGLE');
            
            if iscolumn(w)
                validateattributes(w,{'double'},{'nonempty','finite',...
                    'nrows',getDOF(obj)},'directivity','W');
            else
                validateattributes(w,{'double'},{'2d','nonempty','finite',...
                    'size',[getDOF(obj) numel(freq)]},'directivity','W');
            end
            
            if coder.target('MATLAB')
                myDirectivity = phased.internal.Directivity(...
                    'Sensor',clone(obj),...
                    'PropagationSpeed',c,'WeightsInputPort',true);
            else
                myDirectivity = phased.internal.Directivity(...
                    'Sensor',clonecg(obj),...
                    'PropagationSpeed',c,'WeightsInputPort',true);
            end
            d = step(myDirectivity,freq,angIn,w);
            
        end
        
    end
    
    methods (Access = protected)
        function flag = isInputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function flag = isInputSizeLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function validateInputsImpl(obj,freq,angle)   %#ok<INUSL>

            
            cond = ~isrow(freq) || isempty(freq);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeRowVector','FREQ');
            end
            cond = ~isa(freq,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','FREQ','double');
            end
            cond = ~isreal(freq);
            if cond
                coder.internal.errorIf(cond, ...
                'phased:system:array:expectReal','FREQ');
            end
            
            sz_angle = size(angle);
            cond = ~ismatrix(angle) || isempty(angle);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','ANGLE');
            end
            cond = sz_angle(1) > 2;
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:system:array:NeedTwoRows','ANGLE');
            end
            cond = ~isreal(angle);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:system:array:InvalidAngle','ANGLE');
            end
            cond = ~isa(angle,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','ANGLE','double');
            end
        end
        
        function num = getNumInputsImpl(obj)  %#ok<MANU>
            num = 2;
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.pNumAngles = obj.pNumAngles;
                s.pNumFreqs = obj.pNumFreqs;
                s.pAzimuthOnly = obj.pAzimuthOnly;
                s.pTaper = obj.pTaper;
            end
        end
        
    end
    
    methods (Access = protected, Sealed)
        
        function resp = stepImpl(obj,freq,angArg)
            
            num_angles = obj.pNumAngles;
            if obj.pAzimuthOnly
                ang = [angArg; zeros(1,num_angles)];
            else
                ang = angArg;
            end
            cond = any(ang(1,:)>180) || any(ang(2,:)>90);
            if cond
                coder.internal.errorIf(cond,'phased:step:AzElNotLessEqual');
            end
            cond = any(ang(1,:)<-180) || any(ang(2,:)<-90);
            if cond
                coder.internal.errorIf(cond,'phased:step:AzElNotGreaterEqual');
            end            
            
            if isPolarizationEnabled(obj)
                resp = getPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles);
            else
                resp = getNonPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles);
            end
                
        end
    end
    
    methods (Abstract)
        N = getNumElements(obj);
        %getNumElements Returns the number of elements in the array.
        pos = getElementPosition(obj,eleIdx);
        %getElementPosition Returns the positions of elements in the
        %array.
        w = getTaper(obj);
        %getTaper Return the taper in a column vector
    end
    
    methods (Access = protected, Abstract)
        resp = getPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles)
        %getPolarizationOutputInArrayAzEl Return polarized response
        resp = getNonPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles)
        %getNonPolarizationOutputInArrayAzEl Return nonpolarized response
        arPl = getArrayPlaneOrAxis(obj)
        %getArrayPlaneOrAxis Return array plane or array axis based on the
        %System Object
    end
    
    methods (Access = protected)
        function az_el = convertIncidentToAzEl(obj,inc_ang)
        %convertIncidentToAzEl Converts the incident angle to corresponding
        %normal directions for an element in the array.
        %   AZ_EL = convertIncidentToAzEl(H,INC_ANG) converts the incident
        %   directions specified in INC_ANG (in degrees), in the form of
        %   [azimuth; elevation] to the local azimuth and elevation angles
        %   (in degrees) of sensor elements in the array H. AZ_EL is a
        %   2x(N*M) matrix where M is the number of rows in INC_ANG and N
        %   is the number of elements in the array. Each column of the
        %   AZ_EL represent the translated azimuth and elevation angle, in
        %   the form of [azimuth; elevation], of the corresponding incident
        %   angle in INC_ANG. The first N columns of AZ_EL is the [az; el]
        %   for the elements corresponding to the first angle specified in
        %   INC_ANG. The second N columns of AZ_EL is the [az; el] for the
        %   elements corresponding to the second angle specified in
        %   INC_ANG, and so on.
            
            % normal direction is [0; 0] for ULA and URA
            N = getNumElements(obj);
            az_el = repmat(inc_ang,N,1);
            az_el = reshape(az_el,2,[]);
            
        end
        
    end
    
    methods(Hidden)
        function num = getDOF(obj)
        %getDOF  Degree of freedom of the array
        %   N = getDOF(H) returns the degree of freedom (DOF), N, of the
        %   array H.
        %
        %   % Example:
        %   %   Construct a ULA and then obtain its degree of freedom.
        %
        %   ha = phased.ULA;
        %   dof = getDOF(ha)
            
            num = getNumElements(obj);
        end
    end
    
    
    methods (Hidden, Abstract, Access = {?phased.internal.AbstractArray, ?phased.gpu.internal.AbstractClutterSimulator})
        %Methods used by the GPU ConstantGammaClutter model.
        xscaled = scaleByElemResponse(obj, azin, elin, freq, x, idx);
    end
    
    methods (Access = {?phased.ArrayGain})
        function Rn = getNoiseGainMatrix(obj,freq,c,stang) %#ok<INUSD>
            Rn = eye(getDOF(obj));
        end
    end
    methods (Static,Hidden)
        function a = getAlternateBlock
            a = 'phasedtxrxlib/Narrowband Receive Array';
        end
    end
end

function deleteAssociatedText(src,~)
    annotatebox = getappdata(src,'annotatebox');
    if ishghandle(annotatebox)
        delete(annotatebox);
    end
    titlebox = getappdata(src,'titlebox');
    if ishghandle(titlebox)
        delete(titlebox);
    end
end
