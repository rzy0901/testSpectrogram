classdef (Hidden) AbstractSubarray < matlab.System
%This class is for internal use only. It may be removed in the future.

%AbstractSubarray   Define the AbstractSubarray class.

%   Copyright 2011-2017 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %SubarraySteering   Subarray steering method
        %   Specify the method of subarray steering as one of 'None' |
        %   'Phase' | 'Time' | 'Custom', where the default is 'None'. When
        %   you set this property to 'Phase', the phase shifter is used to
        %   steer the subarray to a steering angle specified as an input.
        %   When you set this property to 'Time', the time delay is used to
        %   steer the subarray to a steering angle specified as an input.
        %   When you set this property to 'Custom', the complex gain
        %   coefficient for each element in the array is specified as an
        %   input.
        SubarraySteering = 'None'
        %PhaseShifterFrequency  Phase shifter frequency (Hz)
        %   Specify the operating frequency (in Hz) of phase shifters used
        %   to perform subarray steering as a positive scalar. This
        %   property is applicable when you set the SubarraySteering
        %   property to 'Phase'. The default value of this property is 3e8.
        PhaseShifterFrequency = 3e8
        %NumPhaseShifterBits    Number of bits in phase shifters
        %   Specify the number of bits used in the phase shifter as a
        %   non-negative integer. This property is applicable when you set
        %   the SubarraySteering property to 'Phase'. The default value of
        %   this property is 0, indicating there is no quantization effect
        %   in the phase shifter.
        NumPhaseShifterBits = 0
    end

    properties (Constant,Hidden)
        SubarraySteeringSet = matlab.system.StringSet(...
            {'None','Phase','Time','Custom'});
    end
    
    properties (Access = protected, Nontunable)
        pSubarrayPosition
        pSubarrayFeed = 'Corporate';
        pFeedWeights
        cArray
        cModuleResponse
        cModuleSteeringVector
    end
    
    properties (Access = protected, PositiveInteger, Nontunable)
        pNumSubarrays
        pNumElements
        pNumAngles
        pNumFreqs
    end
    
    properties (Access = protected, Logical, Nontunable)
        pNoSteeringFlag = true
        pPhaseSteeringFlag = false
        pTimeSteeringFlag = false
        pIsSingleFrequency = false
        pIsSingleWeights = true
    end
    
    methods (Access = protected, Abstract)
        num = calcNumSubarrays(obj)
        num = calcNumElements(obj,subarrayidx)
        pos = calcSubarrayPosition(obj)
        pos = calcElementPosition(obj)
    end
    
    methods (Abstract)
        flag = isPolarizationCapable(obj)
    end
    
    methods (Hidden,Abstract)
        flag = isPolarizationEnabled(obj)
    end
    
    methods (Access=protected,Abstract)
        validateElementWeights(obj,w)
    end
    
    methods
        function set.PhaseShifterFrequency(obj,value)
            sigdatatypes.validateFrequency(value,...
                '','PhaseShifterFrequency',...
                {'scalar','positive'});
            obj.PhaseShifterFrequency = value;
        end
        function set.NumPhaseShifterBits(obj,val)
            validateattributes(val,{'double'},...
                {'scalar','integer','nonnegative'},...
                '','NumPhaseShifterBits');
            obj.NumPhaseShifterBits = val;
        end
    end
    
    methods (Access = protected)
        
        function obj = AbstractSubarray(varargin)
            setProperties(obj,nargin,varargin{:});
        end
        
        function [f,theta] = getSteeringConfiguration(obj,freq,stang_in)
            if isscalar(stang_in)
                stang = [stang_in; 0];
            else
                stang = stang_in;
            end
            cond =  stang(1)>180 || stang(2)>90;
            if cond
                coder.internal.errorIf(cond,'phased:step:AzElNotLessEqual');
            end
            cond = stang(1)<-180 || stang(2)<-90;
            if cond
                coder.internal.errorIf(cond,'phased:step:AzElNotGreaterEqual');
            end
            
            if obj.pPhaseSteeringFlag
                f = obj.PhaseShifterFrequency;
                theta = stang;
            else % obj.pTimeSteeringFlag
                f = freq;
                theta = stang;
            end
        end
        
        function [Na,N,M,L] = calcRespDimensions(obj)
            Na = getNumElements(obj.cArray);
            N = obj.pNumSubarrays;
            M = obj.pNumAngles;
            L = obj.pNumFreqs;
        end
        
    end
    
    methods
        
        function num = getNumSubarrays(obj)
        %getNumSubarrays   Number of subarrays in array
        %   N = getNumSubarrays(H) returns the total number of subarrays,
        %   N, in the array H.
        %
        %   % Example:
        %   %   Construct an array with replicated subarrays and obtain its
        %   %   number of subarrays.
        %
        %   ha = phased.ReplicatedSubarray;
        %   N = getNumSubarrays(ha)
            

            if isempty(coder.target)
                if ~isLocked(obj)
                    validateProperties(obj);
                    num = calcNumSubarrays(obj);
                else
                    num = obj.pNumSubarrays;
                end
            else %cannot use isLocked for codegen since needs to constant fold
                if ~coder.internal.is_defined(obj.pNumSubarrays)
                    validateProperties(obj);
                    num = calcNumSubarrays(obj);
                else
                    num = obj.pNumSubarrays;
                end
            end
        end
        
        function num = getNumElements(obj,subarrayidx)
            %getNumElements   Number of elements in array
            %   N = getNumElements(H) returns the total number of elements, N,
            %   in the array H.
            %
            %   % Example:
            %   %   Construct an array with replicated subarrays and obtain its
            %   %   number of elements.
            %
            %   ha = phased.ReplicatedSubarray;
            %   N = getNumElements(ha)
            
            if nargin < 2
                if isempty(coder.target)
                    if ~isLocked(obj)
                        validateProperties(obj);
                        num = calcNumElements(obj);
                    else
                        num = obj.pNumElements;
                    end
                else%cannot use isLocked for codegen since needs to constant fold
                    if ~coder.internal.is_defined(obj.pNumSubarrays)
                        validateProperties(obj);
                        num = calcNumElements(obj);
                    else
                        num = obj.pNumElements;
                    end
                end
            else
                validateProperties(obj)
                num = calcNumElements(obj,subarrayidx);
            end
        end
        
        function pos = getSubarrayPosition(obj)
        %getSubarrayPosition Subarray positions in array
        %   POS = getSubarrayPosition(H) returns the subarray positions in
        %   the array H. POS is a 3xN matrix where N is the number of
        %   subarrays in H. Each column of POS defines the position, in the
        %   form of [x; y; z] (in meters), of a subarray in the local
        %   coordinate system.
        %
        %   % Example:
        %   %   Construct an array with replicated subarrays and then
        %   %   obtain its subarray positions.
        %
        %   ha = phased.ReplicatedSubarray;
        %   pos = getSubarrayPosition(ha)
            if isempty(coder.target)
                if ~isLocked(obj)
                    validateProperties(obj);
                    pos = calcSubarrayPosition(obj);
                else
                    pos = obj.pSubarrayPosition;
                end
            else%cannot use isLocked for codegen since needs to constant fold
                if ~coder.internal.is_defined(obj.pSubarrayPosition)
                    validateProperties(obj);
                    pos = calcSubarrayPosition(obj);
                else
                    pos = obj.pSubarrayPosition;
                end
            end
        end
        
        function pos = getElementPosition(obj)
        %getElementPosition Element positions in array
        %   POS = getElementPosition(H) returns the element positions in
        %   the array H. POS is a 3xN matrix where N is the number of
        %   elements in H. Each column of POS defines the position, in the
        %   form of [x; y; z] (in meters), of an element in the local
        %   coordinate system.
        %
        %   % Example:
        %   %   Construct an array with replicated subarrays and then
        %   %   obtain its element positions.
        %
        %   ha = phased.ReplicatedSubarray;
        %   pos = getElementPosition(ha)
            
            if ~isLocked(obj)
                validateProperties(obj);
            end
            pos = calcElementPosition(obj);
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
        %                          | 'None' | 'Combined' | 'H' | 'V' |, 
        %                          indicating the combined response, the 
        %                          horizontal polarization response and the 
        %                          vertical polarization response, 
        %                          respectively. The default value is 
        %                          'None'. Setting this parameter to other 
        %                          than 'None' requires the subarray to be 
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
        %                          the array. N is the number of subarrays
        %                          in the array and M is the number of
        %                          frequencies supplied in FREQ if FREQ is 
        %                          a vector. If Weights is a column vector,
        %                          the same weights will be applied to each 
        %                          frequency. If Weights is a matrix,each
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
        %              SteerAngle: A scalar or a length-2 column vector to
        %                          specify the subarray steering angle. If
        %                          SteerAng is a vector, it is in the form
        %                          of [az; el] (in degrees). az must be
        %                          between -180 and 180 while el must be
        %                          between -90 and 90. If SteerAng is a
        %                          scalar, it specifies the azimuth angle
        %                          and the elevation angle is assumed to be
        %                          0. The default value is [0;0]. This
        %                          option is only applicable if the
        %                          SubarraySteering property of H is set to
        %                          'Phase' or 'Time'.
        %          ElementWeights: A matrix or a cell array specifying
        %                          the weights applied to each element
        %                          in the subarray. This parameter is
        %                          applicable when you set the
        %                          SubarraySteering property of H to
        %                          'Custom'. The default value is all ones.
        %                          
        %                          For a ReplicatedSubarray, ElementWeights
        %                          must be an NSExN matrix where NSE is the
        %                          number of elements in each individual 
        %                          subarray and N is the number of 
        %                          subarrays. Each column in ElementWeights
        %                          specifies the weights for the elements 
        %                          in the corresponding subarray.
        %
        %                          For a PartitionedArray, if its
        %                          individual subarrays have the same
        %                          number of elements, ElementWeightsmust
        %                          be an NSExN matrix where NSE is the
        %                          number of elements in each individual
        %                          subarray and N is the number of
        %                          subarrays. Each column in WS specifies
        %                          the weights for the elements in the
        %                          corresponding subarray. 
        %
        %                          If a PartitionedArray's subarrays can
        %                          have different number of elements,
        %                          ElementWeights can be either an NSExN
        %                          matrix, where NSE indicates the number
        %                          of elements in the largest subarray and
        %                          N is the number of subarrays, or a 1xN
        %                          cell array, where N is the number of
        %                          subarrays and each cell contains a
        %                          column vector whose length is the same
        %                          as the number of elements of the
        %                          corresponding subarray.  If WS is a
        %                          matrix, the first K entries in each
        %                          column, where K is the number of
        %                          elements in the corresponding subarray,
        %                          specifies the weights for the elements
        %                          in the corresponding subarray. If WS is
        %                          a cell array, each cell in the array is 
        %                          a column vector specifying the weights 
        %                          for the elements in the corresponding
        %                          subarray. 
        %           AzimuthAngles: A row vector to specify the Azimuthal
        %                          angles and the spacing between these 
        %                          angles over which you can visualize
        %                          the radiation pattern. This parameter is
        %                          only applicable when the RespCut is 'Az'
        %                          or '3D' and the Format is 'Line' or 
        %                          'Polar'. The range of the Azimuthal 
        %                          angles specified should be between -180  
        %                          and 180. The default value is -180:180.
        %         ElevationAngles: A row vector to specify the Elevation
        %                          angles and the spacing between these 
        %                          angles over which you can visualize the
        %                          radiation pattern. This parameter is 
        %                          only applicable when the RespCut is 'El'
        %                          or '3D' and the Format is 'Line' or 
        %                          'Polar'. The range of the Elevation 
        %                          angles specified should be between -90 
        %                          and 90. The default value is -90:90.
        %                   UGrid: A row vector to specify the U 
        %                          coordinates in the U/V Space over which 
        %                          you want to visualize the radiation 
        %                          pattern. This parameter is applicable 
        %                          only when the Format is 'UV'. The range 
        %                          of the U coordinates specified should be
        %                          between -1 and 1. The default value is 
        %                          -1:0.01:1.
        %                   VGrid: A row vector specifying the V 
        %                          coordinates in the U/V Space over which 
        %                          you want to visualize the radiation 
        %                          pattern. This parameter is applicable 
        %                          only when the Format is 'UV' and the 
        %                          RespCut is '3D'. The range of the V 
        %                          coordinates specified should be between 
        %                          -1 and 1. The default value is
        %                          -1:0.01:1.
        %
        %   % Example:
        %   %   Plot the azimuth cut response of a replicated subarray
        %   %   along 0 elevation using a line plot. Assume the operating
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.ReplicatedSubarray;
        %   plotResponse(ha,3e8,physconst('LightSpeed'))
        %
        %   See also phased, phased.ArrayResponse.

            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','plotResponse');
            end
            
            narginchk(3,inf);
            validateattributes(obj,{'phased.internal.AbstractSubarray'},...
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
            
        end
    end
    
    methods
        
        function y = collectPlaneWave(obj,x,ang,freq,c)
        %collectPlaneWave Simulate received plane waves
        %   Y = collectPlaneWave(H,X,ANGLE) returns the received signals at
        %   the sensor array H, formed by subarrays, when the input signals
        %   X arrive at the array from the directions specified in ANGLE.
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
        %   Y is an N column matrix where N is the number of subarrays in
        %   the array. Each column of Y is the received signal at the
        %   corresponding subarray with all incoming signals combined.
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
        %   the array and only models the array factor among subarrays.
        %   Therefore, whether the subarray is steered or not does not
        %   affect the result.
        %
        %   % Example:
        %   %   Simulate the received signal at a 16-element ULA, which is
        %   %   formed by four 4-element ULAs, when two signals are
        %   %   arriving from 10 degrees and 30 degrees azimuth. Both
        %   %   signals have an elevation angle of 0 degrees. Assume the
        %   %   propagation speed is the speed of light and the carrier
        %   %   frequency of the signal is 100 MHz.
        %
        %   ha = phased.ReplicatedSubarray('Subarray',phased.ULA(4),...
        %                       'GridSize',[4 1]);
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
            
            if isempty(coder.target)
                privArray = clone(obj);
                release(privArray);
            else
                privArray = clonecg(obj);
            end
            privHstv = phased.SteeringVector('SensorArray',privArray,...
                'PropagationSpeed',c,'IncludeElementResponse',false);
            sv = step(privHstv,freq,angv);
            y = x*sv.';
        end
        
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
        %     ShowIndex: A vector specifying the element indices to show
        %                in the figure. You can also specify 'All' to show
        %                indices of all elements of the array or 'None' to
        %                suppress indices. The default value is 'None'.
        %
        %  ShowSubarray: A vector specifying the indices of subarrays that
        %                are highlighted in the figure. You can also
        %                specify 'All' to highlight all subarrays in the
        %                array or 'None' to suppress the subarray
        %                highlighting. The default value is 'All'.
        %
        %         Title: A string specifying the title of the plot. The 
        %                default value is 'Array Geometry'. 
        %
        %   % Examples:
        %
        %   % Example 1:
        %   %   Display the geometry of a Replicated Subarray formed by two
        %   %   4x4 URAs.
        %
        %   ha = phased.URA(4);
        %   hsubarray = phased.ReplicatedSubarray('Subarray',ha,...
        %       'Layout','Custom','SubarrayPosition',[0 0 ;-2 2 ;0 0]);
        %   viewArray(hsubarray,'ShowNormals',true);
        %
        %   % Example 2:
        %   %   Display the geometry of an overlapped subarray.
        %
        %   ha = phased.ULA(6);
        %   hsubarray = phased.PartitionedArray('Array',ha,...
        %       'SubarraySelection',[1 1 1 1 0 0;0 0 1 1 1 1]);
        %   viewArray(hsubarray,'ShowNormals',true,'ShowSubarray','all');
        %   
        %   See also phased, phased.ArrayResponse
            
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','viewArray');
            end

            validateattributes(obj,{'phased.ReplicatedSubarray',...
                'phased.PartitionedArray'},...
                {'scalar'},'viewArray','H');
            
            ShowTaper     = false;
            ShowNormals   = false;
            ShowIndex     = 'None';
            ShowSubarray  = 'All';
            Title         = 'Array Geometry';
            elem_pos      = obj.getElementPosition ;
            N             = obj.getNumElements ;
            N_subArr      = obj.getNumSubarrays  ;
            
            if isa(obj,'phased.ReplicatedSubarray')
                subArr = obj.Subarray;
            else
                subArr = obj.Array;
            end
            
            NsubArr_elem = subArr.getNumElements;
            sigutils.pvparse(varargin{:});
            
            defaultPosRGB = [0.04 0.51 0.78];
            defaultNormalRGB = [0.84 0.16 0]; %#ok<NASGU>
                       
            if ~isempty(ShowSubarray)
                if ischar(ShowSubarray)
                    ShowSubarray = validatestring(ShowSubarray,...
                        {'All','None'},'viewArray','ShowSubarray');
                    if strcmp(ShowSubarray,'All')
                        ShowSubarray = 1:N_subArr;
                    else
                        ShowSubarray = [];
                    end
                else
                    sigdatatypes.validateIndex(ShowSubarray,'viewArray',...
                        'ShowSubarray',{'vector','<=',N_subArr});
                end
            end
            
            if ShowNormals
                validateattributes(ShowNormals,{'logical'},{'scalar' },...
                    'viewArray','ShowNormals'); %#ok<UNRCH>
            end

            if ShowTaper
                validateattributes(ShowTaper,{'logical'},{'scalar' },...
                    'viewArray','ShowTaper'); %#ok<UNRCH>
                validatePropertiesImpl(subArr);
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
            % ClnupObj1     = onCleanup(@()hold('off'));
            if ~hold_status
                clf;
                hold on;
            end
            
            % compute element normal direction
            normal_pos_sph = deg2rad(getElementNormal(subArr));
            [normalX, normalY, normalZ] = sph2cart(...
                normal_pos_sph(1,:),normal_pos_sph(2,:),1);
            if isa(obj,'phased.ReplicatedSubarray')  
                % for replicating conformal array
                normalX = repmat(normalX,1,N_subArr);
                normalY = repmat(normalY,1,N_subArr);
                normalZ = repmat(normalZ,1,N_subArr);
            end
            
            % compute subarray normal direction
            if isa(obj,'phased.ReplicatedSubarray')
                norm_subArr  = obj.SubarrayNormal;
                if any(any(norm_subArr))
                    ang = zeros(2,N);
                    for n = 1:N_subArr
                        ang(:,(n-1)*NsubArr_elem+1:n*NsubArr_elem) = ...
                            repmat(norm_subArr(:,n),1,NsubArr_elem);
                    end
                    normal = [normalX' normalY' normalZ'];
                    for k=1:N
                        normal(k,:) = phased.internal.rotazel(...
                            normal(k,:)',ang(:,k))';
                        normalX     = normal(:,1)'; %#ok<NASGU>
                        normalY     = normal(:,2)'; %#ok<NASGU>
                        normalZ     = normal(:,3)'; %#ok<NASGU>
                    end
                end
            end
            
            % compute axes span
            elem_pos_max   = max(elem_pos,[],2);
            elem_pos_min   = min(elem_pos,[],2);
            
            axis_length = max(elem_pos_max - elem_pos_min);
            
            XPos           = axis_length/2;
            YPos           = axis_length/2;
            ZPos           = axis_length/2;
            
            % set marker size and font size
            thresh = (sum((elem_pos(:,1)-elem_pos(:,2)).^2))/(2*YPos);
            
            if thresh<0.012
                sz_marker = 36;
                sz_font = 6;
            elseif thresh>0.21
                sz_marker = 108;
                sz_font = 8;
            else
                sz_marker = 72;
                sz_font   = 8;
            end
            
            % subarray highlight color order
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
            
            
            
            % plot normal directions
            if ShowNormals
                hQuiver = quiver3(...
                    elem_pos(1,:)',elem_pos(2,:)',elem_pos(3,:)',...
                    normalX',normalY',normalZ',0.5,...
                    'Color',[0.84 0.16 0],...
                    'LineWidth',1,...
                    'Tag','Normals',...
                    'MaxHeadSize',0.1); %#ok<UNRCH>
                set(get(get(hQuiver,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off')
            end
            
            % show indices
            if ~isempty(ShowIndex)
                min_xdis = min(sqrt(sum((elem_pos-circshift(...
                    elem_pos,[0 -1])).^2)));
                for h = 1:length(ShowIndex)
                    text(elem_pos(1,ShowIndex(h))+min_xdis/3.5,...
                        elem_pos(2,ShowIndex(h))+min_xdis/3.5,...
                        elem_pos(3,ShowIndex(h))+min_xdis/3.5,...
                    num2str(ShowIndex(h)),'FontSize',sz_font);
                end
            end
            
            % highlight subarray
            
            SubarrayColorMatrix = repmat(defaultPosRGB,...
                numel(elem_pos(1,:)),1);
                        
            if ~isempty(ShowSubarray)
                ShowSubarray = unique(ShowSubarray);
                for s = 1:length(ShowSubarray)
                    subarray_idx= ShowSubarray(s);
                    ColorOrderIdx = mod(s-1,length(ColorOrder))+1;
                                       
                    if isa(obj,'phased.PartitionedArray')
                        select_sub = obj.SubarraySelection;
                        SubarrayColorMatrix(logical(select_sub(subarray_idx,:)),:) = ...
                            repmat(ColorOrder(ColorOrderIdx,:),...
                            numel(find(select_sub(subarray_idx,:))),1);
                    else
                        SubarrayColorMatrix((subarray_idx-1)*NsubArr_elem+1:...
                            subarray_idx*NsubArr_elem,:) = ...
                            repmat(ColorOrder(ColorOrderIdx,:),NsubArr_elem,1);
                    end
                    
                end
                
                % if PartitionedArray, identify overlapped elements and
                % highlight it with white
                ovrLapIdx = [];
                
                if isa(obj,'phased.PartitionedArray')
                    selected_sub = obj.SubarraySelection(ShowSubarray,:);
                    ovrLapIdx = find(sum(logical(selected_sub),1)>1);
                end
                
                if ~isempty(ovrLapIdx)
                    SubarrayColorMatrix(ovrLapIdx,:) = [];
                    ColorMatrixOverlap = repmat([1 1 1],numel(ovrLapIdx),1);
                 
                    ovrLap_elem_pos = elem_pos(:,ovrLapIdx);
                    elem_pos(:,ovrLapIdx) = [];
                    hOverLap = scatter3(ovrLap_elem_pos(1,:)',...
                        ovrLap_elem_pos(2,:)',ovrLap_elem_pos(3,:)',...
                        sz_marker,ColorMatrixOverlap,'filled',...
                        'Tag','Overlapped Elements',...
                        'MarkerEdgeColor','black'); 

                    hoverlaplegend = legend('Overlap','Location','NorthEast','AutoUpdate','off');
                    set(hoverlaplegend,'Tag','Overlap_Legend');
                end
                
            end

            if ShowTaper
                lTaper = abs(getTaper(subArr)); %#ok<UNRCH>
                normTaper = lTaper/max(lTaper);
                if isa(obj,'phased.PartitionedArray') 
                     if ~isempty(ovrLapIdx)
                        normTaper(ovrLapIdx) = [];
                     end
                    SubarrayColorMatrix = SubarrayColorMatrix.*repmat(normTaper,1,3);
                else
                    for subarray_idx = 1:N_subArr
                        SubarrayColorMatrix((subarray_idx-1)*NsubArr_elem+1:...
                                    subarray_idx*NsubArr_elem,:) = ...
                            SubarrayColorMatrix((subarray_idx-1)*NsubArr_elem+1:...
                                    subarray_idx*NsubArr_elem,:) ...
                                    .*repmat(normTaper,1,3);
                    end
                end
            end
            
            % show non-overlapped elements
            hElements = scatter3(elem_pos(1,:)',elem_pos(2,:)',...
                elem_pos(3,:)',sz_marker,SubarrayColorMatrix,'filled',...
                'Tag','Elements','MarkerEdgeColor','black');
            set(hElements,'DeleteFcn',{@deleteAssociatedText});            
            set(get(get(hElements,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
           
            if nargout
                if ~isempty(ovrLapIdx)
                    varargout{1} = [hElements; hOverLap];
                else
                    varargout{1} = hElements;      
                end
            end
            
            origFormat = get(0, 'format');
            format('shortG');
                        
            span = abs(elem_pos_max - elem_pos_min);
            [spanv, spanu] = convert2engstrs(span,3);
            arrayspan_str = {...
                   ['Array Span:         '],...
                   ['    X axis = ' spanv(1,:) ' ' spanu 'm'],...
                   ['    Y axis = ' spanv(2,:) ' ' spanu 'm'],...
                   ['    Z axis = ' spanv(3,:) ' ' spanu 'm']};
            annotatebox = annotation('textbox',[0.7 0 0.2 0.25],'String',arrayspan_str,...
                'FitBoxToText','on','LineStyle','none','FontSize',8,...
                'HorizontalAlignment','right','Tag','Annotation',...
                'Color',[0.501 0.501 0.501]);
            setappdata(hElements,'annotatebox',annotatebox);
            axis_scale = 1.15;
            elem_pos_min = elem_pos_min*axis_scale-axis_length/10;
            
            hxaxis = line([elem_pos_min(1),(XPos/4)+elem_pos_min(1)],...
                [elem_pos_min(2),elem_pos_min(2)],...
                [elem_pos_min(3),elem_pos_min(3)],...
                'Color','k','LineWidth',1,'Tag','X Axis','Clipping','off');
            hyaxis = line([elem_pos_min(1),elem_pos_min(1)],...
                [elem_pos_min(2),(YPos/4)+elem_pos_min(2)],...
                [elem_pos_min(3),elem_pos_min(3)],...
                'Color','k','LineWidth',1,'Tag','Y Axis','Clipping','off');
            hzaxis = line([elem_pos_min(1),elem_pos_min(1)],...
                [elem_pos_min(2),elem_pos_min(2)],...
                [elem_pos_min(3),(ZPos/4)+elem_pos_min(3)],...
                'Color','k','LineWidth',1,'Tag','Z Axis','Clipping','off');
            set(get(get(hxaxis,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
            set(get(get(hyaxis,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
            set(get(get(hzaxis,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
            
            text(1.11*(XPos/4)+elem_pos_min(1),elem_pos_min(2),...
                elem_pos_min(3),'x','Tag','X Label','FontSize',8);
            text(elem_pos_min(1),1.11*(YPos/4)+elem_pos_min(2),...
                elem_pos_min(3),'y','Tag','Y Label','FontSize',8);
            text(elem_pos_min(1),elem_pos_min(2),...
                1.12*(ZPos/4)+elem_pos_min(3),'z','Tag',...
                'Z Label','FontSize',8);
            
            figureRGB = [0.94 0.94 0.94]; 
            axisRGB = [0.49 0.49 0.49]; 
            set(gcf,'Color',figureRGB);
            set(gca,'XColor',axisRGB,...
                'YColor',axisRGB,...
                'ZColor',axisRGB,...
                'Color',figureRGB,...
                'Position',[0 0 1 1]);
            xlabel('x axis (Az 0 El 0) -->','Color',axisRGB);
            ylabel('y axis-->','Color',axisRGB);
            zlabel('z axis (Az 0 El 90) --->','Color',axisRGB);

            view(45,45);
            axis vis3d;
            box off;
            grid off;
            axis off;
            hold off;
            daspect([1 1 1])
            camlight('infinite');

            if ~ShowNormals
                % Show planar arrays in their plane
                if span(1)==0 && span(2)~=0 && span(3)~=0
                    view(90,0),
                elseif span(2)==0 && span(1)~=0 && span(3)~=0
                    view(0,0)
                elseif span(3)==0 && span(1)~=0 && span(2)~=0
                    %view(0,90) will switch back to 2D view
                    % see g1130105
                    view(0,89.9); 
                end
            end
            
            titlebox = annotation('textbox',[0 0.75 1 0.2],...
                'String',Title,'FitBoxToText','off',...
                'FontWeight','demi',...
                'LineStyle','none','HitTest','off',...
                'HorizontalAlignment','center','Tag','Title');
            setappdata(hElements,'titlebox',titlebox);
             
            % clean up
            format(origFormat);
            if hold_status
                hold on;
            else
                hold off;
            end
         end
    end
    
    methods (Hidden)
        function num = getDOF(obj)
        %getDOF  Degree of freedom of the array
        %   N = getDOF(H) returns the degree of freedom (DOF), N, of the
        %   array H.
        %
        %   % Example:
        %   %   Construct an array with replicated subarrays and then
        %   %   obtain its degree of freedom.
        %
        %   ha = phased.ReplicatedSubarray;
        %   dof = getDOF(ha)
            
            num = getNumSubarrays(obj);
        end
    end
    
    methods (Access = protected)
        
        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end
                
        function flag = isInputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.isLocked = isLocked(obj);
            if isLocked(obj)
                s.pSubarrayPosition = obj.pSubarrayPosition;
                % s.pPropagationSpeed = obj.pPropagationSpeed;
                s.pSubarrayFeed = obj.pSubarrayFeed;
                s.pFeedWeights = obj.pFeedWeights;
                s.cArray = saveobj(obj.cArray);
                s.cModuleResponse = saveobj(obj.cModuleResponse);
                s.cModuleSteeringVector = saveobj(obj.cModuleSteeringVector);
                s.pNumSubarrays = obj.pNumSubarrays;
                s.pNumElements = obj.pNumElements;
                s.pNumAngles = obj.pNumAngles;
                s.pNumFreqs = obj.pNumFreqs;
                s.pNoSteeringFlag = obj.pNoSteeringFlag;
                s.pPhaseSteeringFlag = obj.pPhaseSteeringFlag;
                s.pTimeSteeringFlag = obj.pTimeSteeringFlag;
                s.pIsSingleFrequency = obj.pIsSingleFrequency;
                s.pIsSingleWeights = obj.pIsSingleWeights;
            end
        end
        
        function s = loadSubObjects(obj,s)
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cArray = eval(...
                        sprintf('%s.loadobj(s.cArray)',s.cArray.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cArray');
                    obj.cModuleResponse = eval(...
                        sprintf('%s.loadobj(s.cModuleResponse)',s.cModuleResponse.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cModuleResponse');
                    obj.cModuleSteeringVector = eval(...
                        sprintf('%s.loadobj(s.cModuleSteeringVector)',s.cModuleSteeringVector.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cModuleSteeringVector');
                end
            end
        end
        
        function setupImpl(obj,freq,ang,~,~)
            obj.pNumSubarrays = calcNumSubarrays(obj);
            obj.pNumElements = calcNumElements(obj);
            obj.pSubarrayPosition = calcSubarrayPosition(obj);
            
            obj.pNumAngles = size(ang,2);
            obj.pNumFreqs = size(freq,2);
            obj.pIsSingleFrequency = (obj.pNumFreqs==1);
            if ~strncmp(obj.SubarraySteering,'None',1)
                obj.pNoSteeringFlag = false;
                if strncmp(obj.SubarraySteering,'Phase',1)
                    obj.pPhaseSteeringFlag = true;
                    obj.pTimeSteeringFlag = false;
                    obj.pIsSingleWeights = true;
                elseif strncmp(obj.SubarraySteering,'Time',1)
                    obj.pTimeSteeringFlag = true;
                    obj.pPhaseSteeringFlag = false;
                    if ~obj.pIsSingleFrequency
                        obj.pIsSingleWeights = false;
                    end
                else % Custom
                    obj.pPhaseSteeringFlag = false;
                    obj.pTimeSteeringFlag = false;
                    obj.pIsSingleWeights = true; 
                    % only allow one set of weights to specify at element
                    % level because this corresponds to the analog system
                end
            else
                obj.pNoSteeringFlag = true;
                obj.pTimeSteeringFlag = false;
                obj.pPhaseSteeringFlag = false;
                obj.pIsSingleWeights = true;
            end
            
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if (strcmp(prop,'PhaseShifterFrequency') || ...
                    strcmp(prop,'NumPhaseShifterBits')) && ...
                    (~strncmp(obj.SubarraySteering,'Phase',1))
                flag = true;
            end
        end
        
        function releaseImpl(obj)
            %releaseImpl   Free any resources (such as file handles,
            %device drivers, etc) that are used by the object and perform
            %any other end-of-use operations.
            release(obj.cArray);
            release(obj.cModuleResponse);
            if ~obj.pNoSteeringFlag
              release(obj.cModuleSteeringVector);
            end
            
        end
        
        function resetImpl(obj)
        %resetImpl   Initialize/reset any state to its initial values
            reset(obj.cModuleResponse);
            if ~obj.pNoSteeringFlag
              reset(obj.cModuleSteeringVector);
            end
        end
        
        function validateInputsImpl(obj,freq,ang,c,stang)
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
                coder.internal.errorIf(cond,'phased:system:array:expectReal','FREQ');
            end
            
            
            sz_ang = size(ang);
            cond =  ~ismatrix(ang) || isempty(ang);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','ANGLE');
            end
            cond =  sz_ang(1) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:NeedTwoRows','ANGLE');
            end
            cond =  ~isreal(ang);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidAngle','ANGLE');
            end
            cond =  ~isa(ang,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','ANGLE','double');
            end
            
            cond =  ~isscalar(c) || isempty(c);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeScalar','V');
            end
            cond =  ~isa(c,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','V','double');
            end
            cond =  ~isreal(c);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:expectReal','V');
            end
                
            if ~strncmp(obj.SubarraySteering,'None',1)
                if ~strncmp(obj.SubarraySteering,'Custom',1)
                    sz_stang = size(stang);
                    cond =  ~ismatrix(stang) || isempty(stang);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:inputMustBeMatrix','STEERANGLE');
                    end
                    cond =  sz_stang(1) > 2;
                    if cond
                        coder.internal.errorIf(cond,'phased:system:array:NeedTwoRows','STEERANGLE');
                    end
                    cond =  sz_stang(2) > 1;
                    if cond
                        coder.internal.errorIf(cond,'phased:system:array:NeedOneColumn','STEERANGLE');
                    end
                    cond =  ~isreal(stang);
                    if cond
                        coder.internal.errorIf(cond,'phased:system:array:InvalidAngle','STEERANGLE');
                    end
                    cond =  ~isa(stang,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:invalidInputDataType','STEERANGLE','double');
                    end
                else
                    ws = stang;  % weights
                    % does not support multiple weights for multiple
                    % frequency yet because this is still an analog
                    % behavior so at any moment, there is only one set of
                    % weights can be applied.
                    validateElementWeights(obj,ws);
                        
                end

            end
        end
        
        
        function num = getNumInputsImpl(obj)
            %getNumInputsImpl   Give the number of separate input signals
            %that is expected by mUpdate. Can be removed if num is always
            %1.
            
            num = 3;
            if ~strcmp(obj.SubarraySteering,'None')
                num = num+1;
            end
        end
               
    end
    
    methods (Abstract, Hidden, Access = {?phased.gpu.internal.AbstractClutterSimulator, ?phased.internal.AbstractSubarray})
        %Methods used in a GPU ConstantGammaClutter model simulation.
        
        [pos, az, el] = getSubarrayPosAzEl(obj, subArrIdx, azin, elin)     
        sv = getSubArraySteeringVec(obj, subArrIdx, freq, c, stang );
    end
    
    methods 
        
        function varargout = pattern(obj,freq,varargin)
        %pattern Plot subarray response pattern
        %   pattern(H,FREQ) plots the 3D subarray directivity pattern (in
        %   dBi). The operating frequency is specified in FREQ (in Hz) as a
        %   positive scalar.
        %
        %   pattern(H,FREQ,AZ) specifies the azimuth angles (in degrees) as
        %   a vector in AZ and plots the 3D subarray directivity pattern
        %   (in dBi) within the specified azimuth angles.
        %
        %   pattern(H,FREQ,AZ,EL) specifies the elevation angles (in
        %   degrees) as a vector in EL and plots the subarray directivity
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
        %   pattern(H,FREQ,AZ,EL,Name,Value) plots the subarray response
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
        %                          the array. N is the number of subarrays
        %                          in the array and M is the number of
        %                          frequencies supplied in FREQ if FREQ is 
        %                          a vector. If Weights is a column vector,
        %                          the same weights will be applied to each 
        %                          frequency. If Weights is a matrix,each
        %                          column of weight values will be applied
        %                          to the corresponding frequency in FREQ.
        %                          If Weights is a matrix and FREQ is a
        %                          scalar, each column of weight values
        %                          will be applied to the same frequency in
        %                          FREQ.           
        %              SteerAngle: A scalar or a length-2 column vector to
        %                          specify the subarray steering angle. If
        %                          SteerAng is a vector, it is in the form
        %                          of [az; el] (in degrees). az must be
        %                          between -180 and 180 while el must be
        %                          between -90 and 90. If SteerAng is a
        %                          scalar, it specifies the azimuth angle
        %                          and the elevation angle is assumed to be
        %                          0. The default value is [0;0]. This
        %                          option is only applicable if the
        %                          SubarraySteering property of H is set to
        %                          'Phase' or 'Time'.
        %          ElementWeights: A matrix or a cell array specifying
        %                          the weights applied to each element
        %                          in the subarray. This parameter is
        %                          applicable when you set the
        %                          SubarraySteering property of H to
        %                          'Custom'. The default value is all ones.
        %                          
        %                          For a ReplicatedSubarray, ElementWeights
        %                          must be a NSExN matrix where NSE is the 
        %                          number of elements in each individual 
        %                          subarray and N is the number of 
        %                          subarrays. Each column in ElementWeights
        %                          specifies the weights for the elements 
        %                          in the corresponding subarray.
        %
        %                          For a PartitionedArray, if its
        %                          individual subarrays have the same
        %                          number of elements, ElementWeightsmust
        %                          be an NSExN matrix where NSE is the
        %                          number of elements in each individual
        %                          subarray and N is the number of
        %                          subarrays. Each column in WS specifies
        %                          the weights for the elements in the
        %                          corresponding subarray. 
        %
        %                          If a PartitionedArray's subarrays can
        %                          have different number of elements,
        %                          ElementWeights can be either an NSExN
        %                          matrix, where NSE indicates the number
        %                          of elements in the largest subarray and
        %                          N is the number of subarrays, or a 1xN
        %                          cell array, where N is the number of
        %                          subarrays and each cell contains a
        %                          column vector whose length is the same
        %                          as the number of elements of the
        %                          corresponding subarray.  If WS is a
        %                          matrix, the first K entries in each
        %                          column, where K is the number of
        %                          elements in the corresponding subarray,
        %                          specifies the weights for the elements
        %                          in the corresponding subarray. If WS is
        %                          a cell array, each cell in the array is 
        %                          a column vector specifying the weights 
        %                          for the elements in the corresponding
        %                          subarray. 
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
        %   %   Plot the azimuth cut response of a replicated subarray
        %   %   along 0 elevation using a line plot. Assume the operating
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.ReplicatedSubarray;
        %   pattern(ha,3e8,-180:180,0,'CoordinateSystem','rectangular');
        %
        %   See also phased, phased.ArrayResponse.
            
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                    'phased:Waveform:CodegenNotSupported','pattern');
            end
            
            narginchk(2,inf);
            nargoutchk(0,3);
            validateattributes(obj,...
                {'phased.internal.AbstractSubarray'},...
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
        %              SteerAngle: A scalar or a length-2 column vector to
        %                          specify the subarray steering angle. If
        %                          SteerAng is a vector, it is in the form
        %                          of [az; el] (in degrees). az must be
        %                          between -180 and 180 while el must be
        %                          between -90 and 90. If SteerAng is a
        %                          scalar, it specifies the azimuth angle
        %                          and the elevation angle is assumed to be
        %                          0. The default value is [0;0]. This
        %                          option is only applicable if the
        %                          SubarraySteering property of H is set to
        %                          'Phase' or 'Time'.
        %          ElementWeights: A matrix or a cell array specifying
        %                          the weights applied to each element
        %                          in the subarray. This parameter is
        %                          applicable when you set the
        %                          SubarraySteering property of H to
        %                          'Custom'. The default value is all ones.
        %                          
        %                          For a ReplicatedSubarray, ElementWeights
        %                          must be a NSExN matrix where NSE is the 
        %                          number of elements in each individual 
        %                          subarray and N is the number of 
        %                          subarrays. Each column in ElementWeights
        %                          specifies the weights for the elements 
        %                          in the corresponding subarray.
        %
        %                          For a PartitionedArray, if its
        %                          individual subarrays have the same
        %                          number of elements, ElementWeightsmust
        %                          be an NSExN matrix where NSE is the
        %                          number of elements in each individual
        %                          subarray and N is the number of
        %                          subarrays. Each column in WS specifies
        %                          the weights for the elements in the
        %                          corresponding subarray. 
        %
        %                          If a PartitionedArray's subarrays can
        %                          have different number of elements,
        %                          ElementWeights can be either an NSExN
        %                          matrix, where NSE indicates the number
        %                          of elements in the largest subarray and
        %                          N is the number of subarrays, or a 1xN
        %                          cell array, where N is the number of
        %                          subarrays and each cell contains a
        %                          column vector whose length is the same
        %                          as the number of elements of the
        %                          corresponding subarray.  If WS is a
        %                          matrix, the first K entries in each
        %                          column, where K is the number of
        %                          elements in the corresponding subarray,
        %                          specifies the weights for the elements
        %                          in the corresponding subarray. If WS is
        %                          a cell array, each cell in the array is 
        %                          a column vector specifying the weights 
        %                          for the elements in the corresponding
        %                          subarray. 
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
        %   %   Plot the azimuth cut response of a partitioned array
        %   %   along 0 and 10 degrees elevation. Assume the operating 
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.PartitionedArray;
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
                {'phased.internal.AbstractSubarray'},...
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
        %   along elevation at 0 degree azimuth angle (in degrees). The
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
        %              SteerAngle: A scalar or a length-2 column vector to
        %                          specify the subarray steering angle. If
        %                          SteerAng is a vector, it is in the form
        %                          of [az; el] (in degrees). az must be
        %                          between -180 and 180 while el must be
        %                          between -90 and 90. If SteerAng is a
        %                          scalar, it specifies the azimuth angle
        %                          and the elevation angle is assumed to be
        %                          0. The default value is [0;0]. This
        %                          option is only applicable if the
        %                          SubarraySteering property of H is set to
        %                          'Phase' or 'Time'.
        %          ElementWeights: A matrix or a cell array specifying
        %                          the weights applied to each element
        %                          in the subarray. This parameter is
        %                          applicable when you set the
        %                          SubarraySteering property of H to
        %                          'Custom'. The default value is all ones.
        %                          
        %                          For a ReplicatedSubarray, ElementWeights
        %                          must be a NSExN matrix where NSE is the 
        %                          number of elements in each individual 
        %                          subarray and N is the number of 
        %                          subarrays. Each column in ElementWeights
        %                          specifies the weights for the elements 
        %                          in the corresponding subarray.
        %
        %                          For a PartitionedArray, if its
        %                          individual subarrays have the same
        %                          number of elements, ElementWeightsmust
        %                          be an NSExN matrix where NSE is the
        %                          number of elements in each individual
        %                          subarray and N is the number of
        %                          subarrays. Each column in WS specifies
        %                          the weights for the elements in the
        %                          corresponding subarray. 
        %
        %                          If a PartitionedArray's subarrays can
        %                          have different number of elements,
        %                          ElementWeights can be either an NSExN
        %                          matrix, where NSE indicates the number
        %                          of elements in the largest subarray and
        %                          N is the number of subarrays, or a 1xN
        %                          cell array, where N is the number of
        %                          subarrays and each cell contains a
        %                          column vector whose length is the same
        %                          as the number of elements of the
        %                          corresponding subarray.  If WS is a
        %                          matrix, the first K entries in each
        %                          column, where K is the number of
        %                          elements in the corresponding subarray,
        %                          specifies the weights for the elements
        %                          in the corresponding subarray. If WS is
        %                          a cell array, each cell in the array is 
        %                          a column vector specifying the weights 
        %                          for the elements in the corresponding
        %                          subarray. 
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
        %   %   Plot the elevation cut response of a replicated subarray
        %   %   along 0 and 10 degrees azimuth. Assume the operating 
        %   %   frequency is 300 MHz.
        %
        %   ha = phased.ReplicatedSubarray;
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
                {'phased.internal.AbstractSubarray'},...
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
        %             SteerAngle: A scalar or a length-2 column vector 
        %                         specifying the subarray steering angle
        %                         (in degrees). If it is a vector, the
        %                         angle is specified in the form of
        %                         [AzimuthAngle; ElevationAngle]. If it is
        %                         a scalar, it represents the azimuth angle
        %                         and the elevation angle is assumed to be
        %                         0. This parameter is applicable when you
        %                         set the SubarraySteering property of H to
        %                         either 'Phase' or 'Time'. The default
        %                         value is [0;0].
        %         ElementWeights: A matrix or a cell array specifying
        %                         the weights applied to each element
        %                         in the subarray. This parameter is
        %                         applicable when you set the
        %                         SubarraySteering property of H to
        %                         'Custom'. The default value is all ones.
        %                          
        %                         For a ReplicatedSubarray, ElementWeights
        %                         must be a NSExN matrix where NSE is the 
        %                         number of elements in each individual 
        %                         subarray and N is the number of 
        %                         subarrays. Each column in ElementWeights
        %                         specifies the weights for the elements in
        %                         the corresponding subarray.
        %
        %                          For a PartitionedArray, if its
        %                          individual subarrays have the same
        %                          number of elements, ElementWeightsmust
        %                          be an NSExN matrix where NSE is the
        %                          number of elements in each individual
        %                          subarray and N is the number of
        %                          subarrays. Each column in WS specifies
        %                          the weights for the elements in the
        %                          corresponding subarray. 
        %
        %                          If a PartitionedArray's subarrays can
        %                          have different number of elements,
        %                          ElementWeights can be either an NSExN
        %                          matrix, where NSE indicates the number
        %                          of elements in the largest subarray and
        %                          N is the number of subarrays, or a 1xN
        %                          cell array, where N is the number of
        %                          subarrays and each cell contains a
        %                          column vector whose length is the same
        %                          as the number of elements of the
        %                          corresponding subarray.  If WS is a
        %                          matrix, the first K entries in each
        %                          column, where K is the number of
        %                          elements in the corresponding subarray,
        %                          specifies the weights for the elements
        %                          in the corresponding subarray. If WS is
        %                          a cell array, each cell in the array is 
        %                          a column vector specifying the weights 
        %                          for the elements in the corresponding
        %                          subarray. 
        %
        %   % Examples:
        %
        %   % Example 1:
        %   %   Compute the boresight directivity of a subarray formed by 
        %   %   aligning two 10-element, quarter-wavelength spacing ULAs 
        %   %   side by side. The element is assumed to be isotropic.
        %
        %   c = 3e8; fc = 3e8; lambda = c/fc;
        %   myArray = phased.ReplicatedSubarray('Subarray',...
        %       phased.ULA(10,lambda/4));
        %   d = directivity(myArray,fc,0,'PropagationSpeed',c)
        %
        %   % Example 2:
        %   %   Compute the 30 degrees azimuth directivity of a subarray 
        %   %   formed by aligning two 10-element, quarter-wavelength 
        %   %   spacing ULAs side by side. The subarray is phase steered 
        %   %   toward 30 degrees azimuth.
        %
        %   c = 3e8; fc = 3e8; lambda = c/fc; ang = [30;0];
        %   myArray = phased.PartitionedArray('Array',...
        %       phased.ULA(20,lambda/4),'SubarraySelection',...
        %       [ones(1,10) zeros(1,10);zeros(1,10) ones(1,10)],...
        %       'SubarraySteering','Phase','PhaseShifterFrequency',fc);
        %   myStv = phased.SteeringVector('SensorArray',myArray,...
        %       'PropagationSpeed',c);
        %   d = directivity(myArray,fc,ang,'PropagationSpeed',c,...
        %           'Weights',step(myStv,fc,ang),'SteerAngle',ang)
        %
        %   See also phased, phased.ArrayResponse.
        

            phased.internal.narginchk(3,9,nargin);
            
            defaultPropagationSpeed = physconst('lightspeed');
            defaultWeights = ones(getDOF(obj),1);
            if ~strcmp(obj.SubarraySteering,'None')
                if strncmp(obj.SubarraySteering,'Phase',1) || ...
                        strncmp(obj.SubarraySteering,'Time',1)
                    defaultSteerAngle = [0;0];
                    if coder.target('MATLAB')
                        p = inputParser;
                        addParameter(p,'PropagationSpeed',defaultPropagationSpeed);
                        addParameter(p,'Weights',defaultWeights);
                        addParameter(p,'SteerAngle',defaultSteerAngle);

                        parse(p,varargin{:});
                        c = p.Results.PropagationSpeed;
                        w = p.Results.Weights;
                        stang = p.Results.SteerAngle;

                    else
                        parms = struct('PropagationSpeed',uint32(0), ...
                                    'Weights',uint32(0),...
                                    'SteerAngle',uint32(0));

                        pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
                        c = eml_get_parameter_value(pstruct.PropagationSpeed,...
                            defaultPropagationSpeed,varargin{:});
                        w = eml_get_parameter_value(pstruct.Weights,...
                            defaultWeights,varargin{:});
                        stang = eml_get_parameter_value(pstruct.SteerAngle,...
                            defaultSteerAngle,varargin{:});
                    end
                else
                    defaultElementWeights = ones(getNumElements(obj),1);
                    if coder.target('MATLAB')
                        p = inputParser;
                        addParameter(p,'PropagationSpeed',defaultPropagationSpeed);
                        addParameter(p,'Weights',defaultWeights);
                        addParameter(p,'ElementWeights',defaultElementWeights);

                        parse(p,varargin{:});
                        c = p.Results.PropagationSpeed;
                        w = p.Results.Weights;
                        ws = p.Results.ElementWeights;

                    else
                        parms = struct('PropagationSpeed',uint32(0), ...
                                    'Weights',uint32(0),...
                                    'ElementWeights',uint32(0));

                        pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
                        c = eml_get_parameter_value(pstruct.PropagationSpeed,...
                            defaultPropagationSpeed,varargin{:});
                        w = eml_get_parameter_value(pstruct.Weights,...
                            defaultWeights,varargin{:});
                        ws = eml_get_parameter_value(pstruct.ElementWeights,...
                            defaultElementWeights,varargin{:});
                    end
                end
            else
                if coder.target('MATLAB')
                    p = inputParser;
                    addParameter(p,'PropagationSpeed',defaultPropagationSpeed);
                    addParameter(p,'Weights',defaultWeights);

                    parse(p,varargin{:});
                    c = p.Results.PropagationSpeed;
                    w = p.Results.Weights;

                else
                    parms = struct('PropagationSpeed',uint32(0), ...
                                'Weights',uint32(0),...
                                'SteerAngle',uint32(0));

                    pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
                    c = eml_get_parameter_value(pstruct.PropagationSpeed,...
                        defaultPropagationSpeed,varargin{:});
                    w = eml_get_parameter_value(pstruct.Weights,...
                        defaultWeights,varargin{:});
                end
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
        
            if coder.target('MATLAB')
                myDirectivity = phased.internal.Directivity(...
                    'Sensor',clone(obj),...
                    'PropagationSpeed',c,'WeightsInputPort',true);
            else
                myDirectivity = phased.internal.Directivity(...
                    'Sensor',clonecg(obj),...
                    'PropagationSpeed',c,'WeightsInputPort',true);
            end
            
            if iscolumn(w)
                validateattributes(w,{'double'},{'nonempty','finite',...
                    'nrows',getDOF(obj)},'directivity','W');
            else
                validateattributes(w,{'double'},{'2d','nonempty','finite',...
                    'size',[getDOF(obj) numel(freq)]},'directivity','W');
            end
                
            if ~strcmp(obj.SubarraySteering,'None')
                if strncmp(obj.SubarraySteering,'Phase',1) || ...
                        strncmp(obj.SubarraySteering,'Time',1)
                    sigdatatypes.validateAngle(stang,'directivity',...
                        'STEERANGLE',{'column'});
                    if isscalar(stang)
                        stangIn = [stang;0];
                    else
                        stangIn = stang;
                    end
                    sigdatatypes.validateAzElAngle(stangIn,'directivity','ANGLE');

                    d = step(myDirectivity,freq,angIn,w,stangIn);
                else % Custom
                    Ns = getNumSubarrays(obj);
                    validateattributes(ws,{'double','cell'},...
                        {'nonempty','ncols',Ns},...
                        'directivity','ElementWeights');
                    Nse = zeros(1,Ns);
                    for m = 1:Ns
                        Nse(m) = calcNumElements(obj,m);
                    end
                    if iscell(ws)
                        for m = 1:Ns
                            cond = ~iscolumn(ws{m}) || (numel(ws{m})~=Nse(m));
                            if cond
                                coder.internal.errorIf(cond, ...
                                    'phased:system:array:SubarrayElementWeightsSizeMismatch',...
                                    m,'ElementWeights',Nse(m));
                            end
                            cond = ~isa(ws{m},'double');
                            if cond
                                coder.internal.errorIf(cond, ...
                                    'phased:system:array:SubarrayElementWeightsInvalidDataType',...
                                    m,'ElementWeights','double');
                            end
                        end
                    else
                        validateattributes(ws,{'double'},{'finite','size',[max(Nse) Ns]},...
                            'ElementWeights');
                    end
                    d = step(myDirectivity,freq,angIn,w,ws);
                end
            else
                
                d = step(myDirectivity,freq,angIn,w);
            end
            
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

