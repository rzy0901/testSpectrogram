classdef (Sealed,StrictDefaults) ConstantGammaClutter < phased.internal.AbstractConstantGammaClutter 
%ConstantGammaClutter   Constant gamma clutter simulation
%   H = phased.ConstantGammaClutter creates a constant gamma clutter
%   simulation System object, H. This object simulates the clutter return
%   of a monostatic radar system using constant gamma model.
%
%   H = phased.ConstantGammaClutter(Name,Value) creates a constant gamma
%   clutter object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H) returns an NxM matrix, Y, that contains the collected
%   clutter return. M is the number of radiating/collecting subarrays in
%   the radar system if Sensor contains subarrays, or the number of
%   radiating/collecting sensors otherwise. When you set the OutputFormat
%   property to 'Samples', N is specified in the NumSamples property. When
%   you set the OutputFormat property to 'Pulses', N is the total number of
%   samples in the next L pulses where L is specified in the NumPulses
%   property. This syntax is available when you set the
%   TransmitSignalInputPort property to false.
%
%   The clutter simulation makes the following assumptions:
%       * The radar system is monostatic
%       * The propagation is in free space
%       * The terrain is homogeneous
%       * The clutter patch is stationary during coherence time
%       * The signal is narrowband
%       * The radar system maintains a constant height during simulation
%       * The radar system maintains a constant speed during simulation
%
%   Y = step(H,X) specifies the transmit signal in X as a column vector. In
%   general, transmit signal refers to the signal at the output of the
%   transmitter when the transmitter is turned on during each pulse. This
%   syntax is available when you set the TransmitSignalInputPort property
%   to true.
%
%   Y = step(H,STEER) uses STEER as the subarray steering angle (in
%   degrees). STEER can be a scalar or a length-2 column vector. If STEER
%   is a vector, it is in the form of [AzimuthAngle; ElevationAngle]. If
%   STEER is a scalar, it represents the azimuth angle and the elevation
%   angle is assumed to be 0. This syntax is only applicable when you use
%   subarrays in the Sensor property and set the SubarraySteering property
%   in the Sensor to either 'Phase' or 'Time'.
%
%   Y = step(H,WS) uses WS as the weights applied to each element in the
%   subarray. WS can be either a matrix or a cell array. This syntax is
%   only applicable when you use subarrays in the Sensor property and set
%   the SubarraySteering property in the Sensor to 'Custom'.
%   
%   If the Sensor property is a phased.ReplicatedSubarray, WS must be an
%   NSExN matrix where NSE is the number of elements in each individual
%   subarray and N is the number of subarrays. Each column in WS specifies
%   the weights for the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its individual
%   subarrays have same number of elements, WS must be an NSExN matrix
%   where NSE is the number of elements in each individual subarray and N
%   is the number of subarrays. Each column in WS specifies the weights for
%   the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its subarrays
%   can have different number of elements, WS can be either an NSExN
%   matrix, where NSE indicates the number of elements in the largest
%   subarray and N is the number of subarrays, or a 1xN cell array, where N
%   is the number of subarrays and each cell contains a column vector whose
%   length is the same as the number of elements of the corresponding
%   subarray.  If WS is a matrix, the first K entries in each column, where
%   K is the number of elements in the corresponding subarray, specifies
%   the weights for the elements in the corresponding subarray. If WS is a
%   cell array, each cell in the array is a column vector specifying the
%   weights for the elements in the corresponding subarray. 
%
%   Y = step(H,PRFIDX) specifies the index of pulse repetition frequency
%   (PRF), PRFIDX, as a positive integer. The index is used to identify the
%   entries specified in the PRF property. This syntax applies when you set
%   the PRFSelectionInputPort property to true. Use this syntax for the
%   cases where the transmit pulse needs to be dynamically selected. Under
%   such situations, the PRF property includes a list of predetermined
%   choices of PRFs. During the simulation, based on PRF index input, one
%   of the PRFs is selected as the PRF for the next transmission.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,STEER,PRFIDX)
%
%   or
%
%   Y = step(H,X,WS,PRFIDX)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H) and y = H() are
%   equivalent.
%
%   ConstantGammaClutter methods:
%       step     - Simulate clutter using constant gamma model (see above)
%       release  - Allow property value and input characteristics changes
%       clone    - Create a constant gamma clutter object with same
%                  property values
%       isLocked - Locked status (logical)
%       <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset random number and time count for clutter
%                  simulation
%
%   ConstantGammaClutter properties:
%       Sensor                   - Handle of the sensor
%       Gamma                    - Terrain gamma value
%       EarthModel               - Earth model
%     	MaximumRange             - Maximum range
%     	AzimuthCoverage          - Azimuth coverage
%      	PatchAzimuthWidth        - Clutter patch azimuth span
%       CoherenceTime            - Clutter coherence time
%       PropagationSpeed         - Propagation speed
%       SampleRate               - Sample rate
%       PRF                      - Pulse repetition frequency
%       PRFSelectionInputPort    - Enable PRF selection input
%    	OutputFormat             - Output signal format
%      	NumPulses                - Number of pulses in output
%     	NumSamples               - Number of samples in output
%       OperatingFrequency       - Operating frequency
%    	TransmitSignalInputPort  - Enable transmit signal input
%     	TransmitERP              - Effective transmitted power
%       PlatformHeight           - Radar height
%     	PlatformSpeed            - Radar speed
%    	PlatformDirection        - Radar motion direction
%    	BroadsideDepressionAngle - Broadside depression angle
%      	SeedSource               - Source of seed for random number
%                                  generator
%      	Seed                     - Seed for random number generator
%
%   % Examples:
%
%   % Example 1:
%   %   Assume a radar system with a 4 element ULA. The sample rate is 1
%   %   MHz. The PRF is 10 kHz. The propagation speed is 300,000 km/s. The
%   %   operating frequency is 300 MHz. The radar platform is flying 1 km
%   %   above the ground with a path parallel to the ground along the array
%   %   axis. The platform's speed is 2000 m/s. The mainlobe has a
%   %   depression angle of 30 degrees. The terrain has a gamma value of 0
%   %   dB.
%   %
%   %   Simulate the clutter return for 10 pulses, assuming the earth is
%   %   flat. The maximum clutter range of interest is 5 km and the maximum
%   %   azimuth coverage is +/- 60 degrees. The effective transmitted power
%   %   is 5 kw.
%
%   c = 3e8; prf = 10e3; fs = 1e6; fc = 3e8; height = 1000; speed = 2000;
%   tergamma = 0; Rmax = 5000; Npulse = 10; Azcov = 120; lambda = c/fc;
%   tpower = 5000; Nsamp = fs/prf; Nele = 4; depang = 30;
%
%   % create clutter object
%   ha = phased.ULA('NumElements',Nele,'ElementSpacing',lambda/2);
%   clutter = phased.ConstantGammaClutter('Sensor',ha,...
%           'PropagationSpeed',c,'OperatingFrequency',fc,'PRF',prf,...
%           'SampleRate',fs,'Gamma',tergamma,'EarthModel','Flat',...
%           'TransmitERP',tpower,'PlatformHeight',height,...
%           'PlatformSpeed',speed,'PlatformDirection',[90;0],...
%           'BroadsideDepressionAngle',depang,'MaximumRange',Rmax,...
%           'AzimuthCoverage',Azcov);
%
%   % simulate clutter
%   csig = zeros(Nsamp,Nele,Npulse);
%   for m = 1:Npulse
%       csig(:,:,m) = clutter();
%   end
%
%   % plot angle-Doppler response of the clutter at the 20th range bin
%   hresp = phased.AngleDopplerResponse('SensorArray',ha,...
%               'OperatingFrequency',fc,'PropagationSpeed',c,'PRF',prf);
%   plotResponse(hresp,shiftdim(csig(20,:,:)),'NormalizeDoppler',true);
%
%   % Example 2:
%   %   Use the same radar system as the previous example but use a
%   %   rectangular waveform as an additional input. The pulse width of the
%   %   waveform is 2 microseconds.
%
%   c = 3e8; prf = 10e3; fs = 1e6; fc = 3e8; height = 1000; speed = 2000;
%   tergamma = 0; Rmax = 5000; Npulse = 10; Azcov = 120; lambda = c/fc;
%   tpower = 5000; Nsamp = fs/prf; Nele = 4; depang = 30; pw = 2e-6;
%
%   % create clutter object
%   ha = phased.ULA('NumElements',Nele,'ElementSpacing',lambda/2);
%   clutter = phased.ConstantGammaClutter('Sensor',ha,...
%           'PropagationSpeed',c,'OperatingFrequency',fc,'PRF',prf,...
%           'SampleRate',fs,'Gamma',tergamma,'EarthModel','Flat',...
%           'PlatformHeight',height,'TransmitSignalInputPort',true,...
%           'PlatformSpeed',speed,'PlatformDirection',[90;0],...
%           'BroadsideDepressionAngle',depang,'MaximumRange',Rmax,...
%           'AzimuthCoverage',Azcov);
%
%   % simulate clutter
%   tx_sig = tpower*ones(floor(pw*fs),1);
%   csig = zeros(Nsamp,Nele,Npulse);
%   for m = 1:Npulse
%       csig(:,:,m) = clutter(tx_sig);
%   end
%
%   % plot angle-Doppler response of the clutter at the 20th range bin
%   hresp = phased.AngleDopplerResponse('SensorArray',ha,...
%               'OperatingFrequency',fc,'PropagationSpeed',c,'PRF',prf);
%   plotResponse(hresp,shiftdim(csig(20,:,:)),'NormalizeDoppler',true);
%
%   See also phased, phased.BarrageJammer, surfacegamma.


%   Copyright 2010-2017 The MathWorks, Inc.
%

%   Reference
%   [1] Maurice Long, Radar Reflectivity of Land and Sea, 3rd Ed. Artech
%       House, 2001
%   [2] J. R. Guerci, Space Time Adaptive Processing for Radar, Artech
%       House, 2003
%   [3] James Ward, Space-time Adaptive Processing for Airborne Radar,
%       Lincoln Lab Technical Report, 1994
%   [4] David Greig, Time Domain Clutter Model for Airborne Pulse Doppler
%       Radar, IEE Colloquium on Radar System Modeling, 1998

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    methods
        function obj = ConstantGammaClutter(varargin)
            %ConstantGammaClutter   Construct the ConstantGammaClutter
            %class.
            obj@phased.internal.AbstractConstantGammaClutter(varargin{:});
            end
                end
                    
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            
            sensorProp = matlab.system.display.internal.Property('Sensor', ...
                'CustomPresenter', 'phased.internal.SensorArrayDialog', ...
                'CustomPresenterPropertyGroupsArgument', 'antenna', ...
                'Description', 'Sensor array');
            sensorTab = matlab.system.display.SectionGroup('Title', 'Sensor Array', 'PropertyList', {sensorProp});
            
            
            dSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false,'UseClassDefault',false,...
                'Default','Property');
            dSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false,'UseClassDefault',false,...
                'Default','randi(65535,1)');
            
            dTransmitSignalInputPort = matlab.system.display.internal.Property(...
                'TransmitSignalInputPort','IsGraphical',false);
            
            radarProps = {...
                'OperatingFrequency',...
                dTransmitSignalInputPort,...
                'TransmitERP',...
                'PlatformHeight',...
                'PlatformSpeed',...
                'PlatformDirection',...
                'BroadsideDepressionAngle',...
                dSeedSource,...
                dSeed};
            
            radarTab = matlab.system.display.SectionGroup('Title', 'Radar', ...
                'PropertyList', radarProps);
            
            waveformProps = {...
                'SampleRate',...
                'PRF', ...
                'PRFSelectionInputPort',...
                'OutputFormat',...
                'NumPulses',...
                'NumSamples'};
            
            waveformSection = matlab.system.display.Section('Title', 'Reflected signal', ...
                'PropertyList', waveformProps);
            clutterProps = {...
                'Gamma',...
                'EarthModel',...
                'MaximumRange',...
                'AzimuthCoverage',...
                'PatchAzimuthWidth',...
                'CoherenceTime',...
                'PropagationSpeed'};
            clutterSection = matlab.system.display.Section('Title', 'Clutter', ...
                'PropertyList', clutterProps);
            
            mainTab = matlab.system.display.SectionGroup('TitleSource', 'Auto', ...
                'Sections', [clutterSection waveformSection]);
            
            groups = [mainTab,radarTab,sensorTab];
        end
        
        end
        
    methods (Static,Hidden)
        function a = getAlternateBlock
            a = 'phasedenvlib/Constant Gamma Clutter';
        end
    end
end


% [EOF]
