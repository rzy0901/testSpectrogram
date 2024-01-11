classdef (Sealed,StrictDefaults) ScatteringMIMOChannel < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.SampleTime
%ScatteringMIMOChannel  Scattering based MIMO propagation channel
%   H = phased.ScatteringMIMOChannel creates a multipath propagation
%   channel System object, H. This object simulates signal propagation by
%   radiating the signal from the transmit array towards multiple
%   scatterers, reflecting the signal from the scatterers, and propagate
%   back to the receive array. The model applies time delay, gain, phase
%   shift, and atmospheric loss to the input signal.
%
%   H = phased.ScatteringMIMOChannel(Name,Value) returns a line of sight
%   propagation channel object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) propagates the input signal, X, from the transmit array
%   through the channel, H, and outputs the received signal, Y, at the
%   receive array. This syntax is applicable when you set the Polarization
%   property to 'None' or 'Combined'.
%
%   X is an LxNt matrix where L is the number of snapshots and Nt is the
%   number of elements in the transmit array. Each column in X represents
%   the input signal to the corresponding sensor element in the transmit
%   array.
%
%   Y is an LxNr matrix where Nr is the number of elements in the receive
%   array. Each column in Y represents the output signal from the
%   corresponding sensor element in the receive array.
%
%   [YH,YV] = step(H,XH,XV) propagates the input signal XH and XV, to H
%   polarization port and V polarization port, respectively, at the
%   transmit array through the channel, H, and outputs the received signal
%   YH and YV, from H polarization port and V polarization port,
%   respectively, at the receive array. This syntax is applicable when you
%   set the Polarization property to 'Dual'.
%
%   Both XH and XV are LxNt matrices where L is the number of snapshots and
%   Nt is the number of elements in the transmit array. Each column in XH
%   and XV represents the input signal to the H and V polarization ports at
%   corresponding sensor element in the transmit array.
%
%   YH and YV are LxNr matrices where Nr is the number of elements in the
%   receive array. Each column in YH and YV represents the output signal
%   from the H and V polarization ports at corresponding sensor element in
%   the receive array.
%
%   [...] = step(...,TXPOS,TXVEL,TXAX) specifies the transmit array
%   position, TXPOS, as a 3-element column vector in the form of [x;y;z]
%   (in meters); the transmit array velocity, TXVEL, as a 3-element column
%   vector in the form of [Vx;Vy;Vz] (in meters/second); and the transmit
%   array orientation, TXAX, as a 3x3 matrix whose columns represent the x,
%   y, and z axes of the transmit array's local coordinates, respectively.
%   Each axis is specified in the form of [x;y;z]. This syntax applies when
%   you set the TransmitArrayMotionSource property to 'Input port'.
%
%   [...] = step(...,RXPOS,RXVEL,RXAX) specifies the receive array
%   position, RXPOS, as a 3-element column vector in the form of [x;y;z]
%   (in meters); the receive array velocity, RXVEL, as a 3-element column
%   vector in the form of [Vx;Vy;Vz] (in meters/second); and the receive
%   array orientation, RXAX, as a 3x3 matrix whose columns represent the x,
%   y, and z axes of the receive array's local coordinates, respectively.
%   Each axis is specified in the form of [x;y;z]. This syntax applies when
%   you set the ReceiveArrayMotionSource property to 'Input port'.
%
%   [...] = step(...,SCATPOS,SCATVEL,SCATCOEF) specifies the scatterer
%   position, SCATPOS, as a 3xNs matrix where Ns is the number of
%   scatterers. Each column in SCATPOS specifies the position in the form
%   of [x;y;z] (in meters) for a scatterer. SCATVEL is also a 3xNs matrix
%   whose columns specify the scatterer velocity in the form of
%   [Vx;Vy;Vz] (in meters/second) for the corresponding scatterers.
%   SCATCOEF specifies the scattering coefficients as an Ns-element row
%   vector. Each element in the vector represents the scattering
%   coefficient of the corresponding scatterer. This syntax applies when
%   you set the ScattererSpecificationSource property to 'Input port' and
%   the Polarization property to 'None'.
%
%   [...] = step(...,SCATPOS,SCATVEL,SCATMAT,SCATAXES) specifies the
%   scatterer position, SCATPOS, as a 3xNs matrix where Ns is the number of
%   scatterers. Each column in SCATPOS specifies the position in the form
%   of [x;y;z] (in meters) for a scatterer. SCATVEL is also a 3xNs matrix
%   whose columns specify the scatterer velocity in the form of [Vx;Vy;Vz]
%   (in meters/second) for the corresponding scatterers. SCATMAT specifies
%   the scattering matrices as a 2x2xNs array. Each page in the array
%   represents the scattering matrix of the corresponding scatterer in the
%   form of [Shh Shv;Svh Svv] where Shv is the scattering coefficient when
%   the input is V polarized and the output is H polarized. SCATAXES
%   specifies the scatterer orientation axes as a 3x3xNs array. Each page
%   in the array represents the x, y, and z axes of the corresponding
%   scatterer's local coordinates, respectively. Each axis is specified in
%   the form of [x;y;z]. This syntax applies when you set the
%   ScattererSpecificationSource property to 'Input port' and the
%   Polarization property to 'Combined' or 'Dual'.
%
%   [...,CR,TAU] = step(...) also returns the channel response for the
%   channel. CR is either an NtxNrxNs matrix when you set the
%   SimulateDirectPath property to false, or an NtxNrx(Ns+1) matrix when
%   you set the SimulateDirectPath property to true. The pages in CR are
%   the channel matrices for the paths through corresponding scatterers.
%   The channel matrix correspond to the direct line of sight path is the
%   last page in CR when you set the SimulateDirectPath property to true.
%   TAU is either a 1xNs vector when you set the SimulateDirectPath
%   property to false, or a 1x(Ns+1) vector when you set the
%   SimulateDirectPath property to true. The entries in TAU are the path
%   delays for the corresponding paths. The path delay for the direct line
%   of sight path is the last entry in TAU when you set the
%   SimulateDirectPath property to true. This syntax applies you set the
%   ChannelResponseOutputPort property to true and the Polarization
%   property to 'None' or 'Combined'.
%
%   [...,CR_HH,CR_HV,CR_VH,CR_VV,TAU] = step(...) also returns the channel
%   response for the channel. CR_HH is the channel response from input H
%   polarization to output H polarization. CR_HV is the channel response
%   from input H polarization to output V polarization. CR_VH is the
%   channel response from input V polarization to output H polarization.
%   CR_VV is the channel response from input V polarization to output V
%   polarization. Each channel response is either an NtxNrxNs matrix when
%   you set the SimulateDirectPath property to false, or an NtxNrx(Ns+1)
%   matrix when you set the SimulateDirectPath property to true. The pages
%   in each channel response matrix are the channel matrices for the paths
%   through corresponding scatterers. The channel matrix correspond to the
%   direct line of sight path is the last page in the channel response
%   matrix when you set the SimulateDirectPath property to true. TAU is
%   either a 1xNs vector when you set the SimulateDirectPath property to
%   false, or a 1x(Ns+1) vector when you set the SimulateDirectPath
%   property to true. The entries in TAU are the path delays for the
%   corresponding paths. The path delay for the direct line of sight path
%   is the last entry in TAU when you set the SimulateDirectPath property
%   to true. This syntax applies you set the ChannelResponseOutputPort
%   property to true and the Polarization property to 'Dual'.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [Y,CR,TAU] = step(H,X,...
%      TXPOS,TXVEL,TXAX,RXPOS,RXVEL,RXAX,SCATPOS,SCATVEL,SCATCOEF)
%
%   [YH,YV,CRHH,CRHV,CRVH,CRVV,TAU] = step(H,XH,XV,...
%      TXPOS,TXVEL,TXAX,RXPOS,RXVEL,RXAX,SCATPOS,SCATVEL,SCATCOEF,SCATAXES)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ScatteringMIMOChannel methods:
%
%   step     - Propagate signal through MIMO channel (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create MIMO propagation channel object with same property 
%              values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the propagation channel
%
%   ScatteringMIMOChannel properties:
%
%   TransmitArray                - Transmit sensor array
%   ReceiveArray                 - Receive sensor array
%   PropagationSpeed             - Propagation speed 
%   CarrierFrequency             - Signal carrier frequency 
%   Polarization                 - Polarization configuration
%   SpecifyAtmosphere            - Specify atmosphere parameters
%   Temperature                  - Temperature 
%   DryAirPressure               - Dry air pressure 
%   WaterVapourDensity           - Water vapour density 
%   LiquidWaterDensity           - Liquid water density 
%   RainRate                     - Rain rate 
%   SampleRate                   - Sample rate 
%   SimulateDirectPath           - Simulate direct path propagation
%   ChannelResponseOutputPort    - Output channel response
%   MaximumDelaySource           - Source of maximum delay  
%   MaximumDelay                 - Maximum delay 
%   TransmitArrayMotionSource    - Source of transmit array motion
%   TransmitArrayPosition        - Position of the transmit array 
%   TransmitArrayOrientationAxes - Orientation of the transmit array
%   ReceiveArrayMotionSource     - Source of receive array motion
%   ReceiveArrayPosition         - Position of the receive array 
%   ReceiveArrayOrientationAxes  - Orientation of the receive array
%   ScattererSpecificationSource - ScattererSpecification
%   NumScatterers                - Number of scatterers
%   ScattererPositionBoundary    - Boundary of scatterer positions
%   ScattererPosition            - Position of scatterers
%   ScattererCoefficient         - Scattering coefficients
%   ScatteringMatrix             - Scattering matrix 
%   ScattererOrientationAxes     - Orientation of the scatterers
%   SeedSource                   - Source of seed for random number 
%                                  generator
%   Seed                         - Seed for random number generator
%
%   % Examples:
%
%   % Example 1:
%   %   Create a MIMO channel at 30 GHz with an 16-element transmit array 
%   %   and a 16-element receive array. Assume the elements have cosine 
%   %   response shapes and the arrays are standard uniform linear arrays. 
%   %   The channel will have 100 randomly generated static scatterers. 
%   %   The transmit array is located at [0;0;50] meters and the receive 
%   %   array is at [200;0;0] meters. Use the channel to compute the 
%   %   propagated signal. Assume the sample rate for the signal is 10 MHz.
%
%   fc = 30e9;
%   c = 3e8;
%   lambda = c/fc;
%   fs = 10e6;
%   txarray = phased.ULA('Element',phased.CosineAntennaElement,...
%       'NumElements',16,'ElementSpacing',lambda/2);
%   rxarray = phased.ULA('Element',phased.CosineAntennaElement,...
%       'NumElements',16,'ElementSpacing',lambda/2);
%
%   chan = phased.ScatteringMIMOChannel(...
%       'TransmitArray',txarray,...
%       'ReceiveArray',rxarray,...
%       'PropagationSpeed',c,...
%       'CarrierFrequency',fc,...
%       'SampleRate',fs,...
%       'TransmitArrayPosition',[0;0;50],...
%       'ReceiveArrayPosition',[200;0;0],...
%       'NumScatterers',100);
%
%   x = randi(2,[100 16])-1;
%   y = chan(x);
%   
%   % Example 2:
%   %   Create a MIMO channel at 30 GHz with an 16-element transmit array 
%   %   and a 64-element receive array. Assume the elements have cosine 
%   %   response shapes and the arrays are standard uniform linear arrays. 
%   %   The transmit array is located at [0;0;50] meters. 
%   %
%   %   The receive array has an initial position at [200;0;0] meters and 
%   %   is moving at a speed of [10;0;0] meters/second. There are 200 
%   %   static scatterers randomly located on xy plane in a square centered
%   %   at [200;0;0] and with a side length of 100 meters.
%   %   
%   %   Use the channel to compute the propagated signal. Assume the sample
%   %   rate for the signal is 10 MHz and the frame length is 1000 samples.
%   %   Collect 5 frames of received signal.
%    
%   fc = 30e9;
%   c = 3e8;
%   lambda = c/fc;
%   fs = 10e6;
%   txarray = phased.ULA('Element',phased.CosineAntennaElement,...
%       'NumElements',16,'ElementSpacing',lambda/2);
%   rxarray = phased.ULA('Element',phased.CosineAntennaElement,...
%       'NumElements',64,'ElementSpacing',lambda/2);
%
%   Ns = 200;
%   scatpos = [100*rand(1,Ns)+150;100*rand(1,Ns)+150;zeros(1,Ns)];
%   scatcoef = randn(1,Ns)+1i*randn(1,Ns);
%
%   Nframesamp = 1000;
%   Tframe = Nframesamp/fs;
%   rxmobile = phased.Platform('InitialPosition',[200;0;0],...
%       'Velocity',[10;0;0],'OrientationAxesOutputPort',true);
%
%   chan = phased.ScatteringMIMOChannel(...
%       'TransmitArray',txarray,...
%       'ReceiveArray',rxarray,...
%       'PropagationSpeed',c,...
%       'CarrierFrequency',fc,...
%       'SampleRate',fs,...
%       'TransmitArrayPosition',[0;0;50],...
%       'ReceiveArrayMotionSource','Input port',...
%       'ScattererSpecificationSource','Property',...
%       'ScattererPosition',scatpos,...
%       'ScattererCoefficient',scatcoef);
%
%   x = randi(2,[Nframesamp 16])-1;
%   for m = 1:5
%       [rxpos,rxvel,rxax] = rxmobile(Tframe);
%       y = chan(x,rxpos,rxvel,rxax);
%   end
%
%   % Example 3:
%   %   Create a MIMO channel at 30 GHz with an 16-element transmit array 
%   %   and a 64-element receive array. Assume the elements have cosine 
%   %   response shapes and the arrays are standard uniform linear arrays. 
%   %   The transmit array is located at [0;0;50] meters. 
%   %
%   %   The receive array has an initial position at [200;0;0] meters and 
%   %   is moving at a speed of [10;0;0] meters/second. There are 200 
%   %   static scatterers randomly located on xy plane in a square centered
%   %   at [200;0;0] and with a side length of 100 meters.
%   %   
%   %   Use the channel to compute the propagated polarized signal. Assume 
%   %   the sample rate for the signal is 10 MHz and the frame length is 
%   %   1000 samples. Collect 5 frames of received signal.
%
%   fc = 30e9;
%   c = 3e8;
%   lambda = c/fc;
%   fs = 10e6;
%   txarray = phased.ULA('Element',phased.ShortDipoleAntennaElement,...
%       'NumElements',16,'ElementSpacing',lambda/2);
%   rxarray = phased.ULA('Element',phased.ShortDipoleAntennaElement,...
%       'NumElements',64,'ElementSpacing',lambda/2);
%  
%   Ns = 200;
%   scatpos = [100*rand(1,Ns)+150;100*rand(1,Ns)+150;zeros(1,Ns)];
%   temp = randn(1,Ns)+1i*randn(1,Ns);
%   scatcoef = repmat(eye(2),1,1,Ns).*permute(temp,[1 3 2]);
%   scatax = repmat(eye(3),1,1,Ns);
%  
%   Nframesamp = 1000;
%   Tframe = Nframesamp/fs;
%   rxmobile = phased.Platform('InitialPosition',[200;0;0],...
%       'Velocity',[10;0;0],'OrientationAxesOutputPort',true);
%  
%   chan = phased.ScatteringMIMOChannel(...
%       'TransmitArray',txarray,...
%       'ReceiveArray',rxarray,...
%       'PropagationSpeed',c,...
%       'CarrierFrequency',fc,...
%       'SampleRate',fs,...
%       'Polarization','Combined',...
%       'TransmitArrayPosition',[0;0;50],...
%       'ReceiveArrayMotionSource','Input port',...
%       'ScattererSpecificationSource','Property',...
%       'ScattererPosition',scatpos,...
%       'ScatteringMatrix',scatcoef,...
%       'ScattererOrientationAxes',scatax);
%  
%   x = randi(2,[Nframesamp 16])-1;
%   for m = 1:5
%       [rxpos,rxvel,rxax] = rxmobile(Tframe);
%       y = chan(x,rxpos,rxvel,rxax);
%   end
%
%   See also phased, phased.LOSChannel, phased.TwoRayChannel,
%   scatteringchanmtx.

%   Copyright 2016-2017 The MathWorks, Inc.

%   References
%   [1] Robert Health Jr. et al. An Overview of Signal Processing
%   Techniques for Millimeter Wave MIMO Systems
%   [2] David Tse and Pramod Viswanath, Fundamentals of Wireless
%   Communications, Cambridge, 2005
%   [3] Arogyswami Paulraj, Introduction to Space-Time Wireless
%   Communication, Cambridge, 2003

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    % Public, non-tunable properties
    properties(Nontunable)
        %TransmitArray    Transmit sensor array
        %   Specify the transmit sensor array as a handle. The sensor array
        %   must be an array object in the phased package.
        TransmitArray
        %ReceiveArray    Receive sensor array
        %   Specify the receive sensor array as a handle. The sensor array
        %   must be an array object in the phased package.
        ReceiveArray
        %PropagationSpeed Propagation speed (m/s)
        %   Specify the wave propagation speed (in m/s) in free space as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('LightSpeed')
        %CarrierFrequency Signal carrier frequency (Hz)
        %   Specify the carrier frequency (in Hz) of the narrowband signal
        %   as a scalar. The default value of this property is 3e8 (300
        %   MHz).
        CarrierFrequency = 3e8       
        %Polarization   Polarization configuration
        %   Specify the radiator's polarization configuration as one of
        %   'None' | 'Combined' | 'Dual', where the default is 'None'. When
        %   you set this property to 'None', the signal is considered as
        %   polarization independent. When you set this property to
        %   'Combined', the radiated signal contains the polarization
        %   information and the signal reflects the sensor's native
        %   polarization. When you set this property to 'Dual', the H and V
        %   polarization of the sensor can be fed independently.
        Polarization = 'None'
    end
    
    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from the Simulink time engine. Set SampleRateFromInputCheckbox
        %   to false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true
    end
    
    properties (Nontunable)
        %SampleRate Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a scalar. The default value
        %   of this property is 1e6 (1 MHz).
        SampleRate = 1e6    
    end
    
    properties (Nontunable, Logical)
        %SpecifyAtmosphere  Specify atmosphere parameters
        %   Set this property to true to specify atmosphere parameters and
        %   consider the loss due to atmosphere effects. Set this property
        %   to false to ignore atmosphere effects in propagation. The
        %   default value of this property is false.
        SpecifyAtmosphere = false
    end
    
    properties (Nontunable)
        %Temperature    Temperature (degrees Celsius)
        %   Specify the temperature of the propagation channel in degrees
        %   Celsius as a scalar. The default value of this property is 15.
        %   This property only applies when you set the SpecifyAtmosphere
        %   property to true.
        Temperature = 15
        %DryAirPressure     Dry air pressure (Pa)
        %   Specify the dry air pressure of the propagation channel in
        %   Pascal (Pa) as a scalar. The default value of this property is
        %   101325 (1 atmosphere). This property only applies when you set
        %   the SpecifyAtmosphere property to true.
        %
        %   The ITU atmosphere gas model is valid between 1 to 1000 GHz.
        DryAirPressure = 101325
        %WaterVapourDensity Water vapour density (g/m^3)
        %   Specify the water vapour density of the propagation channel in
        %   g/m^3 as a scalar. The default value of this property is 7.5.
        %   This property only applies when you set the SpecifyAtmosphere
        %   property to true.
        %
        %   The ITU atmosphere gas model is valid between 1 to 1000 GHz.
        WaterVapourDensity = 7.5
        %LiquidWaterDensity   Liquid water density (g/m^3)
        %   Specify the liquid water density of the propagation channel in
        %   g/m^3 as a scalar. The default value of this property is 0. The
        %   liquid water density is used to describe the characteristics of
        %   fog or cloud. Typical values for water density are 0.05 for
        %   medium fog and 0.5 for thick fog. This property only applies
        %   when you set the SpecifyAtmosphere property to true.
        %
        %   The ITU fog model is valid between 10 to 1000 GHz.
        LiquidWaterDensity = 0
        %RainRate   Rain rate (mm/h)
        %   Specify the rain rate of the propagation channel in mm/h as a
        %   scalar. The default value of this property is 0, indicating no
        %   rain. This property only applies when you set the
        %   SpecifyAtmosphere property to true.
        %
        %   The ITU rain model is valid between 1 to 1000 GHz.
        RainRate = 0
    end
    
    properties (Nontunable, Logical)
        %SimulateDirectPath     Simulate direct path propagation
        %   Set this property to true to simulate the propagation along the
        %   direct (line of sight) path. Set this property to false to
        %   disable simulating the direct path. The default value of this
        %   property is false.
        SimulateDirectPath = false
        %ChannelResponseOutputPort  Output channel response
        %   Set this property to true to output channel response. Set this
        %   property to false to disable output channel response. The
        %   default value of this property is false.
        ChannelResponseOutputPort = false
    end
    
    properties (Nontunable)
        %TransmitArrayMotionSource  Source of transmit array motion
        %   Specify the source of transmit array motion as one of
        %   'Property' | 'Input port', where the default is 'Property'.
        %   When you set this property to 'Property', the transmit array is
        %   assumed to be stationary and its location and orientation are
        %   specified through properties. When you set this property to
        %   'Input port', the transmit array location, velocity, and
        %   orientation are specified through input arguments.
        TransmitArrayMotionSource = 'Property'
        %TransmitArrayPosition   Position of the transmit array (m)
        %   Specify the position of the transmit array's phase center as a
        %   3x1 vector in the form of [x; y; z] (in meters). The default
        %   value of this property is [0; 0; 0]. This property applies when
        %   you set the TransmitArrayMotionSource property to 'Property'.
        TransmitArrayPosition = [0;0;0]
        %TransmitArrayOrientationAxes   Orientation of the transmit array
        %   Specify the 3 axes that defines the local (x, y, z) coordinate
        %   system for the transmit array as a 3x3 matrix (each column
        %   corresponds to an axis). The 3 axes must be orthonormal. The
        %   default value of this property is [1 0 0;0 1 0;0 0 1]. This
        %   property applies when you set the TransmitArrayMotionSource
        %   property to 'Property'
        TransmitArrayOrientationAxes = eye(3)
        %ReceiveArrayMotionSource  Source of receive array motion
        %   Specify the source of receive array motion as one of 'Property'
        %   | 'Input port', where the default is 'Property'. When you set
        %   this property to 'Property', the receive array is assumed to be
        %   stationary and its location and orientation are specified
        %   through properties. When you set this property to 'Input port',
        %   the receive array location, velocity, and orientation are
        %   specified through input arguments.
        ReceiveArrayMotionSource = 'Property'
        %ReceiveArrayPosition   Position of the receive array (m)
        %   Specify the position of the receive array's phase center as a
        %   3x1 vector in the form of [x; y; z] (in meters). The default
        %   value of this property is [0; 0; 0]. This property applies when
        %   you set the ReceiveArrayMotionSource property to 'Property'.
        ReceiveArrayPosition = [physconst('LightSpeed')/1e5;0;0]
        %ReceiveArrayOrientationAxes   Orientation of the receive array
        %   Specify the 3 axes that defines the local (x, y, z) coordinate
        %   system for the receive array as a 3x3 matrix (each column
        %   corresponds to an axis). The 3 axes must be orthonormal. The
        %   default value of this property is [-1 0 0;0 -1 0;0 0 1] so the
        %   receive array faces the transmit array if the array is in the
        %   yz plane. This property applies when you set the
        %   ReceiveArrayMotionSource property to 'Property'.
        ReceiveArrayOrientationAxes = azelaxes(180,0)
        %ScattererSpecificationSource  Scatterer specification
        %   Specify the source of scatterer specification as one of 'Auto'
        %   | 'Property' | 'Input port', where the default is 'Auto'. When
        %   you set this property to 'Auto', the scatterers are assumed to
        %   be stationary and their locations and scattering coefficients
        %   are randomly generated. When you set this property to
        %   'Property', the scatterers are assumed to be stationary and
        %   their locations and scattering coefficients are specified
        %   through properties. When you set this property to 'Input port',
        %   the scatterer locations, velocities, and scattering
        %   coefficients are specified through input arguments.
        ScattererSpecificationSource = 'Auto'
        %NumScatterers  Number of scatterers
        %   Specify the number of scatterers as a nonnegative integer. The
        %   default value of this property is 1. This property applies when
        %   you set the ScattererSpecificationSource to 'Auto'.
        NumScatterers = 1
        %ScattererPositionBoundary Boundary of scatterer positions
        %   Specify the boundaries in which scatterers are located as a
        %   either a 2-element row vector or a 3x2 matrix. 
        %
        %   If the boundary is a 2-element row vector, the value is
        %   interpreted as [min max] that applies to all 3 dimensions. If
        %   the boundary is a 3x2 matrix, the value specifies the boundary
        %   in all 3 dimensions using the form [x_min x_max;y_min y_max;
        %   z_min z_max].
        %
        %   The default value of this property is [0 1000]. This property
        %   applies when you set the ScattererSpecificationSource to
        %   'Auto'.
        ScattererPositionBoundary = [0 1000]
        %ScattererPosition   Positions of scatterers (m)
        %   Specify the positions of scatterers as a 3xNs matrix in the
        %   form of [x; y; z] (in meters) where Ns is the number of
        %   scatterers. The default value of this property is [0; 0; 0].
        %   This property applies when you set the
        %   ScattererSpecificationSource property to 'Property'.
        ScattererPosition = [physconst('LightSpeed')*5e-6;0;0]
        %ScattererCoefficient    Scattering coefficients
        %   Specify the scatterer scattering coefficients as a 1xNs vector
        %   where Ns is the number of scatterers. Each entry in the vector
        %   represents the scattering coefficient of the corresponding
        %   scatterer. The default value of this property is 1. This
        %   property applies when you set the ScattererSpecificationSource
        %   to 'Property' and the Polarization property to 'None'.
        ScattererCoefficient = 1
        %ScatteringMatrix  Scattering matrix 
        %   Specify the scattering matrices of the scatterers as a 2x2xNs
        %   array where Ns is the number of scatterers. Each page of the
        %   array in the form of [s_hh s_hv;s_vh s_vv], where s_hv
        %   specifies the scattering when the input signal is vertically
        %   polarized and the reflected signal is horizontally polarized.
        %   This property applies when you set the ScatteringMatrixSource
        %   property to 'Property' and the Polarization property to
        %   'Combined' or 'Dual'. The default value of this property is [1
        %   0;0 1]. 
        ScatteringMatrix = eye(2)
        %ScattererOrientationAxes   Orientation of the scatterers
        %   Specify the 3 axes that defines the local (x, y, z) coordinate
        %   system for scatterers as a 3x3xNs matrix where Ns is the number
        %   of scatterers. Within each page of the array, each column
        %   represents an axis of the local coordinate system of the
        %   corresponding scatterer. The 3 axes must be orthonormal. The
        %   default value of this property is [1 0 0;0 1 0;0 0 1]. This
        %   property applies when you set the ScatteringMatrixSource
        %   property to 'Property' and the Polarization property to
        %   'Combined' or 'Dual'.
        ScattererOrientationAxes = eye(3)
    end
   
    properties (Nontunable)
        %MaximumDelaySource  Source of maximum delay 
        %   Specify how the maximum delay is specified as one of 'Auto' |
        %   'Property', where the default is 'Auto'. When you set this
        %   property to 'Auto', the channel automatically allocates the
        %   memory to simulate the propagation delay. When you set this
        %   property to 'Property', the maximum delay is specified via
        %   MaximumDelay property and any signal arrives beyond maximum
        %   delay is ignored.
        %
        %   To use ScatteringMIMOChannel in MATLAB Function Block in
        %   Simulink, set this property to 'Property'.
        MaximumDelaySource = 'Auto'
        %MaximumDelay    Maximum delay (s)
        %   Specify the maximum delay (in seconds) as a positive scalar.
        %   This property indicates the maximum delay. If the delay is
        %   beyond this value, the signal is ignored. This property applies
        %   when you set the MaximumDelaySource property to 'Property'. The
        %   default value of this property is 10e-6;
        MaximumDelay = 10e-6
    end
    
    properties (Nontunable)
        %SeedSource   Source of seed for random number generator
        %   Specify how the random numbers are generated as one of 'Auto' |
        %   'Property', where the default is 'Auto'. When you set this
        %   property to 'Auto', the random numbers are generated using the
        %   default MATLAB random number generator. When you set this
        %   property to 'Property', a private random number generator is
        %   used with a seed specified by the value of the Seed property.
        %   This property applies when you set the
        %   ScattererSpecificationSource to 'Auto'.
        %
        %   To use this object with Parallel Computing Toolbox software,
        %   set this property to 'Auto'.
        SeedSource = 'Auto'
        %Seed     Seed for random number generator
        %   Specify the seed for the random number generator as a
        %   non-negative integer. The integer must be less than 2^32. This
        %   property applies when you set the ScattererSpecificationSource
        %   property to 'Auto' and the SeedSource property to 'Property'.
        %   The default value of this property is 0.
        Seed = 0
    end
    
    properties (Constant, Hidden)
        PolarizationSet = matlab.system.StringSet({'None','Combined','Dual'});
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties(Constant, Hidden)
        TransmitArrayMotionSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        ReceiveArrayMotionSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        ScattererSpecificationSourceSet = matlab.system.StringSet({'Auto','Property','Input port'});
        MaximumDelaySourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        SeedSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    % Pre-computed constants
    properties(Access = private, Nontunable)
        cTxSteeringVector
        cRxSteeringVector
        cBuffer
        cFractionalDelayFilter
        cCoefficientSource
        cPositionSource
        pPolarizationConfigurationIndex
    end
    
    properties(Access = private)
        pCurrentTime 
        pScattererPosition
        pScattererVelocity
        pScattererCoefficient
        pScattererOrientationAxes
        pScatteringMatrix
        pReferenceDistance
    end
    
    properties(Access = private, Nontunable)
        pLambda
        %Sample rate, in MATLAB, specified by property but in Simulink,
        %specified by engine
        pSampleRate
    end
    
    properties (Access = private, Logical)
        pScattererInitialized = false
    end
    
    properties (Access = private, Logical, Nontunable)
        pConstantTransmitArrayMotion
        pConstantReceiveArrayMotion
        pConstantScattererMotion

    end
    
    methods
        function set.TransmitArray(obj,value)
            validateattributes( value, ...
                {'phased.internal.AbstractArray'},...
                {'scalar'}, '', 'Sensor');
            obj.TransmitArray = value;
        end
        function set.ReceiveArray(obj,value)
            validateattributes( value, ...
                {'phased.internal.AbstractArray'},...
                {'scalar'}, '', 'Sensor');
            obj.ReceiveArray = value;
        end
        function set.CarrierFrequency(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'positive'},'ScatteringMIMOChannel','CarrierFrequency');
            obj.CarrierFrequency = value;
        end
        function set.SampleRate(obj,value)
            validateattributes(value,{'double'},{'scalar','positive',...
                'finite'},'ScatteringMIMOChannel','SampleRate');
            obj.SampleRate = value;
        end
        function set.PropagationSpeed(obj,value)
            sigdatatypes.validateSpeed(value,'ScatteringMIMOChannel',...
                'PropagationSpeed',{'scalar','positive'});
            obj.PropagationSpeed = value;
        end
        function set.Temperature(obj,value)
            obj.Temperature = sigdatatypes.validateCFTemperature(value,...
                '','Temperature',{'scalar','>=',-273.15});
        end
        
        function set.DryAirPressure(obj,value)
            obj.DryAirPressure = sigdatatypes.validatePressure(value,...
                '','DryAirPressure',{'scalar','positive'});
        end
        
        function set.WaterVapourDensity(obj,value)
            obj.WaterVapourDensity = sigdatatypes.validateDensity(value,...
                '','WaterVapourDensity',{'scalar'});
        end
        
        function set.LiquidWaterDensity(obj,value)
            obj.LiquidWaterDensity = sigdatatypes.validateDensity(value,...
                '','LiquidWaterDensity',{'scalar'});
        end
        
        function set.RainRate(obj,value)
            obj.RainRate = sigdatatypes.validateSpeed(value,...
                '','RainRate',{'scalar'});
        end
        function set.TransmitArrayPosition(obj,value)
            sigdatatypes.validate3DCartCoord(value,...
                'ScatteringMIMOChannel','TransmitArrayPosition',...
                {'ncols',1});
            obj.TransmitArrayPosition = value;
        end
        function set.TransmitArrayOrientationAxes(obj,val)
            validateattributes(val,{'double'},{'finite','nonnan','nonempty','real',...
                                'size',[3 3]},'ScatteringMIMOChannel',...
                                'TransmitArrayOrientationAxes');

            %Check that all Orientation axis are orthonormal
            cond =  norm(val'*val-eye(3)) > sqrt(eps);
            if cond
                coder.internal.errorIf(cond, ...
                  'phased:phased:expectedOrthonormalAxes','TransmitArrayOrientationAxes');
            end
            obj.TransmitArrayOrientationAxes = val;
        end
        function set.ReceiveArrayPosition(obj,value)
            sigdatatypes.validate3DCartCoord(value,...
                'ScatteringMIMOChannel','ReceiveArrayPosition',...
                {'ncols',1});
            obj.ReceiveArrayPosition = value;
        end
        function set.ReceiveArrayOrientationAxes(obj,val)
            validateattributes(val,{'double'},{'finite','nonnan','nonempty','real',...
                                'size',[3 3]},'ScatteringMIMOChannel',...
                                'ReceiveArrayOrientationAxes');

            %Check that all Orientation axis are orthonormal
            cond =  norm(val'*val-eye(3)) > sqrt(eps);
            if cond
                coder.internal.errorIf(cond, ...
                  'phased:phased:expectedOrthonormalAxes','ReceiveArrayOrientationAxes');
            end
            obj.ReceiveArrayOrientationAxes = val;
        end
        function set.NumScatterers(obj,value)
            validateattributes(value,{'double'},...
                {'nonnan','nonempty','finite','integer','nonnegative','scalar'},...
                'ScatteringMIMOChannel','NumScatterers');
            obj.NumScatterers = value;
        end
        function set.ScattererPositionBoundary(obj,value)
            validateattributes(value,{'double'},...
                {'nonnan','nonempty','finite','2d','real','ncols',2},...
                'ScatteringMIMOChannel','ScattererPositionBoundary');
            cond = ~(size(value,1)==1 || size(value,1)==3);
            if cond
                coder.internal.errorIf(cond, ...
                  'phased:phased:expectedNumRows','ScattererPositionBoundary','1 3');
            end
            cond = any(value(:,1)>value(:,2));
            if cond
                coder.internal.errorIf(cond, ...
                  'phased:phased:expectedNondecreasingRange');
            end            
            obj.ScattererPositionBoundary = value;
        end
        function set.ScattererPosition(obj,value)
            sigdatatypes.validate3DCartCoord(value,...
                'ScatteringMIMOChannel','ScattererPosition');
            obj.ScattererPosition = value;
        end
        function set.ScattererCoefficient(obj,value)
            validateattributes(value,{'double'},...
                {'nonnan','finite','nonempty','row'},...
                'ScatteringMIMOChannel','ScattererCoefficient');
            obj.ScattererCoefficient = value;
        end
        function set.ScatteringMatrix(obj,value)
            validateattributes(value,{'double'},...
                {'nonnan','finite','nonempty','size',[2 2 size(value,3)]},...
                'ScatteringMIMOChannel','ScatteringMatrix');
            obj.ScatteringMatrix = value;
        end
        function set.ScattererOrientationAxes(obj,val)
            validateattributes(val,{'double'},{'finite','nonnan','nonempty','real',...
                                'size',[3 3 size(val,3)]},'ScatteringMIMOChannel',...
                                'ScattererOrientationAxes');

            %Check that all Orientation axis are orthonormal
            for m = 1:size(val,3)
                cond =  norm(val(:,:,m)'*val(:,:,m)-eye(3)) > sqrt(eps);
                if cond
                    coder.internal.errorIf(cond, ...
                      'phased:phased:expectedOrthonormalAxes','ScattererOrientationAxes');
                end
            end
            obj.ScattererOrientationAxes = val;
        end
        function set.MaximumDelay(obj,value)
            sigdatatypes.validateDuration(value,'ScatteringMIMOChannel',...
                'MaximumDelay',{'scalar','positive'});
            obj.MaximumDelay = value;
        end
        function set.Seed(obj,value)
            validateattributes(value,{'double'},{'scalar','nonnegative',...
                'finite','nonnan','nonempty'},'',...
                'Seed');
            obj.Seed = value;
        end
    end

    methods
        % Constructor
        function obj = ScatteringMIMOChannel(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
            if isempty(coder.target)
                if isempty(obj.TransmitArray)
                    obj.TransmitArray = phased.ULA;
                end
                if isempty(obj.ReceiveArray)
                    obj.ReceiveArray = phased.ULA;
                end
            else
                if ~coder.internal.is_defined(obj.TransmitArray)
                    obj.TransmitArray = phased.ULA;
                end
                if ~coder.internal.is_defined(obj.ReceiveArray)
                    obj.ReceiveArray = phased.ULA;
                end
            end

        end
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj,varargin)
            % Perform one-time calculations, such as computing constants
            fs = obj.SampleRate; % property/method duality
            coder.const(fs);
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            obj.pSampleRate = fs;

            if strcmp(obj.Polarization,'None')
                obj.pPolarizationConfigurationIndex = 0;
            elseif strcmp(obj.Polarization,'Combined')
                obj.pPolarizationConfigurationIndex = 1;
            else  % dual
                obj.pPolarizationConfigurationIndex = 2;
            end
            
            polflag = obj.pPolarizationConfigurationIndex~=0;
            obj.cTxSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.TransmitArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'IncludeElementResponse',true,...
                'EnablePolarization',polflag);
            obj.cRxSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.ReceiveArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'IncludeElementResponse',true,...
                'EnablePolarization',polflag);
            obj.pLambda = obj.PropagationSpeed/obj.CarrierFrequency;
            
            if strcmp(obj.MaximumDelaySource,'Auto')
                obj.cBuffer = phased.internal.CircularBuffer(...
                    'BufferLength',1);
            else
                buflen = ceil(obj.MaximumDelay*obj.pSampleRate);
                obj.cBuffer = phased.internal.CircularBuffer(...
                    'FixedLengthBuffer',true,'BufferLength',buflen,...
                    'BufferWidthSource','Auto');
            end
            
            if strcmp(obj.TransmitArrayMotionSource,'Input port')
                txidx = getVararginIndex(obj,'tx');
                txpos = varargin{txidx(1)};
                if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                    rxidx = getVararginIndex(obj,'rx');
                    rxpos = varargin{rxidx(1)};
                else
                    rxpos = obj.ReceiveArrayPosition;
                end
            else
                txpos = obj.TransmitArrayPosition;
                if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                    rxidx = getVararginIndex(obj,'rx');
                    rxpos = varargin{rxidx(1)};
                else
                    rxpos = obj.ReceiveArrayPosition;
                end
            end
            
            obj.pReferenceDistance = norm(rxpos-txpos); % direct path length
            
            obj.cFractionalDelayFilter = dsp.VariableFractionalDelay;

            obj.cCoefficientSource = phased.internal.NoiseSource(...
                'Distribution','Gaussian',...
                'SeedSource',obj.SeedSource,...
                'OutputComplex',true);
            obj.cPositionSource = phased.internal.NoiseSource(...
                'Distribution','Uniform',...
                'SeedSource',obj.SeedSource,...
                'OutputComplex',false);
            if (obj.SeedSource(1) == 'P')
                obj.cCoefficientSource.Seed = obj.Seed;
                obj.cPositionSource.Seed = obj.Seed;
            end
            
            obj.pConstantTransmitArrayMotion = ...
                strcmp(obj.TransmitArrayMotionSource,'Property');
            obj.pConstantReceiveArrayMotion = ...
                strcmp(obj.ReceiveArrayMotionSource,'Property');
            obj.pConstantScattererMotion = ...
                ~strcmp(obj.ScattererSpecificationSource,'Input port');
            
            if ~strcmp(obj.ScattererSpecificationSource,'Input port')
                if strcmp(obj.ScattererSpecificationSource,'Auto')
                    Nscat = obj.NumScatterers;
                    scatpos = bsxfun(@times,diff(obj.ScattererPositionBoundary,1,2),...
                        step(obj.cPositionSource,1,[3,Nscat]))+...
                        obj.ScattererPositionBoundary(:,1);
                    scatvel = zeros(3,Nscat);
                    scatcoef = step(obj.cCoefficientSource,1,[1,Nscat]);
                    obj.pScattererPosition = scatpos;
                    obj.pScattererVelocity = scatvel;
                    if polflag
                        obj.pScatteringMatrix = permute(scatcoef,[1 3 2]).*eye(2);
                        obj.pScattererOrientationAxes = repmat(eye(3),1,1,Nscat);
                    else
                        obj.pScattererCoefficient = scatcoef;
                    end
                    obj.pScattererInitialized = true;
                else % strcmp(obj.ScattererSpecificationSource,'Property')
                    scatpos = obj.ScattererPosition;
                    scatvel = zeros(size(scatpos));
                    obj.pScattererPosition = scatpos;
                    obj.pScattererVelocity = scatvel;
                    if polflag
                        scatcoef = obj.ScatteringMatrix;
                        obj.pScatteringMatrix = scatcoef;
                        obj.pScattererOrientationAxes = obj.ScattererOrientationAxes;
                    else
                        scatcoef = obj.ScattererCoefficient;
                        obj.pScattererCoefficient = scatcoef;
                    end
                    obj.pScattererInitialized = true;
                 end
            end
            
            xidx = getVararginIndex(obj,'sig');
            x = varargin{xidx(1)};
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        end

        function varargout = stepImpl(obj,varargin)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            lambda = obj.pLambda;
            fs = obj.SampleRate;
            
            if ~obj.pConstantTransmitArrayMotion
                txidx = getVararginIndex(obj,'tx');
                txpos = varargin{txidx(1)};
                txvel = varargin{txidx(2)};
                txaxes = varargin{txidx(3)};
            else
                txpos = obj.TransmitArrayPosition;
                txvel = [0;0;0];
                txaxes = obj.TransmitArrayOrientationAxes;
            end

            if ~obj.pConstantReceiveArrayMotion
                rxidx = getVararginIndex(obj,'rx');
                rxpos = varargin{rxidx(1)};
                rxvel = varargin{rxidx(2)};
                rxaxes = varargin{rxidx(3)};
            else
                rxpos = obj.ReceiveArrayPosition;
                rxvel = [0;0;0];
                rxaxes = obj.ReceiveArrayOrientationAxes;
            end

            if ~obj.pConstantScattererMotion
                scatidx = getVararginIndex(obj,'sc');
                scatpos = varargin{scatidx(1)};
                scatvel = varargin{scatidx(2)};
                Nscat = size(scatpos,2);
                switch obj.pPolarizationConfigurationIndex
                    case 0
                        scatcoef = varargin{scatidx(3)};
                    otherwise
                        scatcoef = varargin{scatidx(3)};
                        scatax = varargin{scatidx(4)};
                end
            else  % ConstantScattererMotion and ScattererInitialized
                scatpos = obj.pScattererPosition;
                scatvel = obj.pScattererVelocity;
                Nscat = size(scatpos,2);
                switch obj.pPolarizationConfigurationIndex
                    case 0
                        scatcoef = obj.pScattererCoefficient;
                    otherwise
                        scatcoef = obj.pScatteringMatrix;
                        scatax = obj.pScattererOrientationAxes;
                end
            end

            if Nscat>0
                [dts,txang] = rangeangle(scatpos,txpos,txaxes);
                [dsr,rxang] = rangeangle(scatpos,rxpos,rxaxes);           
                d = dts+dsr;
                vts = radialspeed(scatpos,scatvel,txpos,txvel);
                vsr = radialspeed(scatpos,scatvel,rxpos,rxvel);
                fd = (vts+vsr)./lambda;
                if obj.pPolarizationConfigurationIndex~=0
                    tgtinang = zeros(2,Nscat);
                    tgtoutang = zeros(2,Nscat);
                    for m = 1:Nscat
                        [~,tgtinang(:,m)] = rangeangle(txpos,scatpos(:,m),scatax(:,:,m));
                        [~,tgtoutang(:,m)] = rangeangle(rxpos,scatpos(:,m),scatax(:,:,m));
                    end
                end
            end

            % check direct path
            hasDirectPath = obj.SimulateDirectPath;


            if hasDirectPath
                d0 = norm(rxpos-txpos);  % direct path
                vrt = radialspeed(txpos,txvel,rxpos,rxvel);
                fd0 = vrt/lambda;

                [~,txang0] = rangeangle(rxpos,txpos,txaxes);
                [~,rxang0] = rangeangle(txpos,rxpos,rxaxes);

                if Nscat > 0
                    Nray = Nscat + 1;
                    d_ray = [d(:).' d0];
                    fd_ray = [fd(:).' fd0];
                    txang_ray = [txang txang0];
                    rxang_ray = [rxang rxang0];
                    switch obj.pPolarizationConfigurationIndex
                        case 0
                            scatcoef_ray = [scatcoef 1];
                        otherwise
                            scatcoef_ray = cat(3,scatcoef,eye(2));
                            scatax_ray = cat(3,scatax,eye(3));
                            tgtinang_ray = [tgtinang rxang0];
                            tgtoutang_ray = [tgtoutang txang0];
                    end
                else
                    Nray = 1;
                    d_ray = d0;
                    fd_ray = fd0;
                    txang_ray = txang0;
                    rxang_ray = rxang0;
                    switch obj.pPolarizationConfigurationIndex
                        case 0
                            scatcoef_ray = 1;
                        otherwise
                            scatcoef_ray = eye(2);
                            tgtinang_ray = rxang0;
                            tgtoutang_ray = txang0;
                    end
                end
            else
                Nray = Nscat;
                d_ray = d(:).';
                fd_ray = fd(:).';
                txang_ray = txang;
                rxang_ray = rxang;
                switch obj.pPolarizationConfigurationIndex
                    case 0
                        scatcoef_ray = scatcoef;
                    otherwise
                        scatcoef_ray = scatcoef;
                        scatax_ray = scatax;
                        tgtinang_ray = tgtinang;
                        tgtoutang_ray = tgtoutang;
                end
            end

            if ~obj.SpecifyAtmosphere
                % spreading loss
                sploss = phased.internal.fspl(d_ray,lambda);
                plossfactordb = sploss;
            else
                fc = obj.CarrierFrequency;
                % spreading loss
                sploss = phased.internal.fspl(d_ray,lambda);
                % gas loss
                fc_gas = min(max(fc,1e9),1000e9);
                gasloss = gaspl(d_ray,fc_gas,...
                    obj.Temperature,obj.DryAirPressure,obj.WaterVapourDensity);
                plossfactordb = sploss+gasloss;
                % fog loss
                if obj.LiquidWaterDensity > 0
                    fc_fog = min(max(fc,10e9),1000e9);
                    fogloss = fogpl(d_ray,fc_fog,...
                        obj.Temperature,obj.LiquidWaterDensity);
                    plossfactordb = plossfactordb+fogloss;
                end
                % rain loss
                if obj.RainRate > 0
                    if Nscat > 0
                        elts = phased.LOSChannel.computeElevationAngle(txpos,scatpos,dts);
                        elsr = phased.LOSChannel.computeElevationAngle(scatpos,rxpos,dsr);
                        tau = zeros(1,Nscat);

                        fc_rain = min(max(fc,1e9),1000e9);
                        rainloss_temp = rainpl(dts,fc_rain,obj.RainRate,elts,tau) + ...
                            rainpl(dsr,fc_rain,obj.RainRate,elsr,tau);
                    end
                    if hasDirectPath
                        elts0 = phased.LOSChannel.computeElevationAngle(txpos,rxpos,d0);
                        if Nscat > 0
                            rainloss = [rainloss_temp; rainpl(d0,fc_rain,obj.RainRate,elts0,0)];
                        else
                            rainloss = rainpl(d0,fc_rain,obj.RainRate,elts0,0);
                        end
                    else
                        rainloss = rainloss_temp;
                    end
                    plossfactordb = plossfactordb+rainloss;
                end
            end
            plossfactor_ray = db2mag(plossfactordb(:).');
            switch obj.pPolarizationConfigurationIndex
                case 0
                    amp = scatcoef_ray./plossfactor_ray.*exp(-1i*2*pi*d_ray/lambda);
                otherwise
                    amp = bsxfun(@times,scatcoef_ray,...
                        permute(1./plossfactor_ray.*exp(-1i*2*pi*d_ray/lambda),[1 3 2]));
            end
            
            freq = obj.CarrierFrequency;
            stvts = step(obj.cTxSteeringVector,freq,txang_ray);
            stvsr = step(obj.cRxSteeringVector,freq,rxang_ray);
            
            xidx = getVararginIndex(obj,'sig');
            x = varargin{xidx(1)};
            Nsamp = size(x,1);
            Nrx = getDOF(obj.ReceiveArray);

            t = (0:Nsamp-1).'/fs+obj.pCurrentTime;
            obj.pCurrentTime = obj.pCurrentTime+Nsamp/fs;
            
            % use absolute delay
            delay = d_ray/obj.PropagationSpeed; 
            delay_samp = delay*fs;
            delay_samp_int = fix(delay_samp);
            delay_samp_frac = delay_samp-delay_samp_int;

            delayvec_samp_frac_temp = repmat(delay_samp_frac,Nrx,1);
            delayvec_samp_frac = delayvec_samp_frac_temp(:).';
                    
            switch obj.pPolarizationConfigurationIndex
                case 0
                    tempy = zeros(Nsamp,Nrx*Nray,'like',1+1i);

                    Hchan = zeros([size(stvts,1) size(stvsr,1) Nray],'like',1+1i);
                    for m = 1:Nray
                        tempx = bsxfun(@times,x,exp(1i*2*pi*fd_ray(m)*t)); % doppler
                        tempchanmat = stvsr(:,m)*stvts(:,m).'*amp(m);
                        Hchan(:,:,m) = tempchanmat.';
                        tempy(:,(m-1)*Nrx+(1:Nrx)) = (tempchanmat*tempx.').';
                    end
                    temp_fracdelayed = step(obj.cFractionalDelayFilter,tempy,delayvec_samp_frac);
                    yr = step(obj.cBuffer,temp_fracdelayed,delay_samp_int);
                    y = complex(sum(reshape(yr,Nsamp,Nrx,Nray),3));
                    
                    varargout{1} = y;
                    varargout{2} = Hchan;
                    varargout{3} = delay;
                    
                case 1
                    
                    tempy = zeros(Nsamp,Nrx*Nray,'like',1+1i);

                    Hchan = zeros([size(stvts.H,1) size(stvsr.H,1) Nray],'like',1+1i);
                    
                    for m = 1:Nray
                        tempx = bsxfun(@times,x,exp(1i*2*pi*fd_ray(m)*t)); % doppler
                        
                        % Transmit polarization
                        txsv_h = stvts.H(:,m).';
                        txsv_v = stvts.V(:,m).';
                        
                        % Translate to global coordinates
                        tempsv_global = sph2cartvec([txsv_h;txsv_v;zeros(size(txsv_v))],...
                            txang_ray(1,m),txang_ray(2,m));
                        tempsv_global = phased.internal.local2globalvec(tempsv_global,txaxes);
                        
                        % Translate to target local coordinates
                        lcl_tgtinsv = phased.internal.global2localvec(...
                            tempsv_global,scatax_ray(:,:,m));
                        lcl_tgtinsv = cart2sphvec(lcl_tgtinsv,tgtinang_ray(1,m),tgtinang_ray(2,m));
                        lcl_tgtoutsv_sph = lcl_tgtinsv(1:2,:);
                        
                        % Apply path loss and scattering matrix
                        lcl_tgtoutsv_sph = amp(:,:,m)*lcl_tgtoutsv_sph;
                        
                        % Translate to global coordinates
                        lcl_tgtoutsv = sph2cartvec(...
                            [lcl_tgtoutsv_sph;zeros(1,size(lcl_tgtoutsv_sph,2))],...
                            tgtoutang_ray(1,m),tgtoutang_ray(2,m));
                        tgtsv_global = phased.internal.local2globalvec(...
                            lcl_tgtoutsv,scatax_ray(:,:,m));
                       
                        % Receive polarization
                        rxsv_h = stvsr.H(:,m).';
                        rxsv_v = stvsr.V(:,m).';
                        rxsv_local = sph2cartvec([rxsv_h;rxsv_v;zeros(size(rxsv_h))],...
                            rxang_ray(1,m),rxang_ray(2,m));
                        rxsv_global = phased.internal.local2globalvec(rxsv_local,rxaxes);
                        
                        % Form channel matrix
                        tempchanmat = rxsv_global.'*tgtsv_global;
                        Hchan(:,:,m) = tempchanmat.';
                        
                        % Compute output
                        tempy(:,(m-1)*Nrx+(1:Nrx)) = (tempchanmat*tempx.').';
                        
                    end
                    
                    temp_fracdelayed = step(obj.cFractionalDelayFilter,tempy,delayvec_samp_frac);
                    yr = step(obj.cBuffer,temp_fracdelayed,delay_samp_int);
                    
                    y = complex(sum(reshape(yr,Nsamp,Nrx,Nray),3));
                    
                    varargout{1} = y;
                    varargout{2} = Hchan;
                    varargout{3} = delay;
                    
                case 2
                    
                    tempy_h = zeros(Nsamp,Nrx*Nray,'like',1+1i);
                    tempy_v = zeros(Nsamp,Nrx*Nray,'like',1+1i);

                    Hchan_hh = zeros([size(stvts.H,1) size(stvsr.H,1) Nray],'like',1+1i); % H -> H
                    Hchan_hv = zeros([size(stvts.H,1) size(stvsr.V,1) Nray],'like',1+1i); % H -> V
                    Hchan_vh = zeros([size(stvts.V,1) size(stvsr.H,1) Nray],'like',1+1i); % V -> H
                    Hchan_vv = zeros([size(stvts.V,1) size(stvsr.V,1) Nray],'like',1+1i); % V -> V
                    
                    for m = 1:Nray
                        x_h = x;
                        x_v = varargin{xidx(2)};

                        tempx_h = bsxfun(@times,x_h,exp(1i*2*pi*fd_ray(m)*t)); % doppler
                        tempx_v = bsxfun(@times,x_v,exp(1i*2*pi*fd_ray(m)*t)); % doppler
                        
                        % Transmit polarization
                        txsv_h = stvts.H(:,m).';
                        txsv_v = stvts.V(:,m).';
                        
                        % Translate to global coordinates
                        tempsv_global_h = sph2cartvec([txsv_h;zeros(size(txsv_h));zeros(size(txsv_h))],...
                            txang_ray(1,m),txang_ray(2,m));
                        tempsv_global_h = phased.internal.local2globalvec(tempsv_global_h,txaxes);
                        
                        tempsv_global_v = sph2cartvec([zeros(size(txsv_v));txsv_v;zeros(size(txsv_v))],...
                            txang_ray(1,m),txang_ray(2,m));
                        tempsv_global_v = phased.internal.local2globalvec(tempsv_global_v,txaxes);
                        
                        % Translate to target local coordinates
                        lcl_tgtinsv_h = phased.internal.global2localvec(...
                            tempsv_global_h,scatax_ray(:,:,m));
                        lcl_tgtinsv_h = cart2sphvec(lcl_tgtinsv_h,tgtinang_ray(1,m),tgtinang_ray(2,m));
                        lcl_tgtoutsv_h_sph = lcl_tgtinsv_h(1:2,:);
                        
                        lcl_tgtinsv_v = phased.internal.global2localvec(...
                            tempsv_global_v,scatax_ray(:,:,m));
                        lcl_tgtinsv_v = cart2sphvec(lcl_tgtinsv_v,tgtinang_ray(1,m),tgtinang_ray(2,m));
                        lcl_tgtoutsv_v_sph = lcl_tgtinsv_v(1:2,:);

                        % Apply path loss and scattering matrix
                        % rainpl polarization?
                        lcl_tgtoutsv_h_sph = amp(:,:,m)*lcl_tgtoutsv_h_sph;
                        
                        lcl_tgtoutsv_v_sph = amp(:,:,m)*lcl_tgtoutsv_v_sph;

                        % Translate to global coordinates
                        lcl_tgtoutsv_h = sph2cartvec(...
                            [lcl_tgtoutsv_h_sph;zeros(1,size(lcl_tgtoutsv_h_sph,2))],...
                            tgtoutang_ray(1,m),tgtoutang_ray(2,m));
                        tgtsv_global_h = phased.internal.local2globalvec(...
                            lcl_tgtoutsv_h,scatax_ray(:,:,m));
                       
                        lcl_tgtoutsv_v = sph2cartvec(...
                            [lcl_tgtoutsv_v_sph;zeros(1,size(lcl_tgtoutsv_v_sph,2))],...
                            tgtoutang_ray(1,m),tgtoutang_ray(2,m));
                        tgtsv_global_v = phased.internal.local2globalvec(...
                            lcl_tgtoutsv_v,scatax_ray(:,:,m));

                        % Receive polarization
                        rxsv_h = stvsr.H(:,m).';
                        rxsv_v = stvsr.V(:,m).';
                        
                        rxsv_local_h = sph2cartvec([rxsv_h;zeros(size(rxsv_h));zeros(size(rxsv_h))],...
                            rxang_ray(1,m),rxang_ray(2,m));
                        rxsv_global_h = phased.internal.local2globalvec(rxsv_local_h,rxaxes);
                        
                        rxsv_local_v = sph2cartvec([zeros(size(rxsv_v));rxsv_v;zeros(size(rxsv_v))],...
                            rxang_ray(1,m),rxang_ray(2,m));
                        rxsv_global_v = phased.internal.local2globalvec(rxsv_local_v,rxaxes);
                        
                        % Form channel matrix
                        tempchanmat_hh = rxsv_global_h.'*tgtsv_global_h;
                        Hchan_hh(:,:,m) = tempchanmat_hh.';
                        tempchanmat_hv = rxsv_global_h.'*tgtsv_global_v;
                        Hchan_vh(:,:,m) = tempchanmat_hv.';
                        
                        tempchanmat_vh = rxsv_global_v.'*tgtsv_global_h;
                        Hchan_hv(:,:,m) = tempchanmat_vh.';
                        tempchanmat_vv = rxsv_global_v.'*tgtsv_global_v;
                        Hchan_vv(:,:,m) = tempchanmat_vv.';
                        
                        % Compute output
                        tempy_h(:,(m-1)*Nrx+(1:Nrx)) = (tempchanmat_hh*tempx_h.').'+(tempchanmat_hv*tempx_v.').';
                        tempy_v(:,(m-1)*Nrx+(1:Nrx)) = (tempchanmat_vh*tempx_h.').'+(tempchanmat_vv*tempx_v.').';
                        
                    end
                    
                    temp_fracdelayed = step(obj.cFractionalDelayFilter,...
                        [tempy_h tempy_v],[delayvec_samp_frac delayvec_samp_frac]);
                    yr = step(obj.cBuffer,temp_fracdelayed,[delay_samp_int delay_samp_int]);
                    
                    y_h = complex(sum(reshape(yr(:,1:Nrx*Nray),Nsamp,Nrx,Nray),3));
                    y_v = complex(sum(reshape(yr(:,Nrx*Nray+1:end),Nsamp,Nrx,Nray),3));
                    
                    varargout{1} = y_h;
                    varargout{2} = y_v;
                    varargout{3} = Hchan_hh;
                    varargout{4} = Hchan_hv;
                    varargout{5} = Hchan_vh;
                    varargout{6} = Hchan_vv;
                    varargout{7} = delay;
                    
            end
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            reset(obj.cTxSteeringVector);
            reset(obj.cRxSteeringVector);
            reset(obj.cBuffer);
            reset(obj.cFractionalDelayFilter);
            reset(obj.cCoefficientSource);
            reset(obj.cPositionSource);
            obj.pScattererInitialized = false;
            obj.pCurrentTime = 0;
        end

        function releaseImpl(obj)
            % Release resources, such as file handles
            release(obj.cTxSteeringVector);
            release(obj.cRxSteeringVector);
            release(obj.cBuffer);
            release(obj.cFractionalDelayFilter);
            release(obj.cCoefficientSource);
            release(obj.cPositionSource);
        end

        function validateInputsImpl(obj,varargin)
            coder.extrinsic('num2str');            

            % Validate inputs to the step method at initialization
            xidx = getVararginIndex(obj,'sig');
            x = varargin{xidx};
            cond =  ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','X','double');
            end
            cond =  ~ismatrix(x);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeMatrix','X');
            end
            Nt = getDOF(obj.TransmitArray);
            cond = (size(x,2)~=Nt);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:phased:invalidColumnNumbers','X',num2str(Nt));
            end            
            validateNumChannels(obj,x);
            if strcmp(obj.TransmitArrayMotionSource,'Input port')
                txidx = getVararginIndex(obj,'tx');
                txpos = varargin{txidx(1)};
                txvel = varargin{txidx(2)};
                txaxes = varargin{txidx(3)};
                
                cond =  ~isa(txpos,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','TXPOS','double');
                end
                cond =  ~isreal(txpos);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'TXPOS');
                end
                cond = any(isnan(txpos(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'TXPOS');
                end
                cond = (size(txpos)~=[3 1]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:expectedSize', 'TXPOS', num2str([3 1]));
                end
                
                cond =  ~isa(txvel,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','TXVEL','double');
                end
                cond =  ~isreal(txvel);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'TXVEL');
                end
                cond = any(isnan(txvel(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'TXVEL');
                end
                cond = (size(txvel)~=[3 1]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:expectedSize', 'TXVEL', num2str([3 1]));
                end
                
                cond =  ~isa(txaxes,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','TXAXES','double');
                end
                cond =  ~isreal(txaxes);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'TXAXES');
                end
                cond = any(isnan(txaxes(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'TXAXES');
                end
                cond = (size(txaxes)~=[3 3]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:expectedSize', 'TXAXES', num2str([3 3]));
                end
            end
            if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                rxidx = getVararginIndex(obj,'rx');
                rxpos = varargin{rxidx(1)};
                rxvel = varargin{rxidx(2)};
                rxaxes = varargin{rxidx(3)};
                
                cond =  ~isa(rxpos,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','RXPOS','double');
                end
                cond =  ~isreal(rxpos);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'RXPOS');
                end
                cond = any(isnan(rxpos(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'RXPOS');
                end
                cond = (size(rxpos)~=[3 1]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:expectedSize', 'RXPOS', num2str([3 1]));
                end
                
                cond =  ~isa(rxvel,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','RXVEL','double');
                end
                cond =  ~isreal(rxvel);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'RXVEL');
                end
                cond = any(isnan(rxvel(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'RXVEL');
                end
                cond = (size(rxvel)~=[3 1]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:expectedSize', 'RXVEL', num2str([3 1]));
                end
                
                cond =  ~isa(rxaxes,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','RXAXES','double');
                end
                cond =  ~isreal(rxaxes);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'RXAXES');
                end
                cond = any(isnan(rxaxes(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'RXAXES');
                end
                cond = (size(rxaxes)~=[3 3]);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:expectedSize', 'RXAXES', num2str([3 3]));
                end
            end
            if strcmp(obj.ScattererSpecificationSource,'Input port')
                scidx = getVararginIndex(obj,'sc');
                scpos = varargin{scidx(1)};
                scvel = varargin{scidx(2)};
                
                cond =  ~isa(scpos,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','SCATPOS','double');
                end
                cond =  ~isreal(scpos);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'SCATPOS');
                end
                cond = any(isnan(scpos(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'SCATPOS');
                end
                cond = (size(scpos,1)~=3);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:invalidRowNumbers', 'SCATPOS', 3);
                end
                
                cond =  ~isa(scvel,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','SCATVEL','double');
                end
                cond =  ~isreal(scvel);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:NeedReal', 'SCATVEL');
                end
                cond = any(isnan(scvel(:)));
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:step:expectedNonNaN', 'SCATVEL');
                end
                cond = (size(scvel,1)~=3);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:phased:invalidRowNumbers', 'SCATVEL', 3);
                end
                
                cond = (size(scpos)~=size(scvel));
                if cond
                    coder.internal.errorIf(cond,'phased:phased:sizeMismatch',...
                        'SCATPOS','SCATVEL');
                end
                
                if strcmp(obj.Polarization,'None')
                    sccoef = varargin{scidx(3)};
                    cond =  ~isa(sccoef,'double');
                    
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','SCATCOEF','double');
                    end
                    cond = any(isnan(sccoef));
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:expectedNonNaN', 'SCATCOEF');
                    end
                    cond = any(isinf(sccoef));
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:expectedFinite', 'SCATCOEF');
                    end
                    cond = (size(sccoef,2)~=size(scpos,2));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:NumColumnsMismatch',...
                            'SCATPOS','SCATCOEF');
                    end                
                else
                    scmat = varargin{scidx(3)};
                    scaxes = varargin{scidx(4)};
                    npages = size(scpos,2);
                    
                    cond =  ~isa(scmat,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','SCATMAT','double');
                    end
                    cond = any(isnan(scmat(:)));
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:expectedNonNaN', 'SCATMAT');
                    end
                    cond = any(isinf(scmat(:)));
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:expectedFinite', 'SCATMAT');
                    end
                    cond = (size(scmat,3)~=npages);
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:invalidPageNumbers',...
                            'SCATMAT',npages);
                    end
                    
                    cond =  ~isa(scaxes,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','SCATAXES','double');
                    end
                    cond = any(isnan(scaxes(:)));
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:expectedNonNaN', 'SCATAXES');
                    end
                    cond = any(isinf(scaxes(:)));
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:expectedFinite', 'SCATAXES');
                    end
                    cond = ~all(isreal(scaxes(:)));
                    if cond
                        coder.internal.errorIf(cond, ...
                             'phased:step:NeedReal', 'SCATAXES');
                    end
                    cond = (size(scaxes,3)~=npages);
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:invalidPageNumbers',...
                            'SCATAXES',npages);
                    end     
                    if npages > 1
                        cond = (size(scaxes)~=[3 3 npages]);
                        if cond
                            coder.internal.errorIf(cond, ...
                                'phased:phased:expectedSize', 'SCATAXES', num2str([3 3 npages]));
                        end
                    else
                        cond = (size(scaxes)~=[3 3]);
                        if cond
                            coder.internal.errorIf(cond, ...
                                'phased:phased:expectedSize', 'SCATAXES', num2str([3 3]));
                        end
                    end
                end
                
            end
        end
        
        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
            if strcmp(obj.ScattererSpecificationSource,'Property')
                if strcmp(obj.Polarization,'None')
                    cond = ~isequal(size(obj.ScattererPosition,2),size(obj.ScattererCoefficient,2));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:NumColumnsMismatch',...
                            'ScattererPosition','ScattererCoefficient');
                    end
                else
                    npages = size(obj.ScattererPosition,2);
                    cond = ~isequal(npages,size(obj.ScatteringMatrix,3));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:invalidPageNumbers',...
                            'ScatteringMatrix',npages);
                    end
                    cond = ~isequal(npages,size(obj.ScattererOrientationAxes,3));
                    if cond
                        coder.internal.errorIf(cond,'phased:phased:NumColumnsMismatch',...
                            'ScattererOrientationAxes',npages);
                    end
                end
            end
            if strcmp(obj.ScattererSpecificationSource,'Auto') 
                cond = ((obj.NumScatterers == 0) && ~obj.SimulateDirectPath);
                if cond
                    coder.internal.errorIf(cond,'phased:phased:prop1MustBeWhenProp2Is',...
                        'SimulateDirectPath','true','NumScatterers','0');
                end
            end
        end

        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            if strcmp(prop,'MaximumDelay') && ...
                    strcmp(obj.MaximumDelaySource,'Auto')
                flag = true;
            elseif strcmp(obj.TransmitArrayMotionSource,'Input port') && ...
                    (strcmp(prop,'TransmitArrayPosition') || ...
                    strcmp(prop,'TransmitArrayOrientationAxes'))
                flag = true;
            elseif strcmp(obj.ReceiveArrayMotionSource,'Input port') && ...
                    (strcmp(prop,'ReceiveArrayPosition') || ...
                    strcmp(prop,'ReceiveArrayOrientationAxes'))
                flag = true;
            elseif strcmp(obj.ScattererSpecificationSource,'Input port') && ...
                    (strcmp(prop,'ScattererPosition') || ...
                    strcmp(prop,'ScattererCoefficient') || ...
                    strcmp(prop,'ScattererPositionBoundary') || ...
                    strcmp(prop,'ScatteringMatrix') || ...
                    strcmp(prop,'ScattererOrientationAxes') || ...
                    strcmp(prop,'NumScatterers'))
                flag = true;
            elseif strcmp(obj.ScattererSpecificationSource,'Auto') && ...
                    (strcmp(prop,'ScattererPosition') || ...
                    strcmp(prop,'ScatteringMatrix') || ...
                    strcmp(prop,'ScattererOrientationAxes') || ...
                    strcmp(prop,'ScattererCoefficient'))
                flag = true;
            elseif strcmp(obj.ScattererSpecificationSource,'Property') && ...
                    strcmp(obj.Polarization,'None') && ...
                    (strcmp(prop,'ScattererPositionBoundary') || ...
                    strcmp(prop,'ScatteringMatrix') || ...
                    strcmp(prop,'ScattererOrientationAxes') || ...
                    strcmp(prop,'NumScatterers'))
                flag = true;
            elseif strcmp(obj.ScattererSpecificationSource,'Property') && ...
                    ~strcmp(obj.Polarization,'None') && ...
                    (strcmp(prop,'ScattererPositionBoundary') || ...
                    strcmp(prop,'ScattererCoefficient') || ...
                    strcmp(prop,'NumScatterers'))
                flag = true;
            elseif ~strcmp(obj.ScattererSpecificationSource,'Auto') && ...
                    (strcmp(prop,'Seed') || strcmp(prop,'SeedSource'))
                flag = true;
            elseif strcmp(prop,'Seed') && ~strcmp(obj.SeedSource,'Property')
                flag = true;
            elseif ~obj.SpecifyAtmosphere && ...
                    (strcmp(prop,'DryAirPressure') || ...
                    strcmp(prop,'WaterVapourDensity') || ...
                    strcmp(prop,'LiquidWaterDensity') || ...
                    strcmp(prop,'RainRate') || ...
                    strcmp(prop,'Temperature') )
                flag = true;
            else
                flag = false;
            end
        end

        function idx = getVararginIndex(obj,inpname)
            if strcmp(obj.Polarization,'None')
                if strcmp(inpname,'sig')
                    idx = 1;
                elseif strcmp(inpname,'tx')
                    % pConstantTransmitArrayMotion is already false to access this
                    % code
                    % txpos, txvel, txaxes
                    idx = [2 3 4];
                elseif strcmp(inpname,'rx')
                    % ConstantReceiveArrayMotion is already false to access this
                    % code
                    % rxpos, rxvel, rxaxes
                    if strcmp(obj.TransmitArrayMotionSource,'Input port')
                        idx = [5 6 7];
                    else
                        idx = [2 3 4];
                    end
                elseif strcmp(inpname,'sc')
                    % ConstantScattererMotion is already false to access this
                    % code
                    % scatpos, scatvel
                    if strcmp(obj.TransmitArrayMotionSource,'Input port')
                        if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                            idx = [8 9 10];
                        else
                            idx = [5 6 7];
                        end
                    else
                        if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                            idx = [5 6 7];
                        else
                            idx = [2 3 4];
                        end
                    end

                end
            elseif strcmp(obj.Polarization,'Combined')
                if strcmp(inpname,'sig')
                    idx = 1;
                elseif strcmp(inpname,'tx')
                    % pConstantTransmitArrayMotion is already false to access this
                    % code
                    % txpos, txvel, txaxes
                    idx = [2 3 4];
                elseif strcmp(inpname,'rx')
                    % ConstantReceiveArrayMotion is already false to access this
                    % code
                    % rxpos, rxvel, rxaxes
                    if strcmp(obj.TransmitArrayMotionSource,'Input port')
                        idx = [5 6 7];
                    else
                        idx = [2 3 4];
                    end
                elseif strcmp(inpname,'sc')
                    % ConstantScattererMotion is already false to access this
                    % code
                    % scatpos, scatvel
                    if strcmp(obj.TransmitArrayMotionSource,'Input port')
                        if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                            idx = [8 9 10 11];
                        else
                            idx = [5 6 7 8];
                        end
                    else
                        if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                            idx = [5 6 7 8];
                        else
                            idx = [2 3 4 5];
                        end
                    end

                end
            else % dual
                if strcmp(inpname,'sig')
                    idx = [1 2];
                elseif strcmp(inpname,'tx')
                    % pConstantTransmitArrayMotion is already false to access this
                    % code
                    % txpos, txvel, txaxes
                    idx = [3 4 5];
                elseif strcmp(inpname,'rx')
                    % ConstantReceiveArrayMotion is already false to access this
                    % code
                    % rxpos, rxvel, rxaxes
                    if strcmp(obj.TransmitArrayMotionSource,'Input port')
                        idx = [6 7 8];
                    else
                        idx = [3 4 5];
                    end
                elseif strcmp(inpname,'sc')
                    % ConstantScattererMotion is already false to access this
                    % code
                    % scatpos, scatvel
                    if strcmp(obj.TransmitArrayMotionSource,'Input port')
                        if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                            idx = [9 10 11 12];
                        else
                            idx = [6 7 8 9];
                        end
                    else
                        if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                            idx = [6 7 8 9];
                        else
                            idx = [3 4 5 6];
                        end
                    end

                end
            end
        end
        

        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);

            % Set private and protected properties
            if isLocked(obj)
                s.cTxSteeringVector = saveobj(obj.cTxSteeringVector);
                s.cRxSteeringVector = saveobj(obj.cRxSteeringVector);
                s.cBuffer = saveobj(obj.cBuffer);
                s.cFractionalDelayFilter = saveobj(obj.cFractionalDelayFilter);
                s.cCoefficientSource = saveobj(obj.cCoefficientSource);
                s.cPositionSource = saveobj(obj.cPositionSource);
                s.pLambda = obj.pLambda;
                s.pCurrentTime = obj.pCurrentTime;
                s.pScattererPosition = obj.pScattererPosition;
                s.pScattererVelocity = obj.pScattererVelocity;
                s.pScattererCoefficient = obj.pScattererCoefficient;
                s.pScatteringMatrix = obj.pScatteringMatrix;
                s.pScattererInitialized = obj.pScattererInitialized;
                s.pReferenceDistance = obj.pReferenceDistance;
                s.pConstantTransmitArrayMotion = obj.pConstantTransmitArrayMotion;
                s.pConstantReceiveArrayMotion = obj.pConstantReceiveArrayMotion;
                s.pConstantScattererMotion = obj.pConstantScattererMotion;
                s.pSampleRate = obj.pSampleRate;
                s.pPolarizationConfigurationIndex = obj.pPolarizationConfigurationIndex;
                s.pScattererOrientationAxes = obj.pScattererOrientationAxes;
            end
                
        end

        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked 
                obj.cTxSteeringVector = ...
                    phased.SteeringVector.loadobj(s.cTxSteeringVector);
                s = rmfield(s,'cTxSteeringVector');
                obj.cRxSteeringVector = ...
                    phased.SteeringVector.loadobj(s.cRxSteeringVector);
                s = rmfield(s,'cRxSteeringVector');
                obj.cBuffer = ...
                    phased.internal.CircularBuffer.loadobj(s.cBuffer);
                s = rmfield(s,'cBuffer');
                obj.cFractionalDelayFilter = ...
                    dsp.VariableFractionalDelay.loadobj(s.cFractionalDelayFilter);
                s = rmfield(s,'cFractionalDelayFilter');
                obj.cCoefficientSource = ...
                    phased.internal.NoiseSource.loadobj(s.cCoefficientSource);
                s = rmfield(s,'cCoefficientSource');
                obj.cPositionSource = ...
                    phased.internal.NoiseSource.loadobj(s.cPositionSource);
                s = rmfield(s,'cPositionSource');
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s

            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        %% Simulink functions
        function flag = isInputComplexityLockedImpl(obj,index)  
            if ~strcmp(obj.Polarization,'Dual')
                if index == 1
                    flag = false;
                else 
                    if strcmp(obj.ScattererSpecificationSource,'Input port')
                        scidx = getVararginIndex(obj,'sc');
                        if index == scidx(3)
                            flag = true;
                        else
                            flag = false;
                        end
                    else
                        flag = true;
                    end
                end
            else
                if index == 1 || index == 2
                    flag = false;
                else 
                    if strcmp(obj.ScattererSpecificationSource,'Input port')
                        scidx = getVararginIndex(obj,'sc');
                        if index == scidx(3)
                            flag = true;
                        else
                            flag = false;
                        end
                    else
                        flag = true;
                    end
                end
            end
        end
        
        function flag = isInputSizeLockedImpl(obj,index) 
            % Return true if input size is not allowed to change while
            % system is running
            if ~strcmp(obj.Polarization,'Dual')
                if index == 1
                    flag = false;
                else
                    flag = true;
                end
            else % Dual polarization
                if (index==1)||(index==2)
                    flag = false;
                else
                    flag = true;
                end
            end
        end

        function num = getNumInputsImpl(obj)
            % Define total number of inputs for system with optional inputs
            if ~strcmp(obj.Polarization,'Dual')
                num = 1;
            else
                num = 2;
            end
            if strcmp(obj.TransmitArrayMotionSource,'Input port')
                num = num+3;
            end
            if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                num = num+3;
            end
            if strcmp(obj.ScattererSpecificationSource,'Input port')
                if strcmp(obj.Polarization,'None')
                    num = num+3;
                else  % Combined or Dual
                    num = num+4;
                end
            end
        end

        function num = getNumOutputsImpl(obj)
            % Define total number of outputs for system with optional
            % outputs
            if obj.ChannelResponseOutputPort
                if strcmp(obj.Polarization,'Dual')
                    num = 7;
                else
                    num = 3;
                end
            else
                if strcmp(obj.Polarization,'Dual')
                    num = 2;
                else
                    num = 1;
                end
            end
        end

        function varargout = getOutputSizeImpl(obj)
            % Return size for each output port
            sz = propagatedInputSize(obj,1);
            data_sz = [sz(1) getDOF(obj.ReceiveArray)];
            if obj.ChannelResponseOutputPort
                Nt = getDOF(obj.TransmitArray);
                Nr = getDOF(obj.ReceiveArray);
                if strcmp(obj.ScattererSpecificationSource,'Auto')
                    Ns = obj.NumScatterers;
                elseif strcmp(obj.ScattererSpecificationSource,'Property')
                    Ns = size(obj.ScattererPosition,2);
                else % Input port
                    sidx = getVararginIndex(obj,'sc');
                    sc_sz = propagatedInputSize(obj,sidx(1)+1); % data is the first input
                    Ns = sc_sz(2);
                end
                if obj.SimulateDirectPath
                    Ns = Ns+1;
                end
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {data_sz,[Nt Nr Ns],[1 Ns]}; 
                else
                    varargout = {data_sz,data_sz,[Nt Nr Ns],[Nt Nr Ns],[Nt Nr Ns],[Nt Nr Ns],[1 Ns]}; 
                end
            else
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {data_sz};
                else
                    varargout = {data_sz,data_sz};
                end
            end

        end

        function varargout = getOutputDataTypeImpl(obj)
            % Return data type for each output port
            dtype = propagatedInputDataType(obj,1);
            if obj.ChannelResponseOutputPort
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {dtype,'double','double'};
                else
                    varargout = {dtype,dtype,'double','double','double','double','double'};
                end
            else
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {dtype};
                else
                    varargout = {dtype,dtype};
                end
            end
        end

        function varargout = isOutputComplexImpl(obj) 
            % Return true for each output port with complex data
            if obj.ChannelResponseOutputPort
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {true true true};
                else
                    varargout = {true true true true true true true};
                end
            else
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {true};
                else
                    varargout = {true true};
                end
            end
        end

        function varargout = isOutputFixedSizeImpl(obj) 
            % Return true for each output port with fixed size
            flag = propagatedInputFixedSize(obj,1);
            if obj.ChannelResponseOutputPort
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {flag,true,true};
                else
                    varargout = {flag,flag,true,true,true,true,true};
                end
            else
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {flag};
                else
                    varargout = {flag,flag};
                end
            end
        end
        
        function icon = getIconImpl(obj) %#ok<MANU>
            % Return text as string or cell array of strings for the System
            % block icon
            icon = 'Scattering\nMIMO Channel'; % Use class name
        end

        function varargout = getInputNamesImpl(obj)
            % Return input port names for System block
            if ~strcmp(obj.Polarization,'Dual')
                signame = {'X'};
            else
                signame = {'XH','XV'};
            end
            
            if strcmp(obj.TransmitArrayMotionSource,'Input port')
                if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                    txrxname = {'TxPos','TxVel','TxAxes','RxPos','RxVel','RxAxes'};
                else
                    txrxname = {'TxPos','TxVel','TxAxes'};
                end
            else
                if strcmp(obj.ReceiveArrayMotionSource,'Input port')
                    txrxname = {'RxPos','RxVel','RxAxes'};
                else
                    txrxname = {};
                end
            end
            
            if strcmp(obj.ScattererSpecificationSource,'Input port')
                if strcmp(obj.Polarization,'None')
                    scatname = {'ScatPos','ScatVel','ScatCoef'};
                else  % Combined or Dual
                    scatname = {'ScatPos','ScatVel','ScatMat','ScatAxes'};
                end
            else
                scatname = {};
            end
                        
            varargout = [signame txrxname scatname];
        end

        function varargout = getOutputNamesImpl(obj)
            % Return output port names for System block
            if obj.ChannelResponseOutputPort
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {'Y','CS','Tau'};
                else
                    varargout = {'YH','YV','CSHH','CSHV','CSVH','CSVV','Tau'};
                end
            else
                if ~strcmp(obj.Polarization,'Dual')
                    varargout = {'Y'};
                else
                    varargout = {'YH','YV'};
                end
            end
        end
    end
    
    methods(Access = protected)
        function group = getPropertyGroups(obj)
            % Same display for short and long (no link), but showing long set
            % of properties from getPropertyGroupsImpl.
            dispGroups = getPropertyGroupsLongImpl(obj);
            group = dispGroups(1);
        end

        function group = getPropertyGroupsLongImpl(obj)
            % Move Sensor property to top of first tab to get it to be first in
            % the display. Call inherited version to convert
            % getPropertyGroupsImpl to MATLAB property groups. Sensor starts
            % out by itself on the third tab.
            allGroups = getPropertyGroupsLongImpl@phased.internal.AbstractSampleRateEngine(obj);
            tempDispPropList1 = [allGroups(1:3).PropertyList];
            propNameIdx = strcmp('ReceiveArray',tempDispPropList1);
            tempDispGroup = matlab.mixin.util.PropertyGroup(['ReceiveArray',tempDispPropList1(~propNameIdx)]);
            tempDispPropList2 = [tempDispGroup.PropertyList];
            propNameIdx = strcmp('TransmitArray',tempDispPropList2);
            dispGroup1 = matlab.mixin.util.PropertyGroup(['TransmitArray',tempDispPropList2(~propNameIdx)]);
            
            motionGroup = allGroups(4);
            TxMotionGroup = matlab.mixin.util.PropertyGroup(motionGroup.PropertyList(1:3));
            RxMotionGroup = matlab.mixin.util.PropertyGroup(motionGroup.PropertyList(4:6));
            ScatGroup = matlab.mixin.util.PropertyGroup(motionGroup.PropertyList(7:end));

            group = [dispGroup1 TxMotionGroup RxMotionGroup ScatGroup];
        end
    end

    methods(Static, Access = protected)
        %% Simulink customization functions
        function header = getHeaderImpl
            % Define header panel for System block dialog
            header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:ScatteringMIMOChannelTitle')),...
              'Text',getString(message('phased:library:block:ScatteringMIMOChannelDesc')));
        end

        function group = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
            
            txarrayProp = matlab.system.display.internal.Property('TransmitArray', ...
                'CustomPresenter', 'phased.internal.SensorArrayDialog', ...
                'CustomPresenterPropertyGroupsArgument', 'array', ...
                'Description', 'Transmit array');
            txarrayTab = matlab.system.display.SectionGroup('Title', 'Transmit Array', 'PropertyList', {txarrayProp});
            
            rxarrayProp = matlab.system.display.internal.Property('ReceiveArray', ...
                'CustomPresenter', 'phased.internal.SensorArrayDialog', ...
                'CustomPresenterPropertyGroupsArgument', 'array', ...
                'Description', 'Receive array');
            rxarrayTab = matlab.system.display.SectionGroup('Title', 'Receive Array', 'PropertyList', {rxarrayProp});

            dSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false,'UseClassDefault',false,...
                'Default','Property');
            dSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false,'UseClassDefault',false,...
                'Default','randi(65535,1)');
            
            dMaximumDelaySource = matlab.system.display.internal.Property(...
                'MaximumDelaySource','IsGraphical',false,...
                'UseClassDefault',false,'Default','Property');
            
            mainProps = {...
                'PropagationSpeed', ...
                'CarrierFrequency',...
                'Polarization',...
                'SpecifyAtmosphere',...
                'Temperature',...
                'DryAirPressure',...
                'WaterVapourDensity',...
                'LiquidWaterDensity',...
                'RainRate',...
                'SampleRateFromInputCheckbox',...
                'SampleRate',...
                'SimulateDirectPath',...
                dMaximumDelaySource,...
                'MaximumDelay',...
                'ChannelResponseOutputPort'};
            mainTab = matlab.system.display.SectionGroup('TitleSource', 'Auto', 'PropertyList', mainProps);

            motionProps = {...
                'TransmitArrayMotionSource', ...
                'TransmitArrayPosition',...
                'TransmitArrayOrientationAxes',...
                'ReceiveArrayMotionSource', ...
                'ReceiveArrayPosition',...
                'ReceiveArrayOrientationAxes',...
                'ScattererSpecificationSource', ...
                'NumScatterers',...
                'ScattererPositionBoundary'...
                'ScattererPosition',...
                'ScattererCoefficient',...
                'ScatteringMatrix',...
                'ScattererOrientationAxes',...
                dSeedSource,...
                dSeed};
            motionTab = matlab.system.display.SectionGroup('Title', 'Motion', 'PropertyList', motionProps);
            
            group = [mainTab,txarrayTab,rxarrayTab,motionTab];
        end
    end
end
