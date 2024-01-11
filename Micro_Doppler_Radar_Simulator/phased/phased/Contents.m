% Phased Array System Toolbox
% Version 4.0 (R2018b) 27-Jul-2018 
% 
% Table of Contents (TOC)
% -----------------------
%
% Array Analysis
%   az2broadside                   - Convert azimuth to broadside angle
%   broadside2az                   - Convert broadside angle to azimuth
%   phased.ArrayGain               - Sensor array gain
%   phased.ArrayResponse           - Sensor array response
%   phased.ElementDelay            - Sensor array element delay estimator
%   phased.SteeringVector          - Sensor array steering vector
%   pilotcalib                     - Array calibration using pilot sources
%   steervec                       - Sensor array steering vector
%   taylortaperc                   - Taylor nbar taper for circular aperture
%
% Array Antenna Elements
%   aperture2gain                      - Convert effective aperture to gain
%   gain2aperture                      - Convert gain to effective aperture
%   phased.CosineAntennaElement        - Cosine antenna
%   phased.CrossedDipoleAntennaElement - Crossed dipole antenna
%   phased.CustomAntennaElement        - Custom antenna 
%   phased.IsotropicAntennaElement     - Isotropic antenna
%   phased.ShortDipoleAntennaElement   - Short dipole antenna
% 
% Array Design
%   phased.ConformalArray              - Conformal array
%   phased.HeterogeneousConformalArray - Heterogeneous conformal array
%   phased.HeterogeneousULA            - Heterogeneous uniform linear array
%   phased.HeterogeneousURA            - Heterogeneous uniform rectangular array
%   phased.PartitionedArray            - Phased array partitioned into subarrays
%   phased.ReplicatedSubarray          - Phased array formed by replicated subarrays
%   phased.UCA                         - Uniform circular array
%   phased.ULA                         - Uniform linear array
%   phased.URA                         - Uniform rectangular array
%
% Array Acoustic Elements
%   phased.CustomMicrophoneElement          - Custom microphone
%   phased.OmnidirectionalMicrophoneElement - Omnidirectional microphone
%
% Beamformers
%   cbfweights                         - Narrowband conventional beamformer weights
%   lcmvweights                        - Narrowband LCMV beamformer weights
%   mvdrweights                        - Narrowband MVDR beamformer weights
%   phased.FrostBeamformer             - Frost beamformer
%   phased.GSCBeamformer               - Generalized Sidelobe Canceller
%   phased.LCMVBeamformer              - Narrowband LCMV beamformer
%   phased.MVDRBeamformer              - Narrowband MVDR beamformer
%   phased.PhaseShiftBeamformer        - Narrowband phase shift beamformer
%   phased.SubbandMVDRBeamformer       - Subband MVDR beamformer
%   phased.SubbandPhaseShiftBeamformer - Subband phase shift beamformer
%   phased.TimeDelayBeamformer         - Time delay beamformer
%   phased.TimeDelayLCMVBeamformer     - Time delay LCMV beamformer
%
% Clutter Models
%   billingsleyicm                      - Billingsley's intrinsic clutter motion (ICM) model  
%   depressionang                       - Depression angle of surface target
%   effearthradius                      - Effective earth radius
%   grazingang                          - Grazing angle of surface target
%   horizonrange                        - Horizon range
%   phased.ConstantGammaClutter         - Constant gamma clutter simulation
%   phased.gpu.ConstantGammaClutter     - Constant gamma clutter simulation using a GPU
%   surfacegamma                        - Gamma value for different terrains
%   surfclutterrcs                      - Surface clutter radar cross section
%
% Coordinate System and Motion Modeling
%   dop2speed                   - Convert Doppler shift to speed
%   global2localcoord           - Global to local coordinates conversion
%   local2globalcoord           - Local to global coordinates conversion
%   phased.Platform             - Motion platform
%   phased.ScenarioViewer       - Visualize platform trajectories
%   radialspeed                 - Relative radial speed
%   rangeangle                  - Range and angle calculation
%   speed2dop                   - Convert speed to Doppler shift
%
% Detection
%   albersheim                  - Albersheim's equation
%   npwgnthresh                 - Detection SNR threshold for white Gaussian noise
%   phased.AlphaBetaFilter      - Alpha-beta filter
%   phased.CFARDetector         - Constant false alarm rate (CFAR) detector
%   phased.CFARDetector2D       - Constant false alarm rate detector for 2-D training bands
%   phased.DopplerEstimator     - Doppler estimation for radar detections
%   phased.MatchedFilter        - Matched filter
%   phased.PulseCompressionLibrary - Pulse Compression Library
%   phased.RangeAngleResponse   - Range-angle response
%   phased.RangeDopplerResponse - Range-Doppler response
%   phased.RangeEstimator       - Range estimation for radar detections
%   phased.RangeResponse        - Range response
%   phased.StretchProcessor     - Stretch processor for linear FM waveforms
%   phased.TimeVaryingGain      - Time varying gain control
%   pulsint                     - Pulse integration
%   rocpfa                      - Receiver operating characteristic curves on varying Pfa
%   rocsnr                      - Receiver operating characteristic curves on varying SNR
%   shnidman                    - Shnidman's equation
%   stretchfreq2rng             - Convert frequency offset from stretch processing to range
%
% Direction of Arrival (DOA)
%   espritdoa                       - ESPRIT direction of arrival (DOA)
%   gccphat                         - GCC-PHAT delay estimation
%   musicdoa                        - MUSIC direction of arrival (DOA)
%   phased.BeamscanEstimator        - Beamscan spatial spectrum estimator for ULA
%   phased.BeamscanEstimator2D      - 2-D Beamscan spatial spectrum estimator
%   phased.BeamspaceESPRITEstimator - Beamspace ESPRIT direction of arrival (DOA) estimator
%   phased.ESPRITEstimator          - ESPRIT direction of arrival (DOA) estimator
%   phased.GCCEstimator             - GCC direction of arrival (DOA) estimator
%   phased.MonopulseFeed            - Amplitude monopulse feed
%   phased.MonopulseEstimator       - Amplitude monopulse estimator
%   phased.MUSICEstimator           - MUSIC spatial spectrum estimator for ULA
%   phased.MUSICEstimator2D         - MUSIC spatial spectrum estimator
%   phased.MVDREstimator            - MVDR spatial spectrum estimator for ULA
%   phased.MVDREstimator2D          - 2-D MVDR spatial spectrum estimator
%   phased.RootMUSICEstimator       - Root MUSIC direction of arrival (DOA) estimator
%   phased.RootWSFEstimator         - Root WSF direction of arrival (DOA) estimator
%   phased.SumDifferenceMonopulseTracker   - Sum and difference monopulse for ULA
%   phased.SumDifferenceMonopulseTracker2D - Sum and difference monopulse for URA
%   rootmusicdoa                    - Root MUSIC direction of arrival (DOA)
%   spsmooth                        - Spatial smoothing of a covariance matrix
%
% Environment Models
%   fspl                           - Free space path loss
%   fogpl                          - Path loss due to fog and cloud
%   gaspl                          - Path loss due to atmosphere gaseous absorption
%   phased.FreeSpace               - Free space environment
%   phased.LOSChannel              - Line of sight propagation channel
%   phased.ScatteringMIMOChannel   - Scattering based MIMO propagation channel
%   phased.TwoRayChannel           - Two-ray multipath propagation channel
%   phased.WidebandTwoRayChannel   - Wideband two-ray multipath propagation channel
%   phased.WidebandFreeSpace       - Wideband free space environment
%   phased.WidebandLOSChannel      - Wideband line of sight propagation channel
%   rainpl                         - Path loss due to rain
%
% Jammer Models
%   phased.BarrageJammer        - Barrage jammer
%
% MIMO Communication
%   diagbfweights                - Diagonalization beamforming weights
%   phased.ScatteringMIMOChannel - Scattering based MIMO propagation channel
%   scatteringchanmtx            - Scattering based MIMO channel matrix
%   waterfill                    - Power distribution using water-filling algorithm
%
% Polarization
%   circpol2pol                 - Circular to linear polarization representation conversion
%   pol2circpol                 - Linear to circular polarization representation conversion
%   polellip                    - Polarization ellipse
%   polloss                     - Polarization loss
%   polratio                    - Polarization ratio
%   polsignature                - Polarization signature 
%   stokes                      - Stokes parameters
%
% Radar Analysis
%   blakechart                  - Blake chart for radar
%   radareqpow                  - Radar equation to estimate power
%   radareqrng                  - Radar equation to estimate range
%   radareqsnr                  - Radar equation to estimate SNR
%   radarvcd                    - Vertical coverage diagram for radar
%
% Receiver Models
%   noisepow                    - Noise power at the receiver
%   phased.Collector            - Narrowband signal collector
%   phased.ReceiverPreamp       - Receiver preamp
%   phased.WidebandCollector    - Wideband signal collector
%   sensorcov                   - Covariance matrix of received signal
%   sensorsig                   - Received signal at sensor array
%   systemp                     - System noise temperature
%
% Sonar 
%   phased.IsoSpeedUnderwaterPaths - Constant speed underwater acoustic paths
%   phased.MultipathChannel        - Multipath channel
%   phased.IsotropicHydrophone     - Isotropic hydrophone element
%   phased.IsotropicProjector      - Isotropic projector element
%   phased.BackscatterSonarTarget  - Backscatter sonar point target
%   phased.UnderwaterRadiatedNoise - Underwater radiated noise
%   range2tl                       - Estimate transmission loss from range
%   sonareqsl                      - Sonar equation to estimate source level
%   sonareqsnr                     - Sonar equation to estimate SNR
%   sonareqtl                      - Sonar equation to estimate transmission loss
%   tl2range                       - Estimate range from transmission loss
%
% Space-Time Adaptive Processing (STAP)
%   dopsteeringvec              - Steering vector for Doppler
%   phased.ADPCACanceller       - Adaptive DPCA (ADPCA) pulse canceller
%   phased.AngleDopplerResponse - Angle-Doppler response
%   phased.DPCACanceller        - Displaced phase center array (DPCA) pulse canceller
%   phased.STAPSMIBeamformer    - Sample matrix inversion (SMI) STAP beamformer
%
% Target Models
%   phased.RadarTarget                    - Radar target
%   phased.BackscatterRadarTarget         - Backscatter point radar target
%   phased.WidebandBackscatterRadarTarget - Wideband backscatter point radar target
%
% Transmitter Models
%   phased.Radiator             - Narrowband signal radiator
%   phased.Transmitter          - Pulse transmitter
%   phased.WidebandRadiator     - Wideband signal radiator
%
% Utilities
%   aictest               - Akaike information criterion test
%   azelaxes              - Axes at given azimuth and elevation direction
%   azel2phitheta         - Convert angles from az/el format to phi/theta format
%   azel2phithetapat      - Convert pattern from az/el to phi/theta format
%   azel2uv               - Convert angles from az/el format to u/v format
%   azel2uvpat            - Convert pattern from az/el to u/v format
%   bw2range              - Convert bandwidth to range resolution
%   cart2sphvec           - Convert vector from Cartesian to Spherical representation
%   delayseq              - Delay or advance time sequence
%   mdltest               - Minimum description length test
%   phased.IntensityScope - Display scrolling intensity vectors or arrays
%   physconst             - Physical constants of natural phenomena
%   phitheta2azel         - Convert angles from phi/theta format to az/el format
%   phitheta2azelpat      - Convert pattern from phi/theta to az/el format
%   phitheta2uv           - Convert angles from phi/theta format to u/v format
%   phitheta2uvpat        - Convert pattern from phi/theta to u/v format
%   range2bw              - Convert range resolution to required bandwidth
%   range2time            - Convert propagation distance to propagation time
%   rotx                  - Rotation matrix around x-axis
%   roty                  - Rotation matrix around y-axis
%   rotz                  - Rotation matrix around z-axis
%   sph2cartvec           - Convert vector from Spherical to Cartesian representation
%   time2range            - Convert propagation time to propagation distance
%   unigrid               - Generate uniform grid
%   uv2azel               - Convert angles from u/v format to az/el format
%   uv2azelpat            - Convert pattern from u/v to az/el format
%   uv2phitheta           - Convert angles from u/v format to phi/theta format
%   uv2phithetapat        - Convert pattern from u/v to phi/theta format
%   val2ind               - Convert the grid value to grid index
%
% Waveforms
%   ambgfun                     - Ambiguity function
%   beat2range                  - Convert beat frequency to range
%   dechirp                     - Dechirp FMCW signal
%   pambgfun                    - Periodic ambiguity function
%   phased.FMCWWaveform         - FMCW waveform
%   phased.LinearFMWaveform     - Linear FM pulse waveform
%   phased.MFSKWaveform         - MFSK waveform
%   phased.PhaseCodedWaveform   - Phase-coded pulse waveform
%   phased.PulseWaveformLibrary - Pulse waveform library
%   phased.RectangularWaveform  - Rectangular pulse waveform
%   phased.SteppedFMWaveform    - Stepped FM pulse waveform
%   range2beat                  - Convert range to beat frequency
%   rdcoupling                  - Range Doppler coupling
% 
% Apps
%   radarWaveformAnalyzer       - Analyze performance characteristics of pulsed, 
%                                 frequency modulated, and phase coded waveforms
%   radarEquationCalculator     - Estimate maximum range, peak power and 
%                                 SNR of a radar system
%   sensorArrayAnalyzer         - Analyze beam pattern of linear, planar, 
%                                 and conformal sensor arrays
%   sonarEquationCalculator     - Estimates target range, transmission loss,
%                                 source level and SNR of a sonar system
%
%   <a href="matlab:demo toolbox phased">Examples</a>                    - Phased Array System Toolbox examples
%   <a href="matlab:phasedlib">Simulink library</a>            - Open Phased Array System Toolbox Simulink library
%
% See also SIGNAL, DSP.

%   Copyright 2006-2018 The MathWorks, Inc.
%    

