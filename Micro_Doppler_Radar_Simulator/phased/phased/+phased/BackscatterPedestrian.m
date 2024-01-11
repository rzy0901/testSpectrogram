classdef (Sealed) BackscatterPedestrian < backscatterPedestrian
%BackscatterPedestrian   Backscatter pedestrian model for radar
%
%   phased.BackscatterPedestrian may be removed in a future release. Use
%   backscatterPedestrian instead.
%
%   H = phased.BackscatterPedestrian creates a backscatter pedestrian, H,
%   for radar simulation. This object simulates the signal reflected off
%   the pedestrian while the pedestrian is in motion.
%
%   H = phased.BackscatterPedestrian(Name,Value) returns a backscatter
%   pedestrian model, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   move method syntax:
%
%   [BPPOS,BPVEL,BPAX] = move(H,T,ANGH) returns the pedestrian's body
%   segments position, velocity, and orientation axes at the current time
%   in BPPOS, BPVEL, and BPAX, respectively. It then simulates the walking
%   motion in the next duration, specified in T (in seconds). ANGH
%   specifies the current heading angle (in degrees) as a scalar. The angle
%   is in the xy-plane.
%
%   The pedestrian model includes 16 body segments: left and right feet,
%   left and right lower legs, left and right upper legs, left and right
%   hip, left and right lower arms, left and right upper arms, left and
%   right shoulders, neck, and head. Therefore BPPOS is a 3x16 matrix with
%   each column representing the position of the corresponding body
%   segments in the [x;y;z] format (in meters). BPVEL is also a 3x16 matrix
%   whose columns are velocities of corresponding body segments in the
%   [x;y;z] form (in m/s). BPAX is a 3x3x16 array whose pages are
%   orientation axes of the corresponding body segments. The three columns
%   represents the 3 axes and each column is in [x;y;z] format.
%
%   reflect method syntax:
%
%   Y = reflect(H,X,ANG) returns the reflected signal Y off the pedestrian
%   target due to the input signal X.
%
%   The human model consists of 16 body segments: left and right feet, left
%   and right lower legs, left and right upper legs, left and right hip,
%   left and right lower arms, left and right upper arms, left and right
%   shoulders, neck, and head. Each body segment is represented by a
%   cylinder. Therefore, X is a 16-column matrix whose columns are incident
%   signals to each body segment.
% 
%   ANG is a 2x16 matrix representing the signal's incident direction to
%   each body segment. Each column of ANG specifies the incident direction
%   of the corresponding signal in the form of an [AzimuthAngle;
%   ElevationAngle] pair (in degrees).
%
%   Y is a column vector containing the combined reflected signal from all
%   body segments. The number of rows in Y is the same as the number of
%   rows in X.
%
%   plot method syntax:
%
%   HF = plot(H) returns a plot of the position of the backscatter
%   pedestrian object H at the current time.
%
%   HF is the figure handle.
%
%   BackscatterPedestrian methods:
%
%   move     - Pedestrian motion (see above)
%   reflect  - Signal reflection off pedestrian (see above)
%   plot     - Plot pedestrian's position (see above) 
%   release  - Allow property value and input characteristics changes
%   clone    - Create pedestrian target object with same property values
%   reset    - Reset internal states of the pedestrian target
%
%   BackscatterPedestrian properties:
%
%   Height              - Pedestrian height
%   WalkingSpeed        - Pedestrian walking speed 
%   PropagationSpeed    - Propagation speed
%   OperatingFrequency  - Operating frequency 
%   InitialPosition     - Initial position 
%   InitialHeading      - Initial heading 
%
%   % Example:
%   %   Compute the reflected signal from a pedestrian moving along
%   %   x axis. Assume the radar works at 24 GHz and the signal has
%   %   a 300 MHz bandwidth. The signal is captured at the moment 
%   %   the pedestrian starts to move and 1 second into the
%   %   movement. Assume the radar is at the origin.
%
%   c = 3e8; bw = 3e8; fs = bw; fc = 24e9;
%   wav = phased.LinearFMWaveform('SampleRate',fs,'SweepBandwidth',bw);
%   x = wav();
%
%   chan = phased.FreeSpace('OperatingFrequency',fc,'SampleRate',fs,...
%          'TwoWayPropagation',true);
%   ped = phased.BackscatterPedestrian(...
%          'OperatingFrequency',fc,'InitialPosition',[100;0;0]);
%   rpos = [0;0;0];
%
%   % time 0
%   [bppos,bpvel,bpax] = move(ped,1,0);
%   plot(ped); 
%   xp = chan(repmat(x,1,16),rpos,bppos,[0;0;0],bpvel);
%   [~,ang] = rangeangle(rpos,bppos,bpax);
%   y0 = reflect(ped,xp,ang);
%
%   % time 1
%   [bppos,bpvel,bpax] = move(ped,1,0);
%   plot(ped); 
%   xp = chan(repmat(x,1,16),rpos,bppos,[0;0;0],bpvel);
%   [~,ang] = rangeangle(rpos,bppos,bpax);
%   y1 = reflect(ped,xp,ang);
%
%   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(wav));
%   ymf = mf([y0 y1]);
%   t = (0:size(ymf,1)-1)/fs;
%   figure
%   plot(t,abs(ymf));
%   xlabel('Time (s)'); ylabel('Magnitude'); title('Pedestrian Return')
%
%   See also phased, phased.BackscatterRadarTarget, phased.Platform.
    
%   Copyright 2018-2019 The MathWorks, Inc.

%   Reference
%   [1] Victor Chen, The Micro-Doppler Effect in Radar, Artech House, 2011
%   
%   [2] Boulic, Ronan, et al. A Global Human Walking Model with Real-time
%   Kinematic Personification, The Visual Computer: International Journal
%   of Computer Graphics, Vol. 6, Issue 6, Dec 1990

%#codegen

methods
    function obj = BackscatterPedestrian(varargin)
        obj@backscatterPedestrian(varargin{:});
    end
end

end