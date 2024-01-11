classdef (Sealed,StrictDefaults) WidebandTwoRayChannel < phased.internal.AbstractLOSChannel
%WidebandTwoRayChannel Wideband two-ray multipath propagation channel
%   H = phased.WidebandTwoRayChannel creates a wideband two-ray channel
%   System object, H. This object simulates wideband signal propagation
%   through a two-ray channel over a flat earth. A two-ray channel is a
%   homogeneous, isotropic medium like free space but with a reflecting
%   boundary (often the ground). A two-ray channel propagates a signal from
%   a source to a destination via two paths. The first path, the direct
%   path, is the same as in the free space. The second path, the ground
%   reflected path, propagates from the source to the boundary and then to
%   the destination. The object applies range-dependent time delay, gain
%   and phase shift for both the direct path and the ground reflection path
%   to the input signal. The propagated signal can either be specified as a
%   single signal for both paths, or as two separate signals along each
%   path. The signal is divided into subbands and in each subband the
%   signal is propagated by applying range-dependent time delay, gain and
%   phase shift to the input signal. The resulting signal from all subbands
%   are combined to form the output signal.
%
%   H = phased.WidebandTwoRayChannel(Name,Value) returns a wideband two-ray
%   propagation channel object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,POS1,POS2,VEL1,VEL2) returns the resulting signal Y when
%   the narrowband signal X propagates in free space from either the
%   position of a source POS1 to the positions of one or more destinations
%   POS2 or from the positions of several sources to the position of one
%   destination. POS1 and POS2 have a size of either 3x1 and 3xN
%   respectively or 3xN and 3x1 where each column is in the form of [x; y;
%   z] (in meters). N is the number of positions in POS1 or POS2. VEL1 and
%   VEL2 specify the velocities of the sources and destinations
%   respectively. Similarly VEL1 and VEL2 match the dimensions of POS1 and
%   POS2, respectively. Each column is in the form of [Vx; Vy; Vz] (in
%   meters/second). Note that the reflection boundary (ground) is assumed
%   to be the x-y plane.
%
%   If you set the EnablePolarization property to false, X can be either an
%   MxN matrix or an Mx2N matrix.
%
%   To propagate a single signal along both paths and combine the result at
%   the destination automatically, specify X as an MxN matrix. This happens
%   when the transmit antenna is isotropic or when the two rays are so
%   close that the beam patterns can be considered as identical between the
%   two directions.
%
%   To propagate different signals along two paths, specify X as an Mx2N
%   matrix. In this configuration, X is divided into N 2-column pairs.
%   Within each pair, the first column represents the signal along the
%   direct path and the second column represents the signal along the
%   reflected path.
%
%   If you set the CombinedRaysOutput property to true, Y is an MxN matrix
%   where each column is the combined propagated signal from both the
%   direct path and the reflected path along a propagation path. This is a
%   valid approximation when the receiving antenna is isotropic or when the
%   incoming direction at the receiver is so close that the beam patterns
%   can be considered as identical between the two directions. The
%   propagation paths are defined in the order of the positions specified
%   in POS1 and POS2.
%
%   If you set the CombinedRaysOutput property to false, Y is an Mx2N
%   matrix. The adjacent columns in Y are the propagated signals via the
%   direct path and the reflected path, respectively. 
%
%   If you set the EnablePolarization property to true, X can be either a
%   1xN struct array or a 1x2N struct array. Each element of the struct
%   array contains X, Y, and Z fields. Each field is an Mx1 column vector.
%   The X, Y, and Z fields in the struct specifies the field components
%   along x, y, and z axes in space.
%
%   To propagate a single signal along both paths and combine the result at
%   the destination automatically, specify X as a 1xN struct array. This
%   happens when the transmit antenna is isotropic or when the two rays are
%   so close that the beam patterns can be considered as identical between
%   the two directions.
%
%   To propagate different signals along two paths, specify X as a 1x2N
%   struct array. In this configuration, X is divided into N 2-element
%   pairs. In each pair, the first element represents the signal along the
%   direct path and the second element represents the signal along the
%   reflected path. 
%
%   If you set the CombinedRaysOutput property to true, Y is a 1xN struct
%   array where each element represents the combined propagated signal from
%   both the direct path and the reflected path along a propagation path.
%   This is a valid approximation when the receiving antenna is isotropic
%   or when the incoming direction at the receiver is so close that the
%   beam patterns can be considered as identical between the two
%   directions. The propagation paths are defined in the order of the
%   positions specified in POS1 and POS2.
%
%   If you set the CombinedRaysOutput property to false, Y is a 1x2N
%   struct. The adjacent elements in Y are the propagated signals via the
%   direct path and the reflected path, respectively.
%
%   The output Y represents the signals arriving at the propagation
%   destinations within the current time frame, which is the time occupied
%   by the current input. If it takes longer than the current time frame
%   for the signals to propagate from the origin to the destination, then
%   the output contains no contribution from the input of the current time
%   frame. The output Y can be written as
%
%   Y(t) = X(t-tau)/L
%
%   where tau is the delay and L is the propagation loss. The delay tau can
%   be calculated as R/c where R is the propagation distance and c is the
%   propagation speed. The propagation loss includes free space path loss
%   as well as the loss due to atmosphere and weather. R is the either the
%   direct path propagation distance or the reflected path propagation
%   distance.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   WidebandTwoRayChannel methods:
%
%   step     - Propagate signal from one location to another (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create wideband two-ray channel object with same property 
%              values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the propagation channel
%
%   WidebandTwoRayChannel properties:
%
%   PropagationSpeed            - Propagation speed 
%   OperatingFrequency          - Signal carrier frequency 
%   SpecifyAtmosphere           - Specify atmosphere parameters
%   Temperature                 - Temperature 
%   DryAirPressure              - Dry air pressure 
%   WaterVapourDensity          - Water vapour density 
%   LiquidWaterDensity          - Liquid water density 
%   RainRate                    - Rain rate 
%   SampleRate                  - Sample rate 
%   NumSubbands                 - Number of subbands
%   EnablePolarization          - Enable polarization
%   GroundReflectionCoefficient - Ground reflection coefficient
%   GroundRelativePermittivity  - Ground relative permittivity
%   CombinedRaysOutput          - Combine two rays at output
%   MaximumDistanceSource       - Source of maximum one-way propagation 
%                                 distance
%   MaximumDistance             - Maximum one-way propagation distance
%   MaximumNumInputSamplesSource - Source of maximum number of samples
%                                  of the input signal
%   MaximumNumInputSamples       - Maximum number of samples in input 
%                                  signal
%
%   % Examples:
%
%   % Example 1:
%   %   Calculate the result of propagating a wideband signal in a two-ray 
%   %   environment from a radar at (0, 0, 1) to a target at (300, 200,
%   %   1). Assume both the radar and the target are stationary and the 
%   %   transmitting antenna is isotropic. Combine the signal from two
%   %   paths and compare it to the signal resulting from a free space 
%   %   propagation.
%
%   waveform = phased.LinearFMWaveform;
%   x = waveform();
%   posTx = [0; 0; 1]; posTgt = [300; 200; 1]; 
%   velTx = [0;0;0]; velTgt = [0;0;0];
%
%   % free space propagation
%   env_fs = phased.FreeSpace('SampleRate',waveform.SampleRate);
%   y_fs = env_fs(x,posTx,posTgt,velTx,velTgt);
%
%   % two-ray propagation
%   env_tworay = phased.WidebandTwoRayChannel(...
%           'SampleRate',waveform.SampleRate);
%   y_tworay = env_tworay(x,posTx,posTgt,velTx,velTgt);
%
%   plot(abs([y_tworay y_fs])); legend('Two-ray','Free space');
%   xlabel('Samples'); ylabel('Signal Magnitude');
%
%   % Example 2:
%   %   Calculate the result of propagating a wideband signal in a two-ray 
%   %   environment from a radar at (0, 0, 10) to a target at (300, 200,
%   %   30). Assume both the radar and the target are stationary and the 
%   %   transmitting antenna has a cosine pattern. Combine the signal 
%   %   from two paths and compare it to the signal resulting from a free 
%   %   space propagation. 
%   %
%   %   Note that because the two rays have different gains from the 
%   %   antenna, the two rays need to be considered separately.
%
%   waveform = phased.LinearFMWaveform;
%   x = waveform();
%   posTx = [0; 0; 10]; posTgt = [300; 200; 30]; 
%   velTx = [0;0;0]; velTgt = [0;0;0];
%
%   ant = phased.CosineAntennaElement;
%   antTx = phased.Radiator('Sensor',ant);
%
%   % free space propagation
%
%   % compute the transmit angle toward the target for line of sight (LOS)
%   [~,angLOS] = rangeangle(posTgt,posTx);
%   xLOS = antTx(x,angLOS);         % single column, combined signal
%
%   env_fs = phased.FreeSpace('SampleRate',waveform.SampleRate);
%   y_fs = env_fs(xLOS,posTx,posTgt,velTx,velTgt);
%
%   % two-ray propagation
%   release(antTx);
%
%   % compute the transmit angles toward the target for the two rays
%   [~,angTwoRay] = rangeangle(posTgt,posTx,'two-ray');
%   xTwoRay = antTx(x,angTwoRay);   % two columns, separate rays
%
%   env_tworay = phased.WidebandTwoRayChannel(...
%           'SampleRate',waveform.SampleRate,'CombinedRaysOutput',true);
%   y_tworay = env_tworay(xTwoRay,posTx,posTgt,velTx,velTgt);
%
%   plot(abs([y_tworay y_fs])); legend('Two-ray','Free space');
%   xlabel('Samples'); ylabel('Signal Magnitude');
%
%   % Example 3:
%   %   Calculate the result of propagating a wideband signal in a two-ray 
%   %   environment from a transmitter at (0, 0, 10) to a receiver at (300,
%   %   200,30). Assume both the transmitter and the receiver are 
%   %   stationary. Both the transmit antenna and the receive antenna have
%   %   cosine patterns. The transmit antenna points to the +x direction 
%   %   and the receive antenna points to the -x direction. Simulate the
%   %   received signal after the propagation.
%   %
%   %   Note that because the two rays have different gains from the 
%   %   antenna, the two rays need to be considered separately in both the
%   %   transmitter and the receiver.
%
%   waveform = phased.LinearFMWaveform;
%   x = waveform();
%   posTx = [0; 0; 10]; posRx = [300; 200; 30]; 
%   velTx = [0;0;0]; velRx = [0;0;0];
%   axTx = eye(3); axRx = rotz(180)*eye(3);
%
%   ant = phased.CosineAntennaElement;
%   antTx = phased.Radiator('Sensor',ant);
%   antRx = phased.Collector('Sensor',ant);
%   
%   % two-ray propagation
%
%   % compute the transmit angles toward the receiver for the two rays
%   [~,angT] = rangeangle(posRx,posTx,axTx,'two-ray');
%   xT = antTx(x,angT);          % two columns, separate rays
%
%   env_tworay = phased.WidebandTwoRayChannel(...
%           'SampleRate',waveform.SampleRate,'CombinedRaysOutput',false);
%   y_tworay = env_tworay(xT,posTx,posRx,velTx,velRx);
%
%   % compute the receive angles toward the transmitter for the two rays
%   [~,angR] = rangeangle(posTx,posRx,axRx,'two-ray');
%   yR = antRx(y_tworay,angR);   % two columns, separate rays
%
%   plot(real(yR)); 
%   xlabel('Samples'); ylabel('Signal Magnitude');
%
%   See also phased, phased.LOSChannel, phased.RadarTarget,
%   phased.TwoRayChannel, phased.WidebandLOSChannel, fspl.

%   Copyright 2015-2016 The MathWorks, Inc.

%   Reference
%   [1] Artem Saakian, Radio Wave Propagation Fundamentals, Artech House,
%       2011
%   [2] Constantine Balanis, Advanced Engineering Electromagnetics, John
%       Wiley & Sons, 1989
%   [3] Theodore Rappaport, Wireless Communications: Principles and
%       Practice, 2nd Ed, Prentice Hall, 2002
%   [4] John Seybold, Introduction to RF Propagation, Wiley, 2005
%   [5] Recommendation ITU-R P.838-3, 2005
%   [6] Recommendation ITU-R P.840-3, 2013
%   [7] Recommendation ITU-R P.676-10, 2013


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable, Logical)
        %EnablePolarization  Enable polarization
        %   Set this property to true to enable polarization. Set this
        %   property to false to ignore polarization. The default value of
        %   this property is false. 
        EnablePolarization = false
    end
    
    properties (Nontunable, PositiveInteger)
        %NumSubbands    Number of subbands
        %   Specify the number of subbands used in the subband processing
        %   as a positive integer. The default value of this property is
        %   64.
        NumSubbands = 64
    end
    
    properties (Nontunable)
        %GroundReflectionCoefficient Ground reflection coefficient
        %   Specify the ground reflection coefficient for the field at the
        %   reflection point as a scalar or an N-element row vector where N
        %   is the number of source and destination position pairs. This
        %   property applies when you set the EnablePolarization property
        %   to false. The default value of this property is -1.
        GroundReflectionCoefficient = -1
        %GroundRelativePermittivity Ground relative permittivity
        %   Specify the relative permittivity of the ground at the
        %   reflection point as a positive scalar or an N-element row
        %   vector where N is the number of source and destination position
        %   pairs. This property applies when you set the
        %   EnablePolarization property to true. The default value of this
        %   property is 15.
        GroundRelativePermittivity = 15
    end
    
    properties (Nontunable, Logical)
        %CombinedRaysOutput     Combine two rays at output
        %   Specify whether to combine the propagated signal via the
        %   direct path and the propagated signal via the reflected path
        %   when forming the output signal. The default value is true.
        CombinedRaysOutput = true
    end

    properties (Access = private, Nontunable)
        %Fractional delay filter
        cFractionalDelayFilter
        %subband divider
        cSubbandDivider
        %subband combiner
        cSubbandCombiner
        %subband frequencies
        pSubbandFreq
    end
    
    properties (Access = private, Nontunable)
        %pSeparateInputRays  flag to mention whether the two rays are
        %separately processed
        pSeparateInputRays
    end
    
    methods
        function set.GroundReflectionCoefficient(obj,val)
            validateattributes(val,{'double'},{'row','finite'},...
                'WidebandTwoRayChannel','GroundReflectionCoefficient');
            obj.GroundReflectionCoefficient = val;
        end
        function set.GroundRelativePermittivity(obj,val)
            % the real part should be at least 1, MATLAB's '>' happen to do
            % that for a complex number
            validateattributes(val,{'double'},{'row','finite','>=',1},...
                'WidebandTwoRayChannel','GroundReflectionCoefficient');
            obj.GroundRelativePermittivity = val;
        end
    end
    
    methods
        function obj = WidebandTwoRayChannel(varargin)
            obj@phased.internal.AbstractLOSChannel(varargin{:});
        end
    end

    methods (Access = protected)

        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractLOSChannel(obj, prop);
            if (obj.EnablePolarization && ...
                    strcmp(prop,'GroundReflectionCoefficient')) || ...
                    (~obj.EnablePolarization && ...
                    strcmp(prop,'GroundRelativePermittivity'))
                flag = true;
            end
        end
        
        function validateInputsImpl(obj,x,startLoc,endLoc,baseVel,targetVel)
            validateInputsImpl@phased.internal.AbstractLOSChannel(obj,x,startLoc,endLoc,baseVel,targetVel);
            num_pospair = max(size(startLoc,2),size(endLoc,2));
            if ~obj.EnablePolarization
                cond = ~isscalar(obj.GroundReflectionCoefficient) && ...
                    numel(obj.GroundReflectionCoefficient) ~= num_pospair;
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:phased:invalidColumnDimension','GroundReflectionCoefficient',1,num_pospair);
                end
            else
                cond = ~isscalar(obj.GroundRelativePermittivity) && ...
                    numel(obj.GroundRelativePermittivity) ~= num_pospair;
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:phased:invalidColumnDimension','GroundRelativePermittivity',1,num_pospair);
                end
            end
        end
        
        function validateInputSignal(obj,x) 
            if obj.EnablePolarization
                cond = ~isstruct(x);
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','X','struct');
                end
                flag_hasXYZ = isfield(x(1),'X') && isfield(x(1),'Y') && isfield(x(1),'Z');
                cond = ~flag_hasXYZ;
                if cond
                    coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationXYZStruct');
                end
                % cond = ~isscalar(x);
                % if cond
                %     coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationArrayStruct','X');
                % end
                if flag_hasXYZ
                    xsize = size(x);
                    for m = 1:xsize(2)
                        x_x = x(m).X;
                        x_y = x(m).Y;
                        x_z = x(m).Z;
                        cond =  ~isa(x_x,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).X',m),'double');
                        end
                        cond =  ~iscolumn(x_x) || isempty(x_x);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).X',m));
                        end
                        cond =  ~isa(x_y,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).Y',m),'double');
                        end
                        cond =  ~iscolumn(x_y) || isempty(x_y);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Y',m));
                        end
                        cond =  ~isa(x_z,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType',sprintf('X(%d).Z',m),'double');
                        end
                        cond =  ~iscolumn(x_z) || isempty(x_z);
                        if cond
                            coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeColVector',sprintf('X(%d).Z',m));
                        end
                        cond = numel(x_x)~=numel(x_y) || numel(x_x)~=numel(x_z);
                        if cond
                            coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                                'X,Y,Z',sprintf('X(%d)',m));
                        end
                    end
                end
            else
                cond =  ~isa(x,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                         'MATLAB:system:invalidInputDataType','X','double');
                end
            end
        end
        
        function setupImpl(obj,x,startLoc,endLoc,startVel,endVel)
            setupImpl@phased.internal.AbstractLOSChannel(obj,x,startLoc,endLoc,startVel,endVel);
            obj.cFractionalDelayFilter = dsp.VariableFractionalDelay;
            obj.cSubbandDivider = phased.internal.SubbandDivider(...
                'SampleRate',obj.pSampleRate,...
                'NumSubbands',obj.NumSubbands,...
                'OperatingFrequency',obj.OperatingFrequency);
            obj.cSubbandCombiner = phased.internal.SubbandCombiner(...
                'NumSubbands',obj.NumSubbands,...
                'TimeSignalLengthSource','Inherit',...
                'SubbandFrequencyShiftInputPort',true);
            subbandfreq = phased.internal.subbandCenterFrequency(...
                obj.OperatingFrequency,obj.pSampleRate,obj.NumSubbands);
            obj.pLambda = obj.PropagationSpeed./subbandfreq;
            obj.pSubbandFreq = subbandfreq;
           
            ncol_x = size(x,2);
            ncol_pos = max(size(startLoc,2),size(endLoc,2));
            
            if ncol_x == ncol_pos
                obj.pSeparateInputRays = false;
            else % ncol_x == 2*ncol_pos
                obj.pSeparateInputRays = true;
            end
            
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractLOSChannel(obj);
            release(obj.cFractionalDelayFilter);
            release(obj.cSubbandCombiner);
        end

        function resetImpl(obj)
            resetImpl@phased.internal.AbstractLOSChannel(obj);
            reset(obj.cFractionalDelayFilter);
            reset(obj.cSubbandCombiner);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractLOSChannel(obj);
            if isLocked(obj)
                s.cFractionalDelayFilter = saveobj(obj.cFractionalDelayFilter);
                s.pSeparateInputRays = obj.pSeparateInputRays;
                s.cSubbandDivider = saveobj(obj.cSubbandDivider);
                s.cSubbandCombiner = saveobj(obj.cSubbandCombiner);
                s.pSubbandFreq = obj.pSubbandFreq;
            end
        end

        function s = loadSubObjects(obj,s,wasLocked)
            s = loadSubObjects@phased.internal.AbstractLOSChannel(obj,s,wasLocked);
            if wasLocked 
                obj.cFractionalDelayFilter = ...
                    dsp.VariableFractionalDelay.loadobj(s.cFractionalDelayFilter);
                s = rmfield(s,'cFractionalDelayFilter');
                obj.cSubbandDivider = ...
                    phased.internal.SubbandDivider.loadobj(s.cSubbandDivider);
                s = rmfield(s,'cSubbandDivider');
                obj.cSubbandCombiner = ...
                    phased.internal.SubbandCombiner.loadobj(s.cSubbandCombiner);
                s = rmfield(s,'cSubbandCombiner');
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function [pos1,pos2,vel1,vel2] = preparePositionVelocity(obj,pos1_in,pos2_in,vel1_in,vel2_in)  %#ok<INUSL>
            cond = any(pos1_in(3,:).*pos2_in(3,:) < 0);
            if cond
                coder.internal.errorIf(cond,...
                    'phased:rangeangle:invalidHeight','Pos1','Pos2');
            end
            np1 = size(pos1_in,2);
            np2 = size(pos2_in,2);
            if np1 ~= np2  % one of pos represents multiple positions
                if size(pos1_in,2) > 1 
                    pos1m = pos1_in;
                    pos2m = repmat(pos2_in,1,np1);
                    vel1m = vel1_in;
                    vel2m = repmat(vel2_in,1,np1);
                else
                    pos1m = repmat(pos1_in,1,np2);
                    pos2m = pos2_in;
                    vel1m = repmat(vel1_in,1,np2);
                    vel2m = vel2_in;
                end
            else
                pos1m = pos1_in;
                pos2m = pos2_in;
                vel1m = vel1_in;
                vel2m = vel2_in;
            end
            
            pos1 = reshape([pos1m; pos1m],3,[]); 
            
            pos2temp = [pos2m; pos2m];
            pos2temp(6,:) = -pos2temp(6,:);
            pos2 = reshape(pos2temp,3,[]); 
            
            vel1 = reshape([vel1m; vel1m],3,[]);
            vel2 = reshape([vel2m; vel2m],3,[]);
        end
        
        function [xbuf_in,ndelay] = computeDelayedSignal(obj,x,delay) 
            % For between sample delays, compute nDelay and frac
            % delay the signal
            intd = fix(delay);
            fracd = delay-intd; % in samples
            ndelay = intd;
            xbuf_in = step(obj.cFractionalDelayFilter,x,fracd);
        end
        
        function [y,propdelay] = computeMultiplePropagatedSignal(obj,x,...
                ncol_per_path,numOfPropPaths,startLoc,endLoc,baseVel,targetVel,Fs) 
            y = complex(zeros(size(x)));
            k = obj.pRangeFactor;
            lambda = obj.pLambda;
            [propdelay,propdistance,rspeed] = computePropagationDelayVelocity(...
                obj,startLoc,endLoc,baseVel,targetVel);
            
            if ~obj.SpecifyAtmosphere
                % spreading loss
                sploss = k*phased.internal.fspl(propdistance,lambda); 
                plossfactordb = sploss;
            else
                fcsb = obj.pSubbandFreq;
                % spreading loss
                sploss = k*phased.internal.fspl(propdistance,lambda);
                % gas loss
                fcsb_gas = min(max(fcsb,1e9),1000e9);
                gasloss = k*gaspl(propdistance,fcsb_gas,...
                    obj.Temperature,obj.DryAirPressure,obj.WaterVapourDensity);
                plossfactordb = sploss+gasloss;
                % fog loss
                if obj.LiquidWaterDensity > 0
                    fcsb_fog = min(max(fcsb,10e9),1000e9);
                    fogloss = k*fogpl(propdistance,fcsb_fog,...
                        obj.Temperature,obj.LiquidWaterDensity);
                    plossfactordb = plossfactordb+fogloss;
                end
                if obj.RainRate > 0
                    % rain loss parameters
                    el = obj.computeElevationAngle(startLoc,endLoc,propdistance);
                    if ~obj.pIsInputStruct
                        tau = zeros(1,numOfPropPaths);
                    else
                        switch obj.pValidFields
                            case 'XYZ'
                                tau = computePolarizationTiltAngle(obj,reshape(mean(x),3,[]),startLoc,endLoc);
                            case 'HV'
                                tau = computePolarizationTiltAngle(obj,reshape(mean(x),2,[]),startLoc,endLoc);
                            case 'XYZHV'
                                tau = computePolarizationTiltAngle(obj,reshape(mean(x),5,[]),startLoc,endLoc);
                        end
                    end
                    % rain loss
                    fcsb_rain = min(max(fcsb,1e9),1000e9);
                    rainloss = k*rainpl(propdistance,fcsb_rain,...
                        obj.RainRate,el,tau);
                    plossfactordb = plossfactordb+rainloss;
                end
            end
                       
            if obj.EnablePolarization
                % compute rotation angle
                rot_ang = atan2d(endLoc(2,:)-startLoc(2,:),...
                    endLoc(1,:)-startLoc(1,:))-90;  % deg, y-z plane is reference
                % compute incident angle
                inc_ang = asin(abs(endLoc(3,:)-startLoc(3,:))./propdistance); % rad

                % compute perpendicular and parallel reflection
                % coefficients
                er_temp = obj.GroundRelativePermittivity;
                if isscalar(er_temp)
                    er = repmat(er_temp,1,numOfPropPaths);
                else
                    % must be half of numOfPropPaths
                    er = reshape([er_temp;er_temp],1,[]);
                end
                sa = sin(inc_ang);
                ca = sqrt(er-cos(inc_ang).^2);
                esa = er.*sa;
                refcoef_perp = (sa-ca)./(sa+ca);
                refcoef_para = (-esa+ca)./(esa+ca);
            else
                rc_temp = obj.GroundReflectionCoefficient;
                if isscalar(rc_temp)
                    refcoef = repmat(rc_temp,1,numOfPropPaths);
                else
                    refcoef = reshape([rc_temp;rc_temp],1,[]);
                end
            end
            
            fc = obj.OperatingFrequency;
            for pIdx = 1:numOfPropPaths
                plossfactor = sqrt(db2pow(plossfactordb(pIdx,:)));
                fd = k*rspeed(pIdx)./lambda(:).';
                
                wloss = 1./plossfactor;
                wphase = exp(-1i*2*pi*fc*propdelay(pIdx));
                
                if rem(pIdx,2)
                    colidx = (pIdx-1)*ncol_per_path+(1:ncol_per_path);
                    xsubband = step(obj.cSubbandDivider,x(:,colidx));
                    for m = 1:size(xsubband,3)
                        xsubband(:,:,m) = wloss(m)*xsubband(:,:,m);
                    end
                    ytemp = wphase*step(obj.cSubbandCombiner,xsubband,x(:,colidx),fd/Fs);
                    y(:,colidx) = ytemp(1:size(y,1),1:numel(colidx));
                else
                    colidx = (pIdx-1)*ncol_per_path+(1:ncol_per_path);
                    if obj.EnablePolarization
                        % ncol_per_path is 3 (X,Y,Z component)
                        
                        % two-ray model follows oblique incidence in a
                        % plane along the propagation path. So it may not
                        % align with an axis plane yet the incidence angle
                        % needs to be computed in the incident plane and
                        % the parallel/perpendicular components need to
                        % align with the incident plane as well. Z axis is
                        % already aligned, so a rotation is needed for x
                        % and y components. Once rotated, the new x
                        % component is perpendicular component and the new
                        % y and z components are parallel component
                                                                                                
                        % reflection path

                        % in x, y, z
                        xsubband = step(obj.cSubbandDivider,x(:,colidx));
                        for m = 1:size(xsubband,3)
                            xsubband(:,:,m) = wloss(m)*xsubband(:,:,m);
                        end
                        ytemp = wphase*step(obj.cSubbandCombiner,xsubband,x(:,colidx),fd/Fs);
                        y(:,colidx) = ytemp(1:size(y,1),1:numel(colidx));

                        % convert to new coordinates, apply
                        % corresponding coefficients and convert back
                        % to x, y, z
                        rot_mat = rotz(rot_ang(pIdx));
                        y(:,colidx) = (rot_mat*bsxfun(@times,...
                            [refcoef_perp(pIdx);refcoef_para(pIdx)*[1;-1]],...
                            rot_mat'*y(:,colidx).')).';
                    else
                        xsubband = step(obj.cSubbandDivider,x(:,colidx));
                        for m = 1:size(xsubband,3)
                            xsubband(:,:,m) = wloss(m)*xsubband(:,:,m);
                        end
                        ytemp = wphase*step(obj.cSubbandCombiner,xsubband,x(:,colidx),fd/Fs);
                        y(:,colidx) = ytemp(1:size(y,1),1:numel(colidx));
                        
                        y(:,colidx) = refcoef(pIdx)*y(:,colidx);
                    end
                end
            end
        end

        function setRangeFactor(obj)
            obj.pRangeFactor = 1;
        end
        
        function validateNumberOfPositionPairs(obj,x,startLocSize,endLocSize) %#ok<INUSL>
            coder.extrinsic('mat2str');
            coder.extrinsic('num2str');          
            
            ncol_x = size(x,2);
            ncol_pos = max(startLocSize(2),endLocSize(2));
            
            if ncol_x == ncol_pos
                numPosPairs = ncol_x;
            elseif ncol_x == 2*ncol_pos
                numPosPairs = ncol_x/2;
            else
                error(message('phased:phased:invalidColumnDimension','X',ncol_pos,2*ncol_pos));
            end

            if numPosPairs == 1
                expDim = '[3 ';
            else
                expDim = '[3 1] or [3 ';
            end
            cond =  ~(isequal(startLocSize,[3 numPosPairs]) || isequal(startLocSize,[3 1]));
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDimensions','Pos1',...
                        [expDim coder.const(num2str(numPosPairs)) ']'], ...
                               coder.const(mat2str(startLocSize)));
            end


            cond =  ~(isequal(endLocSize,[3 numPosPairs]) || isequal(endLocSize,[3 1]));
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDimensions','Pos2',...
                        [expDim coder.const(num2str(numPosPairs)) ']'], ...
                               coder.const(mat2str(endLocSize)));
            end
            
            % cannot be vectors at the same time if there are more than one
            % signal, taken care of by error checking at line 622

            % at least one postion needs to be vectors
            cond = ~ (isequal(startLocSize,[3 1]) || isequal(endLocSize,[3 1]));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:phased:FreeSpace:AtLeastOneColumnVect','Pos1','Pos2');
            end
            
        end
        
        function num = getSignalBufferWidth(obj)
            % Used in Simulink, nonpolarized
            sz_x = propagatedInputSize(obj,1);
            ncol_x = sz_x(2);
            sz_pos1 = propagatedInputSize(obj,2);
            sz_pos2 = propagatedInputSize(obj,3);
            ncol_pos = max(sz_pos1(2),sz_pos2(2));
            
            if ncol_x == ncol_pos
                num = 2*ncol_x;
            else % ncol_x == 2*ncol_pos
                num = ncol_x;
            end
        end
        
        function x = preProcess(obj,x_in) 
            if obj.pSeparateInputRays
                x = x_in;
            else
                x = reshape([x_in;x_in],[],2*size(x_in,2));
            end
        end
        
        function y = postProcess(obj,y_in) 
            if ~obj.CombinedRaysOutput
                y = y_in;
            else
                nsamp = size(y_in,1);
                if obj.EnablePolarization
                    y = reshape(sum(reshape(y_in,3*nsamp,2,[]),2),nsamp,[]);
                else
                    y = reshape(sum(reshape(y_in,nsamp,2,[]),2),nsamp,[]);
                end
            end
        end
    end
    
    methods (Access = protected)
        function sz_out = getOutputSizeImpl(obj)
            in_dim = propagatedInputSize(obj,1);
            pos1_dim = propagatedInputSize(obj,2);
            pos2_dim = propagatedInputSize(obj,3);
            npos = max(pos1_dim(2),pos2_dim(2));
            if obj.CombinedRaysOutput
                sz_out = [in_dim(1) npos];
            else
                sz_out = [in_dim(1) 2*npos];
            end
        end
    end
    
    methods (Access = protected, Static, Hidden)
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:WidebandTwoRayChannelTitle')),...
              'Text',getString(message('phased:library:block:WidebandTwoRayChannelDesc')));
        end
        function groups = getPropertyGroupsImpl
            groups = matlab.system.display.Section(...
                'phased.WidebandTwoRayChannel');
            groups.PropertyList = groups.PropertyList([12:13 2 6:11 14 15 1 3:5 16:end]);
            dMaximumDistanceSource = matlab.system.display.internal.Property(...
                'MaximumDistanceSource','IsGraphical',false,...
                'UseClassDefault',false,'Default','Property');
            dMaximumNumInputSamplesSource = matlab.system.display.internal.Property(...
                'MaximumNumInputSamplesSource','IsGraphical',false);
            dMaximumNumInputSamples = matlab.system.display.internal.Property(...
                'MaximumNumInputSamples','IsGraphical',false);
            % dSampleRate = matlab.system.display.internal.Property(...
            %     'SampleRate','IsObjectDisplayOnly',true);
            dEnablePolarization = matlab.system.display.internal.Property(...
                'EnablePolarization','IsObjectDisplayOnly',true);
            for m = 1:numel(groups.PropertyList)
                if strcmp(groups.PropertyList{m},'MaximumDistanceSource')
                    groups.PropertyList{m} = dMaximumDistanceSource;
                % elseif strcmp(groups.PropertyList{m},'SampleRate')
                %     groups.PropertyList{m} = dSampleRate;
                elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamplesSource')
                    groups.PropertyList{m} = dMaximumNumInputSamplesSource;
                elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamples')
                    groups.PropertyList{m} = dMaximumNumInputSamples;
                elseif strcmp(groups.PropertyList{m},'EnablePolarization')
                    groups.PropertyList{m} = dEnablePolarization;
                end
            end
        end
    end
    
    methods (Access = protected)
                
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Wideband Two-Ray\nChannel');
        end
    end
    
end