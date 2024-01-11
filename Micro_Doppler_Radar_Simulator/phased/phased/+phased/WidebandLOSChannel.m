classdef (Sealed,StrictDefaults) WidebandLOSChannel < phased.internal.AbstractLOSChannel
%WidebandLOSChannel Wideband free space environment
%   H = phased.WidebandLOSChannel creates a Wideband line of sight (LOS)
%   propagation channel System object, H. This object simulates wideband
%   signal propagation in a line of sight channel by applying not only
%   range-dependent time delay, gain, and phase shift, but also the weather
%   related loss to the input signal. The signal is divided into subbands
%   and in each subband the signal is propagated by applying
%   range-dependent time delay, gain and phase shift to the input signal.
%   The resulting signal from all subbands are combined to form the output
%   signal.
%
%   H = phased.WidebandLOSChannel(Name,Value) returns a wideband line of
%   sight propagation environment object, H, with the specified property
%   Name set to the specified Value. You can specify additional name-value
%   pair arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,POS1,POS2,VEL1,VEL2) returns the resulting signal Y when
%   the wideband signal X propagates in free space from either the
%   position of a source POS1 to the positions of one or more destinations
%   POS2 or from the positions of several sources to the position of one
%   destination. VEL1 and VEL2 specify the velocities of the sources and
%   destinations respectively. POS1 and POS2 have a size of either 3x1 and
%   3xN respectively or 3xN and 3x1 where each column is in the form of [x;
%   y; z] (in meters). Similarly VEL1 and VEL2 have a size of either 3x1
%   and 3xN respectively or 3xN and 3x1 where each column is in the
%   form of [Vx; Vy; Vz] (in meters/second). N is the number of signals to
%   propagate.
%
%   X can be either an N column matrix or a struct. If X is a matrix, Y is
%   a matrix with same dimensions. Each column of X and Y represents the
%   signal at the source and destination respectively of a propagation
%   path. The propagation paths are defined in the order of the positions
%   specified in POS1 and POS2. If X is a struct, it must contain either X,
%   Y, and Z fields or H and V fields. The X, Y, and Z fields represent the
%   x, y, and z components of the polarized signal, respectively and the H
%   and V fields represent the horizontal and vertical components of the
%   polarized signal, respectively. In this case, the output Y is also a
%   struct containing the same fields as X. Each field in Y contains the
%   resulting signal of the corresponding field in X.
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
%   as well as the loss due to atmosphere and weather.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   WidebandLOSChannel methods:
%
%   step     - Propagate signal from one location to another (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create wideband line of sight propagation channel object 
%              with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the propagation channel
%
%   WidebandLOSChannel properties:
%
%   PropagationSpeed       - Propagation speed 
%   OperatingFrequency     - Signal carrier frequency 
%   SpecifyAtmosphere      - Specify atmosphere parameters
%   Temperature            - Temperature 
%   DryAirPressure         - Dry air pressure 
%   WaterVapourDensity     - Water vapour density 
%   LiquidWaterDensity     - Liquid water density 
%   RainRate               - Rain rate 
%   TwoWayPropagation      - Perform two-way propagation
%   SampleRate             - Sample rate 
%   NumSubbands            - Number of subbands
%   MaximumDistanceSource        - Source of maximum one-way propagation 
%                                  distance
%   MaximumDistance        - Maximum one-way propagation distance
%   MaximumNumInputSamplesSource - Source of maximum number of samples
%                                  of the input signal
%   MaximumNumInputSamples       - Maximum number of samples in input 
%                                  signal
%
%   % Example:
%   %   Propagate a wideband signal with 3 tones in there. The center
%   %   frequency is 10 GHz and the three tones are 9.75 GHz, 10 GHz, and
%   %   10.25 GHz, respectively. Plot the spectrum of original signal and
%   %   the propagated signal to observe the Doppler effect. Assume it is
%   %   raining at a rate of 95 mm/h. Also consider the dry air gas loss.
%
%   c = 3e8; Fc = 10e9; Fs = 2e9; F = [-0.25e9 0 0.25e9];
%   rpos = [0;0;0];        rvel = [0;0;0];     % stationary radar
%   tpos = [30/Fs*c;0;0];  tvel = [45;0;0];    % moving target
%   dop = -tvel(1)./(c./(F+Fc));               % expected Doppler
%   t = (0:199)/Fs;
%   x = sum(exp(1i*2*pi*t.'*F),2);
%
%   wbchan = phased.WidebandLOSChannel(...
%       'PropagationSpeed',c,'OperatingFrequency',Fc,...
%       'SpecifyAtmosphere',true,'RainRate',95,...  
%       'SampleRate',Fs);
%   y = wbchan(x,rpos,tpos,rvel,tvel);
%   periodogram([x y],rectwin(size(x,1)),1024,Fs,'centered');
%   ylim([-200 0]); legend('original','propagated');
%
%   See also phased, phased.LOSChannel, phased.RadarTarget,
%   phased.TwoRayChannel, phased.WidebandFreeSpace, fogpl, gaspl, rainpl.

%   Copyright 2015-2016 The MathWorks, Inc.

%   Reference
%   [1] John Proakis, Digital Communications, 4th Ed., McGraw-Hill, 2001
%   [2] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%       McGraw-Hill, 2001 
%   [3] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [4] John Seybold, Introduction to RF Propagation, Wiley, 2005
%   [5] Recommendation ITU-R P.838-3, 2005
%   [6] Recommendation ITU-R P.840-3, 2013
%   [7] Recommendation ITU-R P.676-10, 2013


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable, Logical) 
        %TwoWayPropagation Perform two-way propagation
        %   Set this property to true to perform two-way propagation. Set
        %   this property to false to perform one-way propagation. The
        %   default value of this property is false.
        TwoWayPropagation = false
    end
    
    properties (Nontunable, PositiveInteger)
        %NumSubbands    Number of subbands
        %   Specify the number of subbands used in the subband processing
        %   as a positive integer. The default value of this property is
        %   64.
        NumSubbands = 64
    end
       
    properties (Access = private, Nontunable)
        cFractionalDelayFilter
        %subband divider
        cSubbandDivider
        %subband combiner
        cSubbandCombiner
        %subband frequencies
        pSubbandFreq
    end
    
    methods
        function obj = WidebandLOSChannel(varargin)
            obj@phased.internal.AbstractLOSChannel(varargin{:});
        end
    end

    methods (Access = protected)

        function setupImpl(obj,x,varargin)
            setupImpl@phased.internal.AbstractLOSChannel(obj,x,varargin{:});
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
        end

        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractLOSChannel(obj);
            release(obj.cSubbandDivider);
            release(obj.cSubbandCombiner);
            release(obj.cFractionalDelayFilter);
        end

        function resetImpl(obj)
            resetImpl@phased.internal.AbstractLOSChannel(obj);
            reset(obj.cSubbandDivider);
            reset(obj.cSubbandCombiner);
            reset(obj.cFractionalDelayFilter);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractLOSChannel(obj);
            if isLocked(obj)
                s.cFractionalDelayFilter = saveobj(obj.cFractionalDelayFilter);
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
            
            fc = obj.OperatingFrequency;            
            for pIdx = 1:numOfPropPaths
                plossfactor = sqrt(db2pow(plossfactordb(pIdx,:)));
                fd = k*rspeed(pIdx)./lambda(:).';
                
                wloss = 1./plossfactor;
                wphase = exp(-1i*2*pi*fc*propdelay(pIdx));
                
                colidx = (pIdx-1)*ncol_per_path+(1:ncol_per_path);
                xsubband = step(obj.cSubbandDivider,x(:,colidx));
                for m = 1:size(xsubband,3)
                    xsubband(:,:,m) = wloss(m)*xsubband(:,:,m);
                end
                ytemp = wphase*step(obj.cSubbandCombiner,xsubband,x(:,colidx),fd/Fs);
                y(:,colidx) = ytemp(1:size(y,1),1:numel(colidx));
            end

        end

        function setRangeFactor(obj)
            if obj.TwoWayPropagation
                obj.pRangeFactor = 2;
            else
                obj.pRangeFactor = 1;
            end
        end
        
        function validateNumberOfPositionPairs(obj,x,startLocSize,endLocSize) %#ok<INUSL>
            coder.extrinsic('mat2str');
            coder.extrinsic('num2str');            

            numPosPairs = size(x,2);
            
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
            
            cond =   ~(isequal(startLocSize,[3 numPosPairs]) || ...
                      isequal(endLocSize,[3 numPosPairs]));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:phased:FreeSpace:AtLeastOneNotColumnVect','Pos1','Pos2',numPosPairs);
            end
            
            cond = ~ (isequal(startLocSize,[3 1]) || isequal(endLocSize,[3 1]));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:phased:FreeSpace:AtLeastOneColumnVect','Pos1','Pos2');
            end
            
        end
    end
    
    methods (Access = protected, Static, Hidden)
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:WidebandLOSChannelTitle')),...
              'Text',getString(message('phased:library:block:WidebandLOSChannelDesc')));
        end
        function groups = getPropertyGroupsImpl
            groups = matlab.system.display.Section(...
                'phased.WidebandLOSChannel');
            
            % reorder for NumSubbands
            proplist = groups.PropertyList;
            proplist = proplist([9:10 2:8 1 11:end]); 
            groups.PropertyList = proplist;
                        
            dMaximumDistanceSource = matlab.system.display.internal.Property(...
                'MaximumDistanceSource','IsGraphical',false,...
                'UseClassDefault',false,'Default','Property');
            dMaximumNumInputSamplesSource = matlab.system.display.internal.Property(...
                'MaximumNumInputSamplesSource','IsGraphical',false);
            dMaximumNumInputSamples = matlab.system.display.internal.Property(...
                'MaximumNumInputSamples','IsGraphical',false);
            % dSampleRate = matlab.system.display.internal.Property(...
            %     'SampleRate','IsObjectDisplayOnly',true);
            for m = 1:numel(groups.PropertyList)
                if strcmp(groups.PropertyList{m},'MaximumDistanceSource')
                    groups.PropertyList{m} = dMaximumDistanceSource;
                % elseif strcmp(groups.PropertyList{m},'SampleRate')
                %     groups.PropertyList{m} = dSampleRate;
                elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamplesSource')
                    groups.PropertyList{m} = dMaximumNumInputSamplesSource;
                elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamples')
                    groups.PropertyList{m} = dMaximumNumInputSamples;
                end
            end
        end
    end
    
    methods (Access = protected)
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Wideband LOS\nChannel');
        end
    end
    
end

