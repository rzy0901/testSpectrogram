classdef (Sealed,StrictDefaults) LOSChannel < phased.internal.AbstractLOSChannel
%LOSChannel Line of sight propagation channel
%   H = phased.LOSChannel creates a line of sight (LOS) propagation channel
%   System object, H. This object simulates narrowband signal propagation
%   in a line of sight channel by applying not only range-dependent time
%   delay, gain, and phase shift, but also the weather related loss to the
%   input signal.
%
%   H = phased.LOSChannel(Name,Value) returns a line of sight propagation
%   channel object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,POS1,POS2,VEL1,VEL2) returns the resulting signal Y when
%   the narrowband signal X propagates in the channel from either the
%   position of a source POS1 to the positions of one or more destinations
%   POS2 or from the positions of several sources to the position of one
%   destination. VEL1 and VEL2 specify the velocities of the sources and
%   destinations respectively. POS1 and POS2 have a size of either 3x1 and
%   3xN respectively or 3xN and 3x1 where each column is in the form of [x;
%   y; z] (in meters). Similarly VEL1 and VEL2 have a size of either 3x1
%   and 3xN respectively or 3xN and 3x1 where each column is in the form of
%   [Vx; Vy; Vz] (in meters/second). N is the number of signals to
%   propagate.
%
%   X can be either an N column matrix or a struct. If X is a matrix, Y is
%   a matrix with the same dimensions. Each column of X and Y represents
%   the signal at the source and destination, respectively, of a
%   propagation path. The propagation paths are defined in the order of the
%   positions specified in POS1 and POS2. If X is a struct, it must contain
%   either X, Y, and Z fields or H and V fields. The X, Y, and Z fields
%   represent the x, y, and z components of the polarized signal,
%   respectively and the H and V fields represent the horizontal and
%   vertical components of the polarized signal, respectively. In this
%   case, the output Y is also a struct containing the same fields as X.
%   Each field in Y contains the resulting signal of the corresponding
%   field in X.
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
%   LOSChannel methods:
%
%   step     - Propagate signal from one location to another (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create line of sight propagation channel object with same 
%              property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the propagation channel
%
%   LOSChannel properties:
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
%   MaximumDistanceSource        - Source of maximum one-way propagation 
%                                  distance
%   MaximumDistance        - Maximum one-way propagation distance
%   MaximumNumInputSamplesSource - Source of maximum number of samples
%                                  of the input signal
%   MaximumNumInputSamples       - Maximum number of samples in input 
%                                  signal
%
%   % Example:
%   %   Calculate the result of propagating a signal in a line of sight  
%   %   channel from a radar at (1000, 0, 0) to a target at (300, 200,
%   %   50) in a medium fog day. Assume both the radar and the target are 
%   %   stationary. The operating frequency is 10 GHz.
%
%   channel = phased.LOSChannel('SampleRate',8e3,...
%           'SpecifyAtmosphere',true,'LiquidWaterDensity',0.05,...
%           'OperatingFrequency',10e9);
%   y = channel(ones(10,1),[1000; 0; 0],[300; 200; 50],[0;0;0],[0;0;0])
%
%   See also phased, phased.RadarTarget, phased.FreeSpace, 
%   phased.TwoRayChannel, phased.WidebandLOSChannel, fogpl, gaspl, rainpl.

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
    
    properties (Access = private, Nontunable)
        %Fractional delay filter
        cFractionalDelayFilter
    end
    
    methods
        function obj = LOSChannel(varargin)
            obj@phased.internal.AbstractLOSChannel(varargin{:});
        end
    end

    methods (Access = protected)

        function setupImpl(obj,varargin)
            setupImpl@phased.internal.AbstractLOSChannel(obj,varargin{:})
            obj.pLambda = obj.PropagationSpeed/obj.OperatingFrequency;
            obj.cFractionalDelayFilter = dsp.VariableFractionalDelay;
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractLOSChannel(obj);
            release(obj.cFractionalDelayFilter);
        end

        function resetImpl(obj)
            resetImpl@phased.internal.AbstractLOSChannel(obj);
            reset(obj.cFractionalDelayFilter);
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractLOSChannel(obj);
            if isLocked(obj)
                s.cFractionalDelayFilter = saveobj(obj.cFractionalDelayFilter);
            end
        end

        function s = loadSubObjects(obj,s,wasLocked)
            s = loadSubObjects@phased.internal.AbstractLOSChannel(obj,s,wasLocked);
            if wasLocked 
                obj.cFractionalDelayFilter = ...
                    dsp.VariableFractionalDelay.loadobj(s.cFractionalDelayFilter);
                s = rmfield(s,'cFractionalDelayFilter');
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
            % flush state, at most 1 sample
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
                fc = obj.OperatingFrequency;
                % spreading loss
                sploss = k*phased.internal.fspl(propdistance,lambda);
                % gas loss
                fc_gas = min(max(fc,1e9),1000e9);
                gasloss = k*gaspl(propdistance,fc_gas,...
                    obj.Temperature,obj.DryAirPressure,obj.WaterVapourDensity);
                plossfactordb = sploss+gasloss;
                % fog loss
                if obj.LiquidWaterDensity > 0
                    fc_fog = min(max(fc,10e9),1000e9);
                    fogloss = k*fogpl(propdistance,fc_fog,...
                        obj.Temperature,obj.LiquidWaterDensity);
                    plossfactordb = plossfactordb+fogloss;
                end
                % rain loss
                if obj.RainRate > 0
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

                    fc_rain = min(max(fc,1e9),1000e9);
                    rainloss = k*rainpl(propdistance,fc_rain,...
                        obj.RainRate,el,tau);
                    plossfactordb = plossfactordb+rainloss;
                end
            end
            
            plossfactor = sqrt(db2pow(plossfactordb));
            
            for pIdx = 1:numOfPropPaths
                colidx = (pIdx-1)*ncol_per_path+(1:ncol_per_path);
                y(:,colidx) = exp(-1i*2*pi*k*propdistance(pIdx)/lambda)/plossfactor(pIdx)*...
                        bsxfun(@times,x(:,colidx),...
                               exp(1i*2*pi*k*rspeed(pIdx)/lambda*(propdelay(pIdx)+(0:size(x,1)-1)'/Fs)));                    
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
              'Title',getString(message('phased:library:block:LOSChannelTitle')),...
              'Text',getString(message('phased:library:block:LOSChannelDesc')));
        end
        function groups = getPropertyGroupsImpl
            groups = matlab.system.display.Section(...
                'phased.LOSChannel');
            groups.PropertyList = groups.PropertyList([8:9 2:7 1 10:end]);
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
            str = sprintf('LOS\nChannel');
        end
    end
    
end

