classdef (Sealed,StrictDefaults) IsoSpeedUnderwaterPaths < matlab.System    
%IsoSpeedUnderwaterPaths  Constant speed underwater acoustic paths
%   H = phased.IsoSpeedUnderwaterPaths creates an underwater acoustic
%   environment System object, H. This object generates multiple
%   propagation paths in an underwater acoustic channel using the method of
%   images. The underwater channel has constant sound speed, a water-to-air
%   interface, and a water-to-ground interface. 
%
%   H = phased.IsoSpeedUnderwaterPaths (Name,Value) returns an underwater
%   acoustic environment object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   PATHS = step(H,POS1,POS2,VEL1,VEL2,T) returns the propagation paths
%   matrix, PATHS, in an underwater acoustic channel from the starting
%   position POS1 to the ending position POS2. VEL1 and VEL2 specify the
%   velocity of the origin and the velocity of the destination,
%   respectively. POS1 and POS2 are 3x1 column vectors in the form of [x;
%   y; z] (in meters). VEL1 and VEL2 are 3x1 column vectors in the form of
%   [Vx; Vy; Vz] (in meters/second). T specifies the elapsed time (in
%   seconds) for the current step.
%
%   PATHS is a 3 by N matrix, where N is the number of paths. The first row
%   of PATHS contains propagation time delays (in seconds), the second
%   contains the total reflection coefficient for each path due to
%   interface reflections, and the third contains the spreading loss for
%   each path in dB.
%
%   [PATHS,DOP,ALOSS,RCVANG,SRCANG] = step(H,POS1,POS2,VEL1,VEL2,T) returns
%   the Doppler factor, DOP, frequency-dependent absorption loss, ALOSS,
%   the receiver arrival angles (in degrees), RCVANG, and the source angles
%   (in degrees), SRCANG.
%
%   DOP is an N element row vector, where N is the number of paths. Each
%   element of DOP contains the factor which multiplies the emitted
%   frequency to produce the observed, doppler-shifted frequency. For
%   one-way propagation, VEL1 corresponds to the source, and VEL2 to the
%   receiver.
%
%   ALOSS is a M by N+1 matrix, where M is the number of elements of the
%   LossFrequencies property and N is the number of paths. The first column
%   of M contains frequencies in Hz and the remaining columns contain
%   absorption loss for each path in dB.
%
%   RCVANG is a 2 by N matrix, where N is the number of paths. Each column
%   contains the azimuth and elevation receiver angles for one path.
%
%   SRCANG is a 2 by N matrix, where N is the number of paths. Each column
%   contains the azimuth and elevation source angles for one path.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   IsoSpeedUnderwaterPaths methods:
%
%   step     - Generate underwater channel paths (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create underwater channel object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the object
%
%   IsoSpeedUnderwaterPaths properties:
%
%   ChannelDepth             - Channel depth 
%   PropagationSpeed         - Propagation speed 
%   NumPathsSource           - Source of number of paths
%   NumPaths                 - Number of paths to return
%   CoherenceTime            - Channel coherence times
%   BottomLoss               - Loss from bottom interface 
%   LossFrequencies          - Frequencies for absorption loss
%   TwoWayPropagation        - Perform two-way propagation
%
%   % Example:
%   %   Compute 10 propagation paths between a stationary source and 
%   %   receiver for a 200 m deep channel. Plot the delay times as a
%   %   function of path index.
%
%   proppaths = phased.IsoSpeedUnderwaterPaths('ChannelDepth',200,...
%     'NumPathsSource','Property','NumPaths',10);
%   [paths,dop,aloss] = proppaths([0;0;-160],[10;0;-50],...
%     [0;0;0],[0;0;0],1);
%   figure;
%   plot(1:10,paths(1,:))
%   xlabel('Path Index')
%   ylabel('Delay Time (s)')
%
%   See also phased.BackscatterSonarTarget, phased.MultipathChannel, 
%   phased.IsotropicHydrophone, phased.IsotropicProjector

%   Copyright 2016 The MathWorks, Inc.

%   Reference
%   [1] Urick, Principles of Underwater Sound, Peninsula Publishing, 1983

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %ChannelDepth       Channel depth (m)
        %   Specify the channel depth as a positive scalar in meters. The
        %   default value of this property is 100.
        ChannelDepth=100
        %PropagationSpeed   Propagation speed (m/s)
        %   Specify the wave propagation speed in meters/second as a
        %   positive scalar. The default value of this property is 1520.
        PropagationSpeed = 1520
        %NumPathsSource Source of the number of paths
        %   Specify how the number of paths is defined as one of 'Auto' |
        %   'Property', where the default is 'Auto'. When you set this
        %   property to 'Auto', the object automatically determines the
        %   number of paths based on spreading and reflection losses. When
        %   you set this property to 'Property', the number of paths is
        %   specified via the NumPaths property.
        NumPathsSource = 'Auto'
        %NumPaths       Number of paths
        %   Specify the number of paths to return as a positive integer
        %   between 1 and 51. The default value is 10.
        NumPaths=10
        %CoherenceTime      Coherence time (s)
        %   Specify the channel coherence time interval as a nonnegative
        %   scalar in seconds. The default value is 0.
        CoherenceTime=0
        %BottomLoss         Loss from channel bottom reflections (dB)
        %   Specify the loss in dB for each bottom reflection as a
        %   nonnegative scalar. The default value is 6.
        BottomLoss=6
        %LossFrequencies    Frequencies for absorption loss (Hz)
        %   Specify a set of frequencies, LossFrequencies, as a vector in
        %   Hz. Absorption loss is computed at these frequencies in the
        %   step output ALOSS. The default value is a vector from 1 to 100
        %   kHz with 1 kHz spacing.
        LossFrequencies = (1:100)*1e3
    end

    properties (Nontunable, Logical)
        %TwoWayPropagation  Perform two-way propagation
        %   Set this property to true to perform two-way propagation. Set
        %   this property to false to perform one-way propagation. The
        %   default value of this property is false.
        TwoWayPropagation = false;        
    end
    
    properties(Access = private,Constant)
        % Maximum number of paths (odd number to include direct path)
        pMaxNumPaths = 51
    end
 
    properties(Access = private, Nontunable)
        % Propagation range factor (one-way only)
        pRangeFactor = 1
    end
    
    properties(Access = private)
        % Current geo-time
        pCurrentTime
        % Channel depth
        pChannelDepth
        % Paths angles at receiver
        prcvAng
        % Path angles at source
        psrcAng
        % Paths
        pPaths
        % Doppler
        pdopShift
        % Logical indicies of paths to return
        pIndPaths     
    end

    properties(Constant, Hidden)
        NumPathsSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    methods
        function obj = IsoSpeedUnderwaterPaths(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
 
    methods
        function set.ChannelDepth(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'positive'},'IsoSpeedUnderwaterPaths','ChannelDepth');
            obj.ChannelDepth = value;
        end
        function set.PropagationSpeed(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'positive'},'IsoSpeedUnderwaterPaths','PropagationSpeed');
            obj.PropagationSpeed = value;
        end
        function set.NumPaths(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'integer','positive','<=',obj.pMaxNumPaths}, ...
                'IsoSpeedUnderwaterPaths','NumPaths');
            obj.NumPaths = value;
        end
        function set.CoherenceTime(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
                'nonnegative'},'IsoSpeedUnderwaterPaths','CoherenceTime');
            obj.CoherenceTime = value;
        end
        function set.BottomLoss(obj,value)
            validateattributes(value,{'double'},{'scalar','finite',...
              'nonnegative'},'IsoSpeedUnderwaterPaths','BottomLoss');
            obj.BottomLoss = value;
        end
        function set.LossFrequencies(obj,value)
            validateattributes(value,{'double'},{'vector','finite',...
              'nonnegative'},'IsoSpeedUnderwaterPaths','LossFrequencies');
            obj.LossFrequencies = value;
        end
    end
    
    methods (Access = protected)
       
        function flag = isInactivePropertyImpl(obj, prop)
            if (obj.NumPathsSource(1) == 'A') && ...
              strcmp(prop,'NumPaths')
                flag = true;
            else
                flag = false;
            end
        end
        
        function validateInputsImpl(obj,startLoc,endLoc,baseVel,targetVel,T) %#ok<INUSL>
            coder.extrinsic('mat2str');
            coder.extrinsic('num2str');            
            
            % Check that positions and velocities are real doubles and
            % non-nan.
            cond =  ~isa(startLoc,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','Pos1','double');
            end
            cond =  ~isreal(startLoc);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:NeedReal', 'Pos1');
            end
            cond = any(any(isnan(startLoc)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Pos1');
            end           
            
            cond =  ~isa(endLoc,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','Pos2','double');
            end
            cond =  ~isreal(endLoc);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:NeedReal', 'Pos2');
            end
            cond = any(any(isnan(endLoc)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Pos2');
            end
            
            cond =   ~isa(baseVel,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','Vel1','double');
            end
            
            cond =   ~isreal(baseVel);
            if cond
                coder.internal.errorIf(cond, ...
                      'phased:step:NeedReal', 'Vel1');
            end
            cond = any(any(isnan(baseVel)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Vel1');
            end           

            cond =   ~isa(targetVel,'double');
            if cond
                coder.internal.errorIf(cond, ...
                      'MATLAB:system:invalidInputDataType','Vel2','double');
            end
            cond =   ~isreal(targetVel);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:NeedReal', 'Vel2');
            end
            cond = any(any(isnan(targetVel)));
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'Vel2');
            end  
            
            % Check time step is a real scalar double and non-nan.
            cond =  ~isscalar(T);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeScalar','T');
            end
            cond =  ~isa(T,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','T','double');
            end
            cond =  ~isreal(T);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:NeedReal', 'T');
            end
            cond = isnan(T);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:step:expectedNonNaN', 'T');
            end
            
            % Check input vector sizes
            startLocSize = size(startLoc);
            cond = ~isequal(startLocSize,[3 1]);
            if cond
                coder.internal.errorIf(cond,...
                    'MATLAB:system:invalidInputDimensions','Pos1',...
                    coder.const(mat2str([3 1])),coder.const(mat2str(startLocSize)));
            end
            endLocSize = size(endLoc);
            cond = ~isequal(endLocSize,[3 1]);
            if cond
                coder.internal.errorIf(cond,...
                    'MATLAB:system:invalidInputDimensions','Pos2',...
                    coder.const(mat2str([3 1])),coder.const(mat2str(endLocSize)));
            end
            baseVelSize = size(baseVel);
            cond = ~isequal(baseVelSize,[3 1]);
            if cond
                coder.internal.errorIf(cond,...
                    'MATLAB:system:invalidInputDimensions','Vel1',...
                    coder.const(mat2str([3 1])),coder.const(mat2str(baseVelSize)));
            end
            targetVelSize = size(targetVel);
            cond = ~isequal(targetVelSize,[3 1]);
            if cond
                coder.internal.errorIf(cond,...
                    'MATLAB:system:invalidInputDimensions','Vel2',...
                    coder.const(mat2str([3 1])),coder.const(mat2str(targetVelSize)));
            end
        end
            
        function setupImpl(obj,startLoc,endLoc,baseVel,targetVel,~)

            obj.pChannelDepth=obj.ChannelDepth;
            obj.pCurrentTime=0;
            obj.pIndPaths = false(1,obj.pMaxNumPaths);
           
            if obj.TwoWayPropagation
                obj.pRangeFactor = 2;
            else
                obj.pRangeFactor = 1;
            end
            
            % Compute initial paths
            [rcvAng,srcAng,pathLength,numRefl,dop,rspeedTx,rspeedTgt] = ...
              computePaths(obj,startLoc,endLoc,baseVel,targetVel);
            
            if obj.NumPathsSource(1) == 'A'
            % Exclude paths that have a spreading plus bottom loss
            % less than 20 dB plus the direct path loss. pRangeRactor is 1
            % for one-way propagation and two for two-way propagation.
            % Units are dB.
                directSpreadLoss = obj.pRangeFactor*20*log10(pathLength(1)/obj.pRangeFactor);
                spreadBottomLoss = obj.pRangeFactor*20*log10(pathLength/obj.pRangeFactor) + ...
                  obj.BottomLoss*numRefl(2,:);
                obj.pIndPaths = spreadBottomLoss <= 20 + directSpreadLoss;
                obj.pPaths = zeros(4,obj.pMaxNumPaths);
            else
                obj.pPaths = zeros(4,obj.NumPaths);
            end
            
            updatePaths(obj,rcvAng,srcAng,pathLength,numRefl,dop,rspeedTx,rspeedTgt);
        end
        
        function [paths, dopShift, aloss, rcvAng, srcAng] = stepImpl(obj,startLoc,endLoc,baseVel,targetVel,T)
            
            validateLocations(obj,startLoc,endLoc,T);
            
            % Update channel information based on coherence time. Update
            % the channel every time current time passes a mulitple of the
            % coherence time. Check that the remainder is less than the
            % period to update only once for each multiple.
            if obj.CoherenceTime == 0 || obj.pCurrentTime == 0 || ...
                    obj.pCurrentTime >= obj.CoherenceTime && mod(obj.pCurrentTime,obj.CoherenceTime) < T      
                % Update paths
                [rcvAng,srcAng,pathLength,numRefl,dop,rspeedTx,rspeedTgt] = computePaths(obj,startLoc,endLoc,baseVel,targetVel);
                updatePaths(obj,rcvAng,srcAng,pathLength,numRefl,dop,rspeedTx,rspeedTgt);
            end
            
            % Update current time
            obj.pCurrentTime = obj.pCurrentTime + T;
            
            % Assign output arguments
            rcvAng = obj.prcvAng;
            srcAng = obj.psrcAng;
            paths = obj.pPaths(2:end,:); % Return delays, reflection coeffs, and spreading loss, but not path length.
            dopShift = obj.pdopShift;
            aloss = [obj.LossFrequencies(:) ...
              bsxfun(@times,thorpeloss(obj.LossFrequencies(:)),obj.pPaths(1,:))]; 
            
            % Assign nans for unused paths when PathNumSource is Auto
            if obj.NumPathsSource(1) == 'A'
                paths(:,~obj.pIndPaths) = nan;
                dopShift(~obj.pIndPaths) = nan;
                aloss(:,[false ~obj.pIndPaths]) = nan;
            end
            
        end
        
        function validateLocations(obj,startLoc,endLoc,T)
            coder.extrinsic('mat2str');
            
            % Check that T is non-negative
            cond = T < 0;
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:expectedNonnegative','T')
            end
            
            if isempty(coder.target)	
                % Check that z position of source and receiver is in the range
                % [-1*obj.ChannelDepth 0]
                cond = startLoc(3) > 0 || startLoc(3) < -1*obj.ChannelDepth;
                if cond
                    coder.internal.errorIf(cond,...
                      'phased:IsoSpeedUnderwaterPaths:DepthOutofRange',coder.const(mat2str([-1*obj.ChannelDepth 0])));
                end

                cond = endLoc(3) > 0 || endLoc(3) < -1*obj.ChannelDepth;
                if cond
                    coder.internal.errorIf(cond,...
                      'phased:IsoSpeedUnderwaterPaths:DepthOutofRange',coder.const(mat2str([-1*obj.ChannelDepth 0])));
                end
            else
                % Check that z position of source and receiver is in the range
                % [-1*obj.ChannelDepth 0]
                cond = startLoc(3) > 0 || startLoc(3) < -1*obj.ChannelDepth;
                if cond
                    coder.internal.errorIf(cond,...
                      'phased:IsoSpeedUnderwaterPaths:DepthOutofRange',coder.const('[-1*obj.ChannelDepth 0]'));
                end

                cond = endLoc(3) > 0 || endLoc(3) < -1*obj.ChannelDepth;
                if cond
                    coder.internal.errorIf(cond,...
                      'phased:IsoSpeedUnderwaterPaths:DepthOutofRange',coder.const('[-1*obj.ChannelDepth 0]'));
                end
            end
        end
    end
    
     methods (Access = protected)
        function resetImpl(obj)
            resetImpl@matlab.System(obj);
            obj.pCurrentTime=0; % This will force the branch in step which populates other private properties.
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.pRangeFactor = obj.pRangeFactor;
                s.pCurrentTime = obj.pCurrentTime;
                s.prcvAng = obj.prcvAng;
                s.psrcAng = obj.psrcAng;
                s.pPaths = obj.pPaths;
                s.pdopShift = obj.pdopShift;
                s.pChannelDepth = obj.pChannelDepth;
                s.pIndPaths = obj.pIndPaths;
                s.pRangeFactor = obj.pRangeFactor;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function num = getNumInputsImpl(obj)   %#ok<MANU>
            num = 5;
        end
        
        function num = getNumOutputsImpl(obj)   %#ok<MANU>
            num = 5;
        end
 
        function flag = isInputComplexityLockedImpl(~,~)
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(~,~)  
            flag = false;
        end
     end
     
     methods (Access = protected)
        
        function [rcvAng,srcAng,pathLength,numRefl,dop,rspeedTx,rspeedTgt] = computePaths(obj,startLoc,endLoc,baseVel,targetVel)
      
                % Compute paths using the method of images.
                if obj.NumPathsSource(1) == 'A'
                    numPaths = obj.pMaxNumPaths;
                else
                    numPaths = obj.NumPaths;
                end
                
                [rcvAng,srcAng,pathLength,numRefl,TxPosImage,TxVelImage] = ...
                  phased.internal.imageMethod(obj.pChannelDepth,startLoc,endLoc,baseVel,numPaths);

                cond =  any(pathLength == 0);
                if cond
                coder.internal.errorIf(cond, ...
                     'phased:phased:FreeSpace:InvalidPropagationDistance');
                end
            
                % Compute doppler shift
                rspeedTx = calcRadialSpeed(endLoc,0*targetVel,TxPosImage,TxVelImage);
                rspeedTgt = calcRadialSpeed(endLoc,targetVel,TxPosImage,0*TxVelImage);
                
                dop = bsxfun(@plus,obj.PropagationSpeed,rspeedTgt)./bsxfun(@minus,obj.PropagationSpeed,rspeedTx);
                
                % Compensate for two-way propagation
                if obj.TwoWayPropagation
                    pathLength = pathLength*obj.pRangeFactor; 
                    numRefl = numRefl*obj.pRangeFactor;
                    dop = dop.*bsxfun(@plus,obj.PropagationSpeed,rspeedTx)./bsxfun(@minus,obj.PropagationSpeed,rspeedTgt);
                end
        end
        
        function updatePaths(obj,rcvAng,srcAng,pathLength,numRefl,dop,rspeedTx,rspeedTgt)
            obj.prcvAng = rcvAng; 
            obj.psrcAng = srcAng;
            obj.pPaths(1,:) = pathLength;
            if ~obj.TwoWayPropagation
                obj.pPaths(2,:) = pathLength./bsxfun(@plus,obj.PropagationSpeed,rspeedTgt);
            else
                % Compute the total delay, which is the one-way delay from
                % transmitter to target plus the one-way delay back. The
                % distance back depends on the relative radial speed of the
                % transmitter and target. Here, pathLength is the two-way
                % path length, so we divide by two for each one-way delay.
                delay1 = pathLength/2./bsxfun(@plus,obj.PropagationSpeed,rspeedTgt);
                obj.pPaths(2,:) = delay1+(pathLength/2-delay1.*(rspeedTx+rspeedTgt))./bsxfun(@plus,obj.PropagationSpeed,rspeedTx);
            end
            % The net reflection coefficient contains top bounces (-1) and
            % bottom bounces, which reduce the amplitude by the bottom
            % loss factor.
            obj.pPaths(3,:) = (-1).^numRefl(1,:).*(sqrt(db2pow(-obj.BottomLoss))).^(numRefl(2,:));
            obj.pPaths(4,:) = obj.pRangeFactor*20*log10(pathLength/obj.pRangeFactor); %spherical spreading
            obj.pdopShift = dop;
        end
        
     end
    methods (Hidden, Static)
        function flag = isAllowedInSystemBlock(obj) %#ok<INUSD>
            flag = false;
        end
    end
end

function y = thorpeloss(f)
%thorpeloss  Frequency-dependent attenuation 
    % Compute frequency-dependent attenuation in dB/m using thorpe's
    % equation (pg. 108 in [1]).
    freqnkHz=f/1e3;         
    y = ((3.3*10^-3)+(.11*freqnkHz.^2)./(1+freqnkHz.^2)+(44*freqnkHz.^2)./(4100+freqnkHz.^2)+(3*10^-4)*freqnkHz.^2);
    y = y/1000; % dB/m
end

function rspeed = calcRadialSpeed(tgtpos,tgtvel,refpos,refvel)
%calcRadialSpeed    Compute radial speed
%   RSPEED = calcRadialSpeed(POS,VEL,REFPOS,REFVEL) compute the relative
%   speed RSPEED (in m/s) for a target at position POS (in meters) with a
%   velocity VEL (in m/s) relative to the reference position REFPOS (in
%   meters) and reference velocity REFVEL (in m/s).

%   This is the same functionality as radialspeed function. However,
%   because we already done the input validation here, we want to skip the
%   validation to improve the performance. In addition, here all position
%   and velocity are always column vectors and the target and reference can
%   never be colocated, so it simplifies the computation too.

tgtdirec = bsxfun(@minus,tgtpos,refpos);
veldirec = bsxfun(@minus,tgtvel,refvel);

%Take the 2-norm of veldirec and tgtdirec
rn = sqrt(sum(tgtdirec.^2));

% negative sign to ensure that incoming relative speed is positive
rspeed = -(sum(veldirec.*tgtdirec)./rn);

end


        

