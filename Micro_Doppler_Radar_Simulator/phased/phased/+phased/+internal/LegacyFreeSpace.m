classdef (Sealed,Hidden, StrictDefaults) LegacyFreeSpace < phased.internal.AbstractLegacyFreeSpace
%This class is for internal use only. It may be removed in the future.

%LegacyFreeSpace preserves the FreeSpace behavior before 15b, where the
%delay is rounded to next sample.

%FreeSpace Free space environment
%   H = phased.internal.LegacyFreeSpace creates a free space environment
%   System object, H. This object simulates narrowband signal propagation
%   in free space, by applying range-dependent time delay, gain and phase
%   shift to the input signal.
%
%   H = phased.internal.LegacyFreeSpace(Name,Value) returns a free space
%   environment object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,POS1,POS2,VEL1,VEL2) returns the resulting signal Y when
%   the narrowband signal X propagates in free space from either the
%   position of a source POS1 to the positions of one or more destinations
%   POS2 or from the positions of several sources to the position of one
%   destination. VEL1 and VEL2 specify the velocities of the sources and
%   destinations respectively. POS1 and POS2 have a size of either 3x1 and
%   3xN respectively or 3xN and 3x1 where each column is in the form of [x;
%   y; z] (in meters). Similarly VEL1 and VEL2 have a size of either 3x1
%   and 3xN respectively or 3xN and 3x1 where each column is in the
%   form of [Vx; Vy; Vz] (in meters/second). N is the number of signals to
%   propagate and can only be equal to 1 when the signal is polarized.
%
%   X can be either an N column matrix or a struct. If X is a matrix, Y is
%   a matrix with same dimensions and each column of X and Y represent the
%   signal at the source and destination respectively of a propagation
%   path. The propagation paths are defined in the order of the positions
%   specified in POS1 and POS2. If X is a struct, it must contain either X,
%   Y, and Z fields or H and V fields. The X, Y, and Z fields represent the
%   x, y, and z components of the polarized signal, respectively and the H
%   and V fields represent the horizontal and vertical components of the
%   polarized signal, respectively. In this case, the output Y is also a
%   struct containing the same fields as X. Each field in Y contains the
%   resulting signal of the corresponding field in X,
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
%   propagation speed. The free space path loss is given by
%
%   L = (4*pi*R/lambda)^2
%
%   where lambda is the signal wavelength.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   FreeSpace methods:
%
%   step     - Propagate signal from one location to another (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create free space object with same property values
%   isLocked - Locked status (logical)
%   reset    - Reset internal states of the propagation channel
%
%   FreeSpace properties:
%
%   PropagationSpeed      - Propagation speed 
%   OperatingFrequency    - Signal carrier frequency 
%   TwoWayPropagation     - Perform two-way propagation
%   SampleRate            - Sample rate 
%   MaximumDistanceSource - Source of maximum one-way propagation distance
%   MaximumDistance       - Maximum one-way propagation distance
%
%   % Example:
%   %   Calculate the result of propagating a signal in a free space 
%   %   environment from a radar at (1000, 0, 0) to a target at (300, 200,
%   %   50). Assume both the radar and the target are stationary.
%
%   henv = phased.FreeSpace('SampleRate',8e3);
%   y = step(henv,ones(10,1),[1000; 0; 0],[300; 200; 50],[0;0;0],[0;0;0])
%
%   See also phased, phased.RadarTarget, fspl.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] John Proakis, Digital Communications, 4th Ed., McGraw-Hill, 2001
%   [2] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%       McGraw-Hill, 2001 
%   [3] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005


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
    
    methods
        function obj = LegacyFreeSpace(varargin)
            obj@phased.internal.AbstractLegacyFreeSpace(varargin{:});
        end
    end

    methods (Access = protected)

        function validateInputSignal(obj,x) %#ok<INUSL>
            if isstruct(x)
                flag_hasXYZ = isfield(x(1),'X') && isfield(x(1),'Y') && isfield(x(1),'Z');
                flag_hasHV = isfield(x(1),'H') && isfield(x(1),'V');
                cond = ~flag_hasXYZ && ~flag_hasHV;
                if cond
                    coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationStruct');
                end
                cond = ~isscalar(x);
                if cond
                    coder.internal.errorIf(cond,'phased:polarization:invalidPolarizationArrayStruct','X');
                end
                if flag_hasXYZ
                    x_x = x.X;
                    x_y = x.Y;
                    x_z = x.Z;
                    cond =  ~isa(x_x,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','X.X','double');
                    end
                    cond =  ~iscolumn(x_x) || isempty(x_x);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector','X.X');
                    end
                    cond =  ~isa(x_y,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','X.Y','double');
                    end
                    cond =  ~iscolumn(x_y) || isempty(x_y);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector','X.Y');
                    end
                    cond =  ~isa(x_z,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','X.Z','double');
                    end
                    cond =  ~iscolumn(x_z) || isempty(x_z);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector','X.Z');
                    end
                    cond = numel(x_x)~=numel(x_y) || numel(x_x)~=numel(x_z);
                    if cond
                        coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                            'X,Y,Z','X');
                    end
                end
                if flag_hasHV
                    x_h = x.H;
                    x_v = x.V;
                    cond =  ~isa(x_h,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','X.H','double');
                    end
                    cond =  ~iscolumn(x_h) || isempty(x_h);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector','X.H');
                    end
                    cond =  ~isa(x_v,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:invalidInputDataType','X.V','double');
                    end
                    cond =  ~iscolumn(x_v) || isempty(x_v);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector','X.V');
                    end
                    cond = numel(x_h)~=numel(x_v) ;
                    if cond
                        coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                            'H,V','X');
                    end
                end
                if flag_hasXYZ && flag_hasHV
                    cond = numel(x_x)~=numel(x_h) ;
                    if cond
                        coder.internal.errorIf(cond,'phased:polarization:polarizationStructDimensionMismatch',...
                            'X,Y,Z,H,V','X');
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
        
        function setupImpl(obj,varargin)
            setupImpl@phased.internal.AbstractLegacyFreeSpace(obj,varargin{:});
            obj.pLambda = obj.PropagationSpeed/obj.OperatingFrequency;
        end
        
        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function [xbuf_in,ndelay] = computeDelayedSignal(obj,x,delay) %#ok<INUSL>
            xbuf_in = x;
            ndelay = ceil(delay);
        end
        
        function [y,propdelay] = computeMultiplePropagatedSignal(obj,x,...
                ncol_per_path,numOfPropPaths,startLoc,endLoc,baseVel,targetVel,Fs) 
            y = complex(zeros(size(x)));
            k = obj.pRangeFactor;
            lambda = obj.pLambda;
            [propdelay,propdistance,rspeed] = computePropagationDelayVelocity(...
                obj,startLoc,endLoc,baseVel,targetVel);
            sploss = k*phased.internal.fspl(propdistance,lambda); % spreading loss
            plossfactor = sqrt(db2pow(sploss));
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
              'Title',getString(message('phased:library:block:FreeSpaceTitle')),...
              'Text',getString(message('phased:library:block:FreeSpaceDesc')));
        end
        function groups = getPropertyGroupsImpl
            groups = matlab.system.display.Section(...
                'phased.internal.LegacyFreeSpace');
            groups.PropertyList = groups.PropertyList([2:3 1 4:end]);
            dMaximumDistanceSource = matlab.system.display.internal.Property(...
                'MaximumDistanceSource','IsGraphical',false,...
                'UseClassDefault',false,'Default','Property');
            dSampleRate = matlab.system.display.internal.Property(...
                'SampleRate','IsObjectDisplayOnly',true);
            for m = 1:numel(groups.PropertyList)
                if strcmp(groups.PropertyList{m},'MaximumDistanceSource')
                    groups.PropertyList{m} = dMaximumDistanceSource;
                elseif strcmp(groups.PropertyList{m},'SampleRate')
                    groups.PropertyList{m} = dSampleRate;
                end
            end
        end
    end
    
    methods (Access = protected)
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Free Space\nChannel');
        end
    end
    
end

