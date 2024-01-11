classdef (Sealed,StrictDefaults) SumDifferenceMonopulseTracker < phased.internal.AbstractSumDifferenceMonopulse & ...
     matlab.system.mixin.Propagates
%SumDifferenceMonopulseTracker    Sum and difference monopulse for ULA
%   H = phased.SumDifferenceMonopulseTracker creates a tracker System
%   object, H, that uses sum and difference monopulse algorithms on a
%   uniform linear array (ULA).
%
%   H = phased.SumDifferenceMonopulseTracker(Name,Value) creates a ULA
%   monopulse tracker object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   ANG = step(H,X,STEER) estimates the incoming direction ANG (in
%   degrees) of the input signal X based on an initial guess on the
%   direction, often the current steering angle, STEER (in degrees).
%
%   The tracker uses a sum and difference monopulse algorithm to estimate
%   the direction and the difference steering vector is obtained by phase
%   reversing the latter half of the sum steering vector. Since monopulse
%   processing only requires one snapshot, X is a row vector whose number
%   of columns corresponds to number of channels.
%
%   Both ANG and STEER are scalars specified as broadside angles (in
%   degrees) defined in the array's local coordinate system. A broadside
%   angle must be between [-90 90].
%
%   For details regarding the local coordinate system of the ULA, type
%   phased.ULA.coordinateSystemInfo.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   SumDifferenceMonopulseTracker methods:
%
%   step     - Perform monopulse tracking using ULA (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a ULA monopulse tracker object with same property 
%              values
%   isLocked - Locked status (logical)
%
%   SumDifferenceMonopulseTracker properties:
%
%   SensorArray         - Sensor array
%   PropagationSpeed    - Signal propagation speed 
%   OperatingFrequency  - Operating frequency
%   NumPhaseShifterBits - Number of bits in phase shifters
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Determine the direction of a target at around 60 degrees broadside
%   %   angle of a ULA
%
%   array = phased.ULA('NumElements',4); 
%   steervector = phased.SteeringVector('SensorArray',array);
%   tracker = phased.SumDifferenceMonopulseTracker('SensorArray',array);
%   x = steervector(tracker.OperatingFrequency,60.1).';
%   est_dir = tracker(x,60)
%
%   See also phased, phased.SumDifferenceMonopulseTracker2D,
%   phased.BeamscanEstimator.

%   Copyright 2008-2018 The MathWorks, Inc.

%   Reference
%   [1] Seliktar, Y. Space-Time Adaptive Monopulse Processing, 
%       PhD dissertation, 1998
%   [2] Rhodes, D. Introduction to Monopulse, Artech House, 1980


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Access = private)
        %privRatioGrid - Private property used to save the monopulse mapping
        %between the ratio and sine offset
        privRatioGrid
        %privRatioGrid - Private property used to save the monopulse mapping
        %between the ratio and sine offset
        privSinOffsetGrid
    end

    methods
        function obj = SumDifferenceMonopulseTracker(varargin)

            obj@phased.internal.AbstractSumDifferenceMonopulse(varargin{:});

        end
    end

    methods (Access = protected)
        function privValidateSensorArray(obj,val)  %#ok<INUSL>
            validateattributes(val,{'phased.ULA','phased.HeterogeneousULA'},{'scalar'},...
                '','SensorArray');
        end
    end  

    methods (Access = protected)
        function initializeSensorArray(obj)
            obj.SensorArray = phased.ULA;
        end
        
        function setupImpl(obj,x,stang)     %#ok<INUSD>
            setupImpl@phased.internal.AbstractSumDifferenceMonopulse(obj);
            %populate the ratio.
            [obj.privSinOffsetGrid, obj.privRatioGrid] = trainingMonopulse(obj);
        end
        
        function validateInputsImpl(obj,x,stang)
            validateInputsImpl@phased.internal.AbstractSumDifferenceMonopulse(obj,x);
            cond = ~isa(stang,'float');
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','Steer','float');
            end
            cond = ~isscalar(stang);
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeScalar','Steer');
            end           
            cond = ~isreal(stang);
            if cond
                coder.internal.errorIf(cond,'phased:SumDifferenceMonopulseTracker:Step:NeedReal','Steer');
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSumDifferenceMonopulse(obj);
            if isLocked(obj)
                s.privRatioGrid = obj.privRatioGrid;
                s.privSinOffsetGrid = obj.privSinOffsetGrid;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function ang = stepImpl(obj,x,stang)

            cond = (stang < -90) || (stang > 90);
            if cond
                coder.internal.errorIf(cond,'phased:SumDifferenceMonopulseTracker:Step:OutOfBoundSteeringAngle','Steer');
            end

            % broadside
            % generate the sum and difference steering vector
            hstv = obj.cSteeringVector;
            stvSum = step(hstv,cast(obj.OperatingFrequency,'double'),cast(stang,'double'));
            stvDiff = phaseReverse(obj,stvSum);
            % apply sum and difference to the data and obtain the angle offset
            ratioErr = obj.applyMonopulse(stvSum,stvDiff,x);
            sinOffset = obj.getSinOffset(...
                ratioErr,obj.privRatioGrid,obj.privSinOffsetGrid);
            % convert the offset in sine domain to an angle
            ang = obj.addSinOffsetToRef(stang,sinOffset);
        end       
    end

    methods (Access = private)
        function stvDiff = phaseReverse(obj,stvSum)
        %phaseReverse Phase reverse the steering vectors of ULA monopulse
        %   stvDiff = phaseReverse(Hmp,stvSum,DIR) returns the difference
        %   steering vector stvDiff by phase reversing the sum steering
        %   vector stvSum of a monopulse object Hmp.
            N = obj.pNumSensorArrayElements;
            % A ULA can only distinguish the angles in azimuth domain.
            % Therefore, second half of the steering vector needs to be
            % phase reversed. If the number of elements is odd, then the
            % middle element is weighted by zero.
            stvDiff = -1i*stvSum;
            if ~rem(N,2)
                stvDiff(N/2+1:end) = -stvDiff(N/2+1:end);
            else
                stvDiff((N+3)/2:end) = -stvDiff((N+3)/2:end);
                stvDiff((N+1)/2) = 0;
            end
        end
        
        function [sinGrid, ratioGrid] = trainingMonopulse(obj)
            % generate the sum and difference steering vector
            hstv1 = phased.SteeringVector(...
                 'SensorArray',obj.SensorArray,...
                 'PropagationSpeed',obj.PropagationSpeed);
            % can be removed once var-D is supported.
            hstv2 = phased.SteeringVector(...
                 'SensorArray',obj.SensorArray,...
                 'PropagationSpeed',obj.PropagationSpeed);

            freq = obj.OperatingFrequency;
            stvSum = step(hstv1,cast(freq,'double'),[0; 0]);
            stvDiff = phaseReverse(obj,stvSum);
            % sample the mainbeam with gridN points in sine domain, which is
            % sinGrid.
            gridN = 1024;
            N = obj.pNumSensorArrayElements;
            lambda = obj.PropagationSpeed/freq;
            sz_grid = obj.SensorArray.ElementSpacing;
            sinGrid = linspace(-1/2,1/2,gridN)*lambda/(N*sz_grid);
            degGrid = asind(sinGrid);
            % simulate incoming signal from those offset in sinGrid
            trainingStvMtx = step(hstv2,cast(freq,'double'),cast([degGrid; zeros(size(degGrid))],'double'));
            trainingStvMtx = trainingStvMtx.';  % channels in columns
            % apply monopulse, calculate the ratio between difference and sum
            % steering vector and save the results in ratioGrid.
            ratioGrid = obj.applyMonopulse(stvSum,stvDiff,trainingStvMtx);
        end
    end
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractSumDifferenceMonopulse('ula');
            props = 'NumPhaseShifterBits';
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:SumDifferenceMonopulseTrackerTitle')),...
                'Text',getString(message('phased:library:block:SumDifferenceMonopulseTrackerDesc')));
        end
    end

    methods (Access=protected)
        function varargout = getOutputNamesImpl(~)
            varargout = {'Ang'};
        end
        
        function varargout = getInputNamesImpl(~)
            varargout = {'X','Steer'};
        end
        
        function varargout = getOutputSizeImpl(~)
            varargout{1} = [1 1];
        end
        function varargout = isOutputFixedSizeImpl(~)
            varargout{1} = true;
        end
        function varargout = getOutputDataTypeImpl(obj)
            varargout{1} = propagatedInputDataType(obj,1);
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = false;
        end        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('ULA\nMonopulse');        
        end
end
    
end

