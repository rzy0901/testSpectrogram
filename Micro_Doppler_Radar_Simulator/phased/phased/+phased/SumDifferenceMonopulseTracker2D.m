classdef (Sealed,StrictDefaults) SumDifferenceMonopulseTracker2D < phased.internal.AbstractSumDifferenceMonopulse  & ...
     matlab.system.mixin.Propagates
%SumDifferenceMonopulseTracker2D    Sum and difference monopulse for URA
%   H = phased.SumDifferenceMonopulseTracker2D creates a tracker System
%   object, H, that uses sum and difference monopulse algorithms on a
%   uniform rectangular array (URA).
%
%   H = phased.SumDifferenceMonopulseTracker2D(Name,Value) creates a URA
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
%   Both ANG and STEER are 2x1 vectors in the format of [AzimuthAngle;
%   ElevationAngle] (in degrees). Azimuth angle must be between [-180 180].
%   Elevation angle must be between [-90 90]. All angles are measured in
%   the local coordinate system of the array.
%
%   For details regarding the local coordinate system of the URA, type
%   phased.URA.coordinateSystemInfo.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   SumDifferenceMonopulseTracker2D methods:
%
%   step     - Perform monopulse tracking using URA (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a URA monopulse tracker object with same property 
%              values
%   isLocked - Locked status (logical)
%
%   SumDifferenceMonopulseTracker2D properties:
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
%   %   Determine the direction of a target at around 60 degrees azimuth
%   %   and 20 degrees elevation of a URA.
%
%   array = phased.URA('Size',4); 
%   steervector = phased.SteeringVector('SensorArray',array);
%   tracker = phased.SumDifferenceMonopulseTracker2D('SensorArray',array);
%   x = steervector(tracker.OperatingFrequency,[60.1; 19.5]).';
%   est_dir = tracker(x,[60; 20])
%
%   See also phased, phased.SumDifferenceMonopulseTracker,
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
        %privAzRatioGrid - Private property which holds the info for monopulse
        %mapping in Azimuth direction
        privAzRatioGrid
        %privAzSinOffsetGrid - Private property which holds the info for
        %monopulse mapping in Azimuth direction
        privAzSinOffsetGrid
        %privElRatioGrid - Private property which holds the info for monopulse
        %mapping in Elevation direction
        privElRatioGrid
        %privElSinOffsetGrid - Private property which holds the info for
        %monopulse mapping in Elevation direction
        privElSinOffsetGrid
    end
    properties (Access = private, Nontunable)
        %pNumElementsInRow - Number of elements in each row
        pNumElementsInRow
        %pNumElementsInColumn - Number of elements in each column
        pNumElementsInColumn
    end

    methods
        function obj = SumDifferenceMonopulseTracker2D(varargin)

            obj@phased.internal.AbstractSumDifferenceMonopulse(varargin{:});

        end
    end

    methods (Access = protected)
        function privValidateSensorArray(obj,val)  %#ok<INUSL>
            validateattributes(val,{'phased.URA','phased.HeterogeneousURA'},{'scalar'},...
                '','SensorArray');
        end
    end  

    methods (Access = protected)
        function initializeSensorArray(obj)
            obj.SensorArray = phased.URA;
        end
        
        function setupImpl(obj,x,stang)                                            %#ok<INUSD>
            setupImpl@phased.internal.AbstractSumDifferenceMonopulse(obj);
            obj.pNumElementsInRow = obj.SensorArray.Size(2);
            obj.pNumElementsInColumn = obj.SensorArray.Size(1);
            %populate the ratio.
            [obj.privAzSinOffsetGrid, obj.privAzRatioGrid] = trainingMonopulse(obj,'Az');
            [obj.privElSinOffsetGrid, obj.privElRatioGrid] = trainingMonopulse(obj,'El');
        end

        function validateInputsImpl(obj,x,stang)
            validateInputsImpl@phased.internal.AbstractSumDifferenceMonopulse(obj,x);
            sz_stang = size(stang);
            cond = ~isa(stang,'float');
            if cond          
                coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','Steer','float');
            end
            cond = ~iscolumn(stang) || isempty(stang);
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeColVector','Steer');
            end
            cond = sz_stang(1) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:system:SumDifferenceMonopulse2D:NeedTwoRows','Steer');
            end
            cond = ~isreal(stang);
            if cond
                coder.internal.errorIf(cond,'phased:system:SumDifferenceMonopulse2D:NeedReal','Steer');
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSumDifferenceMonopulse(obj);
            if isLocked(obj)
                s.privAzRatioGrid = obj.privAzRatioGrid;
                s.privAzSinOffsetGrid = obj.privAzSinOffsetGrid;
                s.privElRatioGrid = obj.privElRatioGrid;
                s.privElSinOffsetGrid = obj.privElSinOffsetGrid;
                s.pNumElementsInRow = obj.pNumElementsInRow;
                s.pNumElementsInColumn = obj.pNumElementsInColumn;
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function ang = stepImpl(obj,x,stangArg)

            if isscalar(stangArg)
                stang = [stangArg; 0];
            else
                stang = stangArg;
            end
            cond = stang(1) < -180 || stang(1) > 180;
            if cond
                coder.internal.errorIf(cond,'phased:SumDifferenceMonopulse2D:step:OutOfBoundAzimuth','Steer');
            end
            cond = stang(2) < -90 || stang(2) > 90;
            if cond
                coder.internal.errorIf(cond,'phased:SumDifferenceMonopulse2D:step:OutOfBoundElevation','Steer');
            end

            % process elevation direction
            hstv = obj.cSteeringVector;
            freq = obj.OperatingFrequency;
            stvElSum = step(hstv,cast(freq,'double'),cast(stang,'double'));
            stvElDiff = phaseReverse(obj,stvElSum,'El');
            ratioErrEl = obj.applyMonopulse(stvElSum,stvElDiff,x);
            sinElOffset = obj.getSinOffset(...
                ratioErrEl,obj.privElRatioGrid,obj.privElSinOffsetGrid);     
            ang = cast(stang,class(x));
            ang(2) = obj.addSinOffsetToRef(ang(2),sinElOffset);
            
            % process azimuth direction
            % use updated elevation angle
            stvAzSum = step(hstv,cast(freq,'double'),cast(ang,'double'));
            stvAzDiff = phaseReverse(obj,stvAzSum,'Az');
            ratioErrAz = obj.applyMonopulse(stvAzSum,stvAzDiff,x);
            sinAzOffset = obj.getSinOffset(...
                ratioErrAz,obj.privAzRatioGrid,obj.privAzSinOffsetGrid);     
            sinAzOffset = sinAzOffset/cosd(ang(2));  % adjust cos elevation
            ang(1) = obj.addSinOffsetToRef(ang(1),sinAzOffset);
        end       
    end

    methods (Access = private)
        function stvDiff = phaseReverse(obj,stvSum,dir)
        %phaseReverse Phase reverse the steering vectors of URA monopulse
        %   stvDiff = phaseReverse(H,stvSum,DIR) returns the difference
        %   steering vector stvDiff by phase reversing the sum steering
        %   vector stvSum of a monopulse object H. The phase reverse can
        %   occur in either azimuth direction or elevation direction
        %   depending on the value of DIR. DIR can be one of the following:
        %   [ {'Az'} | 'El].
            N1 = obj.pNumElementsInRow;
            N2 = obj.pNumElementsInColumn;
            stvDiff = -1i*stvSum; % Common for both directions
            if (dir(1) == 'A') %Az
                % Azimuth. Steering vectors for URA are concatenated in
                % columns, therefore, for Azimuth, reverse the second half
                % of the elements. If the number of columns is odd, then
                % the middle column is weighted by zero.
                if ~rem(N1,2)
                    stvDiff(N1/2*N2+1:end) = -stvDiff(N1/2*N2+1:end);
                else
                    stvDiff((N1+1)/2*N2+1:end) = -stvDiff((N1+1)/2*N2+1:end);
                    stvDiff((N1-1)/2*N2+1:(N1+1)/2*N2) = 0;
                end
            else
                % Elevation. Now the lower half needs to be reversed. If
                % the number of columns is odd, then the middle row is
                % weighted by zero.
                for m = 1:N1
                    if ~rem(N2,2)
                        stvDiff((m-1/2)*N2+1:m*N2) = -stvDiff((m-1/2)*N2+1:m*N2);
                    else
                        stvDiff((m*N2-(N2+1)/2+1:m*N2)) = -stvDiff((m*N2-(N2+1)/2+1:m*N2));
                        stvDiff(m*N2-(N2+1)/2) = 0;
                    end
                end
            end
        end
        function [sinGrid, ratioGrid] = trainingMonopulse(obj,dir)
        %trainingAzMonopulse Training of monopulse in azimuth
        %   [SINGRID, RGRID] = trainingAzMonopulse(Hmp) returns the mapping
        %   between the sine offset SINGRID and the corresponding ratios RGRID
        %   of the monopulse object Hmp.
            % generate the sum and difference steering vector
            hstv1 = phased.SteeringVector(...
                 'SensorArray',obj.SensorArray,...
                 'PropagationSpeed',obj.PropagationSpeed);
            % can be removed once var-D is supported.
            hstv2 = phased.SteeringVector(...
                 'SensorArray',obj.SensorArray,...
                 'PropagationSpeed',obj.PropagationSpeed);
            freq = obj.OperatingFrequency;
            % sample the mainbeam with gridN points in sine domain, which is
            % sinGrid.
            gridN = 1024;
            N = obj.pNumElementsInColumn;
            lambda = obj.PropagationSpeed/freq;
            sz_grid = obj.SensorArray.ElementSpacing;
            % generate the sum and difference steering vector
            stvSum = step(hstv1,cast(freq,'double'),[0; 0]);
            if (dir(1) == 'A') %Az
                sinGrid = linspace(-1/2,1/2,gridN)*lambda/(N*sz_grid(2));
                degGrid = asind(sinGrid);
                stvDiff = phaseReverse(obj,stvSum,'Az');
                gridAngle = [degGrid; zeros(size(degGrid))];
            else
                sinGrid = linspace(-1/2,1/2,gridN)*lambda/(N*sz_grid(1));
                degGrid = asind(sinGrid);
                stvDiff = phaseReverse(obj,stvSum,'El');
                gridAngle = [zeros(size(degGrid)); degGrid];
            end
            % simulate incoming signal from those offset in sinGrid
            trainingStvMtx = step(hstv2,cast(freq,'double'),cast(gridAngle,'double'));
            trainingStvMtx = trainingStvMtx.';  % channels in columns
            % apply monopulse, calculate the ratio between difference and sum
            % steering vector and save the results in ratioGrid.
            ratioGrid = obj.applyMonopulse(stvSum,stvDiff,trainingStvMtx);
        end
    end
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractSumDifferenceMonopulse('ura');
            props = 'NumPhaseShifterBits';
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:SumDifferenceMonopulseTracker2DTitle')),...
                'Text',getString(message('phased:library:block:SumDifferenceMonopulseTracker2DDesc')));
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
            varargout{1} = [2 1];
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
            str = sprintf('URA\nMonopulse');        
        end
end

end

