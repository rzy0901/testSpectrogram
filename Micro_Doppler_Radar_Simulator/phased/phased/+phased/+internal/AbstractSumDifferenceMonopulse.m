classdef (Hidden) AbstractSumDifferenceMonopulse < phased.internal.AbstractNarrowbandArrayProcessing & ...
                matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.

%ABSTRACTSUMDIFFERENCEMONOPULSE Class definition for
%phased.AbstractSumDifferenceMonopulse object

%   Copyright 2010-2015 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties(Nontunable)
        %NumPhaseShifterBits    Number of bits in phase shifters
        %   Specify the number of bits used in the phase shifter as a
        %   non-negative integer. The default value of this property is 0,
        %   indicating there is no quantization effect in the phase
        %   shifter.
        NumPhaseShifterBits = 0
    end
    
    properties (Access = protected, Nontunable)
        cSteeringVector;
        pNumSensorArrayElements;
    end

    methods (Access = private)
        %phaseReverse Phase reverse steering vector
        %   Each subclass needs to define its way to phase reverse the
        %   steering vector.
        stvDiff = phaseReverse(obj,stvSum);
        %trainingMonopulse Train the monopulse processor
        %   Each subclass needs to define its way to train the monopulse
        %   processor.
        [ratioGrid, SinOffset] = trainingMonopulse(obj)
    end

    methods (Access = protected)
        function obj = AbstractSumDifferenceMonopulse(varargin)
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
        end
    end

    methods (Access = protected)
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function validateInputsImpl(obj, x, varargin)
            coder.extrinsic('mat2str');
            cond = ~isa(x,'float');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','float');
            end
            sz_x = size(x);
            sz_sa = [1 getNumElements(obj.SensorArray)];
            cond = ~isequal(sz_x, sz_sa);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDimensions', 'X', ...
                    coder.const(mat2str(sz_sa)), coder.const(mat2str(sz_x)));
            end
        end
        
        function setupImpl(obj)
            obj.pNumSensorArrayElements = getNumElements(obj.SensorArray);
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'NumPhaseShifterBits',obj.NumPhaseShifterBits);
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)  %#ok<INUSL>
            if index == 1
                flag = false;
            else % index == 2
                flag = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function resetImpl(obj)
            reset(obj.cSteeringVector);
        end
        
        function releaseImpl(obj)
            release(obj.cSteeringVector);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            if isLocked(obj)
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.pNumSensorArrayElements = obj.pNumSensorArrayElements;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractNarrowbandArrayProcessing(obj,s);
            if isfield(s,'isLocked') 
                if s.isLocked
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
                end
                s = rmfield(s,'isLocked');
            end
        end
        
    end

    methods (Access = protected,Static)
        function monopulseRatio = applyMonopulse(stvSum,stvDiff,x)
        %applyMonopulse Apply monopulse
        %   R = applyMonopulse(SumSteeringVec,DiffSteeringVec,X) returns the
        %   ratio R of the results when difference steering vector
        %   DiffSteeringVec and sum steering vector SumSteeringVec are applied
        %   to the signal snapshot X.
        %
        %   This is a static method.
            x = x.';  % input signal: # of columns = # of channels
            monopulseRatio = real((stvDiff'*x)./(stvSum'*x));
        end

        function sinOffset = getSinOffset(monopulseRatio,ratioGrid,sinOffsetGrid)
        %getSinOffset Obtain sine offset based on the monopulse ratio
        %   SINOFFSET = getSinOffset(R,RGRID,SINGRID) returns the offset
        %   SINOFFSET in the sin domain corresponding to the ratio R using the
        %   mapping between RGRID and SINGRID.
        %
        %   This is a static method.
            
            cond = abs(monopulseRatio) > abs(max(ratioGrid));
            if cond
                coder.internal.errorIf(cond,'phased:SumDifferenceMonopulse:InvalidRatio','Steer');
            end
            sinOffset = interp1(ratioGrid,sinOffsetGrid,monopulseRatio);
        end

        function ang = addSinOffsetToRef(refAng,sinOffset)
        %addSinOffsetToRef Add sine offset to the reference angle
        %   ANG = addSinOffsetToRef(REFANG,SINOFFSET) returns the ANG (in
        %   degrees) that corresponding to the reference REFANG (in degrees)
        %   added with the offset SINOFFSET. Note that SINOFFSET is specified
        %   in the sine domain.
        %
        %   This is a static method.
            sinang = cast(sind(refAng),class(sinOffset))*ones(1,size(sinOffset,2))+sinOffset;
            if sinang > 1
                sinang = cast(1,class(sinOffset));%codegen purpose
            elseif sinang < -1
                sinang = cast(-1,class(sinOffset));%codegen purpose
            end
            ang = asind(sinang);
        end       
    end
    
end



