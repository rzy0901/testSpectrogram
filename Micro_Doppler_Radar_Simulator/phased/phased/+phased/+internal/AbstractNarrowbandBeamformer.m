classdef (Hidden) AbstractNarrowbandBeamformer < phased.internal.AbstractBeamformer
%This class is for internal use only. It may be removed in the future.

%AbstractNarrowbandBeamformer Abstract class for narrowband beamformers

%   Copyright 2009-2016 The MathWorks, Inc.
%     


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)

        %OperatingFrequency     Operating frequency (Hz)
        %   Specify the operating frequency (in Hz) of the beamformer as a
        %   scalar. The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 3e8
        %NumPhaseShifterBits    Number of bits in phase shifters
        %   Specify the number of bits used in the phase shifter as a
        %   non-negative integer. The default value of this property is 0,
        %   indicating there is no quantization effect in the phase
        %   shifter.
        NumPhaseShifterBits = 0
    end
    
    properties (Nontunable, Logical) 
        %WeightsOutputPort    Enable weights output
        %   Set this property to true to output the weights used in the
        %   beamformer. Set this property to false to not output the
        %   weights. The default value of this property is false.
        WeightsOutputPort = false;
    end
    
    properties(Access = protected, Nontunable)
        cSteeringVector
    end      
    
    methods (Access = protected)
        function obj = AbstractNarrowbandBeamformer(varargin)
            obj@phased.internal.AbstractBeamformer(varargin{:});
        end

    end
    
    methods
        function set.OperatingFrequency(obj,val)
            
            sigdatatypes.validateFrequency(val,'phased.internal',...
                'OperatingFrequency',{'double','single'},{'scalar'});
            obj.OperatingFrequency = val;
        end
    end

    methods (Access = protected)
        
        function setupImpl(obj,x)
            setupImpl@phased.internal.AbstractBeamformer(obj,x);
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'NumPhaseShifterBits',obj.NumPhaseShifterBits);
        end
        
        function num = getNumOutputsImpl(obj)
            num = getNumOutputsImpl@phased.internal.AbstractBeamformer(obj);
            if obj.WeightsOutputPort
                num = 2;
            end
        end
        
        function resetImpl(obj)
            reset(obj.cSteeringVector);        
        end
        
        function releaseImpl(obj)
            release(obj.cSteeringVector);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractBeamformer(obj);
            if isLocked(obj)
                s.cSteeringVector = saveobj(obj.cSteeringVector);
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractBeamformer(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
                end
            end
        end
    end
    
end

