classdef (Hidden) AbstractContinuousPhasePulseWaveform < ...
phased.internal.AbstractPulseWaveform
%This class is for internal use only. It may be removed in the future.

%ABSTRACTCONTINUOUSPHASEPULSEWAVEFORM Define the
%ABSTRACTCONTINUOUSPHASEPULSEWAVEFORM class
% This is an abstract class in support of pulse waveform functionality.
% The waveform is continuous within pulses, i.e., no sub-pulses

%   Copyright 2011 The MathWorks, Inc.     

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %DurationSpecification   Method to specify pulse duration
    %   Specify the method to determine pulse duration as one of 'Pulse
    %   width' | 'Duty cycle', where the default is 'Pulse width'. If you
    %   set this property to 'Pulse width', then the pulse duration is
    %   specified via the PulseWidth property. If you set this property to
    %   'Duty cycle', then the pulse duration is computed using the values
    %   specified in both the PRF property and the DutyCycle property.
    DurationSpecification = 'Pulse width';
    %PulseWidth Pulse width (s)
    %   Specify the length of each pulse (in seconds) as a positive scalar.
    %   The default value of this property is 50e-6. This property is
    %   applicable when you set the DurationSpecification property to
    %   'Pulse width'.
    PulseWidth = 50e-6;
    %DutyCycle  Duty cycle
    %   Specify the duty cycle of the pulse waveform as a positive scalar
    %   between 0 and 1. The default value of this property is 0.5. This
    %   property is applicable when you set the DurationSpecification
    %   property to 'Duty cycle'.
    DutyCycle = 0.5;
end

properties(Constant, Hidden)
    DurationSpecificationSet = matlab.system.StringSet({'Pulse width',...
        'Duty cycle'});
end

methods (Access = protected)

    function obj = AbstractContinuousPhasePulseWaveform(varargin)
        %AbstractContinuousPulseWaveform   Construct the
        %AbstractContinuousPulseWaveform class.
        obj@phased.internal.AbstractPulseWaveform(varargin{:});

    end
    
end

methods
    function set.PulseWidth(obj, value)
        validateattributes(value,{'double'},...
            {'scalar','positive','finite'},...
            '','PulseWidth');
        obj.PulseWidth = value;
    end
    function set.DutyCycle(obj, value)
        validateattributes(value,{'double'},...
            {'scalar','positive','>',0,'<',1},...
            '','PulseWidth');
        obj.DutyCycle = value;
    end
end

methods (Access = protected)
    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractPulseWaveform(obj);        
        cond = any(obj.pPulseWidth.*obj.PRF > 1);
        if cond
            coder.internal.errorIf(cond, ...
                'phased:Waveform:NotLessThanOrEqualTo', 'PulseWidth', sprintf( '%5.2e', 1/max( obj.PRF ) ));
        end
    end
    
    function value = getPulseWidth(obj)
        if strcmp(obj.DurationSpecification,'Pulse width')
            value = obj.PulseWidth*ones(1,numel(obj.PRF));
        else
            value = obj.DutyCycle./obj.PRF;
        end
    end
    
    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractPulseWaveform(obj,prop);
        if strcmp(obj.DurationSpecification,'Pulse width') && ...
                strcmp(prop,'DutyCycle')
            flag = true;
        elseif strcmp(obj.DurationSpecification,'Duty cycle') && ...
                strcmp(prop,'PulseWidth')
            flag = true;
        end
    end
    
end

end

% [EOF]
