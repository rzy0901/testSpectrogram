classdef (Sealed,StrictDefaults) SimulinkPhaseCodedWaveform < phased.internal.AbstractPhaseCodedWaveform & ...
        matlab.system.mixin.CustomIcon 
%This class is for internal use only. It may be removed in the future.

%PhaseCodedWaveform   Phase-coded pulse waveform
%   H = phased.PhaseCodedWaveform creates a phase-coded pulse waveform
%   System object, H. This object generates samples of a phase-coded pulse
%   waveform.
%
%   H = phased.PhaseCodedWaveform(Name,Value) creates a phase-coded pulse
%   waveform object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H) returns samples of the phase-coded pulse in a column vector
%   Y. Y can contain either a certain number of pulses or a certain number
%   of samples.
%
%   [Y,PRF] = step(H) returns an additional output PRF, as the pulse
%   repetition frequency when you set the PRFOutputPort property to true.
%   This property can be used only when the OutputFormat property is
%   'Pulses'. When you set the PRFOutputPort property to true it returns
%   the current PRF used by the system.
%
%   Y = step(H,PRFIDX) specifies the index of pulse repetition frequency
%   (PRF), PRFIDX, as a positive integer. The index is used to identify the
%   entries specified in the PRF property. This syntax applies when you set
%   the PRFSelectionInputPort property to true. Use this syntax for the
%   cases where the transmit pulse needs to be dynamically selected. Under
%   such situations, the PRF property includes a list of predetermined
%   choices of PRFs. During the simulation, based on PRF index input, one
%   of the PRFs is selected as the PRF for the next transmission.
%
%   Note that the transmission always finishes the current pulse before
%   starting the next pulse. Therefore, when you set the OutputFormat
%   property to 'Samples' and then specify the NumSamples property to be
%   shorter than a pulse, it is possible that during a given simulation
%   step, if the entire output is needed to finish the previously
%   transmitted pulse, the specified PRFIDX is ignored.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H) and y = H() are
%   equivalent.
%
%   PhaseCodedWaveform methods:
%   
%   step             - Return samples of the phase-coded pulse waveform
%   release          - Allow property value and input characteristics 
%                      changes
%   clone            - Create a phase-coded waveform object with same 
%                      property values
%   isLocked         - Locked status (logical)
%   reset            - Reset states of phase-coded waveform object
%   plot             - Plot the phase-coded pulse waveform
%   bandwidth        - Bandwidth of the waveform
%   getMatchedFilter - Matched filter coefficients for the waveform
%
%   PhaseCodedWaveform properties:
%   
%   SampleRate            - Sample rate 
%   Code                  - Phase code 
%   ChipWidth             - Chip width
%   NumChips              - Number of chips
%   SequenceIndex         - Zadoff-Chu sequence index
%   PRF                   - Pulse repetition frequency
%   PRFSelectionInputPort - Enable PRF selection input 
%   OutputFormat          - Output signal format
%   NumSamples            - Number of samples in output
%   NumPulses             - Number of pulses in output
%   PRFOutputPort         - Enable PRF output port
%
%   % Examples:
%
%   % Example 1:
%   %   Create and plot a Zadoff-Chu coded pulse waveform.
%
%   waveform = phased.PhaseCodedWaveform('Code','Zadoff-Chu',...
%               'NumChips',16,'ChipWidth',1e-6,'OutputFormat','Pulses',...
%               'NumPulses',2);
%   plot(waveform);
%
%   % Example 2:
%   %   Generate output samples of the above waveform and plot the samples.
%
%   waveform = phased.PhaseCodedWaveform('Code','Zadoff-Chu',...
%               'NumChips',16,'ChipWidth',1e-6,'OutputFormat','Pulses',...
%               'NumPulses',2);
%   x = waveform();
%   plot(real(x)); title('Waveform output, real part');
%   xlabel('Samples'); ylabel('Amplitude (V)');
%
%   See also phased, phased.LinearFMWaveform, phased.SteppedFMWaveform.

%   Copyright 2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %SimulationTimeSource   Source of simulation sample time
    %   Set simulation time source to one of 'Derive from waveform
    %   parameters' | 'Inherit from Simulink engine', where the default is
    %   'Derive from waveform parameters'. This property applies only in
    %   Simulink and when you set the PRFSelectionInputPort to true. If the
    %   simulation time is set to 'Derive from waveform parameters', the
    %   block runs at a variable rate determined by the selected PRF if
    %   there are multiple PRFs involved so the elapsed time is
    %   proportional to the number of output samples; otherwise the block
    %   runs at a fixed rate so the elapsed time is a constant.  When you
    %   set the OutputFormat property to 'Pulses'. Otherwise, the block
    %   runs at a fixed rate.
    SimulationTimeSource = 'Derive from waveform parameters'
end

properties(Constant, Hidden)
    SimulationTimeSourceSet = matlab.system.StringSet(...
        {'Derive from waveform parameters',...
        'Inherit from Simulink engine'});
end

properties (Access = protected, Nontunable, Logical)
    %pInheritSampleTime  Inherit time
    %   Set this property to true to inherit sample time from Simulink
    %   engine. Set this property to false to set the block's own sample
    %   time. This property applies only in Simulink and when you set the
    %   PRFSelectionInputPort to true. When there are multiple PRFs
    %   involved and when you set the OutputFormat property to 'Pulses',
    %   the own sample time is set as a discrete variable sample time.
    pInheritSampleTime = false
end

methods

    function obj = SimulinkPhaseCodedWaveform(varargin)
        obj@phased.internal.AbstractPhaseCodedWaveform(varargin{:});
    end

end

methods (Access = protected)
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractPhaseCodedWaveform(obj,prop);
        if strcmp(prop,'SimulationTimeSource') && ...
                (~obj.PRFSelectionInputPort || ~isComplexityPropagated(obj))
            flag = true;
        end
    end
    
    function setupImpl(obj,~)
        setupImpl@phased.internal.AbstractPhaseCodedWaveform(obj);
        obj.pInheritSampleTime = strcmp(obj.SimulationTimeSource,...
            'Inherit from Simulink engine');
    end
    
    function flag = isInDVRMode(obj)
        flag = isComplexityPropagated(obj) && ...
            numel(obj.PRF)~=1 && ...
            ~(obj.PRFSelectionInputPort && obj.pInheritSampleTime);
    end

    function resetImpl(obj)
        resetImpl@phased.internal.AbstractPhaseCodedWaveform(obj)
        if obj.pOutputByPulse && isInDVRMode(obj)
                setNumTicksUntilNextHit(obj,1);
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractPhaseCodedWaveform(obj);
        if isLocked(obj)
            s.pInheritSampleTime = obj.pInheritSampleTime;
        end
    end
    
    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractPhaseCodedWaveform(obj,s);
        obj.pInheritSampleTime = s.pInheritSampleTime;
        s = rmfield(s,'pInheritSampleTime');
    end
    
    function loadObjectImpl(obj,s,wasLocked)
        loadObjectImpl@phased.internal.AbstractPhaseCodedWaveform(obj,s,wasLocked);
    end
end

methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        props = {...
            'SampleRate',...
            'Code',...
            'ChipWidth',...
            'NumChips',...
            'SequenceIndex',...
            'PRF',...
            'PRFSelectionInputPort',...
            'SimulationTimeSource',...
            'OutputFormat',...
            'NumSamples',...
            'NumPulses',...
            'PRFOutputPort'};
        groups = matlab.system.display.Section('Title', 'Parameters', ...
            'PropertyList', props);
    end
    
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:PhaseCodedWaveformTitle')),...
            'Text',getString(message('phased:library:block:PhaseCodedWaveformDesc')));
    end
end

methods (Access=protected)
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Phase-Coded');
    end
end

methods (Access = protected) %For Simulink propagation and mask
    function sts = getSampleTimeImpl(obj)
        if obj.PRFSelectionInputPort && ...
                strcmp(obj.SimulationTimeSource,'Inherit from Simulink engine')
            sts = createSampleTime(obj,'Type','Inherited');
        else
            if strcmp(obj.OutputFormat, 'Samples') || numel(obj.PRF) == 1
                N = obj.getPulseOutputSize(...
                    obj.SampleRate, obj.PRF,obj.OutputFormat, ...
                    obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
                st = phased.internal.samprate2time(obj.SampleRate,N);
                sts = createSampleTime(obj,'Type','Discrete',...
                    'SampleTime',st);
            else
                sts = createSampleTime(obj,'Type','Controllable',...
                    'TickTime',1/obj.SampleRate);
            end
        end
    end
end

end

% [EOF]