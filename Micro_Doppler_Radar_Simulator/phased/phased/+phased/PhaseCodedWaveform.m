classdef (Sealed,StrictDefaults) PhaseCodedWaveform < phased.internal.AbstractPhaseCodedWaveform
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
%   <a href="matlab:help matlab.System/reset   ">reset</a>            - Reset states of phase-coded waveform object
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

%   Copyright 2011-2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

methods

    function obj = PhaseCodedWaveform(varargin)
        obj@phased.internal.AbstractPhaseCodedWaveform(varargin{:});
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
            'OutputFormat',...
            'NumSamples',...
            'NumPulses',...
            'PRFOutputPort'};
        groups = matlab.system.display.Section('Title', 'Parameters', ...
            'PropertyList', props);
    end
    
end

methods (Static,Hidden)
    function a = getAlternateBlock
        a = 'phasedwavlib/Phase-Coded Waveform';
    end
end
end

% [EOF]
