classdef (Sealed,StrictDefaults) RectangularWaveform < phased.internal.AbstractRectangularWaveform
%RectangularWaveform Rectangular pulse waveform
%   H = phased.RectangularWaveform creates a rectangular pulse waveform
%   System object, H. This object generates samples of a rectangular pulse.
%
%   H = phased.RectangularWaveform(Name,Value) creates a rectangular pulse
%   waveform object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H) returns samples of the rectangular pulse in a column vector
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
%   RectangularWaveform methods:
%
%   step             - Return samples of the rectangular pulse waveform
%   release          - Allow property value and input characteristics 
%                      changes
%   clone            - Create a rectangular waveform object with same 
%                      property values
%   isLocked         - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>            - Reset states of rectangular waveform object
%   plot             - Plot the rectangular pulse waveform
%   bandwidth        - Bandwidth of the waveform
%   getMatchedFilter - Matched filter coefficients for the waveform
%
%   RectangularWaveform properties:
%
%   SampleRate            - Sample rate 
%   DurationSpecification - Method to specify pulse duration
%   PulseWidth            - Pulse width 
%   DutyCycle             - Duty cycle
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
%   %   Create and plot a rectangular pulse waveform object. 
%
%   waveform = phased.RectangularWaveform('PulseWidth',1e-5,...
%           'OutputFormat','Pulses','NumPulses',2);
%   plot(waveform);
%
%   % Example 2:
%   %   Generate output samples of the above waveform and plot the samples.
%
%   waveform = phased.RectangularWaveform('PulseWidth',1e-5,...
%           'OutputFormat','Pulses','NumPulses',2);
%   x = waveform();
%   plot(real(x)); title('Waveform output, real part');
%   xlabel('Samples'); ylabel('Amplitude (V)');
%
%   See also phased, phased.LinearFMWaveform, phased.PhaseCodedWaveform,
%   phased.SteppedFMWaveform.

%   Copyright 2008-2017 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

methods (Access = public)
    function obj = RectangularWaveform(varargin)
        obj@phased.internal.AbstractRectangularWaveform(varargin{:});
    end
    
end

methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        props = {...
            'SampleRate',...
            'DurationSpecification',...
            'PulseWidth',...
            'DutyCycle',...
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
        a = 'phasedwavlib/Rectangular Waveform';
    end
end
end 
