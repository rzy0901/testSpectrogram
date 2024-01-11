classdef(Sealed,StrictDefaults) PulseWaveformLibrary < phased.internal.AbstractPulseWaveformLibrary   
%PulseWaveformLibrary Library of Pulse Waveforms
%   H = phased.PulseWaveformLibrary creates a pulse waveform Library System
%   object, H. This object generates samples of the selected waveform from
%   a predefined waveform list.
%
%   H = phased.PulseWaveformLibrary(Name,Value) creates a pulse waveform
%   library object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,IDX) returns samples of IDX-th pulse waveform in the library
%   in a column vector Y. Y contain one pulse of the IDX-th waveform. IDX
%   must be a positive integer.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H) and y = H() are equivalent.
%
%   PulseWaveformLibrary methods:
%
%   step                - Return samples of IDX-th pulse waveform in the 
%                         library
%   release             - Allow property value and input characteristics
%   clone               - Create a PulseWaveform object with same 
%                         property values
%   isLocked            - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>               - Reset states of PulseWaveform object
%   plot                - Plot the IDX-th pulse waveform 
%   getMatchedFilter    - Matched filter coefficients IDX th pulse waveform
%
%   PulseWaveformLibrary properties:
%
%   SampleRate            - Sample rate
%   WaveformSpecification - Specify the type of pulsed waveform and its
%                           properties
%
%   %Examples
%
%   % Example 1:
%   %   Create and plot a 2-element pulse waveform library.
%   waveform1 = {'Rectangular','PRF',1e4, 'PulseWidth', 50e-6};
%   waveform2 = {'LinearFM','PRF',1e4,'PulseWidth',50e-6,...
%                'SweepBandwidth',1e5,'SweepDirection','Up',...
%                'SweepInterval', 'Positive'};
% 
%   wavlib = phased.PulseWaveformLibrary('SampleRate',1e6,...
%               'WaveformSpecification',{waveform1,waveform2});
%   plot(wavlib,1);
%
%   % Example 2:
%   %   Generate output samples of the waveform library and plot the
%   %   samples.
%
%   waveform1 = {'Rectangular','PRF',1e4, 'PulseWidth', 50e-6};
%   waveform2 = {'LinearFM','PRF',1e4,'PulseWidth',50e-6,...
%                'SweepBandwidth',1e5,'SweepDirection','Up',...
%                'SweepInterval', 'Positive'};
% 
%   wavlib = phased.PulseWaveformLibrary('SampleRate',1e6,...
%               'WaveformSpecification',{waveform1,waveform2});
%   x1 = wavlib(1);
%   x2 = wavlib(2);
%   subplot(2,1,1); plot(wavlib,1,'PlotType','real','PulseIdx',1);
%   xlabel('Samples'); ylabel('Amplitude (V)');
%   subplot(2,1,2); plot(wavlib,2,'PlotType','real','PulseIdx',1);
%   xlabel('Samples'); ylabel('Amplitude (V)');

%   Copyright 2017 The MathWorks, Inc.

%   Reference
%   [1] Cochran et.al., Waveform Libraries, Measures of effectiveness for
%   radar scheduling,IEEE Signal Processing Magazine, 2009.
%
%   [2]Richards, M. A., Fundamentals of Radar Signal Processing. New York:
%   McGraw-Hill, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

methods
    % Constructor
    
    function obj = PulseWaveformLibrary(varargin)
        obj@phased.internal.AbstractPulseWaveformLibrary(varargin{:});
    end
end

methods (Access = protected)
    function validateDuplicateSetting (obj,value)
        validateDuplicateWarningCond(obj,value);
    end
end

methods (Static, Hidden, Access = protected)
    
    function groups = getPropertyGroupsImpl
        props = {'SampleRate',...
            'WaveformSpecification',....
            };
        groups = matlab.system.display.Section('Title',...
            'Pulse Waveform Library', ...
            'PropertyList', props);
        
    end
    
end

methods (Static,Hidden)
    function a = getAlternateBlock
        a = 'phasedwavlib/Pulse Waveform Library';
    end
    function offset = getFreqOffset(inp)
        idx  = find(strcmpi(inp,'FrequencyOffset'),1,'last');
        if ~isempty(idx)
            offset = inp{idx+1};
        else
            offset = 0;
        end
    end
end
end
