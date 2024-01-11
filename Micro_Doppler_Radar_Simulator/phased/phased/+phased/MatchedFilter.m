classdef (Sealed,StrictDefaults) MatchedFilter < phased.internal.AbstractMatchedFilter
%MatchedFilter Matched filter
%   H = phased.MatchedFilter creates a matched filter System object, H.
%   This object performs matched filtering on the input data.
%
%   H = phased.MatchedFilter(Name,Value) creates a matched filter object,
%   H, with the specified property Name set to the specified Value. You
%   can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) applies the matched filtering to the input X and returns
%   the filtered result in Y. The filter is applied along the first
%   dimension. Y and X have the same dimensions.
%
%   [Y,G] = step(H,X) returns additional output G as the gain (in dB) of
%   the matched filter when you set the GainOutputPort property to true.
%
%   Y = step(H,X,COEFF) uses the input COEFF as the matched filter
%   coefficients when you set the CoefficientsSource property to 'Input
%   port'. COEFF must be a column vector.
%   
%   You can combine optional input and output arguments when their enabling
%   properties are set. Optional inputs and outputs must be listed in the
%   same order as the order of the enabling properties. For example,
%
%   [Y,G] = step(H,X,COEFF)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   MatchedFilter methods:
%
%   step     - Perform matched filtering (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a matched filter object with same property values
%   isLocked - Locked status (logical)
%
%   MatchedFilter properties:
%
%   CoefficientsSource           - Source of coefficients
%   Coefficients                 - Coefficients
%   SpectrumWindow               - Spectrum window
%   CustomSpectrumWindow         - Custom spectrum window
%   SpectrumRange                - Spectrum window range
%   SampleRate                   - Sample rate
%   SidelobeAttenuation          - Sidelobe attenuation level
%   Beta                         - Kaiser shape parameter
%   Nbar                         - Number of constant level sidelobes
%   GainOutputPort               - Enable SNR gain output
%   MaximumNumInputSamplesSource - Source of maximum number of samples
%                                  of the input signal
%   MaximumNumInputSamples       - Maximum number of samples in input 
%                                  signal
%
%   % Examples:
%
%   % Example 1:
%   %   Apply the matched filter for a linear FM waveform.
%
%   waveform = phased.LinearFMWaveform('PulseWidth',1e-4,'PRF',5e3);
%   x = waveform();
%   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(waveform));
%   y = mf(x);
%   subplot(211),plot(real(x));
%   xlabel('Samples'),ylabel('Amplitude'),title('Input Signal');
%   subplot(212),plot(real(y));
%   xlabel('Samples'),ylabel('Amplitude'),title('Matched Filter Output');
%
%   % Example 2:
%   %   Apply the matched filter, using a Hamming window to do spectrum 
%   %   weighting.
%
%   waveform = phased.LinearFMWaveform('PulseWidth',1e-4,'PRF',5e3);
%   x = waveform();
%   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(waveform),...
%           'SpectrumWindow','Hamming');
%   y = mf(x);
%   subplot(211),plot(real(x));
%   xlabel('Samples'),ylabel('Amplitude'),title('Input Signal');
%   subplot(212),plot(real(y));
%   xlabel('Samples'),ylabel('Amplitude'),title('Matched Filter Output');
%
%   % Example 3:
%   %   Apply the matched filter, using a custom Gaussian window for 
%   %   spectrum weighting.
%
%   waveform = phased.LinearFMWaveform('PulseWidth',1e-4,'PRF',5e3);
%   x = waveform();
%   mf = phased.MatchedFilter('Coefficients',getMatchedFilter(waveform),...
%           'SpectrumWindow','Custom',...
%           'CustomSpectrumWindow',{@gausswin,2.5});
%   y = mf(x);
%   subplot(211),plot(real(x));
%   xlabel('Samples'),ylabel('Amplitude'),title('Input Signal');
%   subplot(212),plot(real(y));
%   xlabel('Samples'),ylabel('Amplitude'),title('Matched Filter Output');
%
%   See also phased, phased.CFARDetector, phased.TimeVaryingGain,
%   phased.StretchProcessor, pulsint.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [2] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001
%   [3] Alan Oppenheim and Ronald Schafer, Discrete-Time Signal Processing,
%       Prentice Hall, 1989


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
methods
    
    function obj = MatchedFilter(varargin)
        obj@phased.internal.AbstractMatchedFilter(varargin{:});
    end

end

methods (Access = protected)
    
    function Fs = computeSampleRate(obj)
        Fs = obj.SampleRate;
    end

    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

end

methods (Static,Hidden,Access=protected)  
    function groups = getPropertyGroupsImpl
        groups = matlab.system.display.Section(...
            'phased.MatchedFilter');
    end
        
end
    
methods (Static,Hidden)
    function a = getAlternateBlock
        a = 'phaseddetectlib/Matched Filter';
    end
end
end

