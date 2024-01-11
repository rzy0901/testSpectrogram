classdef (Sealed, StrictDefaults) SimulinkMatchedFilter < phased.internal.AbstractMatchedFilter & ...
        matlab.system.mixin.CustomIcon 
%This class is for internal use only. It may be removed in the future.
    
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
%   MatchedFilter methods:
%
%   step     - Perform matched filtering (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a matched filter object with same property values
%   isLocked - Locked status (logical)
%
%   MatchedFilter properties:
%
%   CoefficientsSource   - Source of coefficients
%   Coefficients         - Coefficients
%   SpectrumWindow       - Spectrum window
%   CustomSpectrumWindow - Custom spectrum window
%   SpectrumRange        - Spectrum window range
%   SampleRate           - Sample rate
%   SidelobeAttenuation  - Sidelobe attenuation level
%   Beta                 - Kaiser shape parameter
%   Nbar                 - Number of constant level sidelobes
%   GainOutputPort       - Enable SNR gain output
%
%   % Examples:
%
%   % Example 1:
%   %   Apply the matched filter for a linear FM waveform.
%
%   hw = phased.LinearFMWaveform('PulseWidth',1e-4,'PRF',5e3);
%   x = step(hw);
%   hmf = phased.MatchedFilter('Coefficients',getMatchedFilter(hw));
%   y = step(hmf,x);
%   subplot(211),plot(real(x));
%   xlabel('Samples'),ylabel('Amplitude'),title('Input Signal');
%   subplot(212),plot(real(y));
%   xlabel('Samples'),ylabel('Amplitude'),title('Matched Filter Output');
%
%   % Example 2:
%   %   Apply the matched filter, using a Hamming window to do spectrum 
%   %   weighting.
%
%   hw = phased.LinearFMWaveform('PulseWidth',1e-4,'PRF',5e3);
%   x = step(hw);
%   hmf = phased.MatchedFilter('Coefficients',getMatchedFilter(hw),...
%           'SpectrumWindow','Hamming');
%   y = step(hmf,x);
%   subplot(211),plot(real(x));
%   xlabel('Samples'),ylabel('Amplitude'),title('Input Signal');
%   subplot(212),plot(real(y));
%   xlabel('Samples'),ylabel('Amplitude'),title('Matched Filter Output');
%
%   % Example 3:
%   %   Apply the matched filter, using a custom Gaussian window for 
%   %   spectrum weighting.
%
%   hw = phased.LinearFMWaveform('PulseWidth',1e-4,'PRF',5e3);
%   x = step(hw);
%   hmf = phased.MatchedFilter('Coefficients',getMatchedFilter(hw),...
%           'SpectrumWindow','Custom',...
%           'CustomSpectrumWindow',{@gausswin,2.5});
%   y = step(hmf,x);
%   subplot(211),plot(real(x));
%   xlabel('Samples'),ylabel('Amplitude'),title('Input Signal');
%   subplot(212),plot(real(y));
%   xlabel('Samples'),ylabel('Amplitude'),title('Matched Filter Output');
%
%   See also phased, phased.CFARDetector, phased.TimeVaryingGain,
%   phased.StretchProcessor, pulsint.

%   Copyright 2016 The MathWorks, Inc.

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
    
properties (Nontunable, Logical)
    %SampleRateFromInputCheckbox Inherit sample rate 
    %   Set SampleRateFromInputCheckbox to true to derive sample rate from
    %   Simulink time engine. Set SampleRateFromInputCheckbox to false to
    %   specify the sample rate. This property applies when used in
    %   Simulink.
    SampleRateFromInputCheckbox = true
end
    
methods
    
    function obj = SimulinkMatchedFilter(varargin)
        obj@phased.internal.AbstractMatchedFilter(varargin{:});
    end

end

methods (Access = protected)
    
    function fs = computeSampleRate(obj)
        if obj.SampleRateFromInputCheckbox
            fs = getSampleRateInSimulation(obj);
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
        else
            fs = obj.SampleRate;
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
        flag = isInactivePropertyImpl@phased.internal.AbstractMatchedFilter(obj, prop);
        if strcmp(prop,'SampleRateFromInputCheckbox') && ...
                strcmp(obj.SpectrumWindow,'None')
            flag = true;
        end
        if strcmp(prop,'SampleRate') && ...
                obj.SampleRateFromInputCheckbox
            flag = true;
        end
    end
end

methods (Static,Hidden,Access=protected)  
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:MatchedFilterTitle')),...
          'Text',getString(message('phased:library:block:MatchedFilterDesc')));
    end
    
    function groups = getPropertyGroupsImpl
        groups = matlab.system.display.Section(...
            'phased.internal.SimulinkMatchedFilter');
        groups.PropertyList = groups.PropertyList([2:6 1 7:end]);
        dSpectrumWindow = matlab.system.display.internal.Property(...
            'SpectrumWindow', 'StringSetValues', {'None','Hamming',...
            'Chebyshev','Hann','Kaiser','Taylor'});
        dMaximumNumInputSamplesSource = matlab.system.display.internal.Property(...
            'MaximumNumInputSamplesSource','IsGraphical',false);
        dMaximumNumInputSamples = matlab.system.display.internal.Property(...
            'MaximumNumInputSamples','IsGraphical',false);
        for m = 1:numel(groups.PropertyList)
            if strcmp(groups.PropertyList{m},'SpectrumWindow')
                groups.PropertyList{m} = dSpectrumWindow;
            elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamplesSource')
                groups.PropertyList{m} = dMaximumNumInputSamplesSource;
            elseif strcmp(groups.PropertyList{m},'MaximumNumInputSamples')
                groups.PropertyList{m} = dMaximumNumInputSamples;
            end
        end
    end
        
end
    
methods (Access = protected)
    function varargout = getInputNamesImpl(obj)  
        if strcmp(obj.CoefficientsSource,'Input port')
            varargout = {'X','Coeff'};
        else
            varargout = {'X'};
        end
    end

    function varargout = getOutputNamesImpl(obj) 
        if obj.GainOutputPort
            varargout = {'Y','G'};
        else
            varargout = {'Y'};
        end
    end
        
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Matched\nFilter');
    end
    
end

end

