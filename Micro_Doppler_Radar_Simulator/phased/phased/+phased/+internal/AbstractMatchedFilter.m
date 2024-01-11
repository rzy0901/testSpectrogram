classdef (Hidden) AbstractMatchedFilter < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%   Copyright 2016-2017 The MathWorks, Inc.

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
    
properties (Nontunable) 
    %CoefficientsSource Source of coefficients
    %   Specify the source of matched filter coefficients using one of
    %   'Property' | 'Input port', where the default is 'Property'.
    CoefficientsSource = 'Property';   
end

properties
    %Coefficients Coefficients
    %   Specify matched filter coefficients as a column vector. The default
    %   value is [1;1]. This property applies when you set the
    %   CoefficientsSource property to 'Property'. This property is
    %   tunable.
    Coefficients = [1;1];
end

properties (Nontunable) 
    %SpectrumWindow Spectrum window
    %   Specify the window used for spectrum weighting using one of 'None'
    %   | 'Hamming' | 'Chebyshev' | 'Hann' | 'Kaiser' | 'Taylor' |
    %   'Custom', where the default is 'None'. Spectrum weighting is often
    %   used with linear FM waveform to reduce the sidelobes in the time
    %   domain.
    SpectrumWindow = 'None';
    %CustomSpectrumWindow   Custom spectrum window
    %   Specify the user-defined window for spectrum weighting using a
    %   function handle or a cell array. The default value of this property
    %   is @hamming. This property applies when you set the SpectrumWindow
    %   property to 'Custom'.
    %
    %   If CustomSpectrumWindow is a function handle, the specified
    %   function takes the window length as the input and generates
    %   appropriate window coefficients.
    %
    %   If CustomSpectrumWindow is a cell array, then the first cell must
    %   be a function handle. The specified function takes the window
    %   length as the first input argument, with other additional input
    %   arguments if necessary, and generates appropriate window
    %   coefficients. The remaining entries in the cell array are the
    %   additional input arguments to the function, if any.
    CustomSpectrumWindow = @hamming;
    %SpectrumRange  Spectrum window range (Hz)
    %   Specify the spectrum region on which the spectrum window is applied
    %   as a 1x2 vector in the form of [StartFrequency EndFrequency] (in
    %   Hz). The default value of this property is [0 1e5]. This property
    %   applies when you set the SpectrumWindow property to a value other
    %   than 'None'.
    %   
    %   Note that both StartFrequency and EndFrequency are measured in
    %   baseband, i.e., within [-Fs/2 Fs/2] where Fs is the sample rate
    %   that you specify in the SampleRate property. StartFrequency cannot
    %   be larger than EndFrequency.
    SpectrumRange = [0 1e5];
    %SampleRate     Sample rate (Hz)
    %   Specify the matched filter coefficients sample rate (in Hz) as a
    %   positive scalar. The default value of this property is 1e6. This
    %   property applies when you set SpectrumWindow property to a value
    %   other than 'None'.
    SampleRate = 1e6;
    %SidelobeAttenuation Sidelobe attenuation level (dB)
    %   Specify the sidelobe attenuation level (in dB) of a Chebyshev or
    %   Taylor window as a positive scalar. The default value of this
    %   property is 30. This property applies when you set the
    %   SpectrumWindow property to 'Chebyshev' or 'Taylor'.
    SidelobeAttenuation = 30;
    %Beta   Kaiser shape parameter
    %   Specify the parameter that affects the Kaiser window sidelobe
    %   attenuation as a nonnegative scalar. Please refer to <a href="matlab:help kaiser">kaiser</a> for  
    %   more details. The default value of this property is 0.5. This 
    %   property applies when you set the SpectrumWindow property to 
    %   'Kaiser'.
    Beta = 0.5;       
end

properties (Nontunable,PositiveInteger)
    %Nbar Number of constant level sidelobes
    %   Specify the number of nearly constant level sidelobes adjacent to
    %   the mainlobe in a Taylor window as a positive integer. The default
    %   value of this property is 4. This property applies when you set the
    %   SpectrumWindow property to 'Taylor'.
    Nbar = 4;
end



properties (Nontunable, Logical)  
    %GainOutputPort Enable SNR gain output
    %   Set this property to true to output the matched filter gain. The
    %   default value is false.
    GainOutputPort = false;
end

properties (Nontunable)
    %MaximumNumInputSamplesSource  Source of maximum number of samples
    %                       of the input signal
    %   Specify how the maximum number of samples of the input signal is
    %   specified as one of 'Auto' | 'Property', where the default is
    %   'Auto'. When you set this property to 'Auto', the matched filter 
    %   automatically allocates the memory to buffer the input signal.
    %   When you set this property to 'Property', the maximum number of
    %   samples in the input signal is specified via MaximumNumInputSamples
    %   property and any input signal longer than that value is truncated.
    %   This property applies when you set the MaximumDistanceSource
    %   property to 'Property'. The default value of this property is
    %   'Auto'.
    %
    %   To use the object in MATLAB Function Block in Simulink with
    %   variable-size signal, set this property to 'Property' and set
    %   the MaximumNumInputSamples property.
    MaximumNumInputSamplesSource = 'Auto'
end

properties (Nontunable, PositiveInteger)
    %MaximumNumInputSamples Maximum number of samples in input signal
    %   Specify the maximum number of samples in the input signal as a
    %   positive scalar. The input signal is the first input, X, and the
    %   number of samples is number of rows in X. This property applies
    %   when you set the MaximumNumInputSamplesSource property to
    %   'Property'. The default value of this property is 100.
    MaximumNumInputSamples = 100;
end

properties(Constant, Hidden)
    CoefficientsSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    SpectrumWindowSet = matlab.system.StringSet({'None','Hamming',...
        'Chebyshev','Hann','Kaiser','Taylor','Custom'});
    MaximumNumInputSamplesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
end

properties (Access = protected, Nontunable, Logical) 
    %pCoefficientsViaProp - private flag whether the coefficients are
    %specified via property
    pCoefficientsViaProp;
end

properties (Access = protected)
    %pCubeDim - 3D cube dimension
    pCubeDim = [-1 -1 -1]
end

properties (Access = protected)
    %pMFGain - private property to hold the gain of matched filter. 
    pMFGain;
    %cFFTFilter - private FFTFilter
    cFFTFilter;
    %pWinCoeff - private property to hold coefficients of spectrum window
    pWinCoeff;
end

properties (Access = protected, Logical, Nontunable) 
    %pInput3DFlag - indicate whether the input is a cube or not
    pInput3DFlag 
end

methods 

    function set.Coefficients(obj,value)
        validateattributes( value, { 'double','single' }, { 'vector', 'finite' }, '', 'Coefficients');
        obj.Coefficients = value(:);
    end
    
    function set.SidelobeAttenuation(obj,value)
        validateattributes( value, { 'double','single' }, { 'scalar', 'positive', 'finite' }, '', 'SidelobeAttenuation');
        obj.SidelobeAttenuation = value;
    end
    
    function set.Beta(obj,value)
        validateattributes( value, { 'double','single' }, { 'scalar', 'nonnegative', 'finite' }, '', 'Beta');
        obj.Beta = value;
    end
    
    function set.SampleRate(obj,value)
        validateattributes( value, { 'double','single' }, { 'scalar', 'positive', 'finite' }, '', 'SampleRate');
        obj.SampleRate = value;
    end
    
    function set.SpectrumRange(obj,value)
        validateattributes( value, { 'double','single' }, { 'real', 'finite', 'size', [ 1, 2 ] }, '', 'SpectrumRange');
        cond = value(1) > value(2);
        if cond
            coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidSpectrumRange');
        end
        obj.SpectrumRange = value;
    end
    
    function set.CustomSpectrumWindow(obj,value)
        cond = ~isa(value,'function_handle') && ~isa(value,'cell');
        if cond
            coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidCustomWindow','CustomSpectrumWindow');
        end
        cond = isa(value,'cell') && ~isa(value{1},'function_handle');
        if cond
            coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidCustomWindowCell','CustomSpectrumWindow');
        end
        obj.CustomSpectrumWindow = value;
    end
end

methods(Access = protected, Abstract)
    fs = computeSampleRate(obj)
end
    
methods (Access = protected)
    
    function obj = AbstractMatchedFilter(varargin)
        setProperties(obj, nargin, varargin{:});
    end

end

methods (Access = protected)
    
    function num = getNumInputsImpl(obj)
        num = 1;
        if strncmpi(obj.CoefficientsSource,'Input port',1)
            num = 2;            
        end
    end
    
    function validateInputsImpl(obj,x,w)
        cond = ndims(x) > 3;
        if cond
            coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidDimensions');
        end
        cond = ~isa(x,'float');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','float');
        end
        
        validateNumChannels(obj,x);
        if ndims(x) == 3 && obj.pCubeDim(3) ~= -1
            validateNumPages(obj,x,obj.pCubeDim(3));
        end
        
        if getNumInputs(obj) > 1
            cond = ~isa(w,'float');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Coeff','float');
            end
            cond = ~iscolumn(w) || isempty(w);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeColVector','Coeff');
            end
        end
    end
    
    function setupImpl(obj,x,~)
        obj.pCoefficientsViaProp = ...
            (obj.CoefficientsSource(1) == 'P'); %Property
        
        obj.pNumInputChannels = size(x,2);
        obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        
        if ndims(x) == 3
            obj.pInput3DFlag = true;
            obj.pCubeDim = size(x);
        else
            obj.pInput3DFlag = false;
        end
        
        if isempty(coder.target)
            cSW = obj.CustomSpectrumWindow;
        else
            cond = (obj.SpectrumWindow(2) == 'u'); %'Custom'
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:phased:MatchedFilter:NoCodegenCustom','SpectrumWindow','Custom');
            end
            %'function_handle' as property not supported in codegen
            cSW = [];
        end
        
        %Fs = getSampleRate(obj,size(x,1),1,obj.SampleRate);
        if ~strcmp(obj.SpectrumWindow,'None')
            Fs = computeSampleRate(obj); 
            
            cond = any(obj.SpectrumRange < -Fs/2) || ...
                    any(obj.SpectrumRange > Fs/2);
            if cond
                coder.internal.errorIf(cond,'phased:phased:MatchedFilter:InvalidBaseband',...
                    num2str(-Fs/2), num2str(Fs/2),...
                    num2str(obj.SpectrumRange(1)), num2str(obj.SpectrumRange(2)));
            end
        else
            Fs = 1;
        end
        
        if isInputDataSizePropagated(obj) % Simulink block
            obj.cFFTFilter = phased.internal.FFTFilter(...
                'CoefficientsSource',obj.CoefficientsSource,...
                'WindowLossOutputPort',true,...
                'SpectrumWindow',obj.SpectrumWindow,...
                'Beta',obj.Beta, ...
                'SidelobeAttenuation', obj.SidelobeAttenuation, ...
                'Nbar', obj.Nbar, ...
                'CustomSpectrumWindow',cSW , ...
                'SpectrumRange',obj.SpectrumRange,...
                'SampleRate',Fs,...
                'MaxNumInputSamplesSource','Property',...
                'MaxNumInputSamples',getPropagatedNumInputSamples(obj,x));
        else
            if strcmp(obj.MaximumNumInputSamplesSource,'Auto')
                obj.cFFTFilter = phased.internal.FFTFilter(...
                    'CoefficientsSource',obj.CoefficientsSource,...
                    'WindowLossOutputPort',true,...
                    'SpectrumWindow',obj.SpectrumWindow,...
                    'Beta',obj.Beta, ...
                    'SidelobeAttenuation', obj.SidelobeAttenuation, ...
                    'Nbar', obj.Nbar, ...
                    'CustomSpectrumWindow',cSW , ...
                    'SpectrumRange',obj.SpectrumRange,...
                    'SampleRate',Fs,...
                    'MaxNumInputSamplesSource','Property',...
                    'MaxNumInputSamples',getPropagatedNumInputSamples(obj,x));
            else
                obj.cFFTFilter = phased.internal.FFTFilter(...
                    'CoefficientsSource',obj.CoefficientsSource,...
                    'WindowLossOutputPort',true,...
                    'SpectrumWindow',obj.SpectrumWindow,...
                    'Beta',obj.Beta, ...
                    'SidelobeAttenuation', obj.SidelobeAttenuation, ...
                    'Nbar', obj.Nbar, ...
                    'CustomSpectrumWindow',cSW , ...
                    'SpectrumRange',obj.SpectrumRange,...
                    'SampleRate',Fs,...
                    'MaxNumInputSamplesSource','Property',...
                    'MaxNumInputSamples',obj.MaximumNumInputSamples);
            end
        end
        
        processTunedPropertiesImpl(obj);
    end
    
    function processInputSizeChangeImpl(obj,x,~)
        if obj.pInput3DFlag
            obj.pCubeDim = size(x);
        end
    end

    function flag = isInputComplexityLockedImpl(obj,index)
        flag = false;
        if index == 1
            flag = false;
        end
        if (obj.CoefficientsSource(1) == 'I') ... %Input port
                && (index == 2)
            flag = false;
        end
    end
    
    function flag = isOutputComplexityLockedImpl(obj,index)
        if index == 1
            flag = false;
        end
        if obj.GainOutputPort && (index == 2)
            flag = true;
        end
    end
    
    function processTunedPropertiesImpl(obj)
        if obj.pCoefficientsViaProp
            obj.pMFGain = pow2db(real(obj.Coefficients'*obj.Coefficients));
            obj.cFFTFilter.Coefficients = obj.Coefficients;
        end
    end

    function releaseImpl(obj)
        releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
        release(obj.cFFTFilter);
    end
    
    function resetImpl(obj)
        resetImpl@phased.internal.AbstractSampleRateEngine(obj);
        reset(obj.cFFTFilter);
    end
    
    function num = getNumOutputsImpl(obj)
        num = 1;
        if obj.GainOutputPort
            num = 2;
        end
    end
    
    function flag = isInputSizeLockedImpl(~,~)
        flag = false;
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
        s.isLocked = isLocked(obj);
        if isLocked(obj)
            s.cFFTFilter = saveobj(obj.cFFTFilter);
            s.pInput3DFlag = obj.pInput3DFlag;
            s.pCoefficientsViaProp = obj.pCoefficientsViaProp;
            s.pCubeDim = obj.pCubeDim;
            s.pMFGain = obj.pMFGain;
            s.pWinCoeff = obj.pWinCoeff;
        end
    end

    function s = loadSubObjects(obj,s)
        if isfield(s,'isLocked')
            if s.isLocked
                obj.cFFTFilter = phased.internal.FFTFilter.loadobj(s.cFFTFilter);
                s = rmfield(s,'cFFTFilter');
            end
            s = rmfield(s,'isLocked');
        end
    end

    function [sigout,g] = stepImpl(obj,siginArg,w)

        
        classtouse = class(siginArg);
        if obj.pInput3DFlag
            % input is a cube
            % need to reshape, process and transform it back
            sigin = reshape(siginArg,size(siginArg,1),[]);
        else
            sigin = siginArg;
        end
        
        if obj.pCoefficientsViaProp
            [sig, procloss] = step(obj.cFFTFilter,sigin);
            g = cast(obj.pMFGain,classtouse) - procloss;
        else
            w = cast(w,classtouse);
            cond = ~all(isfinite(w));
            if cond        
                coder.internal.errorIf(cond,'phased:step:expectedFinite','COEFF');
            end
            [sig, procloss] = step(obj.cFFTFilter,sigin,w);
            g = pow2db(norm(w)^2) - procloss;
        end
                   
        if obj.pInput3DFlag
            % input is a cube, so we need to transform the result back
            % to its original shape.
            sigout = reshape(sig,size(siginArg));
        else
            sigout = sig(1:size(siginArg,1),1:size(siginArg,2));
        end

    end

    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        if strcmp(prop,'Coefficients') && ...
                strcmp(obj.CoefficientsSource,'Input port')
            flag = true;
        end
        if strcmp(prop,'CustomSpectrumWindow') && ...
                ~strcmp(obj.SpectrumWindow,'Custom')
            flag = true;
        end
        if strcmp(prop,'SpectrumRange') && ...
                strcmp(obj.SpectrumWindow,'None')
            flag = true;
        end
        if strcmp(prop,'SampleRate') && ...
                strcmp(obj.SpectrumWindow,'None')
            flag = true;
        end
        if strcmp(prop,'SidelobeAttenuation') && ...
                ~(strcmp(obj.SpectrumWindow,'Chebyshev') || ...
                strcmp(obj.SpectrumWindow,'Taylor'))
            flag = true;
        end
        if strcmp(prop,'Beta') && ~strcmp(obj.SpectrumWindow,'Kaiser')
            flag = true;
        end
        if strcmp(prop,'Nbar') && ~strcmp(obj.SpectrumWindow,'Taylor')
            flag = true;
        end
        if strcmp(prop,'MaximumNumInputSamples') && ...
                (obj.MaximumNumInputSamplesSource(1) == 'A')
            flag = true;
        end
    end
end

methods (Access = protected)
    function varargout = getOutputSizeImpl(obj)
        varargout{1} = propagatedInputSize(obj,1);
        if obj.GainOutputPort
            varargout{2} = [1 1];
        end
    end
    function varargout = isOutputFixedSizeImpl(obj)
        varargout{1} = propagatedInputFixedSize(obj, 1);
        if obj.GainOutputPort
            varargout{2} = true;
        end
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj,1);
        if obj.GainOutputPort
            varargout{2} = propagatedInputDataType(obj,1);
        end
    end
    function varargout = isOutputComplexImpl(obj)
        varargout{1} = true;
        if obj.GainOutputPort
            varargout{2} = false;
        end
    end
end

end

