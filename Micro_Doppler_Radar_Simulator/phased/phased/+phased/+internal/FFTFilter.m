classdef (Sealed,Hidden, StrictDefaults) FFTFilter < phased.internal.AbstractVarSizeEngine
%This class is for internal use only. It may be removed in the future.

%FFTFilter   FFT FIR filter
%   H = phased.internal.FFTFilter creates an FFT FIR filter System object,
%   H. This object performs filtering using FFT with overlap-add method
%
%   H = phased.internal.FFTFilter(Name,Value) creates an FFT FIR filter
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) applies the FFT FIR filter object H to the input X and
%   returns the filtered result in Y. The filter is applied along the first
%   dimension. Y and X have the same dimensions. 
%
%   Y = step(H,X,COEFF) uses the input COEFF as the FFT FIR filter
%   coefficients when you set the CoefficientsSource property to 'Input
%   port'. COEFF must be a column vector.
%
%   [Y,PL] = step(H,X) returns the extra output PL as the processing loss
%   (in dB) due to the usage of the spectrum window.
%
%   You can combine optional input and output arguments when their enabling
%   properties are set. Optional inputs and outputs must be listed in the
%   same order as the order of the enabling properties. For example,
%
%   [Y,PL] = step(H,X,COEFF)
%   
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   FFTFilter methods:
%
%   step     - Perform matched filtering (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a matched filter object with same property values
%   isLocked - Locked status (logical)
%
%   FFTFilter properties:
%        
%   CoefficientsSource - Source of matched filter coefficients
%   Coefficients       - Matched filter coefficients
%   SpectrumWindow     - FFT filter spectrum window
%   SpectrumRange      - Spectrum window coverage region
%   SampleRate         - Coefficient sample rate
%   GainOutputPort     - Output gain
%
%   % Example:
%   %   Construct an FFT FIR filter to perform matched filtering for a
%   %   linear FM waveform.
%
%   hw = phased.LinearFMWaveform('PulseWidth',1e-4,'PRF',5e3);
%   x = step(hw);
%   hf = phased.internal.FFTFilter('Coefficients',getMatchedFilter(hw));
%   y = step(hf,x);
%   subplot(211),plot(real(x));
%   xlabel('Samples'),ylabel('Amplitude'),title('Input Signal');
%   subplot(212),plot(real(y));
%   xlabel('Samples'),ylabel('Amplitude'),title('Matched Filter Output');

%   Copyright 2010-2016 The MathWorks, Inc.
%     

%   Reference
%   [1] Alan Oppenheim and Ronald Schafer, Discrete-Time Signal Processing,
%       Prentice Hall, 1989
%   [2] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005, pp 204


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable) 
    %CoefficientsSource Source of matched filter coefficients
    %   Specify the source of FFT filter coefficients using one of
    %   'Property' | 'Input port', where the default is 'Property'.
    CoefficientsSource = 'Property';   
end

properties
    %Coefficients FFT filter coefficients
    %   Specify FFT filter coefficients as a column vector. The default
    %   value is 1. This property applies when you set the
    %   CoefficientsSource property to 'Property'. This property is
    %   tunable.
    Coefficients = 1;
end

properties (Nontunable)
    %SpectrumWindow Window for spectrum weighting
    %   Specify the window used for spectrum weighting using one of 'None'
    %   | 'Hamming' | 'Chebyshev' | 'Hann' | 'Kaiser' | 'Taylor' |
    %   'Custom', where the default is 'None'. Spectrum weighting is often
    %   used with linear FM waveform to reduce the sidelobes in the time
    %   domain.
    SpectrumWindow;
    %SpectrumRange  Spectrum window coverage region
    %   Specify the spectrum region on which the spectrum window is applied
    %   as a 2-element row vector in the form of [StartFrequency
    %   EndFrequency] (in Hz). The default value of this property is [0
    %   1e5]. This property applies when you set SpectrumWindow property to
    %   a value other than 'None'.
    %   
    %   Note that both StartFrequency and EndFrequency is measured in
    %   baseband, i.e., within [-Fs/2 Fs/2] where Fs is the sample rate.
    %   StartFrequency cannot be larger than EndFrequency.
    SpectrumRange;
    %SampleRate     Coefficient sample rate
    %   Specify the FFT filter coefficients sample rate (in Hz) as a
    %   positive scalar. The default value of this property is 1e6. This
    %   property applies when you set SpectrumWindow property to a value
    %   other than 'None'.
    SampleRate;
    %Beta   Kaiser window parameter
    %   Specify the parameter that affects the Kaiser window sidelobe
    %   attenuation as a nonnegative scalar. Please refer to <a href="matlab:help kaiser">kaiser</a> for  
    %   more details. The default value of this property is 0.5. This 
    %   property applies when you set the SpectrumWindow property to 
    %   'Kaiser'.
    Beta;
    %SidelobeAttenuation Window sidelobe attenuation level
    %   Specify the sidelobe attenuation level (in dB) of a Chebyshev or
    %   Taylor window as a positive scalar. The default value of this
    %   property is 30. This property applies when you set the
    %   SpectrumWindow property to 'Chebyshev' or 'Taylor'.
    SidelobeAttenuation;
    %CustomSpectrumWindow   User-defined window for spectrum weighting
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
    CustomSpectrumWindow;
end

properties (Nontunable,PositiveInteger)
    %Nbar Number of nearly constant level sidelobes in Taylor window
    %   Specify the number of nearly constant level sidelobes adjacent to
    %   the mainlobe in a Taylor window as a positive integer. The default
    %   value of this property is 4. This property applies when you set the
    %   SpectrumWindow property to 'Taylor'.
    Nbar = 4;
end

properties (Nontunable)
    MaxNumInputSamplesSource = 'Auto'
end

properties (Nontunable)
    MaxNumInputSamples = 0;  
end

properties (Logical, Nontunable)
    %WindowLossOutputPort Output processing loss from spectrum window
    %   Set this property to true to output the process loss (in dB) caused
    %   by spectrum windowing. Set this property to false to not output the
    %   process loss. The default value of this property is false.
    WindowLossOutputPort = false;
end

properties(Constant, Hidden)
    CoefficientsSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    MaxNumInputSamplesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
end

properties (Access = private, Nontunable, Logical) 
    %pCoefficientsViaProp - private flag whether the coefficients are
    %specified via property
    pCoefficientsViaProp;
end

properties (Access = private, Logical)
    %pIsSignalLonger - private property to indicate whether the signal is
    %longer than the length of the state
    pIsSignalLonger
end

properties (Access = private, Nontunable)
    %pFFTLength - private property to hold FFT length
    pFFTLength
    %pNumTotalChannels - private property to hold total channels. If
    %matrix, number of total channels is the same as the number of columns.
    pNumTotalChannels
    %pWindowCoefficients - window coefficients
    pWindowCoefficients
    %pWindowProcessingLoss - window processing loss
    pWindowProcessingLoss
    %pMaxNumInputSamples - maximum number of input samples per channel
    pMaxNumInputSamples = 0
end

properties (Access = private)
    %pOverlapLength - private property to hold overlap length
    pOverlapLength
    %pSignalLength - private property to hold signal length
    pSignalLength
    %pSignalAddIdx - signal overlap add portion index
    pSignalAddIdx
    %pStateRemainingIdx - state remaining portion index
    pStateRemainingIdx
    %pStateUpdateIdx - state update portion index
    pStateUpdateIdx
    %pSignalIdx - signal index
    pSignalIdx
    %pStateIdx - state index
    pStateIdx
end

properties (Access = private)
    %pCoefficientsFFT - FFT of coefficients
    pCoefficientsFFT
    %pOverlapState - overlap state samples in overlap add
    pOverlapState
    cFFT
    cCoeffFFT
    cIFFT
end

properties (Access = private, Logical)
    pSizeInitialized
end
    
methods

    function obj = FFTFilter(varargin)
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
    
    function num = getNumOutputsImpl(obj)
        num = 1;
        if obj.WindowLossOutputPort
            num = 2;
        end
    end
    
    function validateInputsImpl(obj,x,coeff)        
        if ~ismatrix(x)
            error(message('phased:phased:FFTFilter:InvalidDimensions'));
        end
        if ~isa(x,'float')
          matlab.system.internal.error(...
            'MATLAB:system:invalidInputDataType','X','float');
        end
        validateNumChannels(obj,x)
        
        if getNumInputs(obj) > 1
            if ~iscolumn(coeff) || isempty(coeff)
              matlab.system.internal.error(...
                'MATLAB:system:inputMustBeColVector','COEFF');
            end
            if ~isa(coeff,'float')
              matlab.system.internal.error(...
                'MATLAB:system:invalidInputDataType','COEFF','float');
            end
        end
    end
    
    function setupImpl(obj,x,coeff)
        obj.pNumInputChannels = getNumChannels(obj,x);
        obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        obj.pNumTotalChannels = obj.pValidatedNumInputChannels;
        
        obj.pCoefficientsViaProp = ...
            (obj.CoefficientsSource(1) == 'P');%Property
        if obj.pCoefficientsViaProp
            coefflen = numel(obj.Coefficients);
        else
            coefflen = size(coeff,1);
        end
        
        if isInputDataSizePropagated(obj) || strcmp(obj.MaxNumInputSamplesSource,'Auto')
            obj.pMaxNumInputSamples = getPropagatedNumInputSamples(obj,x);
        else
            obj.pMaxNumInputSamples = double(obj.MaxNumInputSamples);
        end
        nfft = obj.pMaxNumInputSamples+coefflen-1;
        %nfft = 2^nextpow2(nfft);
        obj.pFFTLength = nfft;
        obj.cFFT = dsp.FFT('FFTLengthSource','Property','FFTLength',nfft);
        obj.cCoeffFFT = dsp.FFT('FFTLengthSource','Property','FFTLength',nfft);
        obj.cIFFT = dsp.IFFT;

        coder.extrinsic('phased.internal.FFTFilter.initWindowCoefficientsandLoss');
        initProps =...
               coder.internal.const(...
                   obj.initWindowCoefficientsandLoss ( ...
                       obj.SpectrumWindow, ...
                       obj.pFFTLength, ...
                       obj.SpectrumRange, ...
                       obj.SampleRate, ...
                       obj.Beta , ...
                       obj.SidelobeAttenuation, ...
                       obj.Nbar, ...
                       obj.CustomSpectrumWindow, ...
                       class(x)));
        
        obj.pWindowCoefficients = initProps.windowCoefficients;
        obj.pWindowProcessingLoss = initProps.windowProcessingLoss;
        
        obj.pOverlapState = ...
            complex(zeros(obj.pFFTLength,obj.pNumTotalChannels,class(x)));

        processInputSizeChangeImpl(obj,x);
        
        processTunedPropertiesImpl(obj);
    end
    
    function processInputSizeChangeImpl(obj,x,~)
        nfft = obj.pFFTLength;
        obj.pSignalLength = size(x,1);
        obj.pOverlapLength = nfft-obj.pSignalLength;
        obj.pIsSignalLonger = (obj.pSignalLength >= obj.pOverlapLength);
        obj.pSignalIdx = [1 obj.pSignalLength];
        obj.pStateIdx = [obj.pSignalLength+1 obj.pFFTLength];
        if obj.pIsSignalLonger
            obj.pSignalAddIdx = [1 obj.pOverlapLength];
            obj.pStateUpdateIdx = [1 1];  % not used
            obj.pStateRemainingIdx = [1 1]; % not used
        else
            obj.pSignalAddIdx = [1 obj.pSignalLength];
            obj.pStateUpdateIdx = [1 (obj.pOverlapLength-obj.pSignalLength)];
            obj.pStateRemainingIdx = [(obj.pSignalLength+1) obj.pOverlapLength];
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
        if obj.WindowLossOutputPort && (index == 2)
            flag = true;
        end
    end
    
    function processTunedPropertiesImpl(obj)
        if obj.pCoefficientsViaProp
            if isscalar(obj.Coefficients)
                obj.pCoefficientsFFT = obj.Coefficients;
            else
                if isreal(obj.Coefficients)
                    w_fft_in = complex(obj.Coefficients);
                else
                    w_fft_in = obj.Coefficients;
                end
                obj.pCoefficientsFFT = ...
                    step(obj.cCoeffFFT,w_fft_in) .* ...
                    obj.pWindowCoefficients;
            end
        end
    end

    function resetImpl(obj)
        obj.pOverlapState = zeros(size(obj.pOverlapState),'like',obj.pOverlapState);
        reset(obj.cFFT);
        reset(obj.cCoeffFFT);
        reset(obj.cIFFT);
        obj.pSizeInitialized = false;
    end
    
    function releaseImpl(obj)
        release(obj.cFFT);
        release(obj.cCoeffFFT);
        release(obj.cIFFT);
    end
    
    function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
        if index == 1
            flag = false;
        else
            flag = true;
        end
    end

    function fsz_out = isOutputFixedSizeImpl(obj) 
        fsz_out = propagatedInputFixedSize(obj, 1);
    end

    function sz_out = getOutputSizeImpl(obj)
        sz_out = propagatedInputSize(obj,1);
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
        if isLocked(obj)
            s.pCoefficientsViaProp = obj.pCoefficientsViaProp;
            s.pIsSignalLonger = obj.pIsSignalLonger;
            s.pFFTLength = obj.pFFTLength;
            s.pOverlapLength = obj.pOverlapLength;
            s.pSignalLength = obj.pSignalLength;
            s.pNumTotalChannels = obj.pNumTotalChannels;
            s.pSignalAddIdx = obj.pSignalAddIdx;
            s.pStateRemainingIdx = obj.pStateRemainingIdx;
            s.pStateUpdateIdx = obj.pStateUpdateIdx;
            s.pSignalIdx = obj.pSignalIdx;
            s.pStateIdx = obj.pStateIdx;
            s.pCoefficientsFFT = obj.pCoefficientsFFT;
            s.pOverlapState = obj.pOverlapState;
            s.pWindowCoefficients = obj.pWindowCoefficients;
            s.pWindowProcessingLoss = obj.pWindowProcessingLoss;
            s.pSizeInitialized = obj.pSizeInitialized;
            s.cFFT = saveobj(obj.cFFT);
            s.cIFFT = saveobj(obj.cIFFT);
            s.cCoeffFFT = saveobj(obj.cCoeffFFT);
        end
    end
        
    function s = loadSubObjects(obj,s,wasLocked)
        if wasLocked
            if isfield(s,'cFFT')
                obj.cFFT = dsp.FFT.loadobj(s.cFFT);
                s = rmfield(s,'cFFT');
            end
            if isfield(s,'cCoeffFFT')
                obj.cCoeffFFT = dsp.FFT.loadobj(s.cCoeffFFT);
                s = rmfield(s,'cCoeffFFT');
            end
            if isfield(s,'cIFFT')
                obj.cIFFT = dsp.IFFT.loadobj(s.cIFFT);
                s = rmfield(s,'cIFFT');
            end
        end
    end
        
    function loadObjectImpl(obj,s,wasLocked) 
        s = loadSubObjects(obj,s,wasLocked);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end
    
    function [sigout, wpl] = stepImpl(obj,sigin,w)
        if ~obj.pSizeInitialized
            processInputSizeChangeImpl(obj,sigin);
            obj.pSizeInitialized = true;
        end
        
        NFFT = obj.pFFTLength;
        if obj.pCoefficientsViaProp
            coeff_freq = cast(obj.pCoefficientsFFT,'like',sigin);
            if size(sigin,1) == 1
                y_in = bsxfun(@times,coeff_freq,complex(sigin)); % scalar fft
            else
                x_fft_in = zeros(NFFT,obj.pValidatedNumInputChannels,'like',sigin);
                x_fft_in(1:size(sigin,1),1:obj.pValidatedNumInputChannels) = sigin;
                y_in = bsxfun(@times,coeff_freq,step(obj.cFFT,complex(x_fft_in)));
            end
            if isreal(y_in)
                y_ifft_in = complex(y_in);
            else
                y_ifft_in = y_in;
            end
            if size(y_ifft_in,1) == 1
                fftout = y_ifft_in;
            else
                fftout = step(obj.cIFFT,y_ifft_in);
            end
        else
            win_freq = obj.pWindowCoefficients;
            w = cast(w,'like',sigin);
            
            if size(w,1) == 1
                w_fft_out = complex(w);
            else
                w_fft_in = zeros(NFFT,1,'like',w);
                w_fft_in(1:size(w,1)) = w;
                w_fft_out = step(obj.cCoeffFFT,complex(w_fft_in)).*win_freq;
            end
            
            if size(sigin,1) == 1
                x_fft_out = complex(sigin);
            else
                x_fft_in = zeros(NFFT,obj.pValidatedNumInputChannels,'like',sigin);
                x_fft_in(1:size(sigin,1),1:obj.pValidatedNumInputChannels) = sigin;
                x_fft_out = step(obj.cFFT,complex(x_fft_in));
            end
            
            y_in = bsxfun(@times,w_fft_out,x_fft_out);
            
            if size(y_in,1) == 1
                fftout = complex(y_in);
            else
                fftout = step(obj.cIFFT,complex(y_in));
            end
            
        end
        
        sigout_temp = fftout(obj.pSignalIdx(1):obj.pSignalIdx(2),:);
        if isreal(sigout_temp)
            sigout = complex(sigout_temp);
        else
            sigout = sigout_temp;
        end
        
        addIdx1 = obj.pSignalAddIdx(1);
        addIdx2 = obj.pSignalAddIdx(2);
        assert(addIdx2<=NFFT);
        assert(addIdx1>=1);
        addIdx = addIdx1:addIdx2;
        if ~isempty(addIdx)
            sigout(addIdx,:) = sigout(addIdx,:) + ...
                obj.pOverlapState(addIdx,:);
        end
        if obj.pIsSignalLonger
            obj.pOverlapState(1:obj.pOverlapLength,:) = fftout(obj.pStateIdx(1):obj.pStateIdx(2),:);
        else
            temp = obj.pOverlapState(obj.pStateRemainingIdx(1):obj.pStateRemainingIdx(2),:);
            obj.pOverlapState(1:obj.pOverlapLength,:) = fftout(obj.pStateIdx(1):obj.pStateIdx(2),:);
            updateIdx1 = obj.pStateUpdateIdx(1);
            updateIdx2 = obj.pStateUpdateIdx(2);
            assert(updateIdx2<=obj.pFFTLength);
            assert(updateIdx1>=1);
            updateIdx = updateIdx1:updateIdx2;
            obj.pOverlapState(updateIdx,:) = ...
                obj.pOverlapState(updateIdx,:)+temp;
        end
        if obj.WindowLossOutputPort
            wpl = obj.pWindowProcessingLoss;
        end
            
    end

    function flag = isInactivePropertyImpl(obj, prop)
        if strcmp(prop,'Coefficients') && ...
                strcmp(obj.CoefficientsSource,'Input port')
            flag = true;
        elseif strcmp(obj.MaxNumInputSamplesSource,'Auto') && ...
                strcmp(prop, 'MaximumInputSamples')
            flag = true;
        else
            flag = false;
        end
    end

end
methods(Static, Hidden)

    function initProps = initWindowCoefficientsandLoss ( ...
        spectrumWindow, ...
        FFTLength, ...
        spectrumRange, ...
        sampleRate, ...
        beta , ...
        sidelobeAttenuation, ...
        nbar, ...
        customSpectrumWindow, ...
        dt)
        
        if (spectrumWindow(1) == 'N') %None
            window_length = FFTLength;
            windowCoefficients = ones(window_length,1,dt);
        else
            windowCoefficients = zeros(FFTLength,1,dt);
            window_bw = diff(spectrumRange);
            if window_bw == 0
                window_length = 1;
            else
                window_length = max(round(FFTLength/...
                    (sampleRate/window_bw)),1);
            end

            switch spectrumWindow
              case 'Hamming'
                window_coeff = hamming(window_length,'symmetric');
              case 'Hann'
                window_coeff = hann(window_length,'symmetric');
              case 'Kaiser'
                window_coeff = kaiser(window_length,beta);
              case 'Chebyshev'
                window_coeff = chebwin(window_length,sidelobeAttenuation);
              case 'Taylor'
                window_coeff = taylorwin(window_length,nbar, ...
                                         -sidelobeAttenuation);
              case 'Custom'
                if isa(customSpectrumWindow,'function_handle')
                    window_cell = {customSpectrumWindow};
                else
                    window_cell = customSpectrumWindow;
                end
                window_coeff = window_cell{1}(window_length,...
                                              window_cell{2:end});
            end
            window_coeff = cast(window_coeff,dt);

            freq_grid = (0:FFTLength-1)*sampleRate/FFTLength;

                        
            if rem(FFTLength,2)
                freq_grid = fftshift(freq_grid);
                freq_grid(1:(FFTLength-1)/2) = ...
                    -freq_grid(end:-1:end-(FFTLength-1)/2+1);
            else
                freq_grid = fftshift(freq_grid);
                freq_grid(1:FFTLength/2) = ...
                    freq_grid(FFTLength/2+1:end) - sampleRate/2;
            end

            [~,idx] = min(abs(freq_grid-spectrumRange(1)));
            windowCoefficients(idx:idx+window_length-1) = ...
                window_coeff;
            windowCoefficients = ifftshift(windowCoefficients);
        end
        initProps.windowCoefficients = windowCoefficients;  
        initProps.windowProcessingLoss  = -pow2db(abs(...
            sum(windowCoefficients))^2/...
                     (window_length*...
                     (windowCoefficients'*windowCoefficients)));
    end

    
end
end


        
% [EOF]
