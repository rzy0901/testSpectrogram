classdef(Hidden) AbstractPulseWaveformLibrary < phased.internal.AbstractLibrary
%This class is for internal use only. It may be removed in the future.

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


methods (Access = protected)
    % Constructor    
    function obj = AbstractPulseWaveformLibrary(varargin)
        obj@phased.internal.AbstractLibrary(varargin{:});
    end
end

methods
    function mfcoeff = getMatchedFilter(obj,idx,pidx)
    %getMatchedFilter   Matched filter coefficients for the waveform
    %   COEFF = getMatchedFilter(obj,IDX) returns the matched filter
    %   coefficients for the IDX-th waveform. If the waveform consists of a
    %   single pulse, such as a linear FM waveform, COEFF is a column
    %   vector. If the waveform consists of multiple diverse pulses, such
    %   as a stepped FM waveform, COEFF is a matrix whose columns
    %   correspond to the different pulses.
    %
    %   COEFF = getMatchedFiler(obj,IDX,PIDX) specifies the pulse index,
    %   PIDX, as a positive integer. This syntax is only applicable when
    %   the IDX-th waveform in the waveform library is not a custom
    %   waveform.
    %
    %   % Example:
    %   %   Calculate the matched filter coefficients for pulse waveform
    %   %   library and plot the coefficients.
    %
    %   waveform1 = {'Rectangular','PRF',1e4, 'PulseWidth', 50e-6};
    %   waveform2 = {'LinearFM','PRF',1e4,'PulseWidth',50e-6,...
    %   'SweepBandwidth',1e5,'SweepDirection','Up',...
    %   'SweepInterval', 'Positive'};
    %
    %   H = phased.PulseWaveformLibrary('SampleRate',1e6,...
    %               'WaveformSpecification',{waveform1,waveform2});
    %   coeff1 = getMatchedFilter(H,1,1);
    %   subplot(2,1,1);
    %   stem(real(coeff1));title('Matched filter coefficients, real part');
    %
    %   coeff2 = getMatchedFilter(H,2,1);
    %   subplot(2,1,2);
    %   stem(real(coeff2));title('Matched filter coefficients, real part'); 
        coder.extrinsic('phased.PulseWaveformLibrary.getFreqOffset');
        narginchk(2,3);

        N = numel(obj.WaveformSpecification);
        sigdatatypes.validateIndex(idx,'getMatchedFilter','IDX',{'scalar','<=',N});
            
        if isempty(coder.target)
            wasLocked = isLocked(obj);
            if ~wasLocked
                setup(obj,idx);
            end
        else
            if ~coder.internal.is_defined(obj.cWaveformSpecification)
                setup(obj,idx);
            end
        end
        
        mfcoeffs = cell(1,N);
        FreqOffset = zeros(1,N);
        for m = coder.unroll(1:N)
            if ~isa(obj.WaveformSpecification{m}{1},'function_handle')
                if nargin < 3
                    mfcoeffs{m} = getMatchedFilter(...
                        obj.cWaveformSpecification{m});
                else
                    sigdatatypes.validateIndex(pidx,'getMatchedFilter','PIDX',{'scalar'});
                    mfcoeffs{m} = getMatchedFilter(...
                        obj.cWaveformSpecification{m},pidx);
                end
%                 FreqOffset(m) = obj.pFrequencyOffset(m);
                value = obj.WaveformSpecification;
                FreqOffset(m) = coder.const(phased.PulseWaveformLibrary.getFreqOffset(value{m}));
            else
                narginchk(2,2);
                wav = getCustomWaveformOutput(obj,m);
                z = nonzeros(abs(getCustomWaveformOutput(obj,m)));
                nz = numel(z);
                wavout = wav(1:nz(1));
                if isempty(nz)
                    mfcoeffs{m} = 0;
                else
                    mfcoeffs{m} = flipud(conj(wavout));
                end
                FreqOffset(m) = 0;
            end
        end
        
       assert(idx<=N); 
        
        for i = coder.unroll(1:N)
            if i == idx
                deltaf = FreqOffset(i);
                t = (0:size(mfcoeffs{i},1)-1)'/obj.SampleRate;
                mfcoeff = mfcoeffs{i}.*exp(1j*2*pi*deltaf*t);
                break;
            else
                mfcoeff = zeros(size(mfcoeffs{i}));
            end
        end
        mfcoeff = mfcoeff(:);
    end
    
    function varargout = plot(obj,idx,varargin)
    %plot  Plot the waveform
    %   plot(Hwav,IDX) plots the real part of the IDX-th waveform specified
    %   by Hwav. Hwav must be a single waveform object. Note that the plot
    %   does not include silent samples in each pulse.
    %
    %   plot(Hwav,IDX,'PlotType',Type) specifies the type of plot using one
    %   of the following strings for Type: [ real | imag | complex ].
    %
    %   plot(...,'PulseIdx',PID) specifies the index of the pulse to plot.
    %   PID must be a scalar and its default value is 1.
    %
    %   plot(...,LineSpec) specifies the line color, style, or marker as in
    %   plot function in MATLAB. For the case when both real and imaginary
    %   plots are specified, the LineSpec applies to both subplots.
    %
    %   h = plot(...) returns the line handle in the figure. For the case
    %   when both real and imaginary plots are specified, the vector of
    %   handles h will include handles to the lines in both subplots, in
    %   the form of [RealLineHandle; ImagLineHandle].
    %
    %   % Example:
    %   %   Create and plot a linear FM pulse waveform.
    %
    %   waveform1 = {'Rectangular','PRF',1e4, 'PulseWidth', 50e-6};
    %   waveform2 = {'LinearFM','PRF',1e4,'PulseWidth',50e-6,...
    %                'SweepBandwidth',1e5,'SweepDirection','Up',...
    %                'SweepInterval', 'Positive'};
    %   wavlib = phased.PulseWaveformLibrary('SampleRate',1e6,...
    %               'WaveformSpecification',{waveform1,waveform2});
    %   plot(wavlib,1);
    
        narginchk(2,inf)
        validateattributes(idx,{'double'},...
            {'scalar','positive','finite' },'','IDX');
        
        if ~isLocked(obj)
            setup(obj,idx);
        end
        
        if ~isa(obj.cWaveformSpecification{idx},'function_handle')
            
            h = plot(obj.cWaveformSpecification{idx},varargin{:});
            
        else
            
            if isempty (coder.target) %running MATLAB
                % Parse argument list to separate waveform specific arguments.
                wfArg = {};     % PV pairs for the waveform
                arglist = {'PlotType','PulseIdx'};
                if nargin > 2 && ismember(varargin{1},arglist)
                    if nargin > 4 && ismember(varargin{3},arglist)
                        wfArg = varargin(1:4);
                        lineArg = varargin(5:end);
                    else
                        wfArg = varargin(1:2);
                        lineArg = varargin(3:end);
                    end
                else
                    lineArg = varargin;
                end
                
                % Validate waveform P-V pairs
                PlotType = 'real';
                PulseIdx = 1;
                sigutils.pvparse(wfArg{:});
                PlotType = validatestring(PlotType,{'real','imag',...
                    'complex'},'plot','PlotType');
                validateattributes(PulseIdx,{'double'},{'scalar',...
                    'integer','positive'},'plot','PulseIdx');
                
                wav = getCustomWaveformOutput(obj,idx);
                nz = find(abs(wav),1,'last');
                x = wav(1:nz); 
                t = (0:length(x)-1)/obj.SampleRate;
                name = func2str(obj.WaveformSpecification{idx}{1}); % Get custom waveform name
                
                strxlbl = 'Time (s)'; strylbl = 'Amplitude (v)';
                if strcmp(PlotType,'complex')
                    ha1 = subplot(2,1,1);
                    h1 = plot(t,real(x),lineArg{:}); grid on;
                    set(ha1,'Tag','Real_Part');
                    xlabel(strxlbl), ylabel(strylbl);
                    title(sprintf('%s: real part, pulse %d',name,...
                        PulseIdx));
                    ha2 = subplot(2,1,2);
                    h2 = plot(t,imag(x),lineArg{:}); grid on;
                    set(ha2,'Tag','Imaginary_Part');
                    xlabel(strxlbl), ylabel(strylbl);
                    title(sprintf('%s: imaginary part, pulse %d',name,...
                        PulseIdx));
                    h = [h1; h2];
                else
                    if strcmp(PlotType,'real')
                        h = plot(t,real(x),lineArg{:}); grid on;
                        set(gca,'Tag','Real_Part');
                    else
                        h = plot(t,imag(x),lineArg{:}); grid on;
                        PlotType = 'imaginary';
                        set(gca,'Tag','Imaginary_Part');
                    end
                    xlabel(strxlbl), ylabel(strylbl);
                    title(sprintf('%s: %s part, pulse %d',name,PlotType,...
                        PulseIdx));
                end
            else
                coder.internal.assert(false, ...
                    'phased:Waveform:CodegenNotSupported','plot');
            end
        end
        if nargout > 0
            varargout = {h};
        end
        
    end
end

methods(Access = protected)
    
    function flag = isInDVRMode(obj) %#ok<MANU>
        flag = false;
    end

    function setupImpl(obj,~)
        % Perform one-time calculations, such as computing constants
        setupImpl@phased.internal.AbstractLibrary(obj);
    end
    
    function num = getNumInputsImpl(~)
        num = 1; % Idx
    end
    
    function num = getNumOutputsImpl(~)
        num = 1; % Data out
    end
    
    function validateInputsImpl(~,idx)
        
        cond = ~isa(idx,'double');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','IDX','double');
        end
        cond = ~isscalar(idx);
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeScalar','IDX');
        end
        
    end
    
    function idx = validateInputIdx(obj,value)
        idx = sigdatatypes.validateIndex(value,'','IDX',{'<=',...
            numel(obj.WaveformSpecification)});
    end
    

    function y = stepImpl(obj,idx)
        
        currentIdx = idx;
        cwavspec = obj.cWaveformSpecification;
        N = numel(obj.WaveformSpecification);
        y = 0;
        
        idx = validateInputIdx(obj,idx);
        assert(idx<=N);
        for m = coder.unroll(1:N)
            if m == idx
                deltaf = obj.pFrequencyOffset(m);
                if obj.pIsFunctionHandle(m)
                    y = getCustomWaveformOutput(obj,m);
                else
                    if obj.pPreviousIndex ~= m
                        reset(cwavspec{m});
                    end
                    y = step(cwavspec{m});
                    t = (0:size(y,1)-1)'/obj.SampleRate;
                    y = y.*exp(1j*2*pi*deltaf*t);
                end
                break;
            end
        end
        
        obj.pPreviousIndex = currentIdx;

        if isInDVRMode(obj)
            setNumTicksUntilNextHit(obj,size(y,1));
        end
        
    end
    
    function resetImpl(obj)
       resetImpl@phased.internal.AbstractLibrary(obj);
    end
     
    function releaseImpl(obj)
      releaseImpl@phased.internal.AbstractLibrary(obj);
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractLibrary(obj);
    end
    
    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractLibrary(obj,s);
    end
     
    function loadObjectImpl(obj,s,wasLocked)
       loadObjectImpl@phased.internal.AbstractLibrary(obj,s,wasLocked);
    end
end

methods (Access = protected)
    
    function flag = isInputComplexityLockedImpl(obj,index) %#ok<INUSD>
        flag = true;
    end
    
    function flag = isInputSizeLockedImpl(obj,index)  %#ok<INUSD>
        % Return true if input size is not allowed to change while
        % system is running
        flag = true;
    end
    
    function varargout = getOutputSizeImpl(obj)

        PRF = populatePRF(obj); % Function to obtain the PRFs
        prf = min(PRF); % Minimum PRF for setting maximum output dimension
        pulseLength = round(obj.SampleRate/prf);
        varargout{1} = pulseLength*1; % 1 to signify only one pulse is simulated.
        
    end
    
    function varargout = getOutputNamesImpl(~)
        varargout = {'Y'};
    end
    
    function varargout = getInputNamesImpl(~)
        varargout = {'Idx'};
    end
    
    function varargout = isOutputFixedSizeImpl(~)
        varargout{1} = false;
    end
    
    function varargout = getOutputDataTypeImpl(~)
        varargout{1} = 'double';
    end
    
    function varargout = isOutputComplexImpl(~)
        varargout{1} = true;
    end
end

end



