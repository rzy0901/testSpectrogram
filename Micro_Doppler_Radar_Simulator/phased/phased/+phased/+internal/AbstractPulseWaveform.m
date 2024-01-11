classdef (Hidden) AbstractPulseWaveform < phased.internal.AbstractSampleRateEngine  & ...
     matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime 
%This class is for internal use only. It may be removed in the future.

%ABSTRACTPULSEWAVEFORM Define the ABSTRACTPULSEWAVEFORM class
% This is an abstract class in support of pulse waveform functionality.

%   Copyright 2010-2017 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %SampleRate Sample rate (Hz)
    %   Specify the sample rate (in Hz) as a positive scalar. The default
    %   value of this property is 1e6 (1 MHz).
    SampleRate = 1e6;
    %PRF    Pulse repetition frequency (Hz)
    %   Specify the pulse repetition frequency (in Hz) as a positive scalar
    %   or a row vector. The default value of this property is 1e4 (10
    %   kHz). When PRF is a vector, it represents the case of staggered PRF
    %   where the output pulses use elements in the vector as their PRFs
    %   one after another in cycle.
    PRF = 1e4
end

properties (Nontunable)
    %OutputFormat     Output signal format
    %   Specify the format of the output signal as one of 'Pulses' |
    %   'Samples', where the default is 'Pulses'. When you set the
    %   OutputFormat property to 'Pulses', the output is in the form of
    %   multiple pulses where the number of pulses is determined by the
    %   value of the NumPulses property. When you set the OutputFormat
    %   property to 'Samples', the output is in the form of multiple
    %   samples where the number of samples is determined by the value of
    %   the NumSamples property.
    OutputFormat = 'Pulses'
end

properties (Nontunable, PositiveInteger) 
    %NumSamples     Number of samples in output
    %   Specify the number of samples in each output as a positive integer.
    %   This property only applies when you set the OutputFormat property
    %   to 'Samples'. The default value of this property is 100.
    NumSamples = 100
    %NumPulses  Number of pulses in output
    %   Specify the number of pulses in each output as a positive integer.
    %   This property only applies when you set the OutputFormat property
    %   to 'Pulses'. The default value of this property is 1.
    NumPulses = 1;
end

properties (Nontunable, Logical)
    %PRFSelectionInputPort  Enable PRF selection input
    %   Set this property to true to select which predefined PRF to use
    %   during the simulation via input. Set this property to false to use
    %   the PRF property to define the PRF sequence used in the simulation.
    %   The default value of this property is false.
    PRFSelectionInputPort = false
end

properties(Nontunable, Logical)
    %PRFOutputPort  Enable PRF Output
    %   Set this property to true if you wish to display the PRF used
    %   during the simulation as an output. The default value of this
    %   property is false. This property can be used only when the
    %   OutputFormat property is 'Pulses'.
    PRFOutputPort = false
end

properties(Constant, Hidden)
    OutputFormatSet = matlab.system.StringSet({'Pulses','Samples'});
end

properties (Access = protected)
    %first output pulse index
    pOutputStartPulseIndex;
    %indices for output pulses
    pOutputPulseInterval;
    %first output sample index
    pOutputStartSampleIndex;
    %indices for output samples
    pOutputSampleInterval;
    %internal buffer of waveform samples
    pSamples;
end

properties (Access = protected, Nontunable)
    %private property to hold pulse length for each different PRF
    pPulseLength;
    %private property to hold number of samples for waveform
    pNonZeroLength;
    %private property indicating whether to output by pulse
    pOutputByPulse;
    %number of distinct PRF
    pNumDistinctPRF;
    %number of distinct pulses, when there are multiple waveforms and
    %multiple PRFs, this is the least common multiplier of the two
    pNumDistinctPulses;
    %number of distinct samples, total number of samples in distinct pulses
    pNumDistinctSamples;
    %start sample of each distinct pulse
    pDistinctNonZeroStart;
    %end sample of each distinct pulse, non-zero portion
    pDistinctNonZeroEnd;
    %waveform index for each distinct pulse, non-zero portion
    pWaveformIndexForDistinctPulses;
    %number of samples of each distinct pulse
    pDistinctPulseLength
    %size upper bound
    pOutputUpperBound
end

properties (Access = protected, Dependent)
    %pulse width for the pulse, different between phase coded waveforms and
    %other pulse waveforms where the phase is continuous within the pulse
    pPulseWidth
end

methods
    function set.SampleRate(obj, value)
        validateattributes(value,{'double'}, {'scalar',...
            'positive','finite'},...
            '','SampleRate');
        obj.SampleRate = value;
    end
    function set.PRF(obj,value)
        validateattributes( value, { 'double' }, { 'row', 'positive', 'finite' }, '', 'PRF');
        obj.PRF = value;
    end
    function pw = get.pPulseWidth(obj)
        pw = getPulseWidth(obj);
    end
end

methods (Access = protected)
    function obj = AbstractPulseWaveform(varargin)
        setProperties(obj, nargin, varargin{:});
    end
end

methods (Abstract)
    %bandwidth Return bandwidth
    %   BW = bandwidth(Hwav) returns the bandwidth of the waveform Hwav.
    bw = bandwidth(obj)
end

methods (Abstract, Access = protected)
    pidx = getOutputPulseIndex(obj,prfidx)
    sidx = getOutputSampleIndex(obj,prfidx)
end

methods (Access=protected,Abstract)
    %getWaveformName  Get name of the waveform
    wname = getWaveformName(obj)
    %calcDistinctPulseParameters Calculate distinct pulse parameters for
    %each waveform
    calcDistinctPulseParameters(obj);
    %getMatchingWaveform Return the matching waveform
    s = getMatchingWaveform(obj,pidx);
    %getPulseWidth Returns the pulse width
    pw = getPulseWidth(obj);
    %getDefaultMatchedFilterPulseIndex Return the default pulse index for
    %matched filter
    idx = getDefaultMatchedFilterPulseIndex(obj);
end

methods (Access = protected)
    
    function num = getNumInputsImpl(obj) 
        num = 0;
        if obj.PRFSelectionInputPort
            num = num+1;
        end
    end
    
    function num = getNumOutputsImpl(obj)
        num = 1;
        if obj.PRFOutputPort && (obj.OutputFormat(1) == 'P')
            num = num+1;
        end
    end
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        pulseflag = (obj.OutputFormat(1) == 'P'); %Pulses
        if strcmp(prop,'NumSamples') && pulseflag
            flag = true;
        elseif strcmp(prop,'NumPulses') && ~pulseflag
            flag = true;
        elseif strcmp(prop,'PRFOutputPort') && ~pulseflag
            flag = true;
        end
    end
    
    function validatePropertiesImpl(obj)
        q = obj.SampleRate./obj.PRF;
        cond =  any(abs(q-round(q))>eps(q));
        if cond
            coder.internal.errorIf(cond, ...
                'phased:Waveform:NeedRatioInteger', 'SampleRate', 'PRF');
        end
        
    end
    
    function validateInputsImpl(obj,idx) 
        if obj.PRFSelectionInputPort
            validateattributes(idx,{'double'},{'scalar','integer'},...
                'step','PRFIDX');
        end
    end
    
    function setupImpl(obj,~)
        coder.extrinsic('phased.internal.AbstractPulseWaveform.getPulseOutputSize');
        obj.pNumDistinctPRF = numel(obj.PRF);
        obj.pPulseLength = round(obj.SampleRate./obj.PRF);
       
        obj.pOutputByPulse = (obj.OutputFormat(1) == 'P'); %Pulses
        
        calcDistinctPulseParameters(obj);

        obj.pDistinctNonZeroStart = cumsum(obj.pDistinctPulseLength)-...
            obj.pDistinctPulseLength+1;
        obj.pDistinctNonZeroEnd = obj.pDistinctNonZeroStart+...
            obj.pNonZeroLength-1;
        
        obj.pNumDistinctSamples = sum(obj.pDistinctPulseLength);
        obj.pOutputUpperBound = coder.const(obj.getPulseOutputSize(...
            obj.SampleRate, obj.PRF,obj.OutputFormat, ...
            obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort));

    end
    
    function flag = isOutputComplexityLockedImpl(obj,~) %#ok<INUSD>
        flag = false;  % index == 1 || index == 2
    end
    
    function flag = isInDVRMode(obj) %#ok<MANU>
        flag = false;
    end

    function resetImpl(obj)
        if strcmp(obj.OutputFormat,'Pulses') 
            obj.pOutputStartPulseIndex = 1;
            obj.pOutputPulseInterval = 0:obj.NumPulses;
            if isInDVRMode(obj)
                setNumTicksUntilNextHit(obj,1);
            end
        else
            if obj.PRFSelectionInputPort
                obj.pOutputStartSampleIndex = 0; 
            else
                obj.pOutputStartSampleIndex = 1;
            end
            obj.pOutputSampleInterval = 0:obj.NumSamples;
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        s.isLocked = isLocked(obj);
        if isLocked(obj)
            s.pOutputStartPulseIndex = obj.pOutputStartPulseIndex;
            s.pOutputPulseInterval = obj.pOutputPulseInterval;
            s.pOutputStartSampleIndex = obj.pOutputStartSampleIndex;
            s.pOutputSampleInterval = obj.pOutputSampleInterval;
            s.pSamples = obj.pSamples;
            s.pPulseLength = obj.pPulseLength;
            s.pNonZeroLength = obj.pNonZeroLength;
            s.pOutputByPulse = obj.pOutputByPulse;
            s.pNumDistinctPRF = obj.pNumDistinctPRF;
            s.pNumDistinctSamples = obj.pNumDistinctSamples;
            s.pNumDistinctPulses = obj.pNumDistinctPulses;
            s.pDistinctNonZeroStart = obj.pDistinctNonZeroStart;
            s.pDistinctNonZeroEnd = obj.pDistinctNonZeroEnd;
            s.pWaveformIndexForDistinctPulses = obj.pWaveformIndexForDistinctPulses;
            s.pDistinctPulseLength = obj.pDistinctPulseLength;
            s.pOutputUpperBound = obj.pOutputUpperBound;      
        end
    end
        
    function s = loadSubObjects(obj,s)  %#ok<INUSL>
        if isfield(s,'isLocked')
            s = rmfield(s,'isLocked');
        end
    end
    
    function [y,currentPRF] = stepImpl(obj,prfidx)
        if obj.pOutputByPulse
            if obj.PRFSelectionInputPort
                % validation
                cond = (rem(prfidx,1)~=0) || (prfidx>obj.pNumDistinctPRF);
                if cond
                    coder.internal.errorIf(cond,'phased:step:expectedInteger',...
                        'PRFIdx');
                end
                % determine the output pulses
                OutputPulseIndex = getOutputPulseIndex(obj,prfidx);
                % determine output pulse lengths
                OutputPulseLength = obj.pDistinctPulseLength(OutputPulseIndex(1:end));
                % preallocate, make all zeros because most of time duty cycle
                % is small
                if (obj.pNumDistinctPRF == 1)
                    numely = obj.pDistinctPulseLength(1)*obj.NumPulses;
                else
                    numely = sum(OutputPulseLength);
                    %upper bound of y
                    assert(numely <= obj.pOutputUpperBound);
                    if isInDVRMode(obj)
                        setNumTicksUntilNextHit(obj,numely);
                    end
                end
                y = complex(zeros(numely,1));
                % indices for all nonzero samples
                NonZeroIndex = [ones(1,obj.NumPulses); obj.pNonZeroLength(OutputPulseIndex(1:end))]+...
                    ones(2,1)*(cumsum(OutputPulseLength)-OutputPulseLength);
                % assign waveform samples for each pulse
                for m = 1:obj.NumPulses
                    y(NonZeroIndex(1,m):NonZeroIndex(2,m)) = getMatchingWaveform(obj,...
                        obj.pWaveformIndexForDistinctPulses(OutputPulseIndex(m)));
                end 
                if obj.PRFOutputPort
                    currentPRF = obj.PRF(prfidx);
                end
            else
                % determine the output pulses
                OutputPulseIndex = getOutputPulseIndex(obj);
                % determine output pulse lengths
                OutputPulseLength = obj.pDistinctPulseLength(OutputPulseIndex(1:end-1));
                % preallocate, make all zeros because most of time duty cycle
                % is small
                if (obj.pNumDistinctPRF == 1)
                    numely = obj.pDistinctPulseLength(1)*obj.NumPulses;
                else
                    numely = sum(OutputPulseLength);
                    %upper bound of y
                    assert(numely <= obj.pOutputUpperBound);
                    if isInDVRMode(obj)
                        setNumTicksUntilNextHit(obj,numely);
                    end
                end
                
                y = complex(zeros(numely,1));
                % indices for all nonzero samples
                NonZeroIndex = [ones(1,obj.NumPulses); obj.pNonZeroLength(OutputPulseIndex(1:end-1))]+...
                    ones(2,1)*(cumsum(OutputPulseLength)-OutputPulseLength);
                % assign waveform samples for each pulse
                for m = 1:obj.NumPulses
                    y(NonZeroIndex(1,m):NonZeroIndex(2,m)) = getMatchingWaveform(obj,...
                        obj.pWaveformIndexForDistinctPulses(OutputPulseIndex(m)));
                end
                if obj.PRFOutputPort
                    idx = mod(OutputPulseIndex(end-1),numel(obj.PRF))+1;
                    currentPRF = obj.PRF(idx);
                end
            end
        else  % output by samples
            if obj.PRFSelectionInputPort
                % validation
                cond = (rem(prfidx,1)~=0) || (prfidx>obj.pNumDistinctPRF);
                if cond
                    coder.internal.errorIf(cond,'phased:step:expectedInteger',...
                        'PRFIdx');
                end
                % determine the output sample indices
                OutputSampleIndex = getOutputSampleIndex(obj,prfidx);
                % preallocate, make all zeros because most of time duty cycle
                % is small
                y = complex(zeros(obj.NumSamples,1));
                % assign waveform samples for each distinct pulse
                for m = 1:obj.pNumDistinctPulses
                    currentPulseIdx = find(...
                        (OutputSampleIndex(1:end) >= obj.pDistinctNonZeroStart(m)) & ...
                        (OutputSampleIndex(1:end) <= obj.pDistinctNonZeroEnd(m)));
                    if ~isempty(currentPulseIdx)
                        currentPulse = getMatchingWaveform(obj,...
                            obj.pWaveformIndexForDistinctPulses(m));
                        sampleIdx = OutputSampleIndex(currentPulseIdx)-...
                            obj.pDistinctNonZeroStart(m)+1;
                        y(currentPulseIdx) = currentPulse(sampleIdx);
                    end
                end
            else
                % determine the output sample indices
                OutputSampleIndex = getOutputSampleIndex(obj);
                % preallocate, make all zeros because most of time duty cycle
                % is small
                y = complex(zeros(obj.NumSamples,1));
                % assign waveform samples for each distinct pulse
                for m = 1:obj.pNumDistinctPulses
                    currentPulseIdx = find(...
                        (OutputSampleIndex(1:end-1) >= obj.pDistinctNonZeroStart(m)) & ...
                        (OutputSampleIndex(1:end-1) <= obj.pDistinctNonZeroEnd(m)));
                    if ~isempty(currentPulseIdx)
                        currentPulse = getMatchingWaveform(obj,...
                            obj.pWaveformIndexForDistinctPulses(m));
                        sampleIdx = OutputSampleIndex(currentPulseIdx)-...
                            obj.pDistinctNonZeroStart(m)+1;
                        y(currentPulseIdx) = currentPulse(sampleIdx);
                    end
                end
            end
        end
    end
    
end

methods (Static, Hidden)
    
    function idx = getCircularIndex(idx,N)
        idx = mod(idx-1,N)+1;
    end
    function retSz = getPulseOutputSize(sampleRate,PRF,outputFormat,...
            numPulses,numSamples,flag)
          %num of samples for each pulse in PRF vector
          pulseLength = round(sampleRate./PRF(:));
          if strcmp(outputFormat,'Pulses')
              numPRF = numel(PRF);
              if numPRF == 1
                  retSz = pulseLength*numPulses;
              elseif flag
                  retSz = max(pulseLength)*numPulses;
              else
                  %Get total number of samples of all pulse combinations
                  staggeredIndex = bsxfun(@plus,(1:numPRF)'-1,1:numPulses);
                  staggeredIndex = ...
                     phased.internal.AbstractPulseWaveform.getCircularIndex(staggeredIndex,numPRF);
                  currentSizes = sum(pulseLength(staggeredIndex),2);
                  %Get the upperbound
                  retSz = max(currentSizes);                  
              end
          else
              retSz = numSamples;
          end
    end
end

methods
    
    function mfcoeff = getMatchedFilter(obj,pidx)
    %getMatchedFilter   Matched filter coefficients for the waveform
    %   COEFF = getMatchedFilter(obj) returns the matched filter
    %   coefficients for the waveform. If the waveform consists of a single
    %   pulse, such as a linear FM waveform, COEFF is a column vector. If
    %   the waveform consists of multiple diverse pulses, such as a stepped
    %   FM waveform, COEFF is a matrix whose columns correspond to the
    %   different pulses.
    %
    %   COEFF = getMatchedFiler(obj.PIDX) specifies the pulse index, PIDX,
    %   as a positive integer. The plot assumes the output pulses to loop
    %   through the entries specified in the PRF property so the pulse
    %   index is determined accordingly.
    %
    %   % Examples:
    %   
    %   % Example 1:
    %   %   Calculate the matched filter coefficients for a linear FM
    %   %   waveform.
    %
    %   hwav = phased.LinearFMWaveform;
    %   coeff = getMatchedFilter(hwav);
    %   stem(real(coeff)); title('Matched filter coefficients, real part');
    %   
    %   % Example 2:
    %   %   Calculate the matched filter coefficients for a stepped FM
    %   %   waveform.
    %
    %   hwav = phased.SteppedFMWaveform;
    %   coeff = getMatchedFilter(hwav);
    %   stem(real(coeff(:,3))); 
    %   title('Matched filter coefficients for the 3rd pulse, real part');
        validateProperties(obj);
        if nargin < 2
            pidx = getDefaultMatchedFilterPulseIndex(obj);
        else
            sigdatatypes.validateIndex(pidx,'getMatchedFilter','PIDX',...
                {'scalar'});
        end
        mfcoeff = flipud(conj(getMatchingWaveform(obj,pidx)));
    end
    
    function varargout = plot(obj,varargin)
    %plot   Plot the waveform
    %   plot(Hwav) plots the real part of the pulse waveform specified by
    %   Hwav. Hwav must be a single waveform object. Note that the plot
    %   does not include silent samples in each pulse.
    %
    %   plot(Hwav,'PlotType',Type) specifies the type of plot using one of
    %   the following strings for Type: [ {real} | imag | complex ].
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
    %   hw = phased.LinearFMWaveform;
    %   plot(hw);
    if isempty (coder.target) %running MATLAB
        narginchk(1,inf);
        validateattributes(obj,{'phased.internal.AbstractPulseWaveform'},...
            {'scalar'},'waveform','plot');
        
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
        PlotType = validatestring(PlotType,{'real','imag','complex'},...
            'plot','PlotType');
        validateattributes(PulseIdx,{'double'},{'scalar','integer',...
            'positive'},'plot','PulseIdx');
        
        % Generate samples
        wasLocked = isLocked(obj);
        if ~wasLocked
            setup(obj);
        end
        pidx = obj.getCircularIndex(PulseIdx,obj.pNumDistinctPulses);
        x = getMatchingWaveform(obj,pidx);
        t = (0:length(x)-1)/obj.SampleRate;
        
        strxlbl = 'Time (s)'; strylbl = 'Amplitude (v)';
        if strcmp(PlotType,'complex')
            ha1 = subplot(2,1,1);
            h1 = plot(t,real(x),lineArg{:}); grid on;
            set(ha1,'Tag','Real_Part');
            xlabel(strxlbl), ylabel(strylbl);
            title(sprintf('%s: real part, pulse %d', getWaveformName(obj), PulseIdx));
            ha2 = subplot(2,1,2);
            h2 = plot(t,imag(x),lineArg{:}); grid on;
            set(ha2,'Tag','Imaginary_Part');
            xlabel(strxlbl), ylabel(strylbl);
            title(sprintf('%s: imaginary part, pulse %d', getWaveformName(obj), PulseIdx));
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
            title(sprintf('%s: %s part, pulse %d', getWaveformName(obj), PlotType, PulseIdx));
        end
        
        if nargout > 0
            varargout = {h};
        end
        if ~wasLocked
            release(obj);
        end
    else
        coder.internal.assert(false, ...
            'phased:Waveform:CodegenNotSupported','plot');
    end
    end
    
end
methods (Access = protected) %For Simulink propagation and mask
    function varargout = getOutputNamesImpl(obj)
        varargout = {'Y'};
        if ((strcmp(obj.OutputFormat,'Pulses'))&&(obj.PRFOutputPort))
                varargout{2} = 'PRF';
        end
    end
    
    function varargout = getInputNamesImpl(obj)
        if obj.PRFSelectionInputPort
            varargout = {'PRFIdx'};
        else
            varargout = {};
        end
    end
    
    function flag = isInputSizeLockedImpl(obj,~)  %#ok<INUSD>
        flag = true;
    end
        
    function varargout = getOutputSizeImpl(obj)
        varargout{1} =  obj.getPulseOutputSize(...
            obj.SampleRate, obj.PRF,obj.OutputFormat, ...
            obj.NumPulses,obj.NumSamples,obj.PRFSelectionInputPort);
        if (obj.PRFOutputPort && strcmp(obj.OutputFormat,'Pulses'))
            varargout{2} = 1;
        end
    end
    function varargout = isOutputFixedSizeImpl(obj)
        if strcmp(obj.OutputFormat,'Pulses') && ...
                numel(obj.PRF) ~= 1
            %error(message('phased:Waveform:StaggeredNotSupported','PRF','OutputFormat','Pulses'));
            varargout{1} = false;
        else
            varargout{1} = true;
        end
        if obj.PRFOutputPort
            varargout{2} = true;
        end
    end
    function varargout = getOutputDataTypeImpl(obj)  
        varargout{1} = 'double';
        if obj.PRFOutputPort
            varargout{2} = 'double';
        end
    end
    function varargout = isOutputComplexImpl(obj)
        varargout{1} = true;
        if obj.PRFOutputPort
            varargout{2} = false;
        end
    end
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

