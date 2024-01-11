classdef (Sealed,StrictDefaults) MFSKWaveform < matlab.System & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.SampleTime
%MFSKWaveform   MFSK waveform
%   H = phased.MFSKWaveform creates an MFSK waveform System object, H. This
%   object generates samples of an MFSK waveform. MFSK waveform is a
%   combination of frequency shift keying (FSK) and linear frequency
%   modulated continuous wave (LFMCW) waveforms.
%
%   H = phased.MFSKWaveform(Name,Value) creates an MFSK waveform object, H,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H) returns samples of the MFSK waveform in a column vector Y.
%   Y can contain either a certain number of sweeps or a certain number of
%   samples.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H) and y = H() are
%   equivalent.
%
%   MFSKWaveform methods:
%
%   step                - Return samples of the MFSK waveform
%   release             - Allow property value and input characteristics 
%                         changes
%   clone               - Create an MFSK waveform object with same property
%                         values
%   isLocked            - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>               - Reset states of MFSK waveform object
%   plot                - Plot the MFSK waveform
%
%   MFSKWaveform properties:
%
%   SampleRate      - Sample rate 
%   SweepBandwidth  - Sweep bandwidth 
%   StepTime        - Frequency step burst time
%   StepsPerSweep   - Number of frequency steps in a sweep
%   FrequencyOffset - Chirp offset frequency
%   OutputFormat    - Output signal format
%   NumSamples      - Number of samples in output
%   NumSteps        - Number of frequency steps in output
%   NumSweeps       - Number of sweeps in output
%
%   % Examples:
%
%   % Example 1:
%   %   Create and plot an MFSK waveform.
%
%   waveform = phased.MFSKWaveform('SweepBandwidth',1e5,...
%               'OutputFormat','Steps','NumSteps',2);
%   plot(waveform);
%
%   See also phased, phased.FMCWWaveform, phased.LinearFMWaveform.
    
%   Copyright 2014-2017 The MathWorks, Inc.

%   Reference
%   [1] Marc-Michale Meinecke and Hermann Rohling, Combination of LFMCW and
%       FSK Modulation Principles for Automotive Radar Systems, German
%       Radar Symposium GRS2000, 2000
%   [2] Hermann Rohling and Marc-Michale Meinecke, Waveform Design
%       Principles for Automotive Radar Systems, CIE International
%       Conference on Radar, 2001


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %SampleRate Sample rate (Hz)
    %   Specify the sample rate (in Hz) as a positive scalar. The default
    %   value of this property is 1e6 (1 MHz).
    SampleRate = 1e6;
    %SweepBandwidth     Sweep bandwidth (Hz)
    %   Specify the bandwidth of the MFSK sweeping (in Hz) as a positive
    %   scalar. The default value is 1e5 (100 kHz).
    SweepBandwidth = 1e5
    %StepTime   Frequency step burst time (s)
    %   Specify the burst time (in seconds) for each frequency step as a
    %   positive scalar. The default value is 1e-4 (100 microseconds)
    StepTime = 1e-4
    %StepsPerSweep  Number of steps per sweep
    %   Specify the number of frequency steps in a sweep as a positive even
    %   integer. The default value is 64.
    StepsPerSweep = 64  
    %FrequencyOffset    Chirp offset frequency (Hz)
    %   Specify the frequency offset (in Hz) between the two chirp signals
    %   used in the MFSK sweeping as a real scalar. The default value is
    %   1000 (1 kHz).
    FrequencyOffset = 1000
    %OutputFormat     Output signal format
    %   Specify the format of the output signal as one of 'Steps' |
    %   'Sweeps' | 'Samples', where the default is 'Steps'. When you set
    %   the OutputFormat property to 'Steps', the output is in the form of
    %   multiple frequency steps where the number of steps is determined by
    %   the value of the NumSteps property. When you set the OutputFormat
    %   property to 'Sweeps', the output is in the form of multiple sweeps
    %   where the number of sweeps is determined by the value of the
    %   NumSweeps property. When you set the OutputFormat property to
    %   'Samples', the output is in the form of multiple samples where the
    %   number of samples is determined by the value of the NumSamples
    %   property.
    OutputFormat = 'Steps'
end
    
properties (Nontunable, PositiveInteger) 
    %NumSamples     Number of samples in output
    %   Specify the number of samples in each output as a positive integer.
    %   This property only applies when you set the OutputFormat property
    %   to 'Samples'. The default value of this property is 100.
    NumSamples = 100
    %NumSteps  Number of frequency steps in output
    %   Specify the number of frequency steps in each output as a positive
    %   integer. This property only applies when you set the OutputFormat
    %   property to 'Steps'. The default value of this property is 1.
    NumSteps = 1;
    %NumSweeps  Number of sweeps in output
    %   Specify the number of sweeps in each output as a positive integer.
    %   This property only applies when you set the OutputFormat property
    %   to 'Sweeps'. The default value of this property is 1.
    NumSweeps = 1;
end

properties (Constant, Hidden)
    OutputFormatSet = matlab.system.StringSet({'Steps','Samples','Sweeps'});
end

properties (Access = private)
    % internal buffer to hold waveform samples
    pSamples
end

properties (Access = private, PositiveInteger, Nontunable)
    % flag of whether to output by sweep
    pOtuputFormatIdx
end

properties (Access = private, Nontunable)
    % sweep time 
    pSweepTime
    % sweep bandwidth
    pSweepBandwidth
    % number of distinct sweeping sweep
    pNumSweepTimes
    % number of distinct sweeping bandwidth
    pNumSweepBandwidth
    % number of samples in each sweeping sweep
    pNumSamplesPerSweep
    % number of distinct sweeps
    pNumDistinctSweeps
    % number of distinct samples
    pNumDistinctSamples
    % start of each distinct sweep
    pDistinctSweepStartIndex
    % end of each distinct sweep
    pDistinctSweepEndIndex
    % number of samples per distinct sweep
    pNumSamplesPerDistinctSweep
    % output size upper bound
    pOutputUpperBound
end

properties (Access = private)
    % output start sweep index
    pOutputStartSweepIndex
    % output sweep interval
    pOutputSweepInterval
    % output start sample index
    pOutputStartSampleIndex
    % output sample interval
    pOutputSampleInterval
    % output start step index
    pOutputStartStepIndex
    % output step interval
    pOutputStepInterval
end

methods
    function obj = MFSKWaveform(varargin)
        setProperties(obj, nargin, varargin{:});
    end   
    
    function varargout = plot(obj,varargin)
    %plot   Plot the MFSK waveform
    %   plot(Hwav) plots the real part of the waveform specified by Hwav.
    %   Hwav must be a single waveform object.
    %
    %   plot(Hwav,'PlotType',Type) specifies the type of plot using one of
    %   the following strings for Type: [ {real} | imag | complex ].
    %
    %   plot(...,'StepIdx',SID) specifies the index of the frequency step
    %   to plot. SID must be a scalar and its default value is 1.
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
    %   %   Create and plot an MFSK waveform.
    %
    %   hw = phased.MFSKWaveform;
    %   plot(hw);
        
        if ~isempty(coder.target)
            coder.internal.assert(false, ...
             'phased:Waveform:CodegenNotSupported','plot');
        end

        narginchk(1,inf);
        validateattributes(obj,{'phased.MFSKWaveform'},...
            {'scalar'},'','plot');
        
        % Parse argument list to separate waveform specific arguments.
        wfArg = {};     % PV pairs for the waveform
        arglist = {'PlotType','StepIdx'};
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
        StepIdx = 1;
        sigutils.pvparse(wfArg{:});
        PlotType = validatestring(PlotType,{'real','imag','complex'},...
            'plot','PlotType');
        validateattributes(StepIdx,{'double'},{'scalar','integer',...
            'positive'},'plot','StepIdx');
        
        % Generate samples
        wasLocked = isLocked(obj);
        if ~wasLocked
            setup(obj);
        end
        DistinctSweepIdx = obj.getCircularIndex(StepIdx,obj.pNumDistinctSweeps);
        x = obj.pSamples(obj.pDistinctSweepStartIndex(DistinctSweepIdx):...
            obj.pDistinctSweepEndIndex(DistinctSweepIdx));
        t = (0:length(x)-1)/obj.SampleRate;
        
        strxlbl = 'Time (s)'; strylbl = 'Amplitude (v)';
        if strcmp(PlotType,'complex')
            ha1 = subplot(2,1,1);
            h1 = plot(t,real(x),lineArg{:}); grid on;
            set(ha1,'Tag','Real_Part');
            xlabel(strxlbl), ylabel(strylbl);
            title(sprintf('%s: real part, step %d', getWaveformName(obj), StepIdx));
            ha2 = subplot(2,1,2);
            h2 = plot(t,imag(x),lineArg{:}); grid on;
            set(ha2,'Tag','Imaginary_Part');
            xlabel(strxlbl), ylabel(strylbl);
            title(sprintf('%s: imaginary part, step %d', getWaveformName(obj), StepIdx));
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
            title(sprintf('%s: %s part, step %d', getWaveformName(obj), PlotType, StepIdx));
        end

        if nargout > 0
            varargout = {h};
        end
        if ~wasLocked
            release(obj);
        end
    end
end

methods
    function set.SampleRate(obj, value)
        validateattributes(value,{'double'}, {'scalar',...
            'positive','finite'},...
            '','SampleRate');
        obj.SampleRate = value;
    end
    function set.SweepBandwidth(obj, value)
        validateattributes(value,{'double'},...
            {'scalar','positive','finite'},...
            '','SweepBandwidth');
        obj.SweepBandwidth = value;
    end
    function set.StepsPerSweep(obj, value)
        sigdatatypes.validateIndex(value,'','StepsPerSweep',{'scalar','even'});
        obj.StepsPerSweep = value;
    end
    function set.StepTime(obj, value)
        sigdatatypes.validateDuration(value,'','StepTime',{'scalar'});
        obj.StepTime = value;
    end
    function set.FrequencyOffset(obj, value)
        validateattributes(value,{'double'},{'real','scalar','finite','nonempty'},'','FrequencyOffset');
        obj.FrequencyOffset = value;
    end
end
    
methods (Access = protected)
    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        if strcmp(prop,'NumSamples') && (obj.OutputFormat(2) ~= 'a')
            flag = true;
        elseif strcmp(prop,'NumSteps') && (obj.OutputFormat(2) ~= 't')
            flag = true;
        elseif strcmp(prop,'NumSweeps') && (obj.OutputFormat(2) ~= 'w')
            flag = true;
        end
    end
    
    function num = getNumInputsImpl(obj) %#ok<MANU>
        num = 0;
    end
    
    function num = getNumOutputsImpl(obj) %#ok<MANU>
        num = 1;
    end
    
    function validatePropertiesImpl(obj)
        cond = any(rem(obj.StepTime,1./obj.SampleRate));
        if cond
            coder.internal.errorIf(cond,'phased:Waveform:NeedProductInteger','SampleRate','StepTime');
        end
        cond = obj.SampleRate < obj.SweepBandwidth+abs(obj.FrequencyOffset);
        if cond
            coder.internal.warning(...
                'phased:Waveform:MFSKUnderSample', 'SampleRate',...
                'SweepBandwidth', 'FrequencyOffset');
        end
        
    end
    
    function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
        flag = false;  % index == 1
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        if isLocked(obj)
            s.pOutputStartStepIndex = obj.pOutputStartStepIndex;
            s.pOutputStepInterval = obj.pOutputStepInterval;
            s.pOutputStartSweepIndex = obj.pOutputStartSweepIndex;
            s.pOutputSweepInterval = obj.pOutputSweepInterval;
            s.pOutputStartSampleIndex = obj.pOutputStartSampleIndex;
            s.pOutputSampleInterval = obj.pOutputSampleInterval;
            s.pSamples = obj.pSamples;
            s.pOtuputFormatIdx = obj.pOtuputFormatIdx;
            s.pSweepTime = obj.pSweepTime;
            s.pSweepBandwidth = obj.pSweepBandwidth;
            s.pNumSamplesPerSweep = obj.pNumSamplesPerSweep;
            s.pNumSweepTimes = obj.pNumSweepTimes;
            s.pNumSweepBandwidth = obj.pNumSweepBandwidth;
            s.pNumDistinctSweeps = obj.pNumDistinctSweeps;
            s.pNumDistinctSamples = obj.pNumDistinctSamples;
            s.pDistinctSweepStartIndex = obj.pDistinctSweepStartIndex;
            s.pDistinctSweepEndIndex = obj.pDistinctSweepEndIndex;
            s.pNumSamplesPerDistinctSweep = obj.pNumSamplesPerDistinctSweep;
            s.pOutputUpperBound = obj.pOutputUpperBound;
        end
    end
        
    function loadObjectImpl(obj,s,~)
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function resetImpl(obj)
        if (obj.OutputFormat(2) == 'w') %Sweeps
            obj.pOutputStartSweepIndex = 1;
            obj.pOutputSweepInterval = 0:obj.NumSteps;
        elseif (obj.OutputFormat(2) == 'a') %Samples
            obj.pOutputStartSampleIndex = 1;
            obj.pOutputSampleInterval = 0:obj.NumSamples;
        else %steps
            obj.pOutputStartStepIndex = 1;
            obj.pOutputStepInterval = 0:obj.NumSteps;
        end
    end
    
    function setupImpl(obj)
        coder.extrinsic('phased.MFSKWaveform.getSweepOutputSize');
        
        if obj.OutputFormat(2) == 't' %Steps;
            obj.pOtuputFormatIdx = 1; 
        elseif obj.OutputFormat(2) == 'a' %Samples;
            obj.pOtuputFormatIdx = 2;
        else %Sweeps
            obj.pOtuputFormatIdx = 3;
        end
        
        obj.pSweepTime = obj.StepTime;
        obj.pSweepBandwidth = obj.SweepBandwidth;
        
        obj.pNumSweepTimes = numel(obj.pSweepTime);
        obj.pNumSweepBandwidth = numel(obj.pSweepBandwidth);
        obj.pNumDistinctSweeps = obj.StepsPerSweep;
        
        obj.pNumSamplesPerSweep = round(obj.pSweepTime*obj.SampleRate);
        
        obj.pNumSamplesPerDistinctSweep = repmat(obj.pNumSamplesPerSweep,...
            1,obj.pNumDistinctSweeps);
        
        obj.pNumDistinctSamples = sum(obj.pNumSamplesPerDistinctSweep);
        obj.pDistinctSweepStartIndex = cumsum(obj.pNumSamplesPerDistinctSweep) - ...
            obj.pNumSamplesPerDistinctSweep+1;
        obj.pDistinctSweepEndIndex = cumsum(obj.pNumSamplesPerDistinctSweep);
        obj.pSamples = getMatchingWaveform(obj);
        
        obj.pOutputUpperBound = coder.const(obj.getSweepOutputSize(...
             obj.SampleRate, obj.StepTime,obj.OutputFormat, ...
             obj.NumSteps,obj.NumSamples,obj.NumSweeps,obj.StepsPerSweep));
    end
    
    function y = stepImpl(obj)
        switch obj.pOtuputFormatIdx
            case 1  %steps
                
                InitOutputStepIndex = obj.pOutputStartStepIndex+obj.pOutputStepInterval;
                InitOutputStepIndex = obj.getCircularIndex(InitOutputStepIndex,obj.pNumDistinctSweeps);
                obj.pOutputStartStepIndex = InitOutputStepIndex(end);
                OutputStepIndex = InitOutputStepIndex(1:end-1);
                OutputStepsStartSamples = cumsum(obj.pNumSamplesPerDistinctSweep(OutputStepIndex)) - ...
                    obj.pNumSamplesPerDistinctSweep(OutputStepIndex) + 1;
                OutputStepsEndSamples = cumsum(obj.pNumSamplesPerDistinctSweep(OutputStepIndex));

                % preallocate
                numely = obj.pNumSamplesPerSweep*obj.NumSteps;

                y = complex(zeros(numely,1));
                for m = 1:obj.NumSteps
                    y(OutputStepsStartSamples(m):OutputStepsEndSamples(m)) = ...
                        obj.pSamples(obj.pDistinctSweepStartIndex(OutputStepIndex(m)):...
                        obj.pDistinctSweepEndIndex(OutputStepIndex(m)));
                end
            
            case 2  % samples
                OutputSampleIndex = obj.pOutputStartSampleIndex+obj.pOutputSampleInterval;
                OutputSampleIndex = obj.getCircularIndex(OutputSampleIndex,obj.pNumDistinctSamples);
                obj.pOutputStartSampleIndex = OutputSampleIndex(end);
                y = obj.pSamples(OutputSampleIndex(1:end-1));
            case 3  % sweeps
                OutputSweepIndex = obj.pOutputStartSweepIndex+obj.pOutputSweepInterval;
                OutputSweepIndex = obj.getCircularIndex(OutputSweepIndex,1);
                obj.pOutputStartSweepIndex = OutputSweepIndex(end);
                y = repmat(obj.pSamples(:),obj.NumSweeps,1);
                
        end
    end
    
    function wname = getWaveformName(obj) %#ok<MANU>
        wname = 'MFSK waveform';
    end
   
end

methods (Access = private)
    function s = getMatchingWaveform(obj)
        fs = obj.SampleRate;
        freq = linspace(0,obj.SweepBandwidth,obj.pNumDistinctSweeps/2);
        freqMat = [freq;freq+obj.FrequencyOffset];
        tvec = (0:obj.pNumSamplesPerSweep-1).'/fs;
        tMat = bsxfun(@plus,tvec,(0:obj.pNumDistinctSweeps-1)*obj.StepTime);
        tsigMat = exp(1i*2*pi*bsxfun(@times,tMat,freqMat(:).'));
        % making phase continuous
%         phase_temp = 2*pi*freqvec(:).'.*obj.StepTime;
%         phase_comp = cumsum([0 phase_temp(1:end-1)]);
%         s_temp = bsxfun(@times,exp(1i*phase_comp),s_temp);
        s = tsigMat(:);
    end
end
    
methods (Static, Hidden)
    
    function idx = getCircularIndex(idx,N)
        idx = mod(idx-1,N)+1;
    end
    function retSz = getSweepOutputSize(sampleRate,stepTime,outputFormat,...
         numSteps,numSamples,numSweeps,stepsPerSweep)
     
         stepLength = round(stepTime*sampleRate);
         if strcmp(outputFormat,'Steps')
             retSz = stepLength*numSteps;
         elseif strcmp(outputFormat,'Samples')
             retSz = numSamples;
         else %sweeps
             retSz = numSweeps*stepLength*stepsPerSweep;
         end
    end
end

methods (Static,Hidden,Access=protected)  
  function groups = getPropertyGroupsImpl
    groups = matlab.system.display.SectionGroup('phased.MFSKWaveform');
  end
  
  function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:MFSKWaveformTitle')),...
          'Text',getString(message('phased:library:block:MFSKWaveformDesc')));
  end
end
methods (Access = protected) %For Simulink propagation and mask
    function varargout = getOutputNamesImpl(~)
        varargout = {''};
    end
    
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('MFSK');
    end
    
    function varargout = getOutputSizeImpl(obj)
        varargout{1} =  obj.getSweepOutputSize(...
            obj.SampleRate, obj.StepTime,obj.OutputFormat, ...
            obj.NumSteps,obj.NumSamples,obj.NumSweeps,obj.StepsPerSweep);
    end
    function varargout = isOutputFixedSizeImpl(obj) %#ok<MANU>
        varargout{1} = true;
    end
    function varargout = getOutputDataTypeImpl(obj) %#ok<MANU>
        varargout{1} = 'double';
    end
    function varargout = isOutputComplexImpl(obj) %#ok<MANU>
        varargout{1} = true;
    end
    function sts = getSampleTimeImpl(obj)
        N = obj.getSweepOutputSize(...
            obj.SampleRate, obj.StepTime,obj.OutputFormat, ...
            obj.NumSteps,obj.NumSamples,obj.NumSweeps,obj.StepsPerSweep);
        st = phased.internal.samprate2time(obj.SampleRate,N);
        sts = createSampleTime(obj,'Type','Discrete',...
            'SampleTime',st);
    end
end
end
