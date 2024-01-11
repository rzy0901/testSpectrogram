classdef (Sealed,StrictDefaults) FMCWWaveform < matlab.System & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.SampleTime
%FMCWWaveform   FMCW waveform
%   H = phased.FMCWWaveform creates an FMCW waveform System object, H. This
%   object generates samples of a linear frequency modulated continuous
%   wave (FMCW) waveform.
%
%   H = phased.FMCWWaveform(Name,Value) creates an FMCW waveform object, H,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H) returns samples of the FMCW waveform in a column vector Y.
%   Y can contain either a certain number of sweeps or a certain number of
%   samples.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H) and y = H() are
%   equivalent.
%
%   FMCWWaveform methods:
%
%   step                - Return samples of the FMCW waveform
%   release             - Allow property value and input characteristics 
%                         changes
%   clone               - Create an FMCW waveform object with same property
%                         values
%   isLocked            - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>               - Reset states of FMCW waveform object
%   plot                - Plot the FMCW waveform
%
%   FMCWWaveform properties:
%
%   SampleRate     - Sample rate 
%   SweepTime      - Sweep time
%   SweepBandwidth - Sweep bandwidth 
%   SweepDirection - Sweep direction
%   SweepInterval  - Sweep interval
%   OutputFormat   - Output signal format
%   NumSamples     - Number of samples in output
%   NumSweeps      - Number of sweeps in output
%
%   % Examples:
%
%   % Example 1:
%   %   Create and plot an upsweep FMCW waveform.
%
%   waveform = phased.FMCWWaveform('SweepBandwidth',1e5,...
%               'OutputFormat','Sweeps','NumSweeps',2);
%   plot(waveform);
%
%   % Example 2:
%   %   Generate output samples of a triangle sweep FMCW waveform and
%   %   examine the sweep using spectrogram.
%
%   waveform = phased.FMCWWaveform('SweepBandwidth',1e7,...
%               'SampleRate',2e7,'SweepDirection','Triangle',...
%               'NumSweeps',2);
%   x = waveform();
%   spectrogram(x,32,16,32,waveform.SampleRate,'yaxis');
%
%   See also phased, phased.LinearFMWaveform.
    
%   Copyright 2012-2017 The MathWorks, Inc.

%   Reference
%   [1] Merrill Skolnik, Introduction to Radar Systems, McGraw-Hill, 1962
%   [2] Vadim Issakov, Microwave Circuits for 24 GHz Automotive Radar in 
%       Silicon-based Technologies, Springer, 2010


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %SampleRate Sample rate (Hz)
    %   Specify the sample rate (in Hz) as a positive scalar. The default
    %   value of this property is 1e6 (1 MHz).
    SampleRate = 1e6;
    %SweepTime    Sweep time (s)
    %   Specify the time of linear FM sweeps (in seconds) as a positive
    %   scalar or row vector. The default value is 1e-4 (100 microseconds).
    %   When SweepTime is a vector, the output sweeps use elements in the
    %   vector as their sweep times one after another in cycle.
    %
    %   If both SweepTime and SweepBandwidth are nonscalar, the sizes of
    %   the two vectors must be the same.
    SweepTime = 1e-4
    %SweepBandwidth     Sweep bandwidth (Hz)
    %   Specify the bandwidth of the linear FM sweeping (in Hz) as a
    %   positive scalar or row vector. The default value is 1e5 (100 kHz).
    %   When SweepBandwidth is a vector, the output sweeps use elements in
    %   the vector as their sweep bandwidth one after another in cycle.
    %
    %   If both SweepTime and SweepBandwidth are nonscalar, the sizes of
    %   the two vectors must be the same.
    SweepBandwidth = 1e5
    %SweepDirection     Sweep direction
    %   Specify the direction of the linear FM sweep as one of 'Up' |
    %   'Down' | 'Triangle', where the default is 'Up'. 
    %
    %   If you set SweepDirection to 'Triangle', there are two sweeps, one
    %   up and one down, for each pair of SweepTime and SweepBandwidth.
    SweepDirection = 'Up'
    %SweepInterval  Sweep interval
    %   Specify the linear FM sweeping interval using one of 'Positive' |
    %   'Symmetric', where the default is 'Positive'. If SweepInterval is
    %   'Positive', the waveform sweeps in the interval between 0 and B
    %   where B is the sweeping bandwidth. If SweepInterval is 'Symmetric',
    %   the waveform sweeps in the interval between -B/2 and B/2.
    SweepInterval = 'Positive'
    %OutputFormat     Output signal format
    %   Specify the format of the output signal as one of 'Sweeps' |
    %   'Samples', where the default is 'Sweeps'. When you set the
    %   OutputFormat property to 'Sweeps', the output is in the form of
    %   multiple sweeps where the number of sweeps is determined by the
    %   value of the NumSweeps property. When you set the OutputFormat
    %   property to 'Samples', the output is in the form of multiple
    %   samples where the number of samples is determined by the value of
    %   the NumSamples property.
    OutputFormat = 'Sweeps'
end
    
properties (Nontunable, PositiveInteger) 
    %NumSamples     Number of samples in output
    %   Specify the number of samples in each output as a positive integer.
    %   This property only applies when you set the OutputFormat property
    %   to 'Samples'. The default value of this property is 100.
    NumSamples = 100
    %NumSweeps  Number of sweeps in output
    %   Specify the number of sweeps in each output as a positive integer.
    %   This property only applies when you set the OutputFormat property
    %   to 'Sweeps'. The default value of this property is 1.
    NumSweeps = 1;
end

properties (Constant, Hidden)
    OutputFormatSet = matlab.system.StringSet({'Sweeps','Samples'});
    SweepDirectionSet = matlab.system.StringSet({'Up','Down','Triangle'});
    SweepIntervalSet = matlab.system.StringSet({'Positive','Symmetric'});
end

properties (Access = private)
    % internal buffer to hold waveform samples
    pSamples
end

properties (Access = private, Logical, Nontunable)
    % flag of whether to output by sweep
    pOutputBySweep
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
end

methods
    function obj = FMCWWaveform(varargin)
        setProperties(obj, nargin, varargin{:});
    end   
    
    function varargout = plot(obj,varargin)
    %plot   Plot the FMCW waveform
    %   plot(Hwav) plots the real part of the waveform specified by Hwav.
    %   Hwav must be a single waveform object.
    %
    %   plot(Hwav,'PlotType',Type) specifies the type of plot using one of
    %   the following strings for Type: [ {real} | imag | complex ].
    %
    %   plot(...,'SweepIdx',SID) specifies the index of the sweep to
    %   plot. SID must be a scalar and its default value is 1.
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
    %   %   Create and plot an FMCW waveform.
    %
    %   hw = phased.FMCWWaveform;
    %   plot(hw);
        
        if ~isempty(coder.target)
            coder.internal.assert(false, ...
             'phased:Waveform:CodegenNotSupported','plot');
        end

        narginchk(1,inf);
        validateattributes(obj,{'phased.FMCWWaveform'},...
            {'scalar'},'','plot');
        
        % Parse argument list to separate waveform specific arguments.
        wfArg = {};     % PV pairs for the waveform
        arglist = {'PlotType','SweepIdx'};
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
        SweepIdx = 1;
        sigutils.pvparse(wfArg{:});
        PlotType = validatestring(PlotType,{'real','imag','complex'},...
            'plot','PlotType');
        validateattributes(SweepIdx,{'double'},{'scalar','integer',...
            'positive'},'plot','SweepIdx');
        
        % Generate samples
        wasLocked = isLocked(obj);
        if ~wasLocked
            setup(obj);
        end
        DistinctSweepIdx = obj.getCircularIndex(SweepIdx,obj.pNumDistinctSweeps);
        x = obj.pSamples(obj.pDistinctSweepStartIndex(DistinctSweepIdx):...
            obj.pDistinctSweepEndIndex(DistinctSweepIdx));
        t = (0:length(x)-1)/obj.SampleRate;
        
        strxlbl = 'Time (s)'; strylbl = 'Amplitude (v)';
        if strcmp(PlotType,'complex')
            ha1 = subplot(2,1,1);
            h1 = plot(t,real(x),lineArg{:}); grid on;
            set(ha1,'Tag','Real_Part');
            xlabel(strxlbl), ylabel(strylbl);
            title(sprintf('%s: real part, sweep %d', getWaveformName(obj), SweepIdx));
            ha2 = subplot(2,1,2);
            h2 = plot(t,imag(x),lineArg{:}); grid on;
            set(ha2,'Tag','Imaginary_Part');
            xlabel(strxlbl), ylabel(strylbl);
            title(sprintf('%s: imaginary part, sweep %d', getWaveformName(obj), SweepIdx));
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
            title(sprintf('%s: %s part, sweep %d', getWaveformName(obj), PlotType, SweepIdx));
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
            {'row','positive','finite'},...
            '','SweepBandwidth');
        obj.SweepBandwidth = value;
    end
    function set.SweepTime(obj, value)
        sigdatatypes.validateDuration(value,'','SweepTime',{'row'});
        obj.SweepTime = value;
    end
end
    
methods (Access = protected)
    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        pulseflag = (obj.OutputFormat(2) == 'w'); %Sweeps
        if strcmp(prop,'NumSamples') && pulseflag
            flag = true;
        elseif strcmp(prop,'NumSweeps') && ~pulseflag
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
        cond = ~isscalar(obj.SweepTime) && ~isscalar(obj.SweepBandwidth) && ...
                numel(obj.SweepTime) ~= numel(obj.SweepBandwidth);
        if cond
            coder.internal.errorIf(cond,'phased:Waveform:DimensionMismatchWhenNonScalar','SweepTime','SweepBandwidth');
        end
        cond = any(rem(obj.SweepTime,1./obj.SampleRate));
        if cond
            coder.internal.errorIf(cond,'phased:Waveform:NeedProductInteger','SampleRate','SweepTime');
        end
        
    end
    
    function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
        flag = false;  % index == 1
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        if isLocked(obj)
            s.pOutputStartSweepIndex = obj.pOutputStartSweepIndex;
            s.pOutputSweepInterval = obj.pOutputSweepInterval;
            s.pOutputStartSampleIndex = obj.pOutputStartSampleIndex;
            s.pOutputSampleInterval = obj.pOutputSampleInterval;
            s.pSamples = obj.pSamples;
            s.pOutputBySweep = obj.pOutputBySweep;
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
            obj.pOutputSweepInterval = 0:obj.NumSweeps;
        else
            obj.pOutputStartSampleIndex = 1;
            obj.pOutputSampleInterval = 0:obj.NumSamples;
        end
    end
    
    function setupImpl(obj)
        coder.extrinsic('phased.FMCWWaveform.getSweepOutputSize');
        
        obj.pOutputBySweep = (obj.OutputFormat(2) == 'w'); %Sweeps
        
        if (obj.SweepDirection(1) == 'T') %Triangle
            temp = repmat(obj.SweepTime,2,1);
            obj.pSweepTime = temp(:).';
            temp = repmat(obj.SweepBandwidth,2,1);
            obj.pSweepBandwidth = temp(:).';
        else
            obj.pSweepTime = obj.SweepTime;
            obj.pSweepBandwidth = obj.SweepBandwidth;
        end
        
        obj.pNumSweepTimes = numel(obj.pSweepTime);
        obj.pNumSweepBandwidth = numel(obj.pSweepBandwidth);
        obj.pNumDistinctSweeps = max(obj.pNumSweepTimes,obj.pNumSweepBandwidth);
        
        obj.pNumSamplesPerSweep = round(obj.pSweepTime*obj.SampleRate);
        
        obj.pNumSamplesPerDistinctSweep = repmat(obj.pNumSamplesPerSweep,...
            1,obj.pNumDistinctSweeps/obj.pNumSweepTimes);
        
        obj.pNumDistinctSamples = sum(obj.pNumSamplesPerDistinctSweep);
        obj.pDistinctSweepStartIndex = cumsum(obj.pNumSamplesPerDistinctSweep) - ...
            obj.pNumSamplesPerDistinctSweep+1;
        obj.pDistinctSweepEndIndex = cumsum(obj.pNumSamplesPerDistinctSweep);
        obj.pSamples = getMatchingWaveform(obj);
        
        obj.pOutputUpperBound = coder.const(obj.getSweepOutputSize(...
             obj.SampleRate, obj.SweepTime,obj.SweepDirection,obj.OutputFormat, ...
             obj.NumSweeps,obj.NumSamples));
    end
    
    function y = stepImpl(obj)
        if obj.pOutputBySweep
            InitOutputSweepIndex = obj.pOutputStartSweepIndex+obj.pOutputSweepInterval;
            InitOutputSweepIndex = obj.getCircularIndex(InitOutputSweepIndex,obj.pNumDistinctSweeps);
            obj.pOutputStartSweepIndex = InitOutputSweepIndex(end);
            OutputSweepIndex = InitOutputSweepIndex(1:end-1);
            OutputSweepsStartSamples = cumsum(obj.pNumSamplesPerDistinctSweep(OutputSweepIndex)) - ...
                obj.pNumSamplesPerDistinctSweep(OutputSweepIndex) + 1;
            OutputSweepsEndSamples = cumsum(obj.pNumSamplesPerDistinctSweep(OutputSweepIndex));
            
            % preallocate
            if numel(obj.SweepTime) == 1
 	            numely = obj.pNumSamplesPerSweep(1)*obj.NumSweeps;
            else
                numely = sum(obj.pNumSamplesPerDistinctSweep(OutputSweepIndex));
                %upper bound of y
                assert(numely <= obj.pOutputUpperBound);
            end	                   
            
            y = complex(zeros(numely,1));
            for m = 1:obj.NumSweeps
                y(OutputSweepsStartSamples(m):OutputSweepsEndSamples(m)) = ...
                    obj.pSamples(obj.pDistinctSweepStartIndex(OutputSweepIndex(m)):...
                    obj.pDistinctSweepEndIndex(OutputSweepIndex(m)));
            end
            
        else
            OutputSampleIndex = obj.pOutputStartSampleIndex+obj.pOutputSampleInterval;
            OutputSampleIndex = obj.getCircularIndex(OutputSampleIndex,obj.pNumDistinctSamples);
            obj.pOutputStartSampleIndex = OutputSampleIndex(end);
            y = obj.pSamples(OutputSampleIndex(1:end-1));
        end
    end
    
    function wname = getWaveformName(obj) %#ok<MANU>
        wname = 'FMCW waveform';
    end
   
end

methods (Access = private)
    function s = getMatchingWaveform(obj)
        fs = obj.SampleRate;
        s = complex(zeros(sum(obj.pNumSamplesPerDistinctSweep),1));
        for m = 1:obj.pNumDistinctSweeps
            tau = obj.pSweepTime(obj.getCircularIndex(m,obj.pNumSweepTimes));
            N = obj.pNumSamplesPerDistinctSweep(obj.getCircularIndex(m,obj.pNumSweepTimes));
            t = (0:N-1)'/fs;
            idx = (1:N)+sum(obj.pNumSamplesPerDistinctSweep(1:m-1));
            phaseslope = obj.pSweepBandwidth(obj.getCircularIndex(m,obj.pNumSweepBandwidth))/tau;
            if (obj.SweepDirection(1) == 'U') %Up
                if (obj.SweepInterval(1) == 'P') %Positive
                    s(idx) = exp(1i*pi*phaseslope*t.^2);
                else
                    s(idx) = exp(1i*pi*phaseslope*t.*(t-tau));
                end
            elseif (obj.SweepDirection(1) == 'D') %Down
                if (obj.SweepInterval(1) == 'P') %Positive
                    s(idx) = exp(1i*pi*phaseslope*t.*(2*tau-t));
                else
                    s(idx) = exp(1i*pi*phaseslope*t.*(tau-t));
                end
            else  % strcmp(obj.SweepDirection,'Triangle')
                if rem(m,2)
                    if (obj.SweepInterval(1) == 'P') %Positive
                        s(idx) = exp(1i*pi*phaseslope*t.^2);
                    else
                        s(idx) = exp(1i*pi*phaseslope*t.*(t-tau));
                    end
                else
                    if (obj.SweepInterval(1) == 'P') %Positive
                        s(idx) = exp(1i*pi*phaseslope*t.*(2*tau-t));
                    else
                        s(idx) = exp(1i*pi*phaseslope*t.*(tau-t));
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
    function retSz = getSweepOutputSize(sampleRate,sweepTime,sweepDir,outputFormat,...
         numSweeps,numSamples)
         %num of samples for each pulse in PRF vector
         sweepLength = round(sweepTime*sampleRate);
         if strcmp(outputFormat,'Sweeps')
             numSweepTimes = numel(sweepTime);
             if numSweepTimes == 1
                 retSz = sweepLength*numSweeps;
             else
                 if strcmp(sweepDir,'Triangle')
                     sweepTimeIdx = reshape(repmat(1:numSweepTimes,2,1),[],1);
                 else
                     sweepTimeIdx = (1:numSweepTimes).';
                 end
                 %Get total number of samples of all sweep combinations
                 staggeredIndex = bsxfun(@circshift,sweepTimeIdx,-(0:numSweeps-1));
                 currentSizes = sum(sweepLength(staggeredIndex),2);
                 %Get the upperbound
                 retSz =  max(currentSizes);
             end
         else
             retSz = numSamples;
         end
    end
end

methods (Static,Hidden,Access=protected)  
  function groups = getPropertyGroupsImpl
    props = {...
      'SampleRate',...
      'SweepTime',...
      'SweepBandwidth',...
      'SweepDirection',...
      'SweepInterval',...
      'OutputFormat',...
      'NumSamples',...
      'NumSweeps'};
    groups = matlab.system.display.Section('Title', 'Parameters', ...
                                           'PropertyList', props);
  end
  
  function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:FMCWWaveformTitle')),...
          'Text',getString(message('phased:library:block:FMCWWaveformDesc')));
  end
end
methods (Access = protected) %For Simulink propagation and mask
    function varargout = getOutputNamesImpl(~)
        varargout = {''};
    end
    
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('FMCW');
    end
    
    function varargout = getOutputSizeImpl(obj)
        varargout{1} =  obj.getSweepOutputSize(...
            obj.SampleRate, obj.SweepTime,obj.SweepDirection,obj.OutputFormat, ...
            obj.NumSweeps,obj.NumSamples);
    end
    function varargout = isOutputFixedSizeImpl(obj)
        if strcmp(obj.OutputFormat,'Sweeps') && ...
                numel(obj.SweepTime) ~= 1
            error(message('phased:Waveform:StaggeredNotSupported','SweepTime','OutputFormat','Sweeps'));
            %varargout{1} = false;
        else
            varargout{1} = true;
        end
    end
    function varargout = getOutputDataTypeImpl(obj) %#ok<MANU>
        varargout{1} = 'double';
    end
    function varargout = isOutputComplexImpl(obj) %#ok<MANU>
        varargout{1} = true;
    end
    function sts = getSampleTimeImpl(obj)
        N = obj.getSweepOutputSize(...
            obj.SampleRate, obj.SweepTime,obj.SweepDirection,obj.OutputFormat, ...
            obj.NumSweeps,obj.NumSamples);
        st = phased.internal.samprate2time(obj.SampleRate,N);
        sts = createSampleTime(obj,'Type','Discrete',...
            'SampleTime',st);
    end
end
end
