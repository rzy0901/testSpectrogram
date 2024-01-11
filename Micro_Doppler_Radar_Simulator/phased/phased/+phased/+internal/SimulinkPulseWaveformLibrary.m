classdef(Sealed,StrictDefaults) SimulinkPulseWaveformLibrary < phased.internal.AbstractPulseWaveformLibrary & ...
        matlab.system.mixin.CustomIcon 
%This class is for internal use only. It may be removed in the future.

%   Copyright 2017 The MathWorks, Inc.
    
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
%   reset               - Reset states of PulseWaveform object
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

%   Reference
%   [1] Cochran et.al., Waveform Libraries, Measures of effectiveness for
%   radar scheduling,IEEE Signal Processing Magazine, 2009.
%
%   [2]Richards, M. A., Fundamentals of Radar Signal Processing. New York:
%   McGraw-Hill, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %SimulationTimeSource   Source of simulation sample time
    %   Set simulation time source to one of 'Derive from waveform
    %   parameters' | 'Inherit from Simulink engine', where the default is
    %   'Derive from waveform parameters'. This property applies only in
    %   Simulink. If the simulation time is set to 'Derive from waveform
    %   parameters', the block runs at a variable rate determined by the
    %   PRF of the selected waveform; otherwise the block runs at a fixed
    %   rate so the elapsed time is a constant.
    SimulationTimeSource = 'Derive from waveform parameters'
end

properties(Constant, Hidden)
    SimulationTimeSourceSet = matlab.system.StringSet(...
        {'Derive from waveform parameters',...
        'Inherit from Simulink engine'});
end

properties (Access = protected, Nontunable, Logical)
    %pInheritSampleTime  Inherit time
    %   Set this property to true to inherit sample time from Simulink
    %   engine. Set this property to false to set the block's own sample
    %   time. This property applies only in Simulink.
    pInheritSampleTime = false
end

methods
    % Constructor
    
    function obj = SimulinkPulseWaveformLibrary(varargin)
        obj@phased.internal.AbstractPulseWaveformLibrary(varargin{:});
    end
    
end

methods (Access = protected)
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractPulseWaveformLibrary(obj,prop);
        if strcmp(prop,'SimulationTimeSource') && ~isComplexityPropagated(obj)
            flag = true;
        end
    end
    
    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractPulseWaveformLibrary(obj);
        inp = obj.WaveformSpecification;
        validateDuplicateWarningCond(obj,inp);
    end
    
    function setupImpl(obj,~)
        setupImpl@phased.internal.AbstractPulseWaveformLibrary(obj);
        obj.pInheritSampleTime = strcmp(obj.SimulationTimeSource,...
            'Inherit from Simulink engine');
    end
    
    function flag = isInDVRMode(obj)
        flag = isComplexityPropagated(obj) && ~obj.pInheritSampleTime;
    end

    function resetImpl(obj)
        resetImpl@phased.internal.AbstractPulseWaveformLibrary(obj)
        if isInDVRMode(obj)
            setNumTicksUntilNextHit(obj,1);
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractPulseWaveformLibrary(obj);
        if isLocked(obj)
            s.pInheritSampleTime = obj.pInheritSampleTime;
        end
    end
  
    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractPulseWaveformLibrary(obj,s);        
        obj.pInheritSampleTime = s.pInheritSampleTime;        
        s = rmfield(s,'pInheritSampleTime');
    end 
    
    function loadObjectImpl(obj,s,wasLocked)        
        loadObjectImpl@phased.internal.AbstractPulseWaveformLibrary(obj,s,wasLocked);
    end
            
    function sts = getSampleTimeImpl(obj)
        if isComplexityPropagated(obj) && ...
                ~strcmp(obj.SimulationTimeSource,'Inherit from Simulink engine')
            sts = createSampleTime(obj,'Type','Controllable',...
                'TickTime',1/obj.SampleRate);
        else
            sts = createSampleTime(obj,'Type','Inherited');
        end
    end
end

methods (Static, Hidden, Access = protected)
    
    function groups = getPropertyGroupsImpl
        props = {'SampleRate',...
            'WaveformSpecification',...
            'SimulationTimeSource',...
            };
        groups = matlab.system.display.Section('Title',...
            'Pulse Waveform Library', ...
            'PropertyList', props);
        
    end
    
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:PulseWaveformLibraryTitle')),...
            'Text',getString(message('phased:library:block:PulseWaveformLibraryDesc')));
    end
end

methods (Access=protected)
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Pulse Waveform\nLibrary');
    end
end

end
