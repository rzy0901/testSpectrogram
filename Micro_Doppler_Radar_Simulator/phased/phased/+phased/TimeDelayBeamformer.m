classdef (Sealed,StrictDefaults) TimeDelayBeamformer < phased.internal.AbstractTimeDomainBeamformer
%TimeDelayBeamformer Time delay beamformer
%   H = phased.TimeDelayBeamformer creates a time delay beamformer System
%   object, H. This object performs delay and sum beamforming on the
%   received signal using time delays.
%
%   H = phased.TimeDelayBeamformer(Name,Value) creates a time delay
%   beamformer object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) performs time delay beamforming on the input X, and
%   returns the beamformed output in Y. X is an MxN matrix where N is the
%   number of elements of the sensor array. Y is a length-M column vector.
%
%   [Y,W] = step(H,X) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is a
%   length-N column vector. For a time delay beamformer, the weights are
%   constant because the beamformer simply adds all the channels together
%   and scales the result to preserve the signal power.
%
%   Y = step(H,X,ANG) uses ANG as the beamforming direction, when you
%   set the DirectionSource property to 'Input port'. ANG is a length-2
%   column vector in the form of [AzimuthAngle; ElevationAngle] (in
%   degrees). The azimuth angle must be between [-180 180] and the
%   elevation angle must be between [-90 90].
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [Y,W] = step(H,X,ANG)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   TimeDelayBeamformer methods:
%
%   step     - Perform time delay beamforming (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create time delay beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   TimeDelayBeamformer properties:
%
%   SensorArray       - Sensor array
%   PropagationSpeed  - Signal propagation speed
%   SampleRate        - Sample rate
%   DirectionSource   - Source of beamforming direction
%   Direction         - Beamforming direction
%   WeightsOutputPort - Enable weights output
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Apply a time delay beamformer to an 11-element array. The incident
%   %   angle of the signal is -50 degrees in azimuth and 30 degrees in
%   %   elevation.
%
%   % signal simulation
%   array = phased.ULA('NumElements',11,'ElementSpacing',0.04);
%   array.Element.FrequencyRange = [20 20000];
%   fs = 8e3; t = 0:1/fs:0.3;
%   x = chirp(t,0,1,500);
%   c = 340; % Wave propagation speed (m/s)
%   sigcol = phased.WidebandCollector('Sensor',array,...
%            'PropagationSpeed',c,'SampleRate',fs,'ModulatedInput',false);
%   incidentAngle = [-50; 30];
%   x = sigcol(x.',incidentAngle);
%   noise = 0.2*randn(size(x));
%   rx = x+noise;
%
%   % beamforming
%   beamformer = phased.TimeDelayBeamformer('SensorArray',array,...
%         'SampleRate',fs,'PropagationSpeed',c,'Direction',incidentAngle);
%   y = beamformer(rx);
%   plot(t,rx(:,6),'r:',t,y);
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
%
%   See also phased, phased.FrostBeamformer, phased.PhaseShiftBeamformer,
%   phased.TimeDelayLCMVBeamformer, phased.SubbandPhaseShiftBeamformer.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable, Logical) 
        %WeightsOutputPort    Enable weights output
        %   Set this property to true to output the weights used in the
        %   beamformer. Set this property to false to not output the
        %   weights. The default value of this property is false.
        WeightsOutputPort = false;
    end
    
    properties (Access = private, Nontunable)
        % pre-calculated weights
        pWeights;
    end

    methods
        function obj = TimeDelayBeamformer(varargin)
            obj = obj@phased.internal.AbstractTimeDomainBeamformer(varargin{:});
        end

    end

    methods (Access = 'protected')

        function num = getNumOutputsImpl(obj)
            num = getNumOutputsImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            if obj.WeightsOutputPort
                num = 2;
            end
        end
        
        function setupImpl(obj,x,varargin) 
            classtouse=class(x);
            setupImpl@phased.internal.AbstractTimeDomainBeamformer(obj,x);
            %sv = coder.internal.const(ones(obj.pDOF,1));
            %obj.pWeights = coder.internal.const(sv/real(sv'*sv));
            obj.pWeights = cast((ones(obj.pDOF,1)/obj.pDOF),classtouse);
        end
        
        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index) 
            flag = isInputComplexityLockedImpl@phased.internal.AbstractTimeDomainBeamformer(obj,index);
            if strncmpi(obj.DirectionSource,'i',1) && (index == 2)
                flag = true;
            end
        end
        
        function validateInputsImpl(obj,x,varargin)
            validateInputsImpl@phased.internal.AbstractTimeDomainBeamformer(obj,x);
            if strncmpi(obj.DirectionSource,'i',1)
                ang = varargin{1};
                validateInputAngle(obj,ang);
            end
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            if isLocked(obj)
                s.pWeights = obj.pWeights;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractTimeDomainBeamformer(obj,s);
            if isfield(s,'isLocked')
                s = rmfield(s,'isLocked');
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function [y, w] = stepImpl(obj,x,varargin)
            if (obj.DirectionSource(1) == 'P') %Port
                ang = cast(obj.Direction,class(x));
            else
                ang = cast(varargin{1},class(x));
                validateAngleRange(obj,ang);
            end
            x = steer(obj,x,ang);
            w = obj.pWeights;
            y = x*conj(w);

        end

    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractTimeDomainBeamformer('array');
        props = {...
                 'SampleRateFromInputCheckbox',...
                 'SampleRate',...
                 'DirectionSource',...
                 'Direction',...
                 'WeightsOutputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:TimeDelayBeamformerTitle')),...
              'Text',getString(message('phased:library:block:TimeDelayBeamformerDesc')));
      end    
    end
    methods (Access = protected) %for Simulink
        function varargout = isOutputComplexImpl(obj)
            varargout = { propagatedInputComplexity(obj,1), false};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Time Delay\nBeamformer');
        end
    end

end

