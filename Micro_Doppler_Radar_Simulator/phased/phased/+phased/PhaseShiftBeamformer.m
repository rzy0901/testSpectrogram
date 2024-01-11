classdef (Sealed,StrictDefaults) PhaseShiftBeamformer < phased.internal.AbstractNarrowbandBeamformer
%PhaseShiftBeamformer Narrowband phase shift beamformer
%   H = phased.PhaseShiftBeamformer creates a conventional phase shift
%   beamformer System object, H. This object performs phase shift
%   beamforming on the received signal.
%
%   H = phased.PhaseShiftBeamformer(Name,Value) creates a phase shift
%   beamformer object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The phase shift beamformer uses the conventional delay and sum
%   beamforming algorithm. It assumes that the signal is narrowband, thus
%   using a phase shift to approximate the required delay. The beamformer
%   preserves the incoming signal power.
%
%   Step method syntax:
%
%   Y = step(H,X) performs phase shift beamforming on the input X, and
%   returns the beamformed output in Y. X is an MxN matrix where N is the
%   number of subarrays if SensorArray contains subarrays, or the number of
%   elements otherwise. Y is an MxL matrix where L is the number of
%   beamforming directions.
%
%   [Y,W] = step(H,X) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is an
%   NxL matrix.
%
%   Y = step(H,X,ANG) uses ANG as the beamforming direction, when you
%   set the DirectionSource property to 'Input port'. ANG is a 2-row
%   matrix whose columns are in the form of [AzimuthAngle; ElevationAngle]
%   (in degrees). The azimuth angle must be between [-180 180] and the
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
%   PhaseShiftBeamformer methods:
%
%   step     - Perform phase shift beamforming (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create phase shift beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   PhaseShiftBeamformer properties:
%
%   SensorArray          - Sensor array
%   PropagationSpeed     - Signal propagation speed 
%   OperatingFrequency   - Operating frequency 
%   DirectionSource      - Source of beamforming direction
%   Direction            - Beamforming direction 
%   NumPhaseShifterBits  - Number of bits in phase shifters
%   WeightsNormalization - Weights normalizing method
%   WeightsOutputPort    - Enable weights output
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %    Apply phase shift beamforming to the signal received by a 
%   %    5-element ULA. The beamforming direction is 45 degrees azimuth and
%   %    0 degrees elevation.
%   
%   % simulate signal
%   t = (0:1000)';
%   x = sin(2*pi*0.01*t);
%   c = 3e8; Fc = 3e8;
%   incidentAngle = [45; 0];
%   array = phased.ULA('NumElements',5);
%   x = collectPlaneWave(array,x,incidentAngle,Fc,c);
%   noise = 0.1*(randn(size(x)) + 1j*randn(size(x)));
%   rx = x + noise;
%
%   % beamforming
%   beamformer = phased.PhaseShiftBeamformer('SensorArray',array,...
%         'OperatingFrequency',Fc,'PropagationSpeed',c,...
%         'Direction',incidentAngle,'WeightsOutputPort',true);
%   [y,w] = beamformer(rx);
%   plot(t,real(rx(:,3)),'r:',t,real(y));
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
%
%   % plot response pattern
%   figure; pattern(array,Fc,-180:180,0,'PropagationSpeed',c,'Weights',w);
%
%   See also phased, phased.MVDRBeamformer, phased.LCMVBeamformer,
%   phased.SubbandPhaseShiftBeamformer.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %WeightsNormalization   Weights normalizing method 
        %   Specify the method of normalizing the beamformer weights as one
        %   of 'Distortionless' | 'Preserve power', where the default is
        %   'Distortionless'. If WeightsNormalization is 'Distortionless',
        %   the gain toward the beamforming direction is 0 dB. If
        %   WeightsNormalization is 'Preserve power', then the norm of the
        %   weights is 1.
        WeightsNormalization = 'Distortionless'
    end

    properties(Access = private) 
        %pre-calculated weights
        pWeights;
    end
    
    properties (Constant,Hidden)
        WeightsNormalizationSet = matlab.system.StringSet(...
            {'Distortionless','Preserve power'});       
    end

    methods
        function obj = PhaseShiftBeamformer(varargin)
            obj@phased.internal.AbstractNarrowbandBeamformer(varargin{:});
        end

    end


    methods (Access = 'protected')

        function flag = isMultipleInputAnglesAllowed(obj) %#ok<MANU>
            flag = true;
        end
        
        function flag = isSubarraySupported(obj) %#ok<MANU>
            flag = true;
        end
        
        function setupImpl(obj,x,varargin)
            setupImpl@phased.internal.AbstractNarrowbandBeamformer(obj,x);

            if (obj.DirectionSource(1) == 'P') %Property
                obj.pWeights = calcWeights(obj,obj.Direction);
            end
            
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = isInputComplexityLockedImpl@phased.internal.AbstractNarrowbandBeamformer(obj,index);
            if (obj.DirectionSource(1) == 'I')  && (index == 2) %Input
                flag = true;
            end
        end
        
        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractNarrowbandBeamformer(obj);
            if isLocked(obj)
                s.pWeights = obj.pWeights;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractNarrowbandBeamformer(obj,s);
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

            classtouse = class(x);
            if (obj.DirectionSource(1) == 'P') %Property
                w = cast(obj.pWeights,classtouse);
            else
                ang = varargin{1};
                validateAngleRange(obj,ang);
                w = cast(calcWeights(obj,ang),classtouse);
            end
            y = x*conj(w);

        end

        function validateInputsImpl(obj,x,varargin)
            validateInputsImpl@phased.internal.AbstractNarrowbandBeamformer(obj,x);
            if strncmpi(obj.DirectionSource,'i',1)
                ang = varargin{1};
                validateInputAngle(obj,ang);
            end
        end
    end
    
    methods (Access = private)

        function w = calcWeights(obj,ang)

            sv = step(obj.cSteeringVector,cast(obj.OperatingFrequency,'double'),cast(ang,'double'));
            w = complex(zeros(size(sv)));
            if (obj.WeightsNormalization(1) == 'D') %Distortionless
                for m = size(ang,2):-1:1
                    temp_sv = sv(:,m);
                    w(:,m) = temp_sv/real(temp_sv'*temp_sv);  % distortionless
                end
            else
                for m = size(ang,2):-1:1
                    temp_sv = sv(:,m);
                    w(:,m) = temp_sv/norm(temp_sv);  % preserve power
                end
            end

        end
    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractNarrowbandBeamformer('subarray');
        props = {...
          'OperatingFrequency',...
          'DirectionSource',...
          'Direction',...
          'NumPhaseShifterBits',...
          'WeightsNormalization',...
          'WeightsOutputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:PhaseShiftBeamformerTitle')),...
              'Text',getString(message('phased:library:block:PhaseShiftBeamformerDesc')));
      end        
    end
    methods (Access = protected) %for Simulink
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Phase Shift\nBeamformer');
        end
    end
end

