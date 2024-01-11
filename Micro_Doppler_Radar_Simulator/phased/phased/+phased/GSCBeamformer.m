classdef (Sealed,StrictDefaults) GSCBeamformer < phased.internal.AbstractTimeDomainBeamformer
%GSCBeamformer    Generalized Sidelobe Canceller
%   H = GSCBeamformer creates a generalized sidelobe canceller (GSC) System
%   Object, H. This object performs GSC beamforming on the received signal.
%
%   H = GSCBeamformer(Name,Value) creates a generalized sidelobe canceller
%   (GSC) System Object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The Generalized Sidelobe Canceller (GSC) beamforming algorithm is
%   related to the Frost beamformer, but offers some implementation
%   advantages. The GSC beamformer consists of two paths: a conventional
%   beamformer path and a sidelobe canceling path. The algorithm first
%   pre-steers the array to the beamforming direction and then adaptively
%   chooses filter weights to minimize power at the output of the sidelobe
%   canceling path. The final beamformed signal is the difference between
%   the outputs of the two paths.
%
%   Step method syntax:
%
%   Y = step(H,X) performs GSC beamforming on the input X, and returns the
%   beamformed output in Y. X is an NxM matrix where N is the number of
%   samples of the input signal and M is the number of elements of the
%   sensor array. Y is a length-N column vector. N must be larger than the
%   FIR filter length specified in the FilterLength property.
%
%   Y = step(H,X,ANG) uses ANG as the beamforming direction when the
%   DirectionSource property to set to 'Input port'. ANG is a length-2
%   column vector in the form of [AzimuthAngle; ElevationAngle] (in
%   degrees). The azimuth angle must be between [-180 180] and the
%   elevation angle must be between [-90 90].
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   GSCBeamformer methods:
%
%   step     - Perform GSC beamforming (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create GSC beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of GSC Beamformer object
%
%   GSCBeamformer properties:
%
%   SensorArray               - Sensor array
%   PropagationSpeed          - Signal propagation speed
%   SampleRate                - Sample rate
%   FilterLength              - Fixed target signal filter coefficients
%   LMSStepSize               - LMS adaptive filter step size
%   Direction                 - Beamforming direction
%   DirectionSource           - Source of beamforming direction
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Apply a GSC beamformer to an 11-element array. The incident 
%   %   angle of the signal is -50 degrees in azimuth and 30 degrees in
%   %   elevation. Compare the GSC beamformed signal to the output of a 
%   %   Frost beamformer.
%
%   % Signal simulation
%   ha = phased.ULA('NumElements',11,'ElementSpacing',0.04);
%   ha.Element.FrequencyRange = [20 20000];
%   fs = 8e3; t = 0:1/fs:0.3;
%   x = chirp(t,0,1,500);
%   c = 340; % Wave propagation speed (m/s)
%   collector = phased.WidebandCollector('Sensor',ha,...
%             'PropagationSpeed',c,'SampleRate',fs,'ModulatedInput',false);
%   incidentAngle = [-50; 30];
%   x = collector(x.',incidentAngle);
%   noise = 0.2*randn(size(x));
%   rx = x+noise;
% 
%   % Frost beamforming
%   frost = phased.FrostBeamformer('SensorArray',ha,'PropagationSpeed',...
%       c,'SampleRate',fs,'Direction',incidentAngle,'FilterLength',10);
%   yf = frost(rx);
% 
%   % GSC beamforming 
%    gsc = phased.GSCBeamformer('SensorArray',ha,...
%     'PropagationSpeed',c,'SampleRate',fs,'Direction',incidentAngle,...
%     'FilterLength',10);
%   yg = gsc(rx);
% 
%   plot(t,rx(:,6),'r:',t,yf,t,yg);
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Frost','GSC');
%
%   See also phased, phased.TimeDelayBeamformer, phased.MVDRBeamformer,
%   phased.TimeDelayLCMVBeamformer, phased.SubbandPhaseShiftBeamformer.
%   phased.FrostBeamformer.

%   Copyright 2015-2016 The MathWorks, Inc.

% References
%
% [1] Griffiths, Lloyd J., and Charles W. Jim. "An alternative approach to
%     linearly constrained adaptive beamforming." Antennas and Propagation,
%     IEEE Transactions on 30.1 (1982): 27-34.
% [2] Harry Van Trees, Optimum Array Processing, Wiley, 2002
    
%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    
    properties (Nontunable, PositiveInteger)
        %FilterLength      Signal path FIR filter length
        %   Specify the length of the signal path FIR filter. The FIR
        %   filter is a delta function. The default value of this property
        %   is 1.
        FilterLength = 1;
    end
    
    properties (Nontunable)
        %LMSStepSize        Adaptive filter step size
        %   Specify the adaptive filter step size in the LMS algorithm.
        %   This value determines the LMS step size when divided by the
        %   total power in the sidelobe canceling path. The default value
        %   of this property is 0.1.
        LMSStepSize  = 0.1;
    end
    
    properties(Access = protected)
        pAfilt
        pBlockingMatrix
        pSteeringWeights
        pDataBuffer
        pSignalPathFilter
    end
    
    methods
        function obj = GSCBeamformer(varargin)
            obj@phased.internal.AbstractTimeDomainBeamformer(varargin{:});
        end
        
        function set.LMSStepSize(obj,val)
            validateattributes( val, { 'double','single'},...
              {'scalar', 'nonnegative', 'real', 'finite','nonsparse'}, ...
              '', 'LMSStepSize');
            obj.LMSStepSize = val;
        end
        
    end
    
    methods (Access = protected)
        
        function flag = isMultipleInputAnglesAllowed(obj) %#ok<MANU>
            flag = false;
        end
        
        function initializeSteeringWeights(obj)
            M = getNumElements(obj.SensorArray);
            obj.pSteeringWeights = 1/M*ones(M,1);
        end
        
        function initializeSignalPathFilter(obj)
            F = zeros(obj.FilterLength,1);           
            F(end) = 1;
            obj.pSignalPathFilter = F;
        end
        
        function initializeBlockingMatrix(obj)
            M = getNumElements(obj.SensorArray);
            % Use an orthogonal (hadamard) matrix if possible (M must be an
            % integer and a power of 2)
             M2 = nextpow2(M);
             if ~isequal(2^M2,M)
                Ws = toeplitz([1 zeros(1,M-2)],[1 -1 zeros(1,M-2)]);
             else
                OrthoMatrix = hadamard(2^M2);
                Ws = OrthoMatrix(2:end,:);
             end
            obj.pBlockingMatrix = Ws;
        end
        
        function validateInputsImpl(obj,x,varargin)
            validateInputsImpl@phased.internal.AbstractTimeDomainBeamformer(obj,x);
            if nargin == 3
                validateInputAngle(obj,varargin{1});
            end
        end
        
        function setupImpl(obj,x,varargin)
            setupImpl@phased.internal.AbstractTimeDomainBeamformer(obj,x);
            initializeSteeringWeights(obj);
            initializeBlockingMatrix(obj);
            initializeSignalPathFilter(obj);
            
            obj.pDataBuffer = zeros(obj.FilterLength-1,obj.pDOF,'like',x);
            obj.pAfilt =  zeros(obj.FilterLength,size(obj.pBlockingMatrix,1),'like',x);  
        end
        
        function y = stepImpl(obj,x,varargin)
            cond = size(x,1) < obj.FilterLength;
            if cond
                coder.internal.errorIf(cond,'phased:TimeDomainBeamformer:NotEnoughSamples', ...
                  'X', obj.FilterLength);
            end
            
            % Step Parameters
            Wc = obj.pSteeringWeights;
            Ws = obj.pBlockingMatrix;
            alpha = obj.LMSStepSize;
            L =  obj.FilterLength;
                
            % If the direction is provided in step, use it. Otherwise, use
            % the Direction property of the object.
            if nargin == 2
                ang = obj.Direction;
            else
                ang = varargin{1};
            end
            numSteerAngles = size(ang,2);
            
            % Allocate output matrices
            yc = zeros(size(x,1),numSteerAngles,'like',x);
            ya = zeros(size(x,1),numSteerAngles,'like',x);
            y  = zeros(size(x,1),numSteerAngles,'like',x);

            % Solve for the beamformer output for each direction 
            for steerIndex = 1:numSteerAngles
                % Pre-steering for this direction 
                x_presteered = steer(obj,[obj.pDataBuffer; x],ang(:,steerIndex)); 
                
                % Compute the output of the conventional beamformer
                yct = x_presteered * Wc;
                
                % Loop over each filter block, updating the adaptive filter
                % after every iteration.
                for k = 1:size(x,1) 
  
                    % Define the block
                    x_1block = x_presteered((0:L-1)+k,:);
                   
                    % Compute the output of the conventional beamformer
                    % path filter
                    yc(k,steerIndex) = sum(obj.pSignalPathFilter.*yct(k:k+L-1,steerIndex),1);
                    
                    % Compute the output of the sidelobe cancelling path
                    xnull_1block = x_1block * Ws' ;                                   % equation (26) in Ref [1]
                    
                    % Update sidelobe cancelling path weights 
                    P_xnull = sum(sum(abs(xnull_1block).^2,1),2);
                    epa = (yc(k,steerIndex)-sum(sum(obj.pAfilt.*xnull_1block,1),2));  % equation (7.426 in [2])
                    obj.pAfilt = obj.pAfilt + alpha/P_xnull*conj(epa)*xnull_1block;   % equation (7.425 in [2])
                    ya(k,steerIndex) = sum(sum(obj.pAfilt.*xnull_1block,1),2);     

                    % The beamformer output is the difference between the
                    % two paths
                    y(k,steerIndex) =  yc(k,steerIndex)  - ya(k,steerIndex);          % equation (29) in Ref [1]
                    
                end 
            end 
            if L > 1
                obj.pDataBuffer = x(end-L+2:end,1:obj.pValidatedNumInputChannels);
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractTimeDomainBeamformer(obj,s);
            if isfield(s,'isLocked')
                s = rmfield(s,'isLocked');
            end
        end
        
        function loadObjectImpl(obj,s,~)
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            if isLocked(obj)
                s.pAfilt = obj.pAfilt;
                s.pBlockingMatrix = obj.pBlockingMatrix;
                s.pSteeringWeights = obj.pSteeringWeights;
                s.pDataBuffer = obj.pDataBuffer;
                s.pSignalPathFilter = obj.pSignalPathFilter;
                s.pDataBuffer = obj.pDataBuffer;
            end
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            obj.pDataBuffer(:) = 0;
            obj.pAfilt(:) = 0;
        end
        
    end
    
    methods (Access = protected) %for Simulink
        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end
        
        function varargout = isOutputComplexImpl(obj)
            varargout = {propagatedInputComplexity(obj,1)};
        end
        function flag = isInputComplexityLockedImpl(obj,index) 
            flag = isInputComplexityLockedImpl@phased.internal.AbstractTimeDomainBeamformer(obj,index);
            if strncmpi(obj.DirectionSource,'i',1) && (index == 2)
                   flag = true;
            end
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('GSC\nBeamformer');
        end
        function varargout = getOutputNamesImpl(~)
            varargout = {'Y'};
        end
    end
    
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractTimeDomainBeamformer('array');
            props = {...
                'SampleRateFromInputCheckbox',...
                'SampleRate',...
                'FilterLength',...
                'LMSStepSize',...
                'Direction',...
                'DirectionSource'
                };
            groups(1).PropertyList = [groups(1).PropertyList props];
            
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:GeneralizedSidelobeCancellerTitle')),...
                'Text',getString(message('phased:library:block:GeneralizedSidelobeCancellerDesc')));
        end    
    end
    
end

  
