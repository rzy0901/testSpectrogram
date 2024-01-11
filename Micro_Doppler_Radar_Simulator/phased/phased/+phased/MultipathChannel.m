classdef (Sealed,StrictDefaults) MultipathChannel < phased.internal.AbstractMultipathChannel
%MultipathChannel  Multipath propagation channel
%   H = phased.MultipathChannel creates a multipath channel System object,
%   H. This object simulates signal propagation of multiple paths.
%
%   The object applies path-dependent time delays, losses, and incorporates
%   reflection coefficients for the channel interfaces. Losses include
%   spreading and absorption. 
%
%   H = phased.MultipathChannel (Name,Value) returns a multipath channel,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,PATHS,DOP,ALOSS) returns the resulting signal Y when the
%   narrowband signal X propagates along the paths defined by PATHS with
%   Doppler factor, DOP, and absorption losses, ALOSS. X has N columns,
%   where N is the number of multipath arrivals.
%
%   DOP is an N element row vector, where N is the number of paths. Each
%   element of DOP contains the factor which multiplies the emitted
%   frequency to produce the observed, doppler-shifted frequency.
%
%   ALOSS is an M by N+1 matrix, where M is the number of frequencies and N
%   is the number of paths. The first column of M contains frequencies in
%   Hz and the remaining columns contain absorption loss for each path in
%   dB.
%
%   PATHS is a 3 by N matrix, where N is the number of paths. The first row
%   of PATHS contains propagation time delays (in seconds), the second
%   contains the total reflection coefficient for each path due to
%   interface reflections, and the third contains the spreading loss for
%   each path in dB.
%
%   The output Y represents the signals arriving at the propagation
%   destinations within the current time frame, which is the time occupied
%   by the current input. If it takes longer than the current time frame
%   for the signals to propagate from the origin to the destination, then
%   the output contains no contribution from the input of the current time
%   frame. For a single path, the output Y can be written as
%
%   Y(t) = X(t-tau)/L
%
%   where tau is the time delay and L is the propagation loss. The
%   propagation loss includes channel path spreading loss as well as the
%   loss due to absorption in the medium.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   MultipathChannel methods:
%
%   step     - Propagate signal from one location to another (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create multipath channel object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the object
%
%   MultipathChannel properties:
%
%   OperatingFrequency       - Signal carrier frequency
%   SampleRate               - Channel sample rate
%   MaximumDelaySource       - Source of maximum propagation time delay
%   MaximumDelay             - Maximum propagation time delay
%   InterpolationMethod      - Interpolation method 
%
%   % Example:
%   %   Propagate a signal between a stationary source and 
%   %   receiver for a 200 m deep channel. Plot the sum of the received
%   %   signals.
%
%   proppaths = phased.IsoSpeedUnderwaterPaths('ChannelDepth',200);
%   channel = phased.MultipathChannel;
%   [paths,dop,aloss] = proppaths([0;0;-100],[100;0;-80],...
%     [0;0;0],[0;0;0],1);
%   y = channel(ones(500,51),paths,dop,aloss);   
%   figure;
%   plot((0:500-1)/1e3,abs(sum(y,2)));
%   xlabel('Time (s)')
%   ylabel('Received signal')
%
%   See also phased.IsotropicHydrophone, phased.IsotropicProjector, 
%   phased.BackscatterSonarTarget, phased.IsoSpeedUnderwaterPaths

%   Copyright 2016 The MathWorks, Inc.

%   Reference
%   [1] Urick, Principles of Underwater Sound, Peninsula Publishing, 1983

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Access = private, Nontunable)
        % Fractional delay filter
        cFractionalDelayFilter
        % Noise source delay
        pNoiseSourceDelay
        % Resample factor
        pResampleFactor
    end
    
    properties (Access = private)
        pLastSample
    end
    
    methods
        function obj = MultipathChannel(varargin)
            obj = obj@phased.internal.AbstractMultipathChannel(varargin{:});
        end
    end
 
    methods (Access = protected)
        function setupImpl(obj,x,Paths,~,~,~)
            setupImpl@phased.internal.AbstractMultipathChannel(obj,x);
            obj.cFractionalDelayFilter = dsp.VariableFractionalDelay;
            obj.pUsedPaths = ~isnan(Paths(1,:));
            obj.pLastSample = nan(2,size(x,2));
            if isequal(obj.InterpolationMethod,'Oversample')
              obj.pResampleFactor = 15;
            else
              obj.pResampleFactor = 1;
            end
          end
        
        function y = stepImpl(obj,x_in,paths,dop,aloss)
            validatePaths(obj,paths,dop,aloss);
            
            % Compute propgated signal
            [y_temp,bufferDelay] = propagatedSignal(obj,x_in,paths,dop,aloss);
   
            % Store signal in the circular buffer and retrieve output for
            % the current timestep.
            y = step(obj.cBuffer,y_temp,bufferDelay,size(x_in,1));
        end
    end
    
    methods (Access = protected)
        function num = getNumInputsImpl(obj)   %#ok<MANU>
            num = 4;
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractMultipathChannel(obj);
            release(obj.cFractionalDelayFilter);
        end

        function resetImpl(obj)
            resetImpl@phased.internal.AbstractMultipathChannel(obj);
            reset(obj.cFractionalDelayFilter);
            obj.pLastSample(:) = nan;
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractMultipathChannel(obj);
            if isLocked(obj)
                s.cFractionalDelayFilter = saveobj(obj.cFractionalDelayFilter);
                s.pNoiseSourceDelay = obj.pNoiseSourceDelay;
                s.pLastSample = obj.pLastSample;
                s.pResampleFactor = obj.pResampleFactor;
            end
        end

        function s = loadSubObjects(obj,s,wasLocked)
            s = loadSubObjects@phased.internal.AbstractMultipathChannel(obj,s,wasLocked);
            if wasLocked 
                obj.cFractionalDelayFilter = ...
                    dsp.VariableFractionalDelay.loadobj(s.cFractionalDelayFilter);
                s = rmfield(s,'cFractionalDelayFilter');
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Multipath Channel');
        end
        
        function validatePaths(obj,paths,dop,aloss)
            % Check that PATHS any nan's, all values in that column
            % are nan.
            cond =  ~isequal(any(isnan(paths)),all(isnan(paths)));
            if cond
                coder.internal.errorIf(cond,...
                    'phased:MultipathChannel:AllorNoneNaN','PATHS');
            end
            
            % Check that delays and path loss are positive and reflection coefficients are
            % between -1 and 1.
            cond = any(paths(1,:) < 0) || any(paths(3,:) < 0) || ...
              any(paths(2,:) > 1) || any(paths(2,:) < -1);
            
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:MultipathChannel:InvalidPaths');               
            end
            
            % Check that operating frequency is contained within range of
            % the first column of aloss
            cond = obj.OperatingFrequency > max(aloss(:,1)) | ...
              obj.OperatingFrequency < min(aloss(:,1));
            if cond
                coder.internal.errorIf(cond,...
                    'phased:MultipathChannel:OperatingFreqALOSS','OperatingFrequency','ALOSS');
            end
            
            % Check that loss matrix is nonnegative
            cond = any(any(aloss < 0));
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:expectedNonnegative','ALOSS')
            end
            
            % Check that if ALOSS has any nan's, all values in that column
            % are nan.
            cond =  ~isequal(any(isnan(aloss)),all(isnan(aloss)));
            if cond
                coder.internal.errorIf(cond,...
                    'phased:MultipathChannel:AllorNoneNaN','ALOSS');
            end
            
            % Check that dop is positive
            cond = any(dop <= 0);
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:expectedPositive','DOP')
            end
        end
        
        function [y_temp,bufferDelay] = propagatedSignal(obj,x_in,paths,dop,aloss)
            coder.extrinsic('num2str');  
                 
            sx_in = size(x_in,1);
            numPaths = size(x_in,2);
            pathDelay = paths(1,:);
            
            % Check that radius doesn't go to zero over this prf
            cond = any((dop-1)*sx_in/obj.SampleRate > pathDelay);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:MultipathChannel:PathDelayZero');
            end
            
            % Check that the maximum doppler frequency is less than half
            % the sample rate.
            if isempty(coder.target)	 
                if max(abs((dop-1)*obj.OperatingFrequency)) > obj.SampleRate/2
                    warning(message('phased:MultipathChannel:SampleRateBelowDoppler','SampleRate',...
                      coder.const(num2str(2*max(abs((dop-1)*obj.OperatingFrequency))))));
                end
            end
            
            % Reflection coefficients
            reflCoef = paths(2,:);

            % Spreading
            spreadFactor = sqrt(db2pow(-1*(paths(3,:))));
            
            % Absorption 
            if size(aloss,1) > 1
                absorpFactor = sqrt(db2pow(-1*...
                  interp1(aloss(:,1),aloss(:,2:end),obj.OperatingFrequency,'Linear'))); 
            else
                absorpFactor = sqrt(db2pow(-1*aloss(2:end)));
            end
            
            phaseShift = exp(-1i*2*pi*obj.OperatingFrequency*paths(1,:));
            
            % Perform time dilation/compression on the input signal by
            % defining a dilated/compressed receiver grid.
            if any(isnan(obj.pLastSample(2,:)))
                % We don't have a previous input signal. We will set the
                % value for the grid point immediately prior to this signal
                % to 0.
                rcvGrid = repmat(pathDelay*obj.SampleRate,sx_in+1,1) +...
                  bsxfun(@times,(-1:sx_in-1)',1./dop); 
                obj.pLastSample(2,:) = 0;
            else
                rcvGrid = [obj.pLastSample(1,:);repmat(pathDelay*obj.SampleRate,sx_in,1)+ ...
                  bsxfun(@times,(0:sx_in-1)',1./dop)];
            end
            
            % Interpolate the dilated/compressed input signal onto an
            % (integer) sample grid, rcvSampleGrid. This introduces a
            % fractional delay.
            x_rcv = zeros(max(ceil(rcvGrid(end,:)-rcvGrid(1,:))+1),size(x_in,2));
            sx_in_dop = size(x_rcv,1);
            % Delay for each path in samples for placement in buffer.
            bufferDelay = zeros(1,size(x_in,2));
            % Offset in samples from bufferDelay
            % to the (possibly non-integer) path delay in samples.
            DOPoffset = zeros(1,numPaths);
            for i = 1:numPaths
                if obj.pUsedPaths(i)
                    % Define the integer sample grid for this signal,and
                    % interpolate from dilated/compressed grid onto this
                    % grid. Use floor(x)+1 instead of ceil(x) in defining
                    % the grid in case of an integer to avoid overlapping
                    % endpoints between step calls.
                    rcvSampleGrid = floor(rcvGrid(1,i))+1:floor(rcvGrid(end,i)); 
                    bufferDelay(i) = rcvSampleGrid(1);
                    DOPoffset(i) = rcvSampleGrid(1)-pathDelay(i)*obj.SampleRate; 
                    % Create the interpolation grid and values. If
                    % BandlimitedInput is specified, resample the signal.
                    interpGrid = [rcvGrid(1,i) ...
                      rcvGrid(2,i):1/dop(i)/obj.pResampleFactor:rcvGrid(end,i)+(1-1/obj.pResampleFactor)/dop(i)];
                    interpValues = [obj.pLastSample(2,i); resample(x_in(:,i),obj.pResampleFactor,1,100,0.5)];
                    % Interpolate onto the receiver sample grid.
                    x_rcv(1:length(rcvSampleGrid),i) = ...
                      interp1(interpGrid,interpValues,rcvSampleGrid,'linear');
                end
            end

            % Store the last receiver grid sample for next step call.
            % Translate the receiver grid to the grid for the next step
            % call.
            obj.pLastSample(1,:) = rcvGrid(end,:)-repmat(size(x_in,1),1,numPaths); 
            obj.pLastSample(2,:) = x_in(end,:);

            % Compute the propagated signal
            y_temp = x_rcv.*...
              repmat(phaseShift.*spreadFactor.*absorpFactor.*reflCoef,sx_in_dop,1);

            % Doppler compensation. 
            dopTimeGrid = bsxfun(@plus,(0:sx_in_dop-1)',DOPoffset)/obj.pSampleRate;
            y_temp=y_temp.*...
              exp(1i*2*pi*repmat(bsxfun(@minus,dop,1),sx_in_dop,1)*obj.OperatingFrequency.*dopTimeGrid);

            % Set unused paths to zero.
            y_temp(:,~obj.pUsedPaths) = 0;
        end
    end
end