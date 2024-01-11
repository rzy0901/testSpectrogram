classdef (Sealed,StrictDefaults) SubbandMVDRBeamformer < phased.internal.AbstractSubbandBeamformer
%SubbandMVDRBeamformer  Subband MVDR beamformer
%   H = phased.SubbandMVDRBeamformer creates a subband MVDR beamformer
%   System object, H. This object performs subband MVDR beamforming on the
%   received signal.
%
%   H = phased.SubbandMVDRBeamformer(Name,Value) creates a subband MVDR
%   beamformer object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The subband MVDR beamformer separates the signal into several subbands
%   and applies the narrowband MVDR beamforming to the signal in each
%   subband. The beamformed signals in each subband are added to form the
%   output signal.
%
%   Step method syntax:
%
%   Y = step(H,X) performs subband MVDR beamforming on the input X, and
%   returns the beamformed output in Y. X is an MxN matrix where N is the
%   number of subarrays if SensorArray contains subarrays, or the number of
%   elements otherwise. Y is an MxL matrix where L is the number of
%   beamforming directions.
%
%   If you set the TrainingInputPort to false, then X is used to do the
%   training and M must be larger than N*NB, where NB is the number of
%   subbands. If you set the TrainingInputPort to true, then M can be any
%   positive integer.
%
%   Y = step(H,X,XT) uses XT as the training samples to calculate the
%   beamforming weights when you set the TrainingInputPort property to
%   true. XT is a PxN matrix where P must be larger than N*NB.
%
%   Y = step(H,X,ANG) uses ANG as the beamforming direction, when you
%   set the DirectionSource property to 'Input port'. ANG is a 2-row
%   matrix whose columns are in the form of [AzimuthAngle; ElevationAngle]
%   (in degrees). The azimuth angle must be between [-180 180] and the
%   elevation angle must be between [-90 90].
%
%   [Y,W] = step(...) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is an
%   NxKxL matrix where K is the number of subbands specified in the
%   NumSubbands property. Each column of W specifies the narrowband
%   beamforming weights used in the corresponding subband for the
%   corresponding direction.
%
%   [Y,FREQ] = step(...) returns additional output FREQ as the center
%   frequencies of subbands when you set the SubbandsOutputPort property to
%   true. FREQ is a length-K column vector where K is the number of
%   subbands specified in the NumSubbands property.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [Y,W,FREQ] = step(H,X,XT,ANG)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   SubbandMVDRBeamformer methods:
%
%   step     - Beamforming using subband phase shifting (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create subband MVDR beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   SubbandMVDRBeamformer properties:
%
%   SensorArray           - Sensor array
%   PropagationSpeed      - Signal propagation speed
%   OperatingFrequency    - Operating frequency
%   SampleRate            - Sample rate
%   NumSubbands           - Number of subbands
%   DirectionSource       - Source of beamforming direction
%   Direction             - Beamforming direction
%   DiagonalLoadingFactor - Diagonal loading factor
%   TrainingInputPort     - Enable training data input
%   WeightsOutputPort     - Enable weights output
%   SubbandsOutputPort    - Enable subband center frequencies output 
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Apply subband MVDR beamformer to an 11-element ULA. The incident 
%   %   angle of the signal is 10 degrees in azimuth and 30 degrees in 
%   %   elevation.
%
%   % signal simulation
%   mic = phased.ULA('NumElements',11,'ElementSpacing',0.3);
%   fs = 1e3; carrierFreq = 2e3; t = (0:1/fs:2)';
%   x = chirp(t,0,2,fs);
%   c = 1500; % Wave propagation speed (m/s)
%   sigcol = phased.WidebandCollector('Sensor',mic,'PropagationSpeed',c,...
%       'SampleRate',fs,'ModulatedInput',true,...
%       'CarrierFrequency',carrierFreq);
%   incidentAngle = [10;0];
%   x = sigcol(x,incidentAngle);
%   noise = 0.3*(randn(size(x)) + 1j*randn(size(x)));
%   rx = x+noise;
% 
%   % beamforming
%   bf = phased.SubbandMVDRBeamformer('SensorArray',mic,...
%       'Direction',incidentAngle,'OperatingFrequency',carrierFreq,...
%       'PropagationSpeed',c,'SampleRate',fs,'TrainingInputPort',true, ...
%       'SubbandsOutputPort',true, 'WeightsOutputPort',true);
%   [y,w,subbandfreq] = bf(rx,noise);
%   plot(t(1:300),real(rx(1:300,6)),'r:',t(1:300),real(y(1:300)));
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
% 
%   % plot response pattern for 5 bands
%   figure;
%   pattern(mic,subbandfreq(1:5).',-180:180,0,...
%       'PropagationSpeed',c,'Weights',w(:,1:5));
%
%   See also phased, phased.MVDRBeamformer, phased.PhaseShiftBeamformer,
%   phased.SubbandPhaseShiftBeamformer, phased.WidebandCollector.

%   Copyright 2015-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties
        %DiagonalLoadingFactor  Diagonal loading factor
        %   Specify the diagonal loading factor as a positive scalar. The
        %   default value of this property is 0. Diagonal loading is a
        %   technique used to achieve robust beamforming performance,
        %   especially when the sample support is small. This property is
        %   tunable.
        DiagonalLoadingFactor = 0;
    end

    properties (Nontunable, Logical) 
        %TrainingInputPort  Enable training data input
        %   Set this property to true to add input to specify the
        %   additional training data. Set this property to false to not add
        %   input to specify the training data. The default value of this
        %   property is false. When this property is false, the input
        %   signal itself is used as the training data.
        TrainingInputPort = false;
    end

    properties (Access = private, Nontunable)
        cSteeringVector;
        pNumAngles;
        cIFFT
    end

    methods
        function set.DiagonalLoadingFactor(obj,val)
            validateattributes( val, { 'double','single' }, ...
                { 'scalar', 'nonnegative', 'real', 'finite', 'nonempty' },...
                '', 'DiagonalLoadingFactor');
            obj.DiagonalLoadingFactor = val;
        end
    end

    methods

        function obj = SubbandMVDRBeamformer(varargin)
            obj@phased.internal.AbstractSubbandBeamformer(varargin{:});
        end

    end

    methods (Access = protected)
        
        function flag = isMultipleInputAnglesAllowed(obj) %#ok<MANU>
            flag = true;
        end
        
        function flag = isSubarraySupported(obj) %#ok<MANU>
            flag = true;
        end
        
        function num = getNumInputsImpl(obj)
            num = getNumInputsImpl@phased.internal.AbstractSubbandBeamformer(obj);
            if obj.TrainingInputPort
                num = num + 1;
            end
        end
        
        function validateInputsImpl(obj,x,varargin)
            validateInputsImpl@phased.internal.AbstractSubbandBeamformer(obj,x);
            if obj.TrainingInputPort
                xt = varargin{1};                
                validateInputSignal(obj,xt,'XT');
                
                if strncmpi(obj.DirectionSource,'i',1)
                   ang = varargin{2};
                   validateInputAngle(obj,ang);
                end
                
            else
                if strncmpi(obj.DirectionSource,'i',1)
                    ang = varargin{1};
                    validateInputAngle(obj,ang);
                end
                xt = x;
            end

            sz_xt = size(xt);
            if sz_xt(1) < obj.NumSubbands*sz_xt(2) && isempty(coder.target)
                warning(message('phased:beamformer:SMI:InsufficientTrainingSample'));
            end

            
        end      

        function setupImpl(obj,x,xt,ang)
            setupImpl@phased.internal.AbstractSubbandBeamformer(obj,x);
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',cast(obj.PropagationSpeed,'double'));
            
            if (obj.DirectionSource(1) == 'P') %Property
                obj.pNumAngles = size(obj.Direction,2);
            else
                if obj.TrainingInputPort
                    obj.pNumAngles = size(ang,2);
                else
                    obj.pNumAngles = size(xt,2);
                end
            end
            obj.cIFFT = dsp.IFFT;
        end
        
        function flag = isInputSizeLockedImpl(obj,index) 
            if index == 1
                flag = false;
            elseif index == 2 && obj.TrainingInputPort
                flag = false;
            else
                flag = true;
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = isInputComplexityLockedImpl@phased.internal.AbstractSubbandBeamformer(obj,index);
            if (obj.DirectionSource(1) == 'I')  %Input port
                if obj.TrainingInputPort && (index == 3) 
                    flag = true;
                elseif ~obj.TrainingInputPort && (index == 2)
                    flag = true;
                end
            end
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSubbandBeamformer(obj);
            release(obj.cSteeringVector);
            release(obj.cIFFT);
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractSubbandBeamformer(obj);
            reset(obj.cSteeringVector);
            reset(obj.cIFFT);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSubbandBeamformer(obj);
            if isLocked(obj)
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.cIFFT = saveobj(obj.cIFFT);
                s.pNumAngles = obj.pNumAngles;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSubbandBeamformer(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
                    if isfield(s,'cIFFT')
                        obj.cIFFT = dsp.IFFT.loadobj(s.cIFFT);
                        s = rmfield(s,'cIFFT');
                    end
                end
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
        
        function [y,wOut,f] = stepImpl(obj,x_in,varargin)
            classtouse = class(x_in);
            if ~obj.pSizeInitialized
                computeRemainingSamples(obj,size(x_in,1),obj.NumSubbands);
                obj.pSizeInitialized = true;
            end
            
            if obj.TrainingInputPort
                if (obj.DirectionSource(1) == 'P') %Property
                    ang = obj.Direction;
                else
                    ang = varargin{2};
                end
                x = formatInput(obj,x_in,'S');
                xf = subbandFFT(obj,x,x_in,'S');
                xt_in = cast(varargin{1},classtouse);
                xt = formatInput(obj,xt_in,'T');
                xtf = subbandFFT(obj,xt,xt_in,'T');
            else
                if (obj.DirectionSource(1) == 'P') %Property
                    ang = obj.Direction;
                else
                    ang = varargin{1};
                end
                x = formatInput(obj,x_in,'S');
                xf = subbandFFT(obj,x,x_in,'T');
                xtf = xf;
            end
            w = cast(calcWeights(obj,xtf,ang),classtouse);
            
            f = cast(obj.pSubbandCenterFreq,classtouse);
            Nbands = obj.NumSubbands;
            
            y = complex(zeros(size(x_in,1),obj.pNumAngles,classtouse));
            for m = obj.pNumAngles:-1:1
                outputFreqDomain = complex(zeros(size(xf,1),Nbands));
                for freqIdx = 1:Nbands
                    outputFreqDomain(:,freqIdx) = ...
                        xf(:,:,freqIdx)*conj(w(:,freqIdx,m));
                end

                % convert back to time domain
                
                % ensure input to FFT always complex 
                y_fft_in = complex(outputFreqDomain);
                
                if size(y_fft_in,2) == 1
                    % each band has only one sample
                    ytemp = y_fft_in.';  % scalar ifft is itself
                else
                    ytemp = step(obj.cIFFT,y_fft_in.');
                end

                % restore original length column vector
                y(:,m) = formatOutput(obj,ytemp(:));

            end
            
            if ~obj.WeightsOutputPort && obj.SubbandsOutputPort
                wOut = f;  % second input is subbands center frequencies
            else
                wOut = w;
            end

        end

    end

    methods (Access = private)

        function w_mat = calcWeights(obj,xt,ang)
            % the steering vector matrix
            Nbands = obj.NumSubbands;
            Nang = obj.pNumAngles;
            w_mat = complex(zeros(obj.pDOF,Nbands,Nang));
            subbandfreq = obj.pSubbandCenterFreq;
            hstv = obj.cSteeringVector;
            for freqIdx = 1:Nbands
                sv = step(hstv,cast(subbandfreq(freqIdx),'double'),cast(ang,'double'));
                for m = Nang:-1:1
                    w_mat(:,freqIdx,m) = phased.internal.lcmvweights(...
                        xt(:,:,freqIdx),sv(:,m),1,obj.DiagonalLoadingFactor);
                end
            end
        end
    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractSubbandBeamformer('subarray');
        props = {'OperatingFrequency',...
                 'SampleRateFromInputCheckbox',...
                 'SampleRate',...
                 'NumSubbands',...
                 'DiagonalLoadingFactor',...
                 'TrainingInputPort',...
                 'DirectionSource',...
                 'Direction',...
                 'WeightsOutputPort',...
                 'SubbandsOutputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:SubbandMVDRBeamformerTitle')),...
              'Text',getString(message('phased:library:block:SubbandMVDRBeamformerDesc')));
      end            
    end
    methods (Access = protected) %for Simulink    
        function varargout = getInputNamesImpl(obj)
            %Insert XT
            if obj.TrainingInputPort
                [varargout{1:nargout-1}] = getInputNamesImpl@phased.internal.AbstractSubbandBeamformer(obj);
                varargout = {varargout{1} 'XT' varargout{2:end}};
            else
                [varargout{1:nargout}] = getInputNamesImpl@phased.internal.AbstractSubbandBeamformer(obj);
            end
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Subband MVDR\nBeamformer');
        end        
    end    
end

