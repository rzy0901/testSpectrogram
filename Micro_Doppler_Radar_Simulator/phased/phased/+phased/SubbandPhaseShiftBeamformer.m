classdef (Sealed,StrictDefaults) SubbandPhaseShiftBeamformer < phased.internal.AbstractSubbandBeamformer
%SubbandPhaseShiftBeamformer  Subband phase shift beamformer
%   H = phased.SubbandPhaseShiftBeamformer creates a subband phase shift
%   beamformer System object, H. This object performs subband phase shift
%   beamforming on the received signal.
%
%   H = phased.SubbandPhaseShiftBeamformer(Name,Value) creates a subband
%   phase shift beamformer object, H, with the specified property Name set
%   to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   The subband phase shift beamformer separates the signal into several
%   subbands and applies the narrowband phase shift beamforming to the
%   signal in each subband. The beamformed signals in each subband are
%   added to form the output signal.
%
%   Step method syntax:
%
%   Y = step(H,X) performs subband phase shift beamforming on the input X,
%   and returns the beamformed output in Y. X is an MxN matrix where N is
%   the number of subarrays if SensorArray contains subarrays, or the
%   number of elements otherwise. Y is an MxL matrix where L is the number
%   of beamforming directions.
%
%   [Y,W] = step(H,X) returns additional output W as the beamforming
%   weights when you set the WeightsOutputPort property to true. W is an
%   NxKxL matrix where K is the number of subbands specified in the
%   NumSubbands property. Each column of W specifies the narrowband
%   beamforming weights used in the corresponding subband for the
%   corresponding direction.
%
%   [Y,FREQ] = step(H,X) returns additional output FREQ as the center
%   frequencies of subbands when you set the SubbandsOutputPort property to
%   true. FREQ is a length-K column vector where K is the number of
%   subbands specified in the NumSubbands property.
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
%   [Y,W,FREQ] = step(H,X,ANG)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   SubbandPhaseShiftBeamformer methods:
%
%   step     - Beamforming using subband phase shifting (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create phase shift beamformer object with same property
%              values
%   isLocked - Locked status (logical)
%
%   SubbandPhaseShiftBeamformer properties:
%
%   SensorArray        - Sensor array
%   PropagationSpeed   - Signal propagation speed
%   OperatingFrequency - Operating frequency
%   SampleRate         - Sample rate
%   NumSubbands        - Number of subbands
%   DirectionSource    - Source of beamforming direction
%   Direction          - Beamforming direction
%   WeightsOutputPort  - Enable weights output
%   SubbandsOutputPort - Enable subband center frequencies output 
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Apply subband phase shift beamformer to an 11-element ULA. 
%   %   The incident angle of the signal is 10 degrees in azimuth and
%   %   30 degrees in elevation.
%
%   % signal simulation
%   array = phased.ULA('NumElements',11,'ElementSpacing',0.3);
%   array.Element.FrequencyRange = [20 20000];
%   fs = 1e3; carrierFreq = 2e3; t = (0:1/fs:2)';
%   x = chirp(t,0,2,fs);
%   c = 1500; % Wave propagation speed (m/s)
%   sigcol = phased.WidebandCollector('Sensor',array,...
%               'SampleRate',fs,'ModulatedInput',true,...
%               'PropagationSpeed',c,'CarrierFrequency',carrierFreq);
%   incidentAngle = [10; 30];
%   x = sigcol(x,incidentAngle);
%   noise = 0.3*(randn(size(x)) + 1j*randn(size(x))); 
%   rx = x+noise;
%
%   % beamforming
%   bf = phased.SubbandPhaseShiftBeamformer('SensorArray',array,...
%        'Direction',incidentAngle,'OperatingFrequency',carrierFreq,...
%        'PropagationSpeed',c,'SampleRate',fs, ...
%        'SubbandsOutputPort',true, 'WeightsOutputPort',true);
%   [y,w,subbandfreq] = bf(rx);
%   plot(t(1:300),real(rx(1:300,6)),'r:',t(1:300),real(y(1:300)));
%   xlabel('Time'),ylabel('Amplitude'),legend('Original','Beamformed');
%
%   % plot response pattern for 5 bands
%   figure; 
%   pattern(array,subbandfreq(1:5).',-180:180,0,...
%       'PropagationSpeed',c,'Weights',w(:,1:5));
%
%   See also phased, phased.PhaseShiftBeamformer,
%   phased.TimeDelayBeamformer, phased.Collector, phased.WidebandCollector.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Access = private, Nontunable)
        cSteeringVector;
        pNumAngles;
        cIFFT
    end
    properties(Access = private) 
        %pre-calculated weights
        pWeights;
    end

    methods

        function obj = SubbandPhaseShiftBeamformer(varargin)
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
        
        function validateInputsImpl(obj,x,varargin)
            validateInputsImpl@phased.internal.AbstractSubbandBeamformer(obj,x);
            if (obj.DirectionSource(1) == 'I') %Input port
                ang = varargin{1};
                validateInputAngle(obj,ang);
            end
        end      

        function setupImpl(obj,x,ang)
            setupImpl@phased.internal.AbstractSubbandBeamformer(obj,x);
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed);
            
            if (obj.DirectionSource(1) == 'P') %Property
                obj.pNumAngles = size(obj.Direction,2);
                obj.pWeights = calcWeights(obj,obj.Direction);
            else
                obj.pNumAngles = size(ang,2);
            end
            obj.cIFFT = dsp.IFFT;
        end
        
        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = isInputComplexityLockedImpl@phased.internal.AbstractSubbandBeamformer(obj,index);
            if (obj.DirectionSource(1) == 'I')  && (index == 2) %Input port
                flag = true;
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
                s.pWeights = obj.pWeights;
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
            
            classtouse=class(x_in);
            if ~obj.pSizeInitialized
                computeRemainingSamples(obj,size(x_in,1),obj.NumSubbands);
                obj.pSizeInitialized = true;
            end
            
            if (obj.DirectionSource(1) == 'P') %Property
                w = cast(obj.pWeights,classtouse);
            else
                ang = varargin{1};
                validateAngleRange(obj,ang);
                w = cast(calcWeights(obj,ang),classtouse);
            end
            
            f = cast(obj.pSubbandCenterFreq,classtouse);
            Nbands = obj.NumSubbands;
            
            x = formatInput(obj,x_in,'S');
            y = complex(zeros(size(x_in,1),obj.pNumAngles,classtouse));
            xf = subbandFFT(obj,x,x_in,'S');
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

        function w_mat = calcWeights(obj,ang)
            % the steering vector matrix
            Nbands = obj.NumSubbands;
            Nang = obj.pNumAngles;
            w_mat = complex(zeros(obj.pDOF,Nbands,Nang));
            subbandfreq = obj.pSubbandCenterFreq;
            hstv = obj.cSteeringVector;
            for freqIdx = 1:Nbands
                sv = step(hstv,cast(subbandfreq(freqIdx),'double'),cast(ang,'double'));
                for m = Nang:-1:1
                    svtemp = sv(:,m);
                    w_mat(:,freqIdx,m) = svtemp/real(svtemp'*svtemp);
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
                 'DirectionSource',...
                 'Direction',...
                 'WeightsOutputPort',...
                 'SubbandsOutputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:SubbandPhaseShiftBeamformerTitle')),...
              'Text',getString(message('phased:library:block:SubbandPhaseShiftBeamformerDesc')));
      end            
    end
    methods (Access = protected) %for Simulink    
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Subband\nPhase Shift\nBeamformer');
        end        
    end    
end

