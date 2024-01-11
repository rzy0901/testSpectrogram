classdef (Sealed,StrictDefaults) GCCEstimator < phased.internal.AbstractArrayOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime
%GCCEstimator GCC direction of arrival (DOA) estimator
%   H = phased.GCCEstimator creates a GCC direction of arrival estimator
%   System object, H. This object estimates the direction of arrival or
%   time of arrival between sensor array elements using GCC-PHAT algorithm.
%
%   H = phased.GCCEstimator(Name,Value) returns a GCC direction of arrival
%   estimator object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   ANG = step(H,X) returns the direction of arrival, ANG, of the signal X.
%   X is a matrix specifying the received signal at the
%   array, as specified in the SensorArray property of H. The rows of X
%   represent the number of snapshots and the number of columns in X must
%   match the number of elements in the array if an array is used or the
%   number of subarrays if a subarray is used. It is assumed that X
%   contains the signal from a single source.
%
%   If the array is a uniform linear array, ANG is a scalar representing
%   the direction of arrival in terms of a broadside angle (in degrees).
%   Otherwise, ANG is a 2-element column in the form of [azimuth;
%   elevation] angles (in degrees).
%
%   The computation process assumes the single source is in the far field
%   so that the direction of arrival is the same for all sensors. The
%   computation first estimates the correlations using the GCC-PHAT
%   algorithm between specified sensor pairs and then finds the largest
%   peaks in these correlations to identify the delays between each pair of
%   sensors. A least square estimate is used to compute the direction of
%   arrival, ANG from the estimated delays.
%
%   [ANG,TAU] = step(H,X) returns TAU as the delays estimated from the
%   correlations between each pair of sensors. This output is enabled when
%   you set the DelayOutputPort property to true. TAU is a P-element row
%   vector where P is the number of sensor pairs.
%
%   [ANG,R,LAG] = step(H,X) returns R as the estimated correlations between
%   each pair of sensors. This output is enabled when you set the
%   CorrelationOutputPort property to true. R is a matrix with P columns
%   where P is the number of sensor pairs. Each column in P contains the
%   correlation for the corresponding pair of sensors. LAG is a column
%   vector containing the time lags for each sample in the correlation. The
%   time lags are the same for all correlations.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [ANG,TAU,R,LAG] = step(H,X)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   GCCEstimator methods:
%
%   step     - Estimate direction of arrival using GCC-PHAT (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create GCC estimator object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset internal states of the GCC estimator
%
%   GCCEstimator properties:
%
%   SensorArray           - Sensor array
%   PropagationSpeed      - Signal propagation speed
%   SampleRate            - Sample rate 
%   SensorPairSource      - Source of sensor pairs
%   SensorPair            - Sensor pairs
%   DelayOutputPort       - Enable delay output
%   CorrelationOutputPort - Enable correlation output
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Estimate the direction of arrival using GCC-PHAT with an 8-element
%   %   microphone array.
%
%   load chirp;         % load signal
% 
%   c = 340;
%   N = 8;
%   d = 0.5;
%   mic = phased.ULA(N,d,...
%     'Element',phased.OmnidirectionalMicrophoneElement);
% 
%   % simulate signal
%   incident_ang = 17;   % incoming direction
%   rxmic = phased.WidebandCollector('Sensor',mic,...
%       'PropagationSpeed',c,'SampleRate',Fs,...
%       'ModulatedInput',false);
%   x = rxmic(y,incident_ang);
% 
%   % estimate direction of arrival
%   gxcorr = phased.GCCEstimator('SensorArray',mic,...
%       'PropagationSpeed',c,'SampleRate',Fs,...
%       'SensorPairSource','Property','SensorPair',[1:7;2:8]);
%   ang = gxcorr(x)
%
%   See also phased, phased.BeamscanEstimator, phased.RootMUSICEstimator,
%   gccphat.

%   Copyright 2015-2018 The MathWorks, Inc.

%   Reference
%   [1] Charles Knapp and G. Carter, The Generalized Correlation Method for
%       Estimation of Time Delay, IEEE Trans. ASSP, Vol. 24, No. 4, 1976
%   [2] Jean-Marc Valin et al. Robust Sound Source Localization Using a
%       Microphone Array on a Mobile Robot, Proc. IEEE/RSJ International
%       Conference on Intelligent Robot and Systems, pp. 1228-1233, 2003
%   [3] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from Simulink time engine. Set SampleRateFromInputCheckbox to
        %   false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true
    end
    
    properties (Nontunable)
        %SampleRate     Sample rate (Hz)
        %   Specify the sample rate (in Hz) as a positive scalar. The
        %   default value is 1e6 (1 MHz).
        SampleRate = 1e6
        %SensorPairSource   Source of sensor pairs
        %   Specify how the sensor pairs are specified as one of 'Auto' |
        %   'Property', where the default is 'Auto'. When you set this
        %   property to 'Auto', the correlation is computed between the
        %   first element and rest of elements, where the first element
        %   serves as the reference channel. When you set this property to
        %   'Property', the sensor pairs for computing correlation is
        %   specified in the SensorPair property.
        SensorPairSource = 'Auto'
        %SensorPair     Sensor pairs
        %   Specify the pairs of sensors used in correlation computation as
        %   a 2-row matrix. Each column of this matrix specifies a pair of
        %   sensors between which the correlation is computed. The first
        %   row specifies the target sensor and the second row specifies
        %   the reference sensor. The default value of this property is
        %   [2;1], which means the correlation is computed for the second
        %   element in the array relative to the first element.
        SensorPair = [2;1]
    end
    
    properties (Nontunable, Logical)
        %DelayOutputPort  Enable delay output
        %   Set this property to true to output the delay corresponding to
        %   the angle between each sensor pairs. Set this property to false
        %   to not output the delays. The default value of this property is
        %   false.
        DelayOutputPort = false
        %CorrelationOutputPort  Enable correlation output
        %   Set this property to true to output the correlations computed
        %   using GCC-PHAT algorithm as well as the corresponding lags
        %   between each sensor pairs. Set this property to false to not
        %   output the correlations. The default value of this property is
        %   false.
        CorrelationOutputPort = false
    end
    
    properties (Access = protected, Nontunable, PositiveInteger)
        pNumSources = 1
    end
    
    properties(Constant, Hidden)
        SensorPairSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties (Access = private, Nontunable)
        %Sample rate, in MATLAB, specified by property but in Simulink,
        %specified by engine
        pSampleRate
        %Whether there is only one source
        pIsSingleSource
        %Number of FFT points
        cSigFFT
        %FFT for reference channel
        cRefSigFFT
        %IFFT for correlation
        cIFFT
        %Indices specifying the sensor pairs
        pSensorPair
        %Index corresponding to maximum delay between sensor pairs
        pMaxDelayIndex
    end
    
    properties (Access = private)
        %Delay grid for correlation
        pNFFT
        %FFT for target channel
        pDelayVector
        %Delay grid for correlation search
        pDelaySearchVector
        %The pseudo inverse for least square computation of angle
        pLSInvFactor
        %Indices for delays corresponding to the time that may need to
        %travel among sensor pairs
        pValidDelayIndex
        %Indices for correlations without padding to two's power
        pValidCorrelationIndex
    end
    
    properties (Access = private, Logical)
        pSizeInitialized
    end
    
    methods
        function obj = GCCEstimator(varargin)
            obj@phased.internal.AbstractArrayOperation(varargin{:});
        end
    end
    
    methods
        function set.SampleRate(obj,value)
            validateattributes(value,{'double','single'},{'scalar','positive',...
                'finite'},'GCCEstimator','SampleRate');
            obj.SampleRate = value;
        end
        function set.SensorPair(obj,value)
            validateattributes(value,{'double','single'},{'integer','positive',...
                'finite','2d','nrows',2},...
                'GCCEstimator','SensorPair');
            obj.SensorPair = value;
        end
    end
    
    methods (Access = protected)
        function privValidateSensorArray(obj,val)  %#ok<INUSL>
            validateattributes( val, { 'phased.internal.AbstractArray',...
                'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
        end
        
        function num = getNumOutputsImpl(obj)
            num = 1;
            if obj.DelayOutputPort
                num = num+1;
            end
            if obj.CorrelationOutputPort
                num = num+2;
            end
        end
        
        function validatePropertiesImpl(obj)
            Nele = getDOF(obj.SensorArray);
            if ~(obj.SensorPairSource(1) == 'A')
                sigdatatypes.validateIndex(obj.SensorPair,'',...
                    'SensorPair',{'double','single'},{'<=',Nele});
            end
        end
        
        function validateInputsImpl(obj,x)
            validateattributes(x,{'double','single'},{'nonempty','2d','finite',...
                'ncols',getDOF(obj.SensorArray)},'step','X');
            
            validateNumChannels(obj,x);
        end
        
        function processInputSizeChangeImpl(obj,x)
            sz_x = size(x,1);
            Ncorr = 2*sz_x-1;
            NFFT = 2^nextpow2(Ncorr);
            obj.pNFFT = NFFT;
            
            endidx = NFFT/2+1+getValidDisplayIndex(obj,Ncorr);
            corridx = [-(sz_x-1) (sz_x-1)];
            obj.pValidCorrelationIndex = NFFT/2+1+corridx;
            obj.pValidDelayIndex = endidx;
            % delaygrid = (-NFFT/2:NFFT/2-1).'/obj.pSampleRate;
            % obj.pDelayVector = delaygrid;
            % obj.pDelaySearchVector = delaygrid(endidx(1):endidx(2)).';
        end

        function setupImpl(obj,x)
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);

            %obj.pSampleRate = getSampleRate(obj,sz_x,1,obj.SampleRate);
            fs = obj.SampleRate; % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            obj.pSampleRate = fs;           

            obj.pIsSingleSource = (obj.pNumSources == 1);
            [obj.pSensorPair,posvector,obj.pMaxDelayIndex] = getSensorPair(obj,false);
            [U,S,V] = svd(posvector,0);
            rks = rank(S);
            obj.pLSInvFactor = bsxfun(@times,V(:,1:rks),...
                1./diag(S(1:rks,1:rks)).')*U(:,1:rks)'*obj.PropagationSpeed;
            
            processInputSizeChangeImpl(obj,x)
            % obj.cSigFFT = dsp.FFT('FFTLengthSource','Property','FFTLength',obj.pNFFT);
            % obj.cRefSigFFT = dsp.FFT('FFTLengthSource','Property','FFTLength',obj.pNFFT);
            obj.cSigFFT = dsp.FFT('FFTImplementation','Radix-2');
            obj.cRefSigFFT = dsp.FFT('FFTImplementation','Radix-2');
            obj.cIFFT = dsp.IFFT('FFTImplementation','Radix-2');
            
        end
        
        function varargout = stepImpl(obj,x)
            
            classtouse=class(x);
            if ~obj.pSizeInitialized
                processInputSizeChangeImpl(obj,x);
                obj.pSizeInitialized = true;
            end
            pairidx = obj.pSensorPair;
            sig = x(:,pairidx(1,:));
            refsig = x(:,pairidx(2,:));
            
            % GCC-PHAT
            % zero padding
            numpair = size(pairidx,2);
            corrlen = 2*size(x,1)-1;
            NFFT = 2^nextpow2(corrlen);
            assert(NFFT <= 2^nextpow2(2*getPropagatedNumInputSamples(obj,x)));
            sz_sig = size(sig);
            sigp = zeros(NFFT,numpair,'like',sig);
            sigp(1:sz_sig(1),1:numpair) = sig;
            sz_refsig = size(refsig);
            refsigp = zeros(NFFT,numpair,'like',refsig);
            refsigp(1:sz_refsig(1),1:numpair) = refsig;
            R12 = bsxfun(@times,step(obj.cSigFFT,sigp),...
                conj(step(obj.cRefSigFFT,refsigp)));
            r12_temp = fftshift(step(obj.cIFFT,exp(1i*angle(R12))),1);
            endidx = obj.pValidDelayIndex;
            r12 = r12_temp(endidx(1):endidx(2),:);
            
            r12abs = abs(r12);
            idx = zeros(obj.pNumSources,size(r12,2));
            tau = zeros(obj.pNumSources,size(r12,2),classtouse);
            
            % peak detection, single source
            [~,idx(1,:)] = max(r12abs);
            idxoffset = phased.internal.parabolicFit(r12abs,idx);  
            delaygrid = cast(((0:NFFT-1).'-NFFT/2)/obj.pSampleRate,classtouse);
            delaysearchvec = delaygrid(endidx(1):endidx(2)).';
            tau(1:obj.pNumSources,:) = delaysearchvec(idx)+idxoffset/obj.pSampleRate;
            
            varargout{1} = computeAzElAngle(obj,tau);
            if obj.DelayOutputPort
                varargout{2} = tau;
                if obj.CorrelationOutputPort
                    corridx = obj.pValidCorrelationIndex;
                    y3temp = r12_temp(corridx(1):corridx(2),:);
                    varargout{3} = y3temp(1:corrlen,:);
                    y4temp = delaygrid(corridx(1):corridx(2),:);
                    varargout{4} = y4temp(1:corrlen,:);
                end
            elseif obj.CorrelationOutputPort
                corridx = obj.pValidCorrelationIndex;
                y2temp = r12_temp(corridx(1):corridx(2),:);
                varargout{2} = y2temp(1:corrlen,:);
                y3temp = delaygrid(corridx(1):corridx(2),:);
                varargout{3} = y3temp(1:corrlen,:);
            end
            
        end
        
        function ang = computeAzElAngle(obj,tau)
            % compute Az/El pair from delay
            dirvec = -cast(obj.pLSInvFactor,class(tau))*tau.';
            azel = phased.internal.dirvec2azel(dirvec);
            if isa(obj.SensorArray,'phased.ULA') || ...
                    isa(obj.SensorArray,'phased.HeterogeneousULA')
                ang = azel(1,:);
            else
                ang = azel;
            end
        end
        
        function resetImpl(obj)
            reset(obj.cSigFFT);
            reset(obj.cRefSigFFT);
            reset(obj.cIFFT);
            obj.pSizeInitialized = false;
        end
        
        function releaseImpl(obj)
            release(obj.cSigFFT);
            release(obj.cRefSigFFT);
            release(obj.cIFFT);
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            if (obj.SensorPairSource(1) == 'A') && ...  %Auto
                    strcmp(prop, 'SensorPair')
                flag = true;
            else
                flag = false;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractArrayOperation(obj);
            if isLocked(obj)
                s.pSampleRate = obj.pSampleRate;
                s.pNFFT = obj.pNFFT;
                s.cSigFFT = saveobj(obj.cSigFFT);
                s.cRefSigFFT= saveobj(obj.cRefSigFFT);
                s.cIFFT = saveobj(obj.cIFFT);
                s.pDelayVector = obj.pDelayVector;
                s.pDelaySearchVector = obj.pDelaySearchVector;
                s.pIsSingleSource = obj.pIsSingleSource;
                s.pSensorPair = obj.pSensorPair;
                s.pLSInvFactor = obj.pLSInvFactor;
                s.pNumSources = obj.pNumSources;
                s.pValidDelayIndex = obj.pValidDelayIndex;
                s.pValidCorrelationIndex = obj.pValidCorrelationIndex;
                s.pMaxDelayIndex = obj.pMaxDelayIndex;
                s.pSizeInitialized = obj.pSizeInitialized;
            end
        end

        function s = loadSubObjects(obj,s,wasLocked)
            s = loadSubObjects@phased.internal.AbstractArrayOperation(obj,s);
            if isfield(s,'isLocked')
                s = rmfield(s,'isLocked');
            end
            if wasLocked
                if isfield(s,'cSigFFT')
                    obj.cSigFFT = dsp.FFT.loadobj(s.cSigFFT);
                    s = rmfield(s,'cSigFFT');
                end
                if isfield(s,'cRefSigFFT')
                    obj.cRefSigFFT = dsp.FFT.loadobj(s.cRefSigFFT);
                    s = rmfield(s,'cRefSigFFT');
                end
                if isfield(s,'cIFFT')
                    obj.cIFFT = dsp.IFFT.loadobj(s.cIFFT);
                    s = rmfield(s,'cIFFT');
                end
            end
        end

        function loadObjectImpl(obj,s,wasLocked) 
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function flag = isInputComplexityLockedImpl(obj,~)   %#ok<INUSD>
            flag = false;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end

        function flag = isInputSizeLockedImpl(~,~)
            flag = false;
        end
        
        function varargout = getOutputSizeImpl(obj)
            if isa(obj.SensorArray,'phased.ULA') || ...
                    isa(obj.SensorArray,'phased.HeterogeneousULA')
                varargout{1} = [1 1];
            else
                varargout{1} = [2 1];
            end
            sp = getSesorPairIndex(obj);
            npairs = size(sp,2);
            sz_x = propagatedInputSize(obj,1);
            nlag = 2*sz_x(1)-1;
            if obj.DelayOutputPort 
                varargout{2} = [1 npairs];
            end
            if obj.CorrelationOutputPort
                varargout = [varargout {[nlag npairs],[nlag 1]}];
            end
        end
        
        function varargout = isOutputFixedSizeImpl(obj) 
            varargout = {true};
            if obj.DelayOutputPort 
                varargout{2} = true;
            end
            if obj.CorrelationOutputPort
                varargout = [varargout {false, false}];
            end
        end
        
        function varargout = getOutputDataTypeImpl(obj)
            dt_out = propagatedInputDataType(obj,1);
            varargout = {dt_out dt_out dt_out dt_out};
        end
        
        function varargout = isOutputComplexImpl(obj) 
            if obj.DelayOutputPort
                if obj.CorrelationOutputPort
                    varargout = {false false true false};
                else
                    varargout = {false false};
                end
            else
                if obj.CorrelationOutputPort
                    varargout = {false true false};
                else
                    varargout = {false};
                end
            end
        end
   end
    
    
    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractArrayOperation('subarray');
        % dSampleRate = matlab.system.display.internal.Property(...
        %     'SampleRate','IsObjectDisplayOnly',true);
        props = {...
          'SampleRateFromInputCheckbox',...
          'SampleRate',...
          'SensorPairSource',...
          'SensorPair',...
          'CorrelationOutputPort',...
          'DelayOutputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
    end
    
    methods (Static,Hidden,Access = protected) %for Simulink
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:GCCEstimatorTitle')),...
              'Text',getString(message('phased:library:block:GCCEstimatorDesc')));
        end 
    end
    methods (Access = protected) %for Simulink
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('GCC-PHAT\nDOA');
        end
        function varargout = getInputNamesImpl(obj)   %#ok<MANU>
            varargout = {''};
        end

        function varargout = getOutputNamesImpl(obj)  
            if obj.DelayOutputPort
                if obj.CorrelationOutputPort
                    varargout = {'Ang','Tau','Rxy','Lag'};
                else 
                    varargout = {'Ang','Tau'};
                end
            else
                if obj.CorrelationOutputPort
                    varargout = {'Ang','Rxy','Lag'};
                else
                    varargout = {'Ang'};
                end
            end
        end
    end
    
    methods (Access = private)
        function spidx = getSesorPairIndex(obj)
            if (obj.SensorPairSource(1) == 'A')
                Nele = getDOF(obj.SensorArray);
                spidx = [2:Nele;ones(1,Nele-1)];
            else
                spidx = obj.SensorPair;
            end
        end
        
        function [spidx, sp_posvec, maxindex] = getSensorPair(obj,warnflag)
            spidx = getSesorPairIndex(obj);
            if isa(obj.SensorArray,'phased.internal.AbstractSubarray')
                pos = getSubarrayPosition(obj.SensorArray);
            else
                pos = getElementPosition(obj.SensorArray);
            end
            sp_posvec = ((pos(:,spidx(1,:)) - pos(:,spidx(2,:))).');
            
            posdist = sqrt(sum(abs(sp_posvec).^2,2));
            maxdist = max(posdist);
            mindist = min(posdist);
            sampdist = obj.PropagationSpeed/obj.pSampleRate;
            % If signal arrives at the same sample between channels, the
            % result may be inaccurate.
            if warnflag && ((maxdist <= sampdist) || (mindist <= sampdist))
                warning(message('phased:doa:GCCSampleRateTooSmall'));
            end
            maxindex = ceil(maxdist*obj.pSampleRate/obj.PropagationSpeed);
        end
        
        function endidx = getValidDisplayIndex(obj,Ncorr)
            Ncorr_oneside = (Ncorr-1)/2;
            idxspan_oneside = min(obj.pMaxDelayIndex,Ncorr_oneside-1);
            idxspan_oneside = max(idxspan_oneside,1); % at least 3 samples in correlation
            endidx = [-idxspan_oneside idxspan_oneside];           
        end
    end
end