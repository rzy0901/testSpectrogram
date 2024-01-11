classdef (Hidden, Sealed, StrictDefaults) SubbandDivider < phased.internal.AbstractVarSizeEngine & ...
        matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%SubbandDivider     Divide signal into subbands
%   H = phased.internal.SubbandDivider creates a subband divider System
%   object, H. The object divides signal into subbands.
%
%   H = phased.WidebandRadiator(Name,Value) creates a subband divider
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Y = step(H,X) divides signal X into multiple subbands. X is an MxN matrix
%   whose rows are snapshots in time domain. Y is a PxNxNFFT matrix where
%   NFFT is the number of subbands and P is the number of snapshots in each
%   frequency. 
%
%   P can be computed as ceil(M/NFFT). If M is less than NFFT, zeros are
%   padded when performing NFFT. If M is not an integer multiple of NFFT,
%   then the last conversion to subbands is performed on the last NFFT
%   snapshots of X in time. For example, if NFFT is 4 and M is 5, then P is
%   2. The first conversion is performed on snapshots 1 through 4, and the
%   second conversion is performed on snapshots 2 through 5.
%
%   This object is often used in pair with phased.internal.SubbandCombiner.
%
%   SubbandDivider methods:
%
%   step     - Divide signal to subbands (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a subband divider object with same property values
%   isLocked - Locked status (logical)
%
%   SubbandDivider properties:
%
%   OperatingFrequency     - Operating frequency
%   SampleRate             - Sample rate
%   NumSubbands            - Number of subbands
%   EnableWarning          - Enable warning when the signal is short
%
%   % Examples:
%   %   Divide and combine a wideband signal
%   
%   fc = 3e8; fs = 3e7; 
%   x = randn(128,1)+1i*randn(128,1);
%   sd = phased.internal.SubbandDivider('OperatingFrequency',fc,...
%       'SampleRate',fs);
%   xf = step(sd,x);
%   sc = phased.internal.SubbandCombiner('TimeSignalLength',size(x,1));
%   y = step(sc,xf);
%   max(abs(x-y))
%
%   See also phased, phased.internal.SubbandCombiner,
%   phased.internal.subbandCenterFrequency.

%   Copyright 2015 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %OperatingFrequency     Operating frequency (Hz)
        %   Specify the operating frequency (in Hz) of the beamformer as a
        %   scalar. The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 3e8
    end
    
    properties (Nontunable, Dependent)
        %SampleRate     Sample rate
        %   Specify the signal sampling rate (in Hz) as a positive scalar.
        %   The default value of this property is 1e6.
        SampleRate = 1e6
    end
    
    properties (Nontunable, PositiveInteger)  
        %NumSubbands    Number of subbands
        %   Specify the number of subbands used in the subband processing
        %   as a positive integer. The default value of this property is
        %   64.
        NumSubbands = 64
    end
    
    properties (Nontunable, Logical)
        %EnableWarning  Enable warning when the signal is short
        %   Set this property to true to warn when number of samples is
        %   less than the number of the subbands. Set this property to
        %   false to disable the warning. The default value is true.
        EnableWarning = true
    end
    
    properties (Access = protected, Nontunable)
        pSubbandCenterFreq
        pSampleRate = 1e6
        cFFT
    end
    
    properties (Access = protected)
        pNumSnapshots
        pNumRemainingSamples
        pNsnapshots
        pXSize
    end
    
    methods

        function set.OperatingFrequency(obj,val)
            sigdatatypes.validateFrequency(val,'phased.internal',...
                'OperatingFrequency',{'scalar'});
            obj.OperatingFrequency = val;
        end

        function set.SampleRate(obj,val)
            sigdatatypes.validateFrequency(val,'phased.internal',...
                'SampleRate',{'scalar'});
            obj.pSampleRate = val;
        end
        
        function val = get.SampleRate(obj)
            val = obj.pSampleRate;
        end
    end
    
    methods 
        function obj = SubbandDivider(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods (Access = protected)
        function resetImpl(obj)
            reset(obj.cFFT);
        end
        
        function releaseImpl(obj)
            release(obj.cFFT);
        end

        function setupImpl(obj,x)
            cond = obj.pSampleRate >= 2*obj.OperatingFrequency;
            if cond
                coder.internal.errorIf(cond,'phased:SubbandBeamformer:SampleRateMismatch');
            end
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
            Nbands = obj.NumSubbands;
            obj.pSubbandCenterFreq = phased.internal.subbandCenterFrequency(...
                obj.OperatingFrequency,obj.pSampleRate,Nbands);
            
            if Nbands > 1
                obj.cFFT = dsp.FFT('FFTLengthSource','Property','FFTLength',...
                    Nbands);
            else
                obj.cFFT = dsp.FFT;
            end
            
            processInputSizeChangeImpl(obj,x);
        end
                
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
            s.pSampleRate = obj.pSampleRate;
            if isLocked(obj)
                s.pSubbandCenterFreq = obj.pSubbandCenterFreq;
                s.pNumSnapshots = obj.pNumSnapshots;
                s.pNumRemainingSamples = obj.pNumRemainingSamples;
                s.pNsnapshots = obj.pNsnapshots;
                s.pXSize = obj.pXSize;
                s.cFFT = saveobj(obj.cFFT);
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            obj.pSampleRate = s.pSampleRate;
            s = rmfield(s,'pSampleRate');
            if wasLocked
                obj.cFFT = dsp.FFT.loadobj(s.cFFT);
                s = rmfield(s,'cFFT');
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) 
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function y = formatInput(obj,x)
        % format the input sequence to frames, i.e., making the number of
        % sample equal to integer multiple of Nbands. Use
        % pNumRemainingSamples to store information on input sample number
        % N: 0 means N = k*Nbands; negative pNumRemainingSamples means N <
        % Nbands and N = -pNumRemainingSamples; positive
        % pNumRemainingSamples means N=k*Nbands+pNumRemainingSamples.

            Nbands = obj.NumSubbands;
            Nsnapshots = obj.pNsnapshots;
            % sz_x = obj.pXSize;
            if Nsnapshots == 0
                y = zeros(Nbands,obj.pValidatedNumInputChannels,'like',x);
                y(1:size(x,1),:) = x;
                % y = [x; zeros(Nbands-sz_x(1),sz_x(2))];
            elseif obj.pNumSnapshots == Nsnapshots+1
                % When the number of samples is greater than integer
                % multiple of Nbands, add an artificial snapshot to cover
                % the remaining samples.
                y = [x(1:Nsnapshots*Nbands,:);...
                    x(end-Nbands+1:end,:)];
            else
                y = x;
            end

        end
        
        function [yout,freq] = stepImpl(obj,x_in)
        % convert the time domain signal into frequency domain using FFT, x
        % should be output of formatInput method
        
            Nbands = obj.NumSubbands;
            Nsnapshots = floor(size(x_in,1)/Nbands);
            x = formatInput(obj,x_in);
            if obj.pNumSnapshots > obj.pNsnapshots
                Nsnapshots = Nsnapshots+1;
            end
            y = complex(zeros(Nsnapshots,size(x, 2),Nbands)); % snap x dim2 x freq

            for snapI = 1:Nsnapshots
               snapshotInTime = x((1:Nbands)+Nbands*(snapI-1),:);
               x_fft_in = complex(snapshotInTime);
               if size(x_fft_in,1) == 1
                   % only one band
                   snapshotInFreq = x_fft_in;  % scalar fft is itself
               else
                   snapshotInFreq = step(obj.cFFT,x_fft_in);
               end
               y(snapI, :, :) = snapshotInFreq.';
            end
            freq = obj.pSubbandCenterFreq;

            xlen = size(x_in,1);
            Ns = floor(xlen/Nbands);
            if rem(xlen,Nbands)~=0
                yout = y(1:Ns+1,:,:); % snap x dim2 x freq
            else
                yout = y(1:Ns,:,:); % snap x dim2 x freq
            end
        end
        
        function processInputSizeChangeImpl(obj,x)
            sz_x = size(x);
            obj.pXSize = sz_x;
            Nbands = obj.NumSubbands;
            Nsnapshots = floor(sz_x(1)/Nbands);
            obj.pNsnapshots = Nsnapshots;
            Nremaining = sz_x(1) - Nbands*obj.pNsnapshots;
            if Nsnapshots == 0
                if isempty(coder.target) && obj.EnableWarning
                    warning(message('phased:SubbandBeamformer:NotEnoughSamples', Nremaining, Nbands));
                end
                obj.pNumSnapshots = 1;
                % Use negative pNumRemainingSamples to indicate the total
                % input data length is shorter then Nbands
                obj.pNumRemainingSamples = -Nremaining;
            elseif Nremaining > 0
                % When the number of samples is greater than integer
                % multiple of Nbands, add an artificial snapshot to cover
                % the remaining samples.
                obj.pNumSnapshots = Nsnapshots + 1;
                obj.pNumRemainingSamples = Nremaining;
            else
                obj.pNumSnapshots = Nsnapshots;
                obj.pNumRemainingSamples = Nremaining;
            end
            
        end
        
    end
    
    methods (Access = protected)
        function flag = isInputSizeLockedImpl(~,~)
            flag = false;
        end
        
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj, 1);
        end
         
    end
end
    
