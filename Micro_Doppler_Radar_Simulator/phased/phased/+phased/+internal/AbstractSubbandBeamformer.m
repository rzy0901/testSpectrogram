classdef (Hidden) AbstractSubbandBeamformer < phased.internal.AbstractBeamformer & ...
        matlab.system.mixin.SampleTime
%This class is for internal use only. It may be removed in the future.

%AbstractWideband Abstract subband beamformer

%   Copyright 2009-2017 The MathWorks, Inc.


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
    
    properties (Nontunable, Logical)
        %SampleRateFromInputCheckbox Inherit sample rate 
        %   Set SampleRateFromInputCheckbox to true to derive sample rate
        %   from Simulink time engine. Set SampleRateFromInputCheckbox to
        %   false to specify the sample rate. This property applies when
        %   used in Simulink.
        SampleRateFromInputCheckbox = true
    end
    
    properties (Nontunable)
        %SampleRate     Sample rate
        %   Specify the signal sampling rate (in Hz) as a positive scalar.
        %   The default value of this property is 1e6.
        SampleRate = 1e6
    end
    
    properties (Nontunable, Logical) 
        %WeightsOutputPort    Enable weights output
        %   Set this property to true to output the weights used in the
        %   beamformer. Set this property to false to not output the
        %   weights. The default value of this property is false.
        WeightsOutputPort = false;
    end
    
    properties (Nontunable, PositiveInteger)  
        %NumSubbands    Number of subbands
        %   Specify the number of subbands used in the subband processing
        %   as a positive integer. The default value of this property is
        %   64.
        NumSubbands = 64
    end

    properties (Nontunable, Logical) 
        %SubbandsOutputPort     Enable subband center frequencies output
        %   Set this property to true to output center frequencies of each
        %   subband. Set this property to false to not output center
        %   frequencies. The default value of this property is false.
        SubbandsOutputPort = false
    end
    
    properties (Constant, Hidden)
        SampleRateSet = matlab.system.SourceSet({'PropertyOrMethod',...
            'SystemBlock', 'SampleRateFromInputCheckbox',...
            'getSampleRateInSimulation',false})
    end
    
    properties (Access = protected, Nontunable)
        pSubbandCenterFreq
        pSampleRate
        cFFT
    end
    
    properties (Access = protected)
        pNumSnapshots
        pNumRemainingSamples
        pNsnapshots
        pXSize
    end
    
    properties (Access = protected,Logical)
        pSizeInitialized
    end
    
    methods

        function set.OperatingFrequency(obj,val)
            sigdatatypes.validateFrequency(val,'phased.internal',...
                'OperatingFrequency',{'double','single'},{'scalar'});
            obj.OperatingFrequency = val;
        end

        function set.SampleRate(obj,val)
            sigdatatypes.validateFrequency(val,'phased.internal',...
                'SampleRate',{'double','single'},{'scalar'});
            obj.SampleRate = val;
        end
    end
    
    methods (Access = protected)
        function obj = AbstractSubbandBeamformer(varargin)
            obj@phased.internal.AbstractBeamformer(varargin{:});
        end
    end
    
    methods (Access = protected)
        function num = getNumOutputsImpl(obj)
            num = getNumOutputsImpl@phased.internal.AbstractBeamformer(obj);
            if obj.WeightsOutputPort
                num = 2;
            end
            if obj.SubbandsOutputPort
                num = num+1;
            end
        end

        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractBeamformer(obj);
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractBeamformer(obj);
            reset(obj.cFFT);
            obj.pSizeInitialized = false;
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractBeamformer(obj);
            release(obj.cFFT);
        end

        function setupImpl(obj,x)
            setupImpl@phased.internal.AbstractBeamformer(obj,x);
            
            %obj.pSampleRate = getSampleRate(obj,size(x,1),1,obj.SampleRate);
            fs = cast(obj.SampleRate,class(x)); % property/method duality
            cond = ~isscalar(fs) || (fs<=0);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:phased:invalidSampleTime');
            end
            obj.pSampleRate = fs;
            
            cond = obj.pSampleRate >= 2*obj.OperatingFrequency;
            if cond
                coder.internal.errorIf(cond,'phased:SubbandBeamformer:SampleRateMismatch');
            end
            
            Nbands = cast(obj.NumSubbands,class(x));
            obj.pSubbandCenterFreq = phased.internal.subbandCenterFrequency(...
                obj.OperatingFrequency,obj.pSampleRate,Nbands);
            
            if Nbands > 1
                obj.cFFT = dsp.FFT('FFTLengthSource','Property','FFTLength',...
                    Nbands);
            else
                obj.cFFT = dsp.FFT;
            end
            
            computeRemainingSamples(obj,size(x,1),Nbands);
            
        end
        
        function processInputSizeChangeImpl(obj,x,~,~)
            sz_x = size(x);
            obj.pXSize = sz_x;
            computeRemainingSamples(obj,sz_x(1),obj.NumSubbands);
        end
        
        function computeRemainingSamples(obj,len_x,Nbands)
            [Nsnapshots,NumSnapshots,NumRemainingSamples] = ...
                getRemainingSamples(len_x,cast(Nbands,'double'));%%codegen purpose
            obj.pNsnapshots = Nsnapshots;
            obj.pNumSnapshots = NumSnapshots;
            obj.pNumRemainingSamples = NumRemainingSamples;
            
        end
                
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractBeamformer(obj);
            if isLocked(obj)
                s.pSubbandCenterFreq = obj.pSubbandCenterFreq;
                s.pNumSnapshots = obj.pNumSnapshots;
                s.pNumRemainingSamples = obj.pNumRemainingSamples;
                s.pNsnapshots = obj.pNsnapshots;
                s.pXSize = obj.pXSize;
                s.pSampleRate = obj.pSampleRate;
                s.pSizeInitialized = obj.pSizeInitialized;
                s.cFFT = saveobj(obj.cFFT);
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractBeamformer(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    if isfield(s,'cFFT')
                        obj.cFFT = dsp.FFT.loadobj(s.cFFT);
                        s = rmfield(s,'cFFT');
                    end
                    % recover locked sample rate information
                    if isfield(s,'pSampleRate')
                        obj.pSampleRate = s.pSampleRate;
                        s = rmfield(s,'pSampleRate');
                    else
                        obj.pSampleRate = s.SampleRate;
                    end
                end
            end
        end
        
        function y = formatInput(obj,x,type)
        % format the input sequence to frames, i.e., making the number of
        % sample equal to integer multiple of Nbands. Use
        % pNumRemainingSamples to store information on input sample number
        % N: 0 means N = k*Nbands; negative pNumRemainingSamples means N <
        % Nbands and N = -pNumRemainingSamples; positive
        % pNumRemainingSamples means N=k*Nbands+pNumRemainingSamples.

            Nbands = obj.NumSubbands;
            if type(1) == 'S' % signal
                Nsnapshots = obj.pNsnapshots;
                % sz_x = obj.pXSize;
                NumSnapshots = obj.pNumSnapshots;
            else  % training
                [Nsnapshots,NumSnapshots] = ...
                        getRemainingSamples(size(x,1),Nbands);                
            end
            
            if Nsnapshots == 0
                y = zeros(Nbands,obj.pValidatedNumInputChannels,'like',x);
                y(1:size(x,1),:) = x;
                % y = [x; zeros(Nbands-sz_x(1),sz_x(2))];
            elseif NumSnapshots == Nsnapshots+1
                % When the number of samples is greater than integer
                % multiple of Nbands, add an artificial snapshot to cover
                % the remaining samples.
                y = [x(1:Nsnapshots*Nbands,:);...
                    x(end-Nbands+1:end,:)];
            else
                y = x;
            end

        end
        
        function y = formatOutput(obj,x)
        % restore the output to original length

            Nremaining = obj.pNumRemainingSamples;
            if Nremaining > 0
                y = [x(1:(obj.pNumSnapshots-1)*obj.NumSubbands);...
                    x(end-Nremaining+1:end)];
            elseif Nremaining < 0
                y = x(1:-Nremaining);
            else
                y = x;
            end

        end
        
        function y = subbandFFT(obj,x,x_in,type)
        % convert the time domain signal into frequency domain using FFT, x
        % should be output of formatInput method
            
            classtouse=class(x);        
            Nbands = obj.NumSubbands;
            if type(1) == 'S' % signal
                Nsnapshots = floor(size(x_in,1)/Nbands);
                if obj.pNumSnapshots > obj.pNsnapshots
                    Nsnapshots = Nsnapshots+1;
                end
            else % training
                [Nsnapshots,NumSnapshots] = ...
                        getRemainingSamples(size(x,1),Nbands);                
                if NumSnapshots > Nsnapshots
                    Nsnapshots = Nsnapshots+1;
                end
            end
            y = complex(zeros(Nsnapshots,size(x, 2),Nbands,classtouse));

            for snapI = 1:Nsnapshots
               snapshotInTime = x((1:Nbands)+Nbands*(snapI-1),:);
               if size(snapshotInTime,1) == 1
                   % each band has only one sample
                   snapshotInFreq = complex(snapshotInTime);  % scalar fft is itself
               else
                   snapshotInFreq = step(obj.cFFT,complex(snapshotInTime(1:Nbands,1:obj.pValidatedNumInputChannels)));
               end
               y(snapI, :, :) = snapshotInFreq.';
            end

        end
    end
    methods (Access = protected) %for Simulink
        function varargout = getOutputNamesImpl(obj)
            %Insert Freq
            if obj.SubbandsOutputPort
                [varargout{1:nargout-1}] = getOutputNamesImpl@phased.internal.AbstractBeamformer(obj);
                varargout{end+1} = 'Freq';
            else
                [varargout{1:nargout}] = getOutputNamesImpl@phased.internal.AbstractBeamformer(obj);
            end
        end
        function varargout = getOutputSizeImpl(obj)
            
            if obj.SubbandsOutputPort
                [varargout{1:nargout-1}] = getOutputSizeImpl@phased.internal.AbstractBeamformer(obj);
                varargout{end+1} = [obj.NumSubbands 1];
            else
                [varargout{1:nargout}] = getOutputSizeImpl@phased.internal.AbstractBeamformer(obj);                
            end
            if obj.WeightsOutputPort
                szW = varargout{2};
                varargout{2} = [szW(1) obj.NumSubbands szW(2)];
            end

        end
        function varargout = isOutputFixedSizeImpl(obj)
            if obj.SubbandsOutputPort
                [varargout{1:nargout-1}] = isOutputFixedSizeImpl@phased.internal.AbstractBeamformer(obj);
                varargout{end+1} = true;
            else
                [varargout{1:nargout}] = isOutputFixedSizeImpl@phased.internal.AbstractBeamformer(obj);                
            end
        end
        function varargout = getOutputDataTypeImpl(obj)
            dt = propagatedInputDataType(obj,1);
            varargout = {dt dt dt};
        end
        function varargout = isOutputComplexImpl(obj) %#ok<MANU>
            varargout = {true, true, true};
        end

    end    
     
end

function [Nsnapshots,NumSnapshots,NumRemainingSamples] = ...
        getRemainingSamples(len_x,Nbands)
    Nsnapshots = floor(len_x/Nbands);
    Nremaining = len_x - Nbands*Nsnapshots;
    if Nsnapshots == 0
        if isempty(coder.target)
            warning(message('phased:SubbandBeamformer:NotEnoughSamples', Nremaining, Nbands));
        end
        NumSnapshots = 1;
        % Use negative pNumRemainingSamples to indicate the total
        % input data length is shorter then Nbands
        NumRemainingSamples = -Nremaining;
    elseif Nremaining > 0
        % When the number of samples is greater than integer
        % multiple of Nbands, add an artificial snapshot to cover
        % the remaining samples.
        NumSnapshots = Nsnapshots + 1;
        NumRemainingSamples = Nremaining;
    else
        NumSnapshots = Nsnapshots;
        NumRemainingSamples = Nremaining;
    end
end

