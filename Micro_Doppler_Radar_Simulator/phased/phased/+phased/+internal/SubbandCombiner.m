classdef (Hidden,Sealed,StrictDefaults) SubbandCombiner < phased.internal.AbstractVarSizeEngine & ...
        matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%SubbandCombiner     Combine signal from subbands
%   H = phased.internal.SubbandCombiner creates a subband combiner System
%   object, H. The object combines signal from subbands.
%
%   H = phased.WidebandRadiator(Name,Value) creates a subband combiner
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Y = step(H,X) combines signal X from multiple subbands.  X is a
%   PxNxNFFT matrix where NFFT is the number of subbands and P is the
%   number of snapshots in each frequency. Y is an MxN matrix whose rows
%   are snapshots in time domain. The value of M is determined by the value
%   of the TimeSignalLength property.
%
%   Y = step(H,X,XT) specifies the probe signal XT as an MxN matrix where M
%   and N are the dimensions of the output signal, Y. This syntax only
%   applies when you set the TimeSignalLengthSource property to 'Inherit'.
%   This syntax is useful when the output signal is not fixed size.
%
%   Y = step(H,X,FD) specifies the Doppler shifts in subbands in FD so the
%   combined signal contains the correct Dopple information. FD is a
%   length-NFFT row vector whose elements are the Doppler shift for the
%   correspnding subband. This syntax only applies when you set the
%   SubbandFrequencyShiftInputPort property to true.
%
%   You can combine optional input arguments when their enabling
%   properties are set. Optional inputs must be listed in the same order
%   as the order of the enabling properties. For example,
%
%   Y = step(H,X,XT,FD)
%
%   This object should be used in pair with phased.internal.SubbandDivider
%   so when M is not an integer multiple of NFFT, it knows how to revert
%   the operations performed by phased.internal.SubbandDivider.
%
%   SubbandCombiner methods:
%
%   step     - Divide signal to subbands (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a subband combiner object with same property values
%   isLocked - Locked status (logical)
%
%   SubbandCombiner properties:
%
%   NumSubbands                    - Number of subbands
%   TimeSignalLengthSource         - Souce of time domain signal length
%   TimeSignalLength               - Signal length in time
%   SubbandFrequencyShiftInputPort - Enable subband frequency shift input
%   EnableWarning                  - Enable warning when the signal is
%                                    short
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
%   See also phased, phased.internal.SubbandDivider,
%   phased.internal.subbandCenterFrequency.

%   Copyright 2015-2016 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002  
    
%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %NumSubbands    Number of subbands
        %   Specify the number of subbands used in the subband processing
        %   as a positive integer. The default value of this property is
        %   64.
        NumSubbands = 64
    end
    
    properties (Nontunable)
        %TimeSignalLengthSource    Source of time domain signal length
        %   Specify how to determine the time domain signal length as one
        %   of 'Property' | 'Inherit', where the default is 'Property'.
        %   When you set this property to 'Property', the signal length is
        %   determined by the value of the TimeSignalLength property. When
        %   you set this property to 'Inherit', the signal length is the
        %   same as the signal length of the probe signal.
        TimeSignalLengthSource = 'Property'
        %TimeSignalLength    Signal length in time
        %   Specify the length of the signal in time before it is broken
        %   into subbands as a positive integer. The default value of this
        %   property is 64.
        TimeSignalLength = 64
        %SubbandFrequencyShiftInputPort  Enable subband frequency shift
        %input
        %   Set this property to true to allow specifying frequency shifts
        %   in each subband before combining subband signals back to time
        %   domain. This is useful when you want to simulate the Doppler at
        %   each subband. Set this property to false to not allow
        %   specifying frequency shifts in each subband. The default
        %   value of this property is false.
        SubbandFrequencyShiftInputPort = false
    end

    properties (Nontunable, Logical)
        %EnableWarning  Enable warning when the signal is short
        %   Set this property to true to warn when number of samples is
        %   less than the number of the subbands. Set this property to
        %   false to disable the warning. The default value is true.
        EnableWarning = true
    end
    
    properties(Constant, Hidden)
        TimeSignalLengthSourceSet = matlab.system.StringSet(...
            {'Property','Inherit'});
    end
    
    properties (Access = protected, Nontunable)
        cIFFT
    end
    
    properties (Access = protected)
        pNumSnapshots
        pNumRemainingSamples
    end
    
    properties (Access = protected)
        pIFFTMatrix
    end
    
    properties (Access = protected, Logical, Nontunable)
        pIsTimeSignalLengthViaProperty 
    end
    
    methods
        function obj = SubbandCombiner(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods (Access = protected)
        function resetImpl(obj)
            reset(obj.cIFFT);
        end
        
        function releaseImpl(obj)
            release(obj.cIFFT);
        end

        function num = getNumInputsImpl(obj)
            num = 1;
            if obj.SubbandFrequencyShiftInputPort
                num = num+1;
            end
            if strcmp(obj.TimeSignalLengthSource,'Inherit')
                num = num+1;
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if strcmp(obj.TimeSignalLengthSource,'Inherit') && ...
                    strcmp(prop,'TimeSignalLength')
                flag = true;
            end
        end
        
        function processInputSizeChangeImpl(obj,~,xt,~)
            if ~strcmp(obj.TimeSignalLengthSource,'Property')
                len_x = size(xt,1);
                computeRemainingSamples(obj,len_x,obj.NumSubbands);
            end

        end
    
        function setupImpl(obj,x,xt,~) 
            Nbands = obj.NumSubbands;
            obj.pIsTimeSignalLengthViaProperty = ...
                strcmp(obj.TimeSignalLengthSource,'Property');
            if strcmp(obj.TimeSignalLengthSource,'Property')
                len_x = obj.TimeSignalLength;
            else
                len_x = size(xt,1);
            end
            computeRemainingSamples(obj,len_x,Nbands);
            obj.cIFFT = dsp.IFFT;
            if obj.SubbandFrequencyShiftInputPort
                obj.pIFFTMatrix = step(obj.cIFFT,eye(Nbands));
            end
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
        end
                
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
            if isLocked(obj)
                s.cIFFT = saveobj(obj.cIFFT);
                s.pNumSnapshots = obj.pNumSnapshots;
                s.pNumRemainingSamples = obj.pNumRemainingSamples;
                s.pIFFTMatrix = obj.pIFFTMatrix;
                s.pIsTimeSignalLengthViaProperty = obj.pIsTimeSignalLengthViaProperty;
            end
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                obj.cIFFT = dsp.IFFT.loadobj(s.cIFFT);
                s = rmfield(s,'cIFFT');
                if ~isfield(s,'pIsTimeSignalLengthViaProperty')
                    obj.pIsTimeSignalLengthViaProperty = true;
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
        
        function y = formatOutput(obj,x)
        % restore the output to original length

            Nremaining = obj.pNumRemainingSamples;
            Nsubbands = obj.NumSubbands;
            xp = reshape(permute(x,[1 3 2]),[],obj.pValidatedNumInputChannels);
            if Nremaining > 0
                y = [xp(1:(obj.pNumSnapshots-1)*Nsubbands,:);...
                    xp(end-Nremaining+1:end,:)];
            elseif Nremaining < 0
                y = xp(1:-Nremaining,:);
            else
                y = xp;
            end

        end
        
        function y = stepImpl(obj,x_in,xt,fd_in)
        % convert the time domain signal into frequency domain using FFT, x
        % should be output of formatInput method

            % x_in is numsnapshots x channel x bands
            [sz_x1,sz_x2,sz_x3] = size(x_in);
            if ismatrix(x_in) || sz_x3 == 1
                % each band has only one sample
                ytemp = x_in;  % scalar ifft is itself
            else
                if obj.SubbandFrequencyShiftInputPort
                    % fd normalized to fs, row vector
                    if obj.pIsTimeSignalLengthViaProperty
                        fd = xt;
                    else
                        fd = fd_in;
                    end
                    Nbands = obj.NumSubbands;
                    if any(fd~=0)
                        % ifftmtx = obj.pIFFTMatrix.*exp(1i*2*pi*(0:Nbands-1).'*fd);
                        % ytemp = reshape(ifftmtx*reshape(permute(x_in,[3 1 2]),sz_x3,[]),[],sz_x2);
                        ytemp = complex(zeros(sz_x1*sz_x3,sz_x2));
                        x_in_temp = permute(x_in,[3 2 1]);

                        % Anti-aliasing: remove frequencies shifted outside
                        % of the range [-fs/2,fs/2)
                        f0 = phased.internal.subbandCenterFrequency(0,1,Nbands)*Nbands;
                        fbounds = [min(f0) max(f0)];
                        fshifted = fd'*Nbands+f0;
                        isaliased = fshifted>fbounds(2) | fshifted<fbounds(1);
                        x_in_temp(isaliased,:) = 0;
                        
                        for m = 1:sz_x1
                            ifftmtx = obj.pIFFTMatrix.*exp(1i*2*pi*((0:Nbands-1).'+(m-1)*Nbands)*fd);
                            ytemp((m-1)*Nbands+1:m*Nbands,:) = ifftmtx*x_in_temp(:,:,m);
                        end
                    else
                        ifftmtx = obj.pIFFTMatrix;
                        ytemp = reshape(ifftmtx*reshape(permute(x_in,[3 1 2]),sz_x3,[]),[],sz_x2);
                    end
                else
                    ytemp = reshape(step(obj.cIFFT,permute(x_in,[3 1 2])),[],sz_x2);
                end
            end

            % restore original length column vector
            yt = formatOutput(obj,ytemp);
            if obj.pIsTimeSignalLengthViaProperty
                y = yt(1:obj.TimeSignalLength,:);
            else
                y = yt(1:size(xt,1),:);
            end

        end
        
    end
    
    methods (Access = private)
        function computeRemainingSamples(obj,len_x,Nbands)
            Nsnapshots = floor(len_x/Nbands);
            Nremaining = len_x - Nbands*Nsnapshots;
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

        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            if strcmp(obj.TimeSignalLengthSource, 'Property')
                x_len = obj.TimeSignalLength;
            else
                szXT = propagatedInputSize(obj,2);
                x_len = szXT(1);
            end
            varargout{1} = [x_len szX(2)];
        end
    end
end
