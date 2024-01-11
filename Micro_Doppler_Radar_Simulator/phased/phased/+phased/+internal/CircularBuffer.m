classdef (Hidden, Sealed, StrictDefaults) CircularBuffer < phased.internal.AbstractVarSizeEngine
%This class is for internal use only. It may be removed in the future.

%CircularBuffer   Circular buffer
%   H = phased.internal.CircularBuffer(N) creates a circular buffer System
%   object, H, with the buffer length set to N. The object simulates the
%   circular buffering behavior.
%
%   Step syntax:
%
%   Y = step(H,X,M) pushes X into the buffer's (M+1)th position and then
%   outputs Y, whose dimension is the same as X, from the beginning of the
%   buffer. The buffer will dynamically grow if the current buffer cannot
%   accommodate X when it is inserted at (M+1)th position.
%
%   X can be a matrix. Each column of X is considered to be a separate
%   signal.
%
%   CircularBuffer methods:
%
%   step  - See above description for use of this method
%   reset - Reset the buffer
%

%   Copyright 2010-2015 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable, Logical)
        FixedLengthBuffer = false
    end
    
    properties (Nontunable, PositiveInteger)
        BufferLength = 1
    end
    
    properties (Nontunable)
        BufferWidthSource = 'Auto'
    end
    
    properties (Nontunable, PositiveInteger)
        BufferWidth = 1
    end
    
    properties (Nontunable)
        MaxNumInputSamplesSource = 'Auto'
    end
    
    properties (Nontunable)
        MaxNumInputSamples = 0
    end

    properties(Constant, Hidden)
        MaxNumInputSamplesSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
        BufferWidthSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    properties (Access = private)
        % Private property that holds the count of occupied cells in the
        % buffer
        pCount
        % Private property that holds the starting position of the buffer
        pStart
        % Private property that holds the buffer
        pBuffer
    end
    
    properties (Access = private,Nontunable)    
        % Private property that holds number of channels that share the
        % same delay
        pNumSharedChannels = 1
    end
    
    properties (Access = private)
        % Private property to hold the signal length
        pXSize
    end
    
    properties (Access = private, PositiveInteger) 
        % Private property that holds the buffer length
        pBufferLength
        pInitialBufferLength
    end
    
    properties (Access = private, PositiveInteger, Nontunable)
        pFixedBufferLength
    end
    
    properties (Access = private, Nontunable)
        pMaxNumInputSamples 
    end
    
    methods

        function obj = CircularBuffer(varargin)
            %CircularBuffer   Construct the CircularBuffer class.
            
            setProperties(obj, nargin, varargin{:},'BufferLength');
        end
        
    end
    
    methods(Hidden)
        
        function len = getBufferLength(obj)
            len = obj.pBufferLength;
        end
        
    end
    
    methods(Access = protected)
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            if ~obj.FixedLengthBuffer && ...  %Auto
                    (strcmp(prop, 'MaxNumInputSamplesSource') || ...
                    strcmp(prop, 'MaxNumInputSamples'))
                flag = true;
            elseif obj.FixedLengthBuffer && ...
                    strcmp(obj.MaxNumInputSamplesSource,'Auto') && ...
                    strcmp(prop, 'MaxNumInputSamples')
                flag = true;
            elseif strcmp(obj.BufferWidthSource,'Auto') && ...
                    strcmp(prop, 'BufferWidth')
                flag = true;
            else
                flag = false;
            end
        end
        
        function setupImpl(obj,x,M)
            
            setupImpl@phased.internal.AbstractVarSizeEngine(obj,x);
            
            x_sz = getPropagatedNumInputSamples(obj,x);
            if ~obj.FixedLengthBuffer
                obj.pInitialBufferLength = obj.BufferLength;
                N = obj.pInitialBufferLength;
                buflen = N+x_sz;
            else
                obj.pFixedBufferLength = obj.BufferLength;
                N = obj.pFixedBufferLength;
                if isInputDataSizePropagated(obj) || strcmp(obj.MaxNumInputSamplesSource,'Auto')
                    obj.pMaxNumInputSamples = x_sz;
                    buflen = N+obj.pMaxNumInputSamples;
                else
                    obj.pMaxNumInputSamples = obj.MaxNumInputSamples;
                    buflen = N+obj.pMaxNumInputSamples;
                end
            end            
            if strcmp(obj.BufferWidthSource,'Auto')
                obj.pNumInputChannels = getNumChannels(obj,x);
                obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            else
                obj.pNumInputChannels = obj.BufferWidth;
                obj.pValidatedNumInputChannels = obj.BufferWidth;
            end
            
            if ~isempty(coder.target) && ~obj.FixedLengthBuffer
                % Tell codegen that pBuffer is a variable-sized vector.
                if isreal(x)
                    obj.pBuffer = zeros([1 obj.pValidatedNumInputChannels]);
                else
                    obj.pBuffer = complex(zeros([1 obj.pValidatedNumInputChannels]));
                end
            end
            
            if isreal(x)
                obj.pBuffer = zeros([buflen obj.pValidatedNumInputChannels]);
            else
                obj.pBuffer = complex(zeros([buflen obj.pValidatedNumInputChannels]));
            end
            obj.pStart = ones(size(M));
            obj.pCount = zeros(size(M));
            
            if ~isscalar(M)
                M_sz = size(M);
                obj.pNumSharedChannels = obj.pValidatedNumInputChannels/M_sz(2);
            end

            obj.pBufferLength = buflen;
            obj.pXSize = [-1 getNumChannels(obj,x)];
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)  %#ok<INUSL>
            if index == 1
                flag = false;
            else % (index == 2) 
                flag = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function resetImpl(obj)
            obj.pBuffer(:) = 0;
            obj.pStart(:) = 1;
            obj.pCount(:) = 0;
        end
        
        function validateInputsImpl(obj,x,m)
            cond = ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:invalidInputDataType','X','double');
            end
            cond = ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:inputMustBeMatrix','X');
            end
            
            validateNumChannels(obj,x);
                        
            cond =  ~isa(m,'double');
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:invalidInputDataType','M','double');
            end
            cond =  ~isrow(m);
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:inputMustBeRowVector','M');
            end
            cond =  ~isreal(m);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:step:NeedReal', 'M');        
            end
            numDelays = numel(m);
            cond =  numDelays > 1 && ...
                numDelays ~= size(x,2) && rem(size(x,2),numDelays)~=0;
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:target:mustBeScalarOrSameColumns','M','X');
            end
        end
      
        function processInputSizeChangeImpl(obj,x,~)
            obj.pXSize(1) = size(x,1);
        end
    
        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSL>
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end

        function fsz_out = isOutputFixedSizeImpl(obj) 
            fsz_out = propagatedInputFixedSize(obj, 1);
        end
        
        function sz_out = getOutputSizeImpl(obj)
            sz_out = propagatedInputSize(obj,1);
        end        
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
            if isLocked(obj)
                s.pBufferLength = obj.pBufferLength;
                s.pBuffer = obj.pBuffer;
                s.pXSize = obj.pXSize;
                s.pStart = obj.pStart;
                s.pCount = obj.pCount;
                s.pMaxNumInputSamples = obj.pMaxNumInputSamples;
                s.pNumSharedChannels = obj.pNumSharedChannels;
                if ~obj.FixedLengthBuffer
                    s.pInitialBufferLength = obj.pInitialBufferLength;
                else
                    s.pFixedBufferLength = obj.pFixedBufferLength;
                end
            end
        end

        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            if isfield(s,'pXLength')
                obj.pXSize = [s.pXLength 1];
                s = rmfield(s,'pXLength');
            end
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function y = stepImpl(obj,x_in,Narg)
            
            % x can be a matrix with each column a channel. the delay Narg
            % can be either a scalar or a row vector matching the number of
            % columns in x
            x_in_length = size(x_in,1);
            if obj.FixedLengthBuffer 
                if x_in_length > obj.pMaxNumInputSamples
                    x = x_in(1:obj.pMaxNumInputSamples,:);
                else
                    x = x_in;
                end
            else
                x = x_in;
            end
            x_length = size(x,1);
            numDelays = numel(Narg);
            % preallocate            
            % if isreal(x)
            %     y = zeros(x_in_length,obj.pValidatedNumInputChannels);
            % else
            %     y = complex(zeros(x_in_length,obj.pValidatedNumInputChannels));
            % end
            y = zeros(x_in_length,obj.pValidatedNumInputChannels,'like',x);
         
            if numDelays == 1
                if ~isEmptyBuffer(obj)
                    popnum = min(obj.pCount,x_in_length);
                    popidx = translateIndex(obj,...
                        obj.pStart+(0:popnum-1));
                    
                    y(1:popnum,:) = obj.pBuffer(popidx,:);
                    obj.pBuffer(popidx,:) = 0;
                    obj.pStart = translateIndex(obj,popidx(end)+1);
                    obj.pCount = obj.pCount - popnum;
                end
                N = Narg;
                if ~obj.FixedLengthBuffer || N <= obj.pFixedBufferLength
                    
                    % working on the input signal, see if anything needs to be
                    % popped from the current input
                    overlap = min(x_length,x_in_length-N);
                    overlapidx = (1:overlap)+N;
                    y(overlapidx,:) = y(overlapidx,:) + x(1:overlap,:);
                    
                    % push the remaining content
                    if N > 0
                        % When N > 0, there are things needs to be pushed into the
                        % buffer.
                        
                        if overlap >= 0 % N is less than x_length
                            % calculate indices for where to push input
                            insertnode = 1:x_length-overlap;
                            xinsertidx = overlap+(1:x_length-overlap);
                            pushidx = translateIndex(obj,...
                                obj.pStart+insertnode-1);
                            % push, note that it combines with original content
                            obj.pBuffer(pushidx,:) = obj.pBuffer(pushidx,:) + ...
                                x(xinsertidx,:);
                            
                        else
                            if ~obj.FixedLengthBuffer && ...
                                    N > getBufferLength(obj)-x_length
                                % grow the buffer if necessary. Note that the
                                % output is the same length as input, so only
                                % the delay portion will occupy the buffer.
                                resize(obj,min(2*N+x_length,N+1024+x_length),numDelays);
                            end
                            insertnode = (1:x_length)-overlap;
                            xinsertidx = 1:x_length;
                            pushidx = translateIndex(obj,...
                                obj.pStart+insertnode-1);
                            % push, note that it combines with original content
                            obj.pBuffer(pushidx,:) = obj.pBuffer(pushidx,:) + ...
                                x(xinsertidx,:);
                            
                        end
                        
                        % update Count
                        % total content is always x_length + N, so for each
                        % input, the remaining content after output occupies N
                        % samples.
                        if N > obj.pCount % N == x_length - overlap (+ or -)
                            obj.pCount = N;
                        end
                    end
                    
                end

            else  % multiple delays
                K = obj.pNumSharedChannels;
                for m = 1:numDelays
                    % pop, clear the popped buffer, then update start and count
                    if ~isEmptyBuffer(obj,m)
                        % calculate the indices for what to pop
                        popnum = min(obj.pCount(m),x_in_length);
                        popidx = translateIndex(obj,...
                            obj.pStart(m)+(0:popnum-1));
                        popchan = (m-1)*K+(1:K);

                        y(1:popnum,popchan) = obj.pBuffer(popidx,popchan);
                        obj.pBuffer(popidx,popchan) = 0;
                        obj.pStart(m) = translateIndex(obj,popidx(end)+1);
                        obj.pCount(m) = obj.pCount(m) - popnum;
                    end
                end
                
                for m = 1:numDelays
                    N = Narg(m);
                    assert(N>=0);
                    if ~obj.FixedLengthBuffer || N <= obj.pFixedBufferLength

                        % working on the input signal, see if anything needs to be
                        % popped from the current input
                        overlap = min(x_length,x_in_length-N);
                        overlapidx = N+1:N+overlap;
                        overlapchan = (m-1)*K+(1:K);
                        y(overlapidx,overlapchan) = y(overlapidx,overlapchan) + x(1:overlap,overlapchan);

                        % push the remaining content
                        if N > 0
                            % When N > 0, there are things needs to be pushed into the
                            % buffer.

                            if overlap >= 0 % N is less than x_length
                                % calculate indices for where to push input
                                insertnode = 1:x_length-overlap;
                                xinsertidx = overlap+(1:x_length-overlap);
                                pushidx = translateIndex(obj,...
                                    obj.pStart(m)+insertnode-1);
                                % push, note that it combines with original content
                                obj.pBuffer(pushidx,overlapchan) = obj.pBuffer(pushidx,overlapchan) + ...
                                    x(xinsertidx,overlapchan);

                            else
                                if ~obj.FixedLengthBuffer && ...
                                        N > getBufferLength(obj)-x_length
                                    % grow the buffer if necessary. Note that the
                                    % output is the same length as input, so only
                                    % the delay portion will occupy the buffer.
                                    resize(obj,min(2*max(Narg)+x_length,max(Narg)+1024+x_length),numDelays);
                                end
                                insertnode = (1:x_length)-overlap;
                                xinsertidx = 1:x_length;
                                pushidx = translateIndex(obj,...
                                    obj.pStart(m)+insertnode-1);
                                % push, note that it combines with original content
                                obj.pBuffer(pushidx,overlapchan) = obj.pBuffer(pushidx,overlapchan) + ...
                                    x(xinsertidx,overlapchan);

                            end

                            % update Count
                            % total content is always x_length + N, so for each
                            % input, the remaining content after output occupies N
                            % samples.
                            if N > obj.pCount(m) % N == x_length - overlap (+ or -)
                                obj.pCount(m) = N;
                            end
                        end

                    end
                end
                % push the remaining content

                % shrink buffer if only small portion of buffer is used
                % if obj.pCount < floor(getBufferLength(obj)/4) && obj.pCount ~= 0
                %     resize(obj,2*obj.pCount);
                % end
            end
        end
    end
    
    methods(Access = private)
        
        function resize(obj,bufferlen,numDelays)
            
            % the size can only be changed along first dimension, not the
            % second dimension
            
            % create the buffer with the new size N
            xcols = obj.pValidatedNumInputChannels;
            if isreal(obj.pBuffer)
                tempBuffer = zeros(bufferlen,xcols);
            else
                tempBuffer = complex(zeros(bufferlen,xcols));
            end
            
            if numDelays == 1
                if ~isEmptyBuffer(obj)
                    % save current content
                    temp = readBuffer(obj);
                    % copy over the saved content
                    tempBuffer(1:obj.pCount,:) = temp;
                end
            else
                K = obj.pNumSharedChannels;
                for m = 1:numDelays
                    if ~isEmptyBuffer(obj,m)
                        chanidx = (m-1)*K+(1:K);
                        % save current content
                        temp = readBuffer(obj,m,chanidx);
                        % copy over the saved content
                        tempBuffer(1:obj.pCount(m),chanidx) = temp;
                    end
                end
            end
            
            % Resave buffer
            if isreal(obj.pBuffer)
                obj.pBuffer = tempBuffer;
            else
                obj.pBuffer = complex(tempBuffer);
            end
            
            obj.pBufferLength = bufferlen;
                        
            % reset the starting position, note that the count stays
            % the same
            obj.pStart(:) = 1;
            
        end
               
        function c = readBuffer(obj,idx,chanidx)
            if nargin < 2
                c = obj.pBuffer(translateIndex(...
                    obj,obj.pStart+(0:obj.pCount-1)),:);
            else
                c = obj.pBuffer(translateIndex(...
                    obj,obj.pStart(idx)+(0:obj.pCount(idx)-1)),chanidx);
            end
        end           
        
        function b = isEmptyBuffer(obj,idx)
            if nargin < 2
                b = all(obj.pCount == 0);
            else
                b = (obj.pCount(idx) == 0);
            end
        end

        function cidx = translateIndex(obj,idx)
            cidx = mod((idx-1),getBufferLength(obj))+1;
        end
        
    end
    
end


