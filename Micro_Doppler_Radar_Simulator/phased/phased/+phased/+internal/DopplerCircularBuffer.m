classdef (Hidden, Sealed, StrictDefaults) DopplerCircularBuffer < phased.internal.AbstractVarSizeEngine
%This class is for internal use only. It may be removed in the future.

%DopplerCircularBuffer   Circular Buffer for Doppler
%   H = phased.internal.DopplerCircularBuffer(N) creates a circular buffer System
%   object, H, with the buffer length set to N. The object simulates the
%   circular buffering behavior.
%
%   Step syntax:
%
%   Y = step(H,X,M,N) pushes X into the buffer's (M+1)th position and then
%   outputs Y, whose dimension N, from the beginning of the buffer. The
%   buffer will dynamically grow if the current buffer cannot accommodate X
%   when it is inserted at (M+1)th position.
%
%   X can be a matrix. Each column of X is considered to be a separate
%   signal.
%
%   CircularBuffer methods:
%
%   step  - See above description for use of this method
%   reset - Reset the buffer
%

%   Copyright 2016 The MathWorks, Inc.


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
    
    properties(Constant, Hidden)
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
    
    methods

        function obj = DopplerCircularBuffer(varargin)
            %CircularBuffer   Construct the CircularBuffer class.
            setProperties(obj, nargin, varargin{:},'BufferLength');
        end
        
        function displayBuffer(obj)
            disp(obj.pBuffer);
        end
        
    end
    
    methods(Hidden)
        
        function len = getBufferLength(obj)
            len = obj.pBufferLength;
        end
        
    end
    
    methods(Access = protected)
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 3;
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
                flag = false;
        end
        
        function setupImpl(obj,x,M,~)
            
            setupImpl@phased.internal.AbstractVarSizeEngine(obj,x);
            
            if ~obj.FixedLengthBuffer
                obj.pInitialBufferLength = obj.BufferLength;
                N = obj.pInitialBufferLength;
                buflen = N;
            else
                obj.pFixedBufferLength = obj.BufferLength;
                buflen = obj.pFixedBufferLength;
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
        
        function validateInputsImpl(obj,x,m,n)
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
            cond =  numDelays > 1 && numDelays ~= size(x,2);
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:target:mustBeScalarOrSameColumns','M','X');
            end
            
            cond =  ~isa(n,'double');
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:invalidInputDataType','N','double');
            end
            cond =  ~isscalar(n);
            if cond
                coder.internal.errorIf(cond,...
                     'MATLAB:system:inputMustBeScalar','N');
            end
            cond =  ~isreal(n);
            if cond
                coder.internal.errorIf(cond,...
                     'phased:step:NeedReal', 'N');        
            end
            

        end
      
        function processInputSizeChangeImpl(obj,x,~,~)
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

        function y = stepImpl(obj,x_in,Narg,Nout)
            
            % Check that Narg is nonnegative
            cond = any(Narg < 0);
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:expectedNonnegative','Narg')
            end
            
            % Check that Nout is nonnegative
            cond = Nout < 0;
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:expectedNonnegative','Nout')
            end
            
            numDelays = numel(Narg);
            x_in_length = size(x_in,1);
            y = zeros(Nout,obj.pValidatedNumInputChannels,'like',x_in); 
                            
            if numDelays == 1
              
                if ~isEmptyBuffer(obj)
                    popnum = min(obj.pCount,Nout);
                    popidx = translateIndex(obj,...
                        obj.pStart+(0:popnum-1));
                    
                    y(1:popnum,:) = obj.pBuffer(popidx,:);
                    obj.pBuffer(popidx,:) = 0;
                    obj.pStart = translateIndex(obj,popidx(end)+1);
                    obj.pCount = obj.pCount - popnum;
                end
                
                if obj.FixedLengthBuffer 
                    if x_in_length-Nout+Narg > obj.pFixedBufferLength
                        % Delay is too large. Don't propagate signal.
                        x = [];
                    else
                        x = x_in;
                    end
                else
                    x = x_in;
                end

                x_length = size(x,1);
                
                N = Narg;
                M = x_length - Nout + N; % Number of buffer samples to be occupied (including loading zeros). Push onto the buffer when M > 0. 
                
                % working on the input signal, see if anything needs to be
                % popped from the current input
                overlap = min(x_length,Nout-N);
                overlapidx = N+1:N+overlap;
                y(overlapidx,:) = y(overlapidx,:) + x(1:overlap,:);

                % push the remaining content
                if M > 0
                    % There are things needs to be pushed into the
                    % buffer. First, check if we need to resize the buffer.
                    % x_length can be arbitrarily large.

                    if  ~obj.FixedLengthBuffer && M > getBufferLength(obj)
                        % grow the buffer if necessary. Note M samples will
                        % occupy the buffer.
                        resize(obj,min(2*M+x_length,M+1024+x_length),numDelays);
                    end
                        
                    if overlap >= 0 
                        % calculate indices for where to push input
                        insertnode = 1:x_length-overlap;
                        xinsertidx = overlap+(1:x_length-overlap);
                        pushidx = translateIndex(obj,...
                            obj.pStart+insertnode-1);
                        % push, note that it combines with original content
                        obj.pBuffer(pushidx,:) = obj.pBuffer(pushidx,:) + ...
                            x(xinsertidx,:);

                    else
                        % negative overlap means all samples are stored
                        insertnode = (1:x_length)-overlap;
                        xinsertidx = 1:x_length;
                        pushidx = translateIndex(obj,...
                            obj.pStart+insertnode-1);
                        % push, note that it combines with original content
                        obj.pBuffer(pushidx,:) = obj.pBuffer(pushidx,:) + ...
                            x(xinsertidx,:);

                    end

                    % update Count - M is the number of occupied cells
                    % in the current step call (x_length is not fixed).
                    if M > obj.pCount 
                        obj.pCount = M; 
                    end
                end
                    
            else  % multiple delays
              
                for m = 1:numDelays
                  
                    % pop, clear the popped buffer, then update start and count
                    if ~isEmptyBuffer(obj,m)
                        % calculate the indices for what to pop
                        popnum = min(obj.pCount(m),Nout);
                        popidx = translateIndex(obj,...
                            obj.pStart(m)+(0:popnum-1));

                        y(1:popnum,m) = obj.pBuffer(popidx,m);
                        obj.pBuffer(popidx,m) = 0;
                        obj.pStart(m) = translateIndex(obj,popidx(end)+1);
                        obj.pCount(m) = obj.pCount(m) - popnum;
                    end
                end
                
                for m = 1:numDelays
                    
                    N = Narg(m);
                             
                    if obj.FixedLengthBuffer 
                        if x_in_length-Nout+N > obj.pFixedBufferLength
                            x = [];
                        else
                            x = x_in(:,m);
                        end
                    else
                        x = x_in(:,m);
                    end                 
                    
                    x_length = size(x,1);
                    
                    M = x_length-Nout+N; % Number of buffer samples to be occupied (including leading zeros). Push onto the buffer when M > 0.
                    
                    % working on the input signal, see if anything needs to be
                    % popped from the current input
                    overlap = min(x_length,Nout-N);
                    overlapidx = N+1:N+overlap;
                    y(overlapidx,m) = y(overlapidx,m) + x(1:overlap);

                    % push the remaining content
                    if M > 0
                        % There are things needs to be pushed into the
                        % buffer.

                        if  ~obj.FixedLengthBuffer && M > getBufferLength(obj)
                            % grow the buffer if necessary. Note that the
                            % output is the same length as input, so only
                            % the delay portion will occupy the buffer.
                            resize(obj,min(2*M+x_length,M+1024+x_length),numDelays); 
                        end
                            
                        if overlap >= 0 
                            % calculate indices for where to push input
                            insertnode = 1:x_length-overlap;
                            xinsertidx = overlap+(1:x_length-overlap);
                            pushidx = translateIndex(obj,...
                                obj.pStart(m)+insertnode-1);
                            % push, note that it combines with original content
                            obj.pBuffer(pushidx,m) = obj.pBuffer(pushidx,m) + ...
                                x(xinsertidx);
                        else
                            insertnode = (1:x_length)-overlap;
                            xinsertidx = 1:x_length;
                            pushidx = translateIndex(obj,...
                                obj.pStart(m)+insertnode-1);
                            % push, note that it combines with original content
                            obj.pBuffer(pushidx,m) = obj.pBuffer(pushidx,m) + ...
                                x(xinsertidx);
                        end

                        % update Count
                        if M > obj.pCount(m) 
                            obj.pCount(m) = M;
                        end
                   end
               end
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
                for m = 1:numDelays
                    if ~isEmptyBuffer(obj,m)
                        % save current content
                        temp = readBuffer(obj,m);
                        % copy over the saved content
                        tempBuffer(1:obj.pCount(m),m) = temp;
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
               
        function c = readBuffer(obj,idx)
            if nargin < 2
                c = obj.pBuffer(translateIndex(...
                    obj,obj.pStart+(0:obj.pCount-1)),:);
            else
                c = obj.pBuffer(translateIndex(...
                    obj,obj.pStart(idx)+(0:obj.pCount(idx)-1)),idx);
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
