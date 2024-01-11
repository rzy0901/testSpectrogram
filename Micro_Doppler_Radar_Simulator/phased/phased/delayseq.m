function y=delayseq(x,dtIn,fs)
%delayseq Delay or advance time sequence
%   Y = delayseq(X, DELAY) returns the delayed or advanced sequence Y by
%   applying DELAY to the input sequence X. DELAY (in samples) can be
%   integer or non-integer values. When it is negative, the sequence X is
%   advanced. X can be a vector or a matrix. DELAY is a scalar or a vector.
%
%   When X is a column vector, X is delayed by each element of DELAY and
%   the resulting sequence will be stored in corresponding column of Y.
%
%   When X has multiple columns, each column is delayed by corresponding
%   element of DELAY. If DELAY is a scalar, it will be applied to each
%   column of X.
%
%   The output sequence Y always has the same length as input with
%   appropriate truncations or zero padding.
%
%   Y = delayseq(X, DELAY, Fs) specifies DELAY in seconds. Fs is the
%   sampling frequency (in Hz).
%
%   This function supports single and double precision for input data and
%   arguments. If the input data X is single precision, the output data 
%   is single precision. If the input data X is double precision, the 
%   output data is double precision. The precision of the output is
%   independent of the precision of the arguments.
%
%   % Example:
%   %   Delay a sinusoid by 0.001 second.
%
%   fs = 8e3; t = 0:1/fs:0.005; x = sin(2*pi*100*t).';
%   y = delayseq(x,0.001,fs);
%   plot(t,x,'r',t,y,'b'); legend('Original','Delayed');
%
%   See also phased, phased.TimeDelayBeamformer.

%   Copyright 2010-2012 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, Wiley, 2002

%#codegen
%#ok<*EMCA>
    
phased.internal.narginchk(2,3,nargin);
% validate inputs
validateattributes(x, {'double','single'}, {'2d','nonempty','finite','nonnan'},...
    'delayseq', 'X');
validateattributes(dtIn, {'double','single'},...
    {'vector', 'real', 'finite', 'nonnan'}, 'delayseq', 'DELAY');

if nargin > 2
    validateattributes(fs, {'double','single'},...
        {'scalar', 'real', 'positive', 'finite', 'nonnan'},...
        'delayseq', 'Fs');
    % convert dt from time to samples
    dtIn = dtIn*cast(fs,'double'); %%codegen purpose
end


bSingleVecX = true;
if size(x,2)~=1         % not a column vector
    bSingleVecX = false;
    if (length(dtIn) == 1)    % scalar expansion of dt
        dt = dtIn*ones(size(x,2),1);
    else
        dt = dtIn;
    end
else
    dt = dtIn;
end

cond = size(x,2)~=1 && (size(x,2)~=length(dt));
if cond
    coder.internal.errorIf(cond,'phased:delayseq:MismatchedDelay');
end


% initialization
inputLength = size(x,1);
delayInt = round(dt);    % Integer delays in samples
delayFrac = dt - delayInt;
maxLength = inputLength+max(0, max(delayInt)); % maximum sequence length
% Define upperbound
if maxLength > 2*inputLength
    maxLengthLimit = 2*inputLength;
else
    maxLengthLimit = maxLength;
end
nfftLimit = 2^nextpow2(2*inputLength);
if isreal(x)
   output = zeros(maxLengthLimit,length(dt),class(x));
else
   output = complex(zeros(maxLengthLimit,length(dt),class(x)));
end

% Perform delay operation
for colI=1:length(dt)

    endIdx = inputLength+delayInt(colI);
    % no operation if delayed or advanced out of scope
    if (endIdx <= 0) || (endIdx > 2*inputLength)
        continue;
    end

    if bSingleVecX
        tmpx = x(:);
    else
        tmpx = x(:,colI);
    end
    
    if delayFrac(colI)    
        nfft = 2^nextpow2(inputLength + max(0, delayInt(colI)));
        assert(nfft <= nfftLimit);
        binStart = floor(nfft/2);
        % Notice the FFT bins must belong to [-pi, pi].
        fftBin = 2*pi*ifftshift(((0:nfft-1)-binStart).')/nfft;
        if isscalar(tmpx) && nfft==1
            tmpxd = tmpx;
        else
            tmpxd = fft(tmpx,nfft);
            tmpxd = ifft(tmpxd(:).*exp(-1i*dt(colI)*fftBin));
        end
        
        if delayInt(colI) >= 0
            orgStart = delayInt(colI) + 1;
            newStart = delayInt(colI) + 1;
        else
            orgStart = 1;
            newStart = 1;
        end
        orgEnd = endIdx;
        if isreal(x)
           output(newStart:endIdx,colI) = real(tmpxd(orgStart:orgEnd));
        else
           output(newStart:endIdx,colI) = tmpxd(orgStart:orgEnd);
        end
    else
        % Integer sampling shifting
        if delayInt(colI) >= 0
            orgStart = 1;
            newStart = delayInt(colI) + 1;
        else
            orgStart = 1-delayInt(colI);
            newStart = 1;
        end
        orgEnd = inputLength;
        output(newStart:endIdx,colI) = tmpx(orgStart:orgEnd);
    end
    
    
    
end

% output has the same length as input
y = output(1:inputLength, :);

end
