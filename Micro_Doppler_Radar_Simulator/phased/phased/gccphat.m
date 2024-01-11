function [tau,r12,lags] = gccphat(x,varargin)
%GCCPHAT    GCC-PHAT delay estimation
%   TAU = gccphat(X,XREF) computes the time delay TAU between the signal X
%   and the reference signal XREF. The computation assumes the signal comes
%   from a single source. The delay is estimated from the maximum peak
%   location in the cross-correlation between X and XREF. The cross
%   correlation is computed using GCC-PHAT algorithm.
%
%   X can be an N-element column vector or an NxM matrix. N is the number
%   of time samples and M is the number of channels. If X is a matrix, each
%   column is a channel and the cross correlation is computed for each
%   channel. XREF can be an N-element column vector or an NxM matrix. If
%   XREF is a column vector, then all channels in X uses XREF as the
%   reference signal when computing the cross correlation. If XREF is a
%   matrix, its size must be identical to the size of X and the cross
%   correlation is computed between corresponding channels in X and XREF.
%   TAU is an M-element row vector having the same number of columns in X.
%   Each entry in TAU specifies the estimated delay for the corresponding
%   signal pairs in X and XREF.
%
%   TAU = gccphat(X,XREF,Fs) specifies the sampling frequency of the signal
%   as a positive scalar in Fs. The default value of Fs is 1. X and XREF
%   have the same sampling frequency.
%
%   [TAU,R,LAG] = gccphat(...) also returns the cross-correlation in R and
%   the lags for cross-correlation in LAG (in seconds). R is a (2*N-1)xM
%   matrix whose columns are cross-correlations for the corresponding
%   signal pairs in X and XREF. LAG is a (2*N-1)-element column vector
%   containing the corresponding lags for each sample in the
%   cross-correlation. All cross correlations have the same lags.
%
%   [...] = gccphat(X) or [...] = gccphat(X,Fs) returns the estimated
%   delays and cross correlations of all pairs of channels in X. If X has M
%   columns, the resulting TAU and R each has M^2 columns. In these M^2
%   columns, the first M columns contain the delays and cross correlations
%   using the first channel as the reference, the second M columns contain
%   the delays and cross correlations using the second channel as the
%   reference, and so on.
%
%   This function supports single and double precision for input data and
%   arguments. If the input data X is single precision, the output data
%   is single precision. If the input data X is double precision, the
%   output data is double precision. The precision of the output is
%   independent of the precision of the properties and other arguments.
%
%   % Example:
%   %   Estimate the delay between two channels.
%
%   load gong;
%   tau = 0.005;
%   xref = y;
%   x = delayseq(y,tau,Fs);
%
%   tau_est = gccphat(x,xref,Fs)
%
%   See also phased, phased.GCCEstimator, delayseq.

%   Copyright 2018 The MathWorks, Inc.

%#ok<*EMCA>
%#codegen

narginchk(1,3)

classtouse = class(x);

if nargin == 3
    xref = varargin{1};
    fs = varargin{2};
    refflag = true;
elseif nargin == 2
    if isscalar(varargin{1})
        fs = varargin{1};
        refflag = false;
        xref = 0; % define for codegen
    else
        xref = varargin{1};
        fs = 1;
        refflag = true;
    end
else %nargin==1
    refflag = false;
    fs = 1;
    xref = 0; % define for codegen
end

validateattributes(x,{'double','single'},{'2d','nonempty','finite'},'gccphat','A');
[N,M] = size(x);
sigdatatypes.validateFrequency(fs,'gccphat','Fs',{'double','single'},{'scalar'});
if refflag
    if iscolumn(xref)
        validateattributes(xref,{'double','single'},{'nonempty','finite','nrows',N},'gccphat','B');
    else
        validateattributes(xref,{'double','single'},{'nonempty','finite','size',[N M]},'gccphat','B');
    end
end

Ncorr = 2*N-1;
NFFT = 2^nextpow2(Ncorr);
lags = (-(Ncorr-1)/2:(Ncorr-1)/2).';
lags = cast((lags/fs),classtouse);
if refflag
    r12 = privgccphat(x,cast(xref,classtouse),NFFT,Ncorr);
else
    r12 = complex(zeros(Ncorr,M^2,classtouse));
    for m = 1:M
        r12(:,(m-1)*M+(1:M))=privgccphat(x,x(:,m),NFFT,Ncorr);
    end
end

% Assume single source
Npeaks = 1;
tau = zeros(Npeaks,size(r12,2),classtouse);
if Npeaks == 1
    [~,idx] = max(abs(r12));
    tau(1,:) = lags(idx);
end

function r12 = privgccphat(x,xref,NFFT,N)
R12 = bsxfun(@times,fft(x,NFFT),conj(fft(xref,NFFT)));
r12_temp = fftshift(ifft(exp(1i*angle(R12))),1);
r12 = r12_temp(NFFT/2+1-(N-1)/2:NFFT/2+1+(N-1)/2,:);


