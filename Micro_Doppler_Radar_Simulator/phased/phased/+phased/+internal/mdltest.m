function D = mdltest(eigenvals,K,fb)
%This function is for internal use only. It may be removed in the future.

% mdltest Estimate signal subspace dimension using MDL algorithm
% D = mdltest(EIGENVALS,K,FB) returns an estimate of the signal subspace
% dimension using the MDL or MDL-FB algorithm. EIGENVALS are the
% eigenvalues of the signal covariance matrix. The eigenvalues must be in
% descending order. K is the number of snapshots. FB specifies whether
% forward-backward averaging is used when obtaining EIGENVALS. The
% estimated signal subspace dimension D could be any value from 0 to
% numel(EIGENVALS).

%   Copyright 2010 The MathWorks, Inc.

% Reference:
% [1] Harry V. Trees, Optimal Array Processing, 2002

%#ok<*EMCA>
%#codegen

if any(eigenvals==0)
    coder.internal.errorIf(any(eigenvals==0), ...
        'phased:phased:internal:mdltest:ZeroEigvalues');
end


M = numel(eigenvals);

% detect the Signal Subspace Dimension
aux = zeros(M,1);
if fb
    for d=0:M-1,
        aux(d+1) = Ld(d,eigenvals,K,M) +  ...
            d*(2*M-d+1)*log(K)/4;  % Eq (7.511) in [1]
    end
else
    for d=0:M-1,
        aux(d+1) = Ld(d,eigenvals,K,M) + ...
            (d*(2*M-d)+1)*log(K)/2;% Eq (7.508) in [1]
    end
end
[~, idx] = min(aux); %  Eq (7.507) in [1]
D = idx-1;

end


function out = Ld(d,eigenvals,K,M)
% Compute Anderson's sufficient statistic - Eq. (7.502) in [1] Note
% the exponent is done before the product to avoid the numerical
% limitation. When the DOFs is large and many eigenvalues are
% small, the product term approach zero, which makes the log-term
% infinite. In such the case, the statistic breaks down. Thus, the
% exponent is done before the product to prevent the log of
% infinite.
out = ...
    K*(M-d)*log(mean(eigenvals(d+1:M))/...
    prod(eigenvals(d+1:M).^(1/(M-d))));
end
