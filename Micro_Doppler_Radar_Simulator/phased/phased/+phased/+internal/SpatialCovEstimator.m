classdef (Sealed, Hidden, StrictDefaults) SpatialCovEstimator < phased.internal.AbstractVarSizeEngine
%This class is for internal use only. It may be removed in the future.

%ULACxEstimator Covariance matrix estimator for ULA

%   Copyright 2009-2018 The MathWorks, Inc.
%     
    
%   References:
%   [1] Optimum Array Processing, Part IV of Detection, Estimation, and
%       Modulation Theory, Harry Van Trees 2002, p.723
%
%   [2] D.A. Linebarger, R.D DeGroat and E.M. Dowling. Efficient
%       direction-finding methods employing forward/backward averaging.
%       IEEE Trans. Signal Process. vol.SP-42, pp.2136-2145, August 1994
%
%   [3] M Haardt and J.A. Nossek. Unitary Esprit: How to Obtain Increased
%       Estimation accuracy with a Reduced Computational Burden. IEEE 
%       Trans. Signal Process.,vol.SP-43, pp.1232-1242, May 1995


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable, PositiveInteger) 
    %NumSubarrays Number of subarrays for spatial smoothing
    %   Specify the number of subarrays used to average the covariance
    %   matrix for a uniform linear array as a positive scalar integer. For
    %   each additional subarray the method can handle one extra coherent
    %   source signal. The default value is 1, which is equivalent to no
    %   subarray averaging.
    NumSubarrays = 1;
end

properties (Nontunable, Logical) 
    %ForwardBackwardAveraging Forward-backward averaging
    %   Set this property to true to use forward-backward averaging to
    %   estimate the covariance matrix for sensor arrays with conjugate
    %   symmetric array manifold. The default value is true.
    ForwardBackwardAveraging = true;
    %UnitaryTransformed Output is unitary transformed
    %   Set this property to true to perform unitary transform on the 
    %   forward-backward averaged covariance matrix. The default value is
    %   false. This property applies when you set the
    %   ForwardBackwardAveraging property to true.
    UnitaryTransformed = false;
end

methods
    function obj = SpatialCovEstimator(varargin)
        setProperties(obj, nargin, varargin{:});
    end
end

methods (Access = protected)
    
    function validateInputsImpl(obj,x)
        % verify NumSubarrays no greater than channel number
        sz = size(x);
        cond = sz(2) < obj.NumSubarrays;
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:SpatialCovEstimator:InvalidSpatialSmoothing');
        end
    end
    
    function flag = isInputComplexityLockedImpl(~,~) 
        flag = false;
    end

    function flag = isInputSizeLockedImpl(~,~)
        flag = false;
    end

    function Sx = stepImpl(obj,X)
        
        X = cast(X,'double');%g1812952
        L = cast(obj.NumSubarrays,class(X));
        if obj.ForwardBackwardAveraging
            if obj.UnitaryTransformed
                Sx = privComputeSq(X,L);
            else
                Sx = UnitaryCovMtx(X,L);
            end
        else
            Sx = MLCovMtx(X,L);
        end
        
    end

    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        if ~obj.ForwardBackwardAveraging && strcmp(prop,'UnitaryTransformed')
            flag = true;
        end
    end
    
end

end

function Sx = MLCovMtx(X,L)
    % Maximum Likelihood
    K = size(X,1);
    M = size(X,2)-L+1;
    Sx = complex(zeros(M,M));
    for n=1:L
        % Forward-Only Spatial Smoothing
        subOut = X(:,n:n+M-1).';
        Sx =  Sx + 1/K*(subOut*subOut');   % Eq (7.3) and (9.288) in [1]
    end
    Sx = 1/L*Sx;
end
    
function Sx = UnitaryCovMtx(X,L)
    if isreal(X)
        % Requires only 1/2 of the computations
        Sx = local_computesx_realdata(X,L);
    else
        % Compute real matrix Sq
        Sq = privComputeSq(X,L);
        % Unitary transform
        Q = phased.internal.unitarymat(size(Sq,1));
        Sx = Q*Sq*Q';                             % Eq (7.66) in [1]
    end

end
    
function Sq = privComputeSq(X,L)
%privComputeSq Real Symmetric Matrix Sq
%   Sq = privComputeSq(X,L) compute the data matrix Sq as an average of L
%   subarrays. Sq is an always real, symmetric forward-backward  matrix Sq
%   from the snapshots X. If X is real, Sq is block-diagonal.

    [K, N] = size(X);
    M = N-L+1;

    Sq = zeros(M);

    % Sq must be real
    if rem(M,2)
        half = (M-1)/2; % M odd
        for n=1:L      % Spatial Smoothing
            X1 = X(:,n:n+half-1).';                % half x K
            X2 = X(:,n+half+1:n+2*half).';         % half x K
            Xmid = X(:,n+half).';                  % 1    x K
            Z1 = X1+flipud(X2); Z2 = X1-flipud(X2);        % half x K
            % Eq 7 in [3]
            Zfb = [real(Z1) -imag(Z1); ...
                sqrt(2)*real(Xmid) -sqrt(2)*imag(Xmid) ; ...
                imag(Z2) real(Z2)];                        % M x 2K
            Sq = Sq + 1/(2*K)*(Zfb*Zfb.');                   % M x M
        end
    else
        half = M/2; % M even
        for n=1:L  % Spatial Smoothing
            X1 = X(:,n:n+half-1).';                % half x K
            X2 = X(:,n+half:n+2*half-1).';         % half x K
            Z1 = X1+flipud(X2); Z2 = X1-flipud(X2);        % half x K
            % Eq (27) in [2]
            Zfb = [real(Z1) -imag(Z1); imag(Z2) real(Z2)]; % M x 2K
            Sq = Sq + 1/(2*K)*(Zfb*Zfb.');                   % M x M
        end
    end

    Sq = 1/L*Sq;

end
    
function Sx = local_computesx_realdata(X,L)
% Sq must be real and block diagonal.
% Sq construction requires only half the computations compared to
% the complex case

[K,M] = size(X);
M = M-L+1;
Sx = zeros(M);
if rem(M,2)
    half = (M-1)/2; % M odd
    for n=1:L      % Spatial Smoothing
        X1 = X(:,n:n+half-1).';                % half x K
        X2 = X(:,n+half+1:n+2*half).';         % half x K
        Xmid = X(:,n+half).';                  % 1    x K
        Z11 = X1+flipud(X2); Z22 = X1-flipud(X2);      % half x K
        B11 = Z11*Z11.'; B22 = Z22*Z22.';              % half x half
        Sq1 = B11+B22; Sq2 = B11-B22;                  % half x half
        C = 2*Xmid*Z11.';
        D = 4*(Xmid*Xmid.');
        Sx = Sx + 1/(4*K)*[Sq1         C.'         fliplr(Sq2); ...
                           C           D           fliplr(C); ...
                           flipud(Sq2) flipud(C.') rot90(Sq1,2)];
    end
else
    half = M/2; % M even
    for n=1:L  % Spatial Smoothing
        X1 = X(:,n:n+half-1).';                % half x K
        X2 = X(:,n+half:n+2*half-1).';         % half x K
        Z11 = X1+flipud(X2); Z22 = X1-flipud(X2);      % half x K
        B11 = Z11*Z11.'; B22 = Z22*Z22.';              % half x half
        Sq1 = B11+B22; Sq2 = B11-B22;                  % half x half
        % Eq (30) in [2]
        Sx = Sx + 1/(4*K)*[Sq1 fliplr(Sq2); flipud(Sq2) rot90(Sq1,2)];
    end
end
Sx = 1/L*Sx;

end
