classdef (Sealed,StrictDefaults) RootWSFEstimator < phased.internal.AbstractULASubspaceDOA
%RootWSFEstimator Root WSF direction of arrival (DOA) estimator
%   H = phased.RootWSFEstimator creates a root WSF DOA estimator System
%   object, H. The object estimates the signal's direction of arrival using
%   the root weighted subspace fitting (WSF) algorithm with a uniform
%   linear array (ULA).
%
%   H = phased.RootWSFEstimator(Name,Value) creates a root WSF DOA
%   estimator object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   ANG = step(H,X) estimates the DOAs from X using the DOA estimator H. X
%   is a matrix whose columns correspond to channels. ANG is a row vector
%   of the estimated broadside angles (in degrees).
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   RootWSFEstimator methods:
%
%   step     - Perform DOA estimation (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a root WSF DOA estimator object with same property 
%              values
%   isLocked - Locked status (logical)
%
%   RootWSFEstimator properties:
%
%   SensorArray           - Sensor array
%   PropagationSpeed      - Signal propagation speed
%   OperatingFrequency    - Operating frequency
%   NumSignalsSource      - Source of number of signals
%   NumSignalsMethod      - Method to estimate number of signals
%   NumSignals            - Number of signals
%   Method                - Iterative method
%   MaximumIterationCount - Maximum number of iterations
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   % Estimate the DOAs of two signals received by a standard 10-element 
%   % ULA with element spacing 1 meter. The antenna operating frequency is
%   % 150 MHz. The actual direction of the first signal is 10 degrees in
%   % azimuth and 20 degrees in elevation. The direction of the second 
%   % signal is 45 degrees in azimuth and 60 degrees in elevation.
%
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   array = phased.ULA('NumElements',10,'ElementSpacing',1);
%   array.Element.FrequencyRange = [100e6 300e6];
%   fc = 150e6;
%   x = collectPlaneWave(array,[x1 x2],[10 20;45 60]',fc);
%   rs = RandStream('mt19937ar','Seed',0);
%   noise = 0.1/sqrt(2)*(randn(rs,size(x))+1i*randn(rs,size(x)));
%   estimator = phased.RootWSFEstimator('SensorArray',array,...
%                 'OperatingFrequency',fc,...
%                 'NumSignalsSource','Property','NumSignals',2);
%   doas = estimator(x+noise);
%   az = broadside2az(sort(doas),[20 60])
%
%   See also phased, phased.RootMUSICEstimator, broadside2az.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties
    %MaximumIterationCount Maximum number of iterations
    %   Specify the maximum number of iterations as a positive integer
    %   scalar or 'Inf'. The default is 'Inf'. This property is tunable.
    MaximumIterationCount = Inf;
end

properties (Nontunable)
    %Method Iterative method
    %   Specify the iterative method as one of 'IMODE' | 'IQML'. The
    %   default is 'IMODE'.
    Method = 'IMODE';
end
    
properties(Constant = true, Hidden = true)
    MethodSet = matlab.system.StringSet({'IMODE','IQML'});
end

methods

    function obj = RootWSFEstimator(varargin)
        obj = obj@phased.internal.AbstractULASubspaceDOA(varargin{:});
    end
            
    function obj = set.MaximumIterationCount(obj,var)
        if isfinite(var)
            validateattributes( var, { 'double','single' }, { 'scalar', 'integer', 'positive' }, '', 'MaximumIterationCount');
        elseif var ~= inf
            coder.internal.errorIf(var ~= inf,'phased:phased:RootWSFEstimator:InvalidMaxIteration');
        end
        obj.MaximumIterationCount = var; 
    end
 
end

methods (Access = protected)
    
    function setupImpl(obj, X) 
        setupImpl@phased.internal.AbstractULASubspaceDOA(obj,X);
        % Always compute the maximum likelihood cov matrix.
        % No need for spatial smoothing
        obj.cCovEstimator.ForwardBackwardAveraging = false;
    end
    
    function loadObjectImpl(obj,s,wasLocked)
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function doas = stepImpl(obj, X)

        K = obj.pNumSnapshots;
        N = obj.pEffChannels;
        
        % Compute the eigenvectors/eigenvalues of the spatial
        % covariance matrix
        Cx = step(obj.cCovEstimator,X);
        [eigenvals eigenvects] = privEig(obj,Cx);
        
        % Obtain Signal subspace dimension
        D = getNumSignals(obj,eigenvals,K,false);
        if D==0 && isempty(coder.target)
           warning(message('phased:phased:doa:ZeroSourceNumber'));
           doas = zeros(1,0,class(X));
           return;
        end
                            
        % compute constants        
        if (obj.Method(2) == 'M') % 'IMODE'
            % Define handle to the function that will compute Qd            
            Gamasv = eigenvals(1:D);
            sigmaw2 = sum(eigenvals(D+1:N))/(N-D);  % Eq (8.525) in [1]
            Waov = (Gamasv-sigmaw2).^2./Gamasv;  % Eq (8.523) in [1]
            Us = bsxfun(@times,eigenvects(:,1:D),sqrt(Waov(:).'));             
            IN = Us.';
            K = D;
        else
            IN = X;
        end
%        switch lower(obj.Method),
%            case 'imode',
%                % Define handle to the function that will compute Qd
%                Gamas = diag(eigenvals(1:D));
%                sigmaw2 = sum(eigenvals(D+1:N))/(N-D);  % Eq (8.525) in [1]
%                Wao = (Gamas-sigmaw2*eye(D)).^2/Gamas;  % Eq (8.523) in [1] NsigxNsig
%                Us = eigenvects(:,1:D)*sqrtm(Wao);
%                f = @(B)local_compute_Qd(B,Us,N,D); 
%
%            case 'iqml',
%                % Define handle to the function that will compute Qx
%                f = @(B)local_compute_Qx(B,X,N,D,K); 
%        end        
        T = transformation_matrix(D); % Transformation matrix
        
        % Step 1: Initialization
        [m,test,B,b] = local_init(N,D);

        maxiter = obj.MaximumIterationCount;
        while test && m<maxiter
            % Step 2: Compute Q
            f = local_compute_Qx(B,IN,N,D,K); 
            Q = T'*f*T; % Eq (8.539) IMODE or Eq (8.514) IQML in [1]

            % Step 3: Compute the Eigenvectors of Real(Q)
            [V,E] = eig(real(Q));
            [~,idx] = min(diag(E));
            c = V(:,idx);                        % Eq (8.540) in [1]

            % Step 4: Find roots of polynomial b
            aux = b;
            b = T*c;                             % Eq (8.541) in [1]
            b_conj = conj(flipud(b));
            B = toeplitz([b_conj;zeros(N-D-1,1)],[b_conj(1) zeros(1,N-D-1)]);
            m=m+1;                               % Number of iterations
            test = (norm(aux-b)>sqrt(eps));      % Eq (8.543) in [1]
        end

        doas = angle(roots(cast(b,class(X))));                  % Eq (8.542) in [1]

        % Convert angles
        doas = privConvertAngles(obj,doas);
        
    end
    
    end

    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractULASubspaceDOA;
            props = {'Method','MaximumIterationCount'};
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:RootWSFEstimatorTitle')),...
                'Text',getString(message('phased:library:block:RootWSFEstimatorDesc')));
        end    
    end
    methods (Access = protected) %for Simulink
      function str = getIconImpl(obj) %#ok<MANU>
          str = sprintf('Root WSF\nDOA');        
      end    
    end

end

function [m,test,B,b] = local_init(N,D)
% Initialization 

    m = 0;                       % Counter
    test = true;
    B = complex([eye(N-D);zeros(D,N-D)]); % Eq (8.538) in [1]
    b = complex([zeros(D,1);1]);

end

% function Qd = local_compute_Qd(B,Us,N,D)
% % Us,N,D are constants
%
%    R = B'*B;
%    Qd = zeros(D+1,D+1);
%    for d = 1:D, % Loop over the dimension of the signal subspace
%        Ad = zeros(N-D,D+1);
%        for n = 1:N-D,
%            idx = n+D:-1:n;
%            Ad(n,:) = Us(idx,d);       % Eq (8.533) in [1]
%        end
%         Qd = Qd + Ad'/R*Ad;
%    end
%
%end

function Qx = local_compute_Qx(B,X,N,D,K)
% X,N,D are constants

    R = B'*B;
    Qx = complex(zeros(D+1,D+1));
    for k = 1:K, % Loop over the number of snapshots
        Ak = complex(zeros(N-D,D+1));
        for n = 1:N-D,
            idx = n+D:-1:n;
            Ak(n,:) = X(k,idx);       % Eq (8.498) in [1]
        end
        Qx = Qx + Ak'/R*Ak;
    end

end

function T = transformation_matrix(D)
% Transformation Matrix T is (D+1) x (D+1)

    if rem(D,2),
        % D odd
        I = eye((D+1)/2);
        J = fliplr(I);   % Exchange matrix
        T = 1/sqrt(2)*[kron(I,[1 1i]);kron(J,[1 -1i])]; % Eq (8.506) in [1]
    else
        % D even
        I = eye(D/2);
        J = fliplr(I);
        T = 1/sqrt(2)*[...
            [kron(I,[1 1i]) zeros(D/2,1)];...
            [zeros(1,D) 1];...
            [kron(J,[1 -1i]) zeros(D/2,1)]];
    end

end
