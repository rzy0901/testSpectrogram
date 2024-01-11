classdef (Sealed,StrictDefaults) RootMUSICEstimator < phased.internal.AbstractULASubspaceDOA
%RootMUSICEstimator Root MUSIC direction of arrival (DOA) estimator
%   H = phased.RootMUSICEstimator creates a root MUSIC DOA estimator System
%   object, H. The object estimates the signal's direction of arrival using
%   the root MUSIC algorithm with a uniform linear array (ULA). When a
%   uniform circular array (UCA) is used, the algorithm transforms the
%   input to ULA like structure using the phase mode excitation technique.
%
%   H = phased.RootMUSICEstimator(Name,Value) creates a root MUSIC DOA
%   estimator object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   ANG = step(H,X) estimates the DOAs from X using the DOA estimator H. X
%   is a matrix whose columns correspond to channels. ANG is a row vector
%   of the estimated broadside angles (in degrees). This syntax is only
%   applicable when the SensorArray property is a ULA.
%
%   ANG = step(H,X,ElAng) specifies the elevation angle (in degrees) of the
%   signals. ElAng is a scalar between -90 and 90. ANG is a row vector of
%   the estimated azimuth angles (in degrees). This syntax is only
%   applicable when the SensorArray property is a UCA. Phase mode
%   excitation assumes same and known elevation for all signals of
%   interest.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   RootMUSICEstimator methods:
%
%   step     - Perform DOA estimation (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a root MUSIC DOA estimator object with same property 
%              values
%   isLocked - Locked status (logical)
%
%   RootMUSICEstimator properties:
%
%   SensorArray              - Sensor array
%   PropagationSpeed         - Signal propagation speed
%   OperatingFrequency       - Operating frequency
%   NumSignalsSource         - Source of number of signals
%   NumSignalsMethod         - Method to estimate number of signals
%   NumSignals               - Number of signals
%   ForwardBackwardAveraging - Forward-backward averaging
%   SpatialSmoothing         - Spatial smoothing
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Examples:
%
%   % Example 1:
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
%   estimator = phased.RootMUSICEstimator('SensorArray',array,...
%                 'OperatingFrequency',fc,...
%                 'NumSignalsSource','Property','NumSignals',2);
%   doas = estimator(x+noise);
%   az = broadside2az(sort(doas),[20 60])
%
%   % Example 2: 
%   % Estimate the azimuth of two signals received by a 15-element UCA with
%   % a 1 meter radius. The signals are assumed to be at 0 degrees
%   % elevation. The antenna operating frequency is 150 MHz. The actual
%   % direction of the first signal is 10 degrees in azimuth and 4 degrees
%   % in elevation. The direction of the second signal is 45 degrees in
%   % azimuth and -2 degrees in elevation.
% 
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   array = phased.UCA('NumElements',15,'Radius',1);
%   fc = 150e6;
%   x = collectPlaneWave(array,[x1 x2],[10 4;45 -2]',fc);
%   rs = RandStream('mt19937ar','Seed',0);
%   noise = 0.1/sqrt(2)*(randn(rs,size(x))+1i*randn(rs,size(x)));
%   estimator = phased.RootMUSICEstimator('SensorArray',array,...
%                  'OperatingFrequency',fc,...
%                  'NumSignalsSource','Property','NumSignals',2);
%   doas = estimator(x+noise,0);
%   az = sort(doas)
%
%   See also phased, phased.RootWSFEstimator, broadside2az.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] Mathews, C.P., Zoltowski, M.D., "Eigenstructure techniques for 2-D
%   angle estimation with uniform circular arrays". IEEE Transaction on
%   Signal Processing, Vol. 42 No.9, 1994, pp. 2395-2407


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Dependent, Nontunable)
    %SpatialSmoothing Spatial smoothing
    %   Specify the effective reduction in the size of the sensor array due
    %   to spatial smoothing as a nonnegative integer. If the ULA consists
    %   of M elements and the value of this property is L, maximally
    %   overlapped subarrays of M-L elements are formed. For a UCA, M is
    %   the number of the internal ULA like structure.
    %
    %   Covariance matrices are estimated for each subarray of M-L elements
    %   and averaged together to produce a covariance matrix of size
    %   (M-L)x(M-L). The maximum value of L is M-2 resulting in subarrays
    %   consisting of two elements.
    %   
    %   Note that each additional increment in this property (decorrelates)
    %   handles one additional coherent source, but reduces the effective
    %   size of the array aperture. The default value of this property is 0
    %   resulting in no spatial smoothing.
    SpatialSmoothing 
end

properties (Nontunable, Logical) 
    %ForwardBackwardAveraging Forward-backward averaging
    %   Set this property to true to use forward-backward averaging to
    %   estimate the covariance matrix for sensor arrays with conjugate
    %   symmetric array manifold. The default value is false.
    ForwardBackwardAveraging = false;
end
    

% Following properties are only applicable to UCA
properties (Access = protected, Nontunable)
 pPhasedModeTransform = false;
 % K * Radius
 pKR;
 %beamspace or ula space
 pIsBeamspace = true;
end
properties (Access = protected)
 % phase mode beamformer
 pFeH;   
 % Holds either phase mode beamformer
 % or the elevation dependent ULA space
 % beamformer weights
 pConjugatedBFWeights
 % Bessel function diagonal matrix
 pJ
 % last used elevation angle
 pLastElAng; 
end

methods
    function set.SpatialSmoothing(obj,value)
        validateattributes( value, { 'double','single' }, { 'scalar', 'integer', 'nonnegative', 'finite' }, '', 'SpatialSmoothing');
        obj.pSpatialSmoothing = value;
    end
    
    function value = get.SpatialSmoothing(obj)
        value = obj.pSpatialSmoothing;
    end
    
end

methods
    function obj = RootMUSICEstimator(varargin)
        obj = obj@phased.internal.AbstractULASubspaceDOA(varargin{:});
    end
end
methods (Access = private)
    function maxPhaseMode = getMaxPhaseMode(obj)
           wavelength = obj.PropagationSpeed/obj.OperatingFrequency;
           k = 2*pi/wavelength;
           maxPhaseMode = floor(k*obj.SensorArray.Radius); %Eq 3 in [2]
    end
    function updateBeamformer(obj,elAng)
    % update beamformer if elevation changed since last call
        if elAng ~= obj.pLastElAng
            if elAng >= 90 || elAng <= -90
                sigdatatypes.validateAngle(elAng,'','ElAng',...
                                           {'double','single'},{'scalar','<=',90,'>=',-90});
            end
            %our el is starts from x instead of z hence cos instead of sin
            Nuca = getNumElements(obj.SensorArray);
            kr = obj.pKR;
            Jarg = kr*cosd(elAng);
            maxPhaseMode = floor(kr); %Eq 3 in [2]
            halfdiag = real(besselj((1:maxPhaseMode),Jarg));
            besselj0 = real(besselj(0,Jarg));
            obj.pJ = ... %Eq 14 in [2]
                sqrt(Nuca)*diag([fliplr(halfdiag) besselj0 halfdiag]);
            if ~obj.pIsBeamspace 
                %convert from beamspace to ula like space
                %this is needed for spatial smoothing 
                %all sources are assumed to be confined to a given elevation angle.
                obj.pConjugatedBFWeights = obj.pFeH/obj.pJ; %Eq 12,23 in [2]
            end
            obj.pLastElAng = elAng;
        end
    end
end
methods (Access = protected)
    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractULASubspaceDOA(obj);
        if isa(obj.SensorArray,'phased.UCA')         
            cond = obj.SpatialSmoothing && strcmp(obj.NumSignalsSource,'Auto');
            if cond
                coder.internal.errorIf(cond,'phased:phased:RootMUSICEstimator:NoUcaAICSmoothing','aictest','mdltest','X');
            end            
        end
    end 
    function validateInputsImpl(obj,x,varargin)
        validateInputsImpl@phased.internal.AbstractULASubspaceDOA(obj,x);
         if isa(obj.SensorArray,'phased.UCA')
             elAng = varargin{1};
             sigdatatypes.validateAngle(elAng,'','ElAng',...
                    {'double','single'},{'scalar','<=',90,'>=',-90});
         end
    end
    
    function processInputSizeChangeImpl(obj,x, varargin)
        processInputSizeChangeImpl@phased.internal.AbstractULASubspaceDOA(obj,x);
    end
    
    function setupImpl(obj,x,varargin)
       setupImpl@phased.internal.AbstractULASubspaceDOA(obj,x);
       obj.cCovEstimator.ForwardBackwardAveraging = ...
           obj.ForwardBackwardAveraging;
       % apply unitary transform to the covariance matrix when
       % forward-backward averaging is true.
       if obj.cCovEstimator.ForwardBackwardAveraging
           obj.cCovEstimator.UnitaryTransformed  = true;
       end
       sensorArray = obj.SensorArray;
       if isa(sensorArray,'phased.UCA')
           obj.pPhasedModeTransform = true;
           wavelength = obj.PropagationSpeed/obj.OperatingFrequency;
           k = 2*pi/wavelength;
           obj.pKR = k*sensorArray.Radius;
           Nuca = getNumElements(obj.SensorArray);
           maxPhaseMode = getMaxPhaseMode(obj);      
           V = getPhasedModeW(Nuca,-maxPhaseMode:maxPhaseMode);%Eq 11 in [2]
           Cj = diag(1i.^[-maxPhaseMode:0 -1:-1:-maxPhaseMode]); %Eq 10 in [2] 
           %phase mode beamformer
           obj.pFeH = (Cj*V').'; %Eq 9 in [2]  
           obj.pConjugatedBFWeights = obj.pFeH;
           if obj.SpatialSmoothing
               obj.pIsBeamspace = false;
           end
           %Initialize to an invalid value to indicate nothing cached yet
           obj.pLastElAng = 91;
           obj.pJ = zeros(2*maxPhaseMode+1,2*maxPhaseMode+1);
       end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractULASubspaceDOA(obj);
        if isLocked(obj)
            s.pPhasedModeTransform = obj.pPhasedModeTransform;
            s.pFeH = obj.pFeH;
            s.pKR = obj.pKR;
            s.pIsBeamspace = obj.pIsBeamspace;
            s.pConjugatedBFWeights = obj.pConjugatedBFWeights;
            s.pJ = obj.pJ;
            s.pLastElAng = obj.pLastElAng;
        end
    end
    
    function loadObjectImpl(obj,s,~)
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function doas = stepImpl(obj, Xarg,varargin)
    
        classtouse=class(Xarg);    
        fb = obj.ForwardBackwardAveraging;
        phasedModeTransform = obj.pPhasedModeTransform;
        if phasedModeTransform
            elAng = varargin{1};
            updateBeamformer(obj,elAng);
            % Transform input from element space to beamspace
            X=  Xarg*obj.pConjugatedBFWeights;
        else
            X = Xarg;
        end
        Cx = cast(step(obj.cCovEstimator,X),classtouse);
        [eigenvals, eigenvects] = privEig(obj,Cx);
        eigenvals = cast(eigenvals,'double');
        eigenvects = cast(eigenvects,'double');

        % Note when ForwardBackwardAveraging is true, the returned
        % eigenvects are the eigenvectors of Sq instead of the covariance
        % matrix. We have to further deduct eigenvects of the covariance
        % matrix.
        if fb
            eigenvects = local_convert(eigenvects);
        end

        % Obtain signal subspace dimension
        Nsig = getNumSignals(obj,eigenvals,obj.pNumSnapshots,fb);
        if Nsig==0
            if isempty(coder.target)
                warning(message('phased:phased:doa:ZeroSourceNumber'));
            end
            doas = zeros(1,0,classtouse);
            return;
        end
        
        % Core Root-MUSIC algorithm
        N = obj.pEffChannels;
        % Separate the noise eigenvectors
        noise_eigenvects = eigenvects(:,Nsig+1:end); %Qn
        if phasedModeTransform && obj.pIsBeamspace % when in beamspace
            % polynomial D = Sbeamspace'*Qn*Qn'*Sbeamspace
            %               Sula'*J'*Qn*Qn'*J*Sula 
             Qn = obj.pJ*noise_eigenvects;
        else %when in ULA or ULA like space
            % polynomial D = Sula'QnQn'Sula
             Qn = noise_eigenvects;
        end
        doasCM = cast(privCoreRootMUSIC(Qn,Nsig,N),classtouse);
        %doasCM = privCoreRootMUSIC(noise_eigenvects,Nsig,N);
        if phasedModeTransform
            % only need to convert from rad to degrees here since
            % Z = e^j*theta, angle(Z) = theta (root music polynomial in
            % terms of z) Eq 13 in [2]          
            doasRad = doasCM*180/pi;
            % convert to row vector
            doasCol = doasRad(:);
            doas = doasCol.';
        else
            % Convert angles when ULA
            doas = privConvertAngles(obj,doasCM);
        end
    end
    function privValidateSensorArray(~,val)
        %privValidateSensorArray
        %   Certain array operation is limited to certain array geometries.
        %   Each operation can then overload this method to do its own
        %   validation. By default, any array geometry is ok.
        %The array must be a ULA or UCA in phased mode excitation
        cond = ~(isa(val,'phased.ULA') ...
            || isa(val,'phased.HeterogeneousULA') ...
            || isa(val,'phased.UCA'));
        if cond
            coder.internal.errorIf(cond,'phased:phased:RootMUSICEstimator:InvalidArray');
        end
        
        if isa(val,'phased.UCA')
            elem = val.Element;
            cond = ~(isa(elem,'phased.IsotropicAntennaElement') || ...
                   isa(elem,'phased.OmnidirectionalMicrophoneElement') || ...
                   isa(elem,'phased.ShortDipoleAntennaElement') || ...
                   isa(elem,'phased.CustomAntennaElement') || ...
                   isa(elem,'phased.CustomMicrophoneElement'));                   
            if cond
                coder.internal.errorIf(cond,'phased:phased:RootMUSICEstimator:InvalidElement',class(elem));
            end                   
        end
        
    end
    function N = getNumEffectiveElements(obj)
        if isa(obj.SensorArray,'phased.UCA')
            maxPhaseMode = getMaxPhaseMode(obj);
            N = 2*maxPhaseMode+1;
        else
            N = getNumElements(obj.SensorArray);
        end
    end
     function num = getNumInputsImpl(obj)
        if isa(obj.SensorArray,'phased.UCA')
            num = 2;
        else
            num = 1;
        end
     end
end
   
methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl
    groups = getPropertyGroupsImpl@phased.internal.AbstractULASubspaceDOA('ulauca');
    props = {...
      'ForwardBackwardAveraging'...
      'SpatialSmoothing'};
    groups(1).PropertyList = [groups(1).PropertyList props];
    % SpatialSmoothing as dependent on private-only
    groups(1).DependOnPrivatePropertyList =  {'SpatialSmoothing'};
  end
  function header = getHeaderImpl
      header = matlab.system.display.Header(...
          'Title',getString(message('phased:library:block:RootMUSICEstimatorTitle')),...
          'Text',getString(message('phased:library:block:RootMUSICEstimatorDesc')));
  end    
end
methods (Access = protected) %for Simulink
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Root MUSIC\nDOA');        
    end    
    function varargout = getInputNamesImpl(obj)
         if isa(obj.SensorArray,'phased.UCA')
          varargout = {'','ElAng'};   
         else
          varargout = {''};
         end
    end
end

end

function doas = privCoreRootMUSIC(Qn,Nsig,N)
%privCoreRootMUSIC Core Root-MUSIC algorithm
%   DOA = privCoreRootMUSIC(Hdoa,EIGENVECTS,NSIG) separates the signal and
%   noise eigenvectors, forms a polynomial D consisting of a sum of
%   polynomials given by the product of the noise subspace eigenvectors and
%   the reversed and conjugated version and computes the NSIG roots of D
%   closest to the unit circle.


    % Form a polynomial D
    % D consists of a sum of polynomials given by the product of the noise
    % subspace eigenvectors and the reversed and conjugated version.
    % D = S'QnQn'S 
    D = complex(zeros(2*N-1,1));
    for i = 1:N-Nsig
        D = D + conv(Qn(:,i),conj(flipud(Qn(:,i))));
    end
 % Another way   
%      Q = noise_eigenvects; % numElx(numEl-numSig)
%      C = Q*Q';   % numEl x numEl
%      for idx = -(N-1):(N-1); DD(idx+N) = sum(diag(C,idx)); end
%      D = DD';


    % Take the angle of the NSIG roots of D closest to the unit circle.
    roots_D = roots(D);
    roots_D1 = roots_D(abs(roots_D) < 1);
    [~,indx] = sort(abs(abs(roots_D1)-1));
    sorted_roots = roots_D1(indx);
    doas = angle(sorted_roots(1:Nsig)); %get phase of the roots (z)

end

function eigenvects = local_convert(eigenvects)
% Deduct eigenvects of the covariance matrix from eigenvects of Sq

    M = size(eigenvects,1);

    half = floor(M/2);
    Uc1 = eigenvects(1:half,:);

    if rem(M,2)
        Uc2 = eigenvects(half+2:M,:);
        Umid = sqrt(2)*eigenvects(half+1,:);
    else
        Uc2 = eigenvects(half+1:M,:);
        Umid = [];
    end
    % Note: J*X2 with J = fliplr(eye(half)) is equivalent to flipud(X2)
    % eigenvects = 1/sqrt(2)*[Uc1+i*Uc2;Umid;J*(Uc1-i*Uc2)];   % Eq 7.68 in [1]
    eigenvects = 1/sqrt(2)*[Uc1+1i*Uc2;Umid;flipud(Uc1-1i*Uc2)];

end

function w = getPhasedModeW(N,phaseModes)

w = complex(zeros(N,numel(phaseModes)));
colIdx = 1;
elAzimuth = (-(N-1)/2:(N-1)/2)*2*pi/N; %in radians
for phaseMode  = phaseModes
    %Eq 4 in [2], note that the element positions on the circle are ordered
    %differently in our toolbox than in [2]
    w(:,colIdx) = (exp(1i*phaseMode*elAzimuth)/sqrt(N))'; 
    colIdx = colIdx+1;
end
end

