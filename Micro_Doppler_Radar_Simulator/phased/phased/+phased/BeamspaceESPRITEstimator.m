classdef (Sealed,StrictDefaults) BeamspaceESPRITEstimator < phased.internal.AbstractESPRIT
%BeamspaceESPRITEstimator Beamspace ESPRIT direction of arrival (DOA)
%estimator
%   H = phased.BeamspaceESPRITEstimator creates a beamspace ESPRIT DOA
%   estimator System object, H. The object estimates the signal's direction
%   of arrival using the beamspace ESPRIT algorithm with a uniform linear
%   array (ULA).
%
%   H = phased.BeamspaceESPRITEstimator(Name,Value) creates a beamspace 
%   ESPRIT DOA estimator object, H, with the specified property Name set
%   to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
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
%   BeamspaceESPRITEstimator methods:
%
%   step     - Perform DOA estimation (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a beamspace ESPRIT DOA estimator object with same
%              property values
%   isLocked - Locked status (logical)
%
%   BeamspaceESPRITEstimator properties:
%
%   SensorArray        - Sensor array
%   PropagationSpeed   - Signal propagation speed
%   OperatingFrequency - Operating frequency
%   NumSignalsSource   - Source of number of signals
%   NumSignalsMethod   - Method to estimate number of signals
%   NumSignals         - Number of signals
%   SpatialSmoothing   - Spatial smoothing
%   Method             - Type of least square method
%   BeamFanCenter      - Beam fan center direction
%   NumBeamsSource     - Source of number of beams
%   NumBeams           - Number of beams
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision, regardless of the precision of the
%   properties and arguments. If the input data X is double precision,
%   regardless of the precision of the properties and arguments.
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
%   ha = phased.ULA('NumElements',10,'ElementSpacing',1);
%   ha.Element.FrequencyRange = [100e6 300e6];
%   fc = 150e6;
%   x = collectPlaneWave(ha,[x1 x2],[10 20;45 60]',fc);
%   rs = RandStream('mt19937ar','Seed',0);
%   noise = 0.1/sqrt(2)*(randn(rs,size(x))+1i*randn(rs,size(x)));
%   beamspace = phased.BeamspaceESPRITEstimator('SensorArray',ha,...
%                 'OperatingFrequency',fc,...
%                 'NumSignalsSource','Property','NumSignals',2);
%   doas = beamspace(x+noise);
%   az = broadside2az(sort(doas),[20 60])
%
%   See also phased, phased.ESPRITEstimator, broadside2az.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties
    %BeamFanCenter Beam fan center direction (deg)
    %   Specify the direction of the center of the beam fan (in degrees) as
    %   a real scalar value between -90 and 90. The default value is 0.
    %   This property is tunable.
    BeamFanCenter = 0;
end

properties (Nontunable)
    %NumBeamsSource Source of number of beams
    %   Specify the source of the number of beams as one of 'Auto' |
    %   'Property'. The default is 'Auto'. If you set this property to
    %   'Auto', the number of beams equals to N-L, where N is the number of
    %   array elements and L is the value of the SpatialSmoothing property.
    NumBeamsSource = 'Auto';
end

properties (Nontunable, PositiveInteger) 
    %NumBeams Number of beams
    %   Specify the number of beams as a positive scalar integer. The
    %   default value is 2. The lower the number of beams, the greater the
    %   reduction in computational cost. This property applies when you set
    %   the NumBeamsSource to 'Property'. 
    NumBeams = 2;
end

properties(Constant, Hidden)
    NumBeamsSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
end

methods

    function set.BeamFanCenter(obj,var)
        validateattributes( var, { 'double','single' }, { 'scalar', 'real', '<=', 90, '>=', -90 }, '', 'BeamFanCenter');
        obj.BeamFanCenter = var;
    end

    function obj = BeamspaceESPRITEstimator(varargin)
        obj = obj@phased.internal.AbstractESPRIT(varargin{:});
    end
    
end

methods (Access = protected)
    
    
    function effChan = getEffectiveChannel(obj)
        NChannels = getNumElements(obj.SensorArray) - obj.pSpatialSmoothing;
        if strcmp(obj.NumBeamsSource,'Auto')
            effChan = NChannels;
        else
            effChan = obj.NumBeams;
            cond =  effChan > NChannels;
            if cond
                coder.internal.errorIf(cond,'phased:phased:BeamspaceESPRITEstimator:TooManyBeams', NChannels);
            end
        end
    end
    
    function maxNumSig = getMaxNumSignal(obj)
    % NumSignals should be no greater than possible signal subspace dim
        maxNumSig = getEffectiveChannel(obj);
    end

    
    function setupImpl(obj, X)
        setupImpl@phased.internal.AbstractESPRIT(obj, X);
        % always use forward-backward averaging to compute Cx.
        obj.cCovEstimator.ForwardBackwardAveraging = true;
        % don't compute the unitary transformation of Cx
        obj.cCovEstimator.UnitaryTransformed = false;
    end
    
    function loadObjectImpl(obj,s,wasLocked)
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function doasOut = stepImpl(obj, X)

        % Estimate the forward-backward averaged spatial covariance matrix
        % in element space.
        classtouse=class(X);
        Sx = step(obj.cCovEstimator,X);
        Nbs = obj.pEffChannels;
        
        M = size(Sx,1); % effective array element number
            
        uc = round(sin(obj.BeamFanCenter*pi/180)*M/2);
                
        % Compute the real beamspace covariance matrix
        [W, beamfan] = privDFTMatrix(M,Nbs,uc);
        Sxbs = real(W'*Sx*W);
        [eigenvals, eigenvects] = privEig(obj,Sxbs);
        eigenvals = cast(eigenvals,classtouse);% g1812952
        eigenvects = cast(eigenvects,classtouse);

        % Obtain signal subspace dimension
        D = getNumSignals(obj,eigenvals,obj.pNumSnapshots,true);
        if D==0
            if isempty(coder.target)
                warning(message('phased:phased:doa:ZeroSourceNumber'));
            end
            doasOut = zeros(1,0,class(X));            
            return;
        end
        
        % Selection matrices
        [Gamma1, Gamma2] = local_selection_matrices(M,cast(Nbs,classtouse),cast(beamfan,classtouse));
        
        % Core ESPRIT algorithm
        psieig = privCoreESPRIT(obj,eigenvects,D,Gamma1,Gamma2,obj.Method);
        psieigReal = real(psieig);
        % Convert angles
        doasOut = privConvertAngles(obj,2*atan(psieigReal));
        
    end
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractESPRIT(obj, prop);
        if strcmp(obj.NumBeamsSource,'Auto') && strcmp(prop,'NumBeams')
            flag = true;
        end
    end

end

methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl
    groups = getPropertyGroupsImpl@phased.internal.AbstractESPRIT;
    props = {...
      'BeamFanCenter',...
      'NumBeamsSource',...
      'NumBeams'};
    groups(1).PropertyList = [groups(1).PropertyList props];
  end
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:BeamspaceESPRITEstimatorTitle')),...
            'Text',getString(message('phased:library:block:BeamspaceESPRITEstimatorDesc')));
    end    
end
methods (Access = protected) %for Simulink
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Beamspace\nESPRIT DOA');        
    end    
end
end

function [W, beamfan] = privDFTMatrix(M,Nbs,uc)
%privDFTMatrix DFT (a.k.a. Butler) Matrix
%   M=N-L+1 number of effective sensors
%   Nbs number of beams
%   uc center beam (for reduced dimension beamspace i.e. Nbs<M)

    % Define Nbs beam indices
    m = floor(-Nbs/2)+1:floor(Nbs/2); % Eq 3.322 & 3.323 in [1]
    beamfan = m+uc; % uc*2/M increments

    % Wrap if beamfan>floor(M/2)
    idx = beamfan>floor(M/2);
    if any(idx),
        beamfan(idx) = beamfan(idx)-M;
    end

    % Wrap if beamfan<floor(-M/2)+1
    idx = beamfan<floor(-M/2)+1;
    if any(idx),
        beamfan(idx) = beamfan(idx)+M;
    end

    n = (-(M-1)/2:(M-1)/2).';
    W = 1/M*exp(1j*2*pi/M*n*beamfan);

end

function [Gamma1, Gamma2] = local_selection_matrices(M,Nbs,beamfan)
% Generate selection matrix

    aux = cos(pi/M*beamfan);
    lastCol1 = [zeros(Nbs-2,1);aux(Nbs)];
    Gamma1temp = diag(aux(1:Nbs-1))+diag(aux(2:Nbs-1),1); % Eq 9.313
    Gamma1 = [Gamma1temp lastCol1];

    aux = sin(pi/M*beamfan);
    lastCol2 = [zeros(Nbs-2,1);aux(Nbs)];
    Gamma2temp = diag(aux(1:Nbs-1))+diag(aux(2:Nbs-1),1); % Eq 9.314
    Gamma2 = [Gamma2temp lastCol2];

end

