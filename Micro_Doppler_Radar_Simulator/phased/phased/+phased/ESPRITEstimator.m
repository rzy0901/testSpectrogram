classdef (Sealed,StrictDefaults) ESPRITEstimator < phased.internal.AbstractESPRIT
%ESPRITEstimator ESPRIT direction of arrival (DOA) estimator
%   H = phased.ESPRITEstimator creates an ESPRIT DOA estimator System
%   object, H. The object estimates the signal's direction-of-arrival (DOA)
%   using the ESPRIT algorithm with a uniform linear array (ULA).
%
%   H = phased.ESPRITEstimator(Name,Value) creates an ESPRIT DOA estimator 
%   object, H, with the specified property Name set to the specified
%   Value. You can specify additional name-value pair arguments in any
%   order as (Name1,Value1,...,NameN,ValueN).
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
%   ESPRITEstimator methods:
%
%   step     - Perform DOA estimation (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create an ESPRIT DOA estimator object with same property 
%              values
%   isLocked - Locked status (logical)
%
%   ESPRITEstimator properties:
%
%   SensorArray              - Sensor array
%   PropagationSpeed         - Signal propagation speed
%   OperatingFrequency       - Operating frequency
%   NumSignalsSource         - Source of number of signals
%   NumSignalsMethod         - Method to estimate number of signals
%   NumSignals               - Number of signals
%   SpatialSmoothing         - Spatial smoothing
%   Method                   - Type of least square method
%   ForwardBackwardAveraging - Forward-backward averaging
%   RowWeighting             - Row weighting factor
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
%   estimator = phased.ESPRITEstimator('SensorArray',array,...
%               'OperatingFrequency',fc);
%   doas = estimator(x+noise);
%   az = broadside2az(sort(doas),[20 60])
%
%   See also phased, phased.BeamspaceESPRITEstimator, broadside2az.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002
%   [2] Martin Haardt and Josef A. Nossek. Unitary ESPRIT: How to Obtain
%   increased Estimation Accuracy with a Reduced Computational Burden, lEEE
%   TRANSACTIONS ON SIGNAL PROCESSING. VOL. 43, NO. 5, MAY 1995

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable, PositiveInteger) 
    %RowWeighting Row weighting factor
    %   Specify the row weighting factor for signal subspace eigenvectors
    %   as a positive integer scalar. The default is 1. This property
    %   controls the weights applied to the selection matrices. In most
    %   cases the higher value the better. However it can never be greater
    %   than (N-SS)/2 where N is the number of elements of the array and
    %   SS is the spatial smoothing factor.
    RowWeighting = 1;
end

properties (Nontunable,Hidden)
    %VisibleRegion Visible region
    %   Specify the DOA search limits (in degrees) as a real 2-element row
    %   vector. The default is [-90 90]. The vector must be symmetric
    %   around broadside (0 degree). This property applies when you set the
    %   NumSignalsSource property to 'Property'.
    VisibleRegion = [-90 90];
end

properties (Nontunable, Logical) 
    %ForwardBackwardAveraging Forward-backward averaging
    %   Set this property to true to use forward-backward averaging to
    %   estimate the covariance matrix for sensor arrays with conjugate
    %   symmetric array manifold. The default value is false.
    ForwardBackwardAveraging = false;
end
    
properties(Access = private, Nontunable)
    % pDS is the number of element spacing between two subarrays.
    pDS;
end

methods

    function set.VisibleRegion(obj,var)
        if ~isempty(coder.target)
            coder.internal.errorIf(true,'phased:phased:ESPRITEstimator:ObsoleteVisibleRegion');
        end
        warning(message('phased:phased:ESPRITEstimator:ObsoleteVisibleRegion'));
        validateattributes( var, { 'double','single' }, { 'size', [ 1, 2 ], 'real', 'finite' }, '', 'VisibleRegion');
        % check symmetric around 0 and ascending
        if abs(var(1))~=abs(var(2)) || var(1)>=var(2)
            error(message('phased:phased:ESPRITEstimator:AngleNotSymmetric'));
        end
        % check if the values are over 90
        if var(1)<-90 || var(2)>90
            error(message('phased:phased:ESPRITEstimator:AngleInvalid'));
        end
        obj.VisibleRegion = var;
    end

    function obj = ESPRITEstimator(varargin)
        obj = obj@phased.internal.AbstractESPRIT(varargin{:});
    end
end

methods (Access = protected)
    
    function effChan = getEffectiveChannel(obj)
    %Size of each ESPRIT subarray when DS = 1
        DS = 1;
        M = getNumElements(obj.SensorArray) - DS;
        effChan = M - obj.pSpatialSmoothing;   
    end
    
    function maxNumSig = getMaxNumSignal(obj)
    % NumSignals should be no greater than possible signal subspace dim
        maxNumSig = getEffectiveChannel(obj);
    end

    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractESPRIT(obj);
        % verify RowWeighting
        effChannels = getEffectiveChannel(obj);
        cond = (obj.RowWeighting > (effChannels+1)/2);
        if cond
            coder.internal.errorIf(cond,'phased:phased:ESPRITEstimator:InvalidRowWeighting', fix((effChannels+1)/2));
        end
        % Calculate the maximum subarray spacing for given visible region.
        wavelength = obj.PropagationSpeed/obj.OperatingFrequency;
        d_lambda = obj.SensorArray.ElementSpacing/wavelength;
        v = obj.VisibleRegion(2);
        if v == 90 || strcmp(obj.NumSignalsSource,'Auto')
            obj.pDS = 1;
        else
            obj.pDS = floor(1/(2*d_lambda*sin(v*pi/180)));
            % Maximum subarray spacing is also limited by the number of 
            % array elements.
            if obj.pDS > floor(effChannels/2) && isempty(coder.target)
                warning(message('phased:phased:ESPRITEstimator:VisibleRegionTooSmall'));
                obj.pDS = floor(effChannels/2);
            end
        end
    end
           
    function setupImpl(obj, X) 
       setupImpl@phased.internal.AbstractESPRIT(obj, X);
       obj.cCovEstimator.ForwardBackwardAveraging = ...
           obj.ForwardBackwardAveraging;
       % When ForwardBackwardAveraging is true, return Sq instead of the 
       % covariance matrix.
       if obj.cCovEstimator.ForwardBackwardAveraging
           obj.cCovEstimator.UnitaryTransformed = true;
       end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractESPRIT(obj);
        %Being obsoleted, do not save.
        if isfield(s,'VisibleRegion') && isequal(s.VisibleRegion,[-90 90])
            s = rmfield(s,'VisibleRegion');
        end
        if isLocked(obj)
            s.pDS = obj.pDS;
        end
    end
    
    function loadObjectImpl(obj,s,~)
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            % Do not load visible region if default
            % since it will warn
            if ~(strcmp(fn{m},'VisibleRegion') && ...
                    isequal(s.VisibleRegion,[-90 90]))
                obj.(fn{m}) = s.(fn{m});
            end
        end
    end

    function doasOut = stepImpl(obj, X)
    % Element-space ESPRIT 
    
        classtouse = class(X);
        ds = obj.pDS;
        fb = obj.ForwardBackwardAveraging;

        Cx = step(obj.cCovEstimator,X);
        [eigenvals, eigenvects] = privEig(obj,Cx);
        eigenvals = cast(eigenvals,classtouse);%g1812952
        eigenvects = cast(eigenvects,classtouse);

        % Signal subspace dimension 
        D = getNumSignals(obj,eigenvals,obj.pNumSnapshots,fb);
        if D==0
            if isempty(coder.target)
                warning(message('phased:phased:doa:ZeroSourceNumber'));
            end
            doasOut = zeros(1,0,classtouse);
            return;
        end
                            
        % check eigenvalues against source dimension. ESPRIT does not work
        % when there are less then D non-zero eigenvalues.
        D_act = sum(eigenvals>eps(max(abs(eigenvals))));
        cond = D_act < D;
        if cond
            coder.internal.errorIf(cond,'phased:phased:ESPRITEstimator:NotEnoughRank', D, D_act);
        end
        
        % Selection matrices
        [Js1, Js2] = local_selection_matrices(cast(obj.pEffChannels,classtouse)+1,...
            cast(obj.RowWeighting,classtouse),ds,fb);

        % Core ESPRIT algorithm
        psieig = privCoreESPRIT(obj,eigenvects,D,Js1,Js2,obj.Method);

        % Extract angle
        doas = local_extract_angles(psieig,ds,fb);

        % Convert angles
        doasOut = privConvertAngles(obj,doas);

    end
           
end

methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractESPRIT;
        props = {...
            'ForwardBackwardAveraging',...
            'RowWeighting'};
        groups(1).PropertyList = [groups(1).PropertyList props];
    end
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:ESPRITEstimatorTitle')),...
            'Text',getString(message('phased:library:block:ESPRITEstimatorDesc')));
    end    
end
methods (Access = protected) %for Simulink
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('ESPRIT\nDOA');        
    end    
end

end

function [Js1, Js2] = local_selection_matrices(M,ms,ds,fb)
%   Generate selection matrices Js1 and Js2. These two selection matrices
%   are used to create subarray information for ESPRIT algorithms.

    % Row weighting
    Ns = M-ds;
    w = min(ceil(ms),Ns-ceil(ms)+1);                             % Eq 9.133 in [1]
    weights = diag(sqrt([1:w-1 w*ones(1,Ns-2*(w-1)) w-1:-1:1])); % Eq 9.132 in [1]
    O = zeros(Ns,ds);

    if fb
        % Constants
        Qn = phased.internal.unitarymat(M);
        Qns = phased.internal.unitarymat(Ns);
        % Selection matrices
        Js2 = [O weights];
        K = Qns'*Js2*Qn;
        Js1 = real(K);     % Eq 9.147 in [1]
        Js2 = imag(K);     % Eq 9.148 in [1]
    else
        % Selection Matrices
        Js1 = [weights O]; % Eq 9.134 in [1]
        Js2 = [O weights]; % Eq 9.135 in [1]    
    end

end

function doas = local_extract_angles(psieig,ds,fb)
%   Extract angle information estimated from two subarrays based on the
%   distance of the phase center between subarrays.

    if fb
        if isreal(psieig)
            psieigReal = psieig;
        else
            %reliability test for unitary ESPRIT
            %failed, return real part only
            psieigReal = real(psieig);        
        end
        doas = 2*atan(1/ds*psieigReal);
    else
        doas = 1/ds*angle(psieig);
    end

end


