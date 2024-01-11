classdef (Sealed,StrictDefaults) MUSICEstimator < phased.internal.AbstractULASpectralDOA
%MUSICEstimator MUSIC spatial spectrum estimator for ULA
%   H = phased.MUSICEstimator creates a MUSIC spatial spectrum
%   estimator System object, H. The object estimates the incoming
%   signal's spatial spectrum using a narrowband MUSIC beamformer for a
%   uniform linear array (ULA).
%
%   H = phased.MUSICEstimator(Name,Value) creates a MUSIC spatial
%   spectrum estimator object, H, with the specified property Name set
%   to the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) estimates the spatial spectrum from X using the
%   estimator H. X is a matrix whose columns correspond to channels. Y
%   is a column vector representing the magnitude of the estimated
%   spatial spectrum.
%
%   [Y,ANG] = step(H,X) returns additional output ANG as the signal's
%   direction of arrival (DOA) when the DOAOutputPort property is true.
%   ANG is a row vector of the estimated broadside angles (in degrees).
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   MUSICEstimator methods:
%
%   step         - Perform spatial spectrum estimation (see above)
%   release      - Allow property value and input characteristics changes
%   clone        - Create a MUSIC spatial spectrum estimator object with
%                  the same property values
%   isLocked     - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>        - Reset states of MUSIC spatial spectrum estimator object
%   plotSpectrum - Plot spatial spectrum
%
%   MUSICEstimator properties:
%
%   SensorArray              - Sensor array
%   PropagationSpeed         - Signal propagation speed
%   OperatingFrequency       - Operating frequency
%   ForwardBackwardAveraging - Forward-backward averaging
%   ScanAngles               - Scan angles
%   DOAOutputPort            - Enable DOA output
%   NumSignalsSource         - Source of number of signals
%   NumSignalsMethod         - Method to estimate number of signals
%   NumSignals               - Number of signals
%   SpatialSmoothing         - Spatial smoothing for coherent sources 
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
%   %   Estimate the DOAs of two signals received by a standard 10-element
%   %   ULA with element spacing 1 meter. The antenna operating frequency
%   %   is 150 MHz. The actual direction of the first signal is 10 degrees
%   %   in azimuth and 20 degrees in elevation. The direction of the second
%   %   signal is 60 degrees in azimuth and -5 degrees in elevation.
%
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   array = phased.ULA('NumElements',10,'ElementSpacing',1);
%   fc = 150e6;
%   x = collectPlaneWave(array,[x1 x2],[10 20;60 -5]',fc);
%   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
%   estimator = phased.MUSICEstimator('SensorArray',array,...
%                  'OperatingFrequency',fc,...
%                  'DOAOutputPort',true,'NumSignalsSource','Property',...
%                  'NumSignals',2);
%   [y,doas] = estimator(x+noise);
%   doas = broadside2az(sort(doas),[20 -5]);
%   % Plot the spatial spectrum
%   plot(estimator.ScanAngles,y);
%   xlabel('Broadside Angles (degrees)'); ylabel('Magnitude Spectrum');
%
%   % Example 2:
%   %   Repeat the previous example and use plotSpectrum to plot the
%   %   spatial spectrum.
%
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   array = phased.ULA('NumElements',10,'ElementSpacing',1);
%   fc = 150e6;
%   x = collectPlaneWave(array,[x1 x2],[10 20;60 -5]',fc);
%   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
%   estimator = phased.MUSICEstimator('SensorArray',array,...
%                  'OperatingFrequency',fc);
%   y = estimator(x+noise);
%   plotSpectrum(estimator);
%
%   See also phased, phased.MUSICEstimator2D, musicdoa, broadside2az.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %NumSignalsSource Source of number of signals
    %   Specify the source of the number of signals as one of 'Auto' |
    %   'Property'. The default is 'Auto'. If you set this property to
    %   'Auto', the number of signals is estimated by the method
    %   specified by the NumSignalsMethod property.
    NumSignalsSource = 'Auto';
    %NumSignalsMethod Method to estimate number of signals
    %   Specify the method to estimate the number of signals as one of
    %   'AIC' | 'MDL'. The default is 'AIC'. The 'AIC' uses the Akaike
    %   Information Criterion and the 'MDL' uses Minimum Description
    %   Length Criterion. This property applies when you set the
    %   NumSignalsSource property to 'Auto'.
    NumSignalsMethod = 'AIC';
end

properties (Nontunable, PositiveInteger) 
    %NumSignals Number of signals
    %   Specify the number of signals as a positive scalar integer. The
    %   default value is 1. This property applies when you set the
    %   NumSignalsSource property to 'Property'.
    NumSignals = 1;
end

properties(Constant, Hidden)
    NumSignalsSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    NumSignalsMethodSet = matlab.system.StringSet({ 'AIC','MDL' });
end

methods

    function obj = MUSICEstimator(varargin)
        obj = obj@phased.internal.AbstractULASpectralDOA(varargin{:});
    end

end

methods (Access = protected)

    function [scanpattern,doasOut] = stepImpl(obj,X)
    % This method is modified from the superclass method to return the
    % number of signals (either from the property or computed) from
    % privDOASpectrum

        classtouse=class(X);
        % generate spatial magnitude spectrum
        numSignals = privDOASpectrum(obj,X); 
        scanpattern = sqrt(obj.pPattern);  

        if nargout > 1
            [~,locs] = findpeaks(scanpattern,'SortStr','descend');
            D = min(numSignals,length(locs));
            assert(D <= numSignals);            

            if D>0
                doasOut = NaN(1,numSignals,classtouse);%invalid angles
                doasOut(1:D) = obj.ScanAngles(locs(1:D).');
                % Warn if the number of peaks found is less than the number
                % of signals requested.
                if strcmp(obj.NumSignalsSource,'Property') && ....
                    D < obj.NumSignals && isempty(coder.target)
                    warning(message('phased:phased:internal:AbstractULASpectralDOA:LessNumSources', D));
                end
            else %if D==0
                if isempty(coder.target)
                    warning(message('phased:phased:doa:ZeroSourceNumber'));
                end
                doasOut = NaN(1,numSignals,classtouse);%invalid angles
            end
        end

    end

    function Nsig = privDOASpectrum(obj,X)

        fb = obj.ForwardBackwardAveraging;
        classtouse = class(X);

        % Estimate the covariance matrix Cx and compute eigenvalues and
        % eigenvectors of Cx.
        % Estimate the covariance matrix
        Cx = step(obj.cCovEstimator,X);

        % Compute eigenvectors of the covariance matrix
        [eigenvals, eigenvects] = privEig(obj,Cx);
        eigenvals = cast(eigenvals,classtouse);
        eigenvects = cast(eigenvects,classtouse);

        % Obtain signal subspace dimension
        Nsig = getNumSignals(obj,eigenvals,obj.pNumSnapshots,fb);

        if Nsig==0
            if isempty(coder.target)
                warning(message('phased:phased:doa:ZeroSourceNumber'));
            end
            return;
        end
        
        % Form MUSIC denominator matrix from noise subspace
        % eigenvectors
        noise_eigenvects = eigenvects(:, Nsig + 1:end); 
      
        % Compute spatial spectrum. Add a small positive constant to prevent
        % division by zero.
        sv = obj.pSteeringVector;
        D = sum(abs((sv'*noise_eigenvects)).^2,2)+eps(1); % 9.44 in [1]
        obj.pPattern = 1./D; 

    end

    function titlestr = getSpectrumPlotTitle(obj) % #ok<MANU>
        coder.extrinsic('phased.MUSICEstimator.getSpectrumPlotTitleExtrinsic');
        titlestr = obj.getSpectrumPlotTitleExtrinsic();
    end

    function D = getNumSignals(obj,eigenvals,K,fb)
        % Obtain signal subspace dimension
        if strcmp(obj.NumSignalsSource,'Auto')
            if strcmp(obj.NumSignalsMethod,'AIC')
                D = phased.internal.aictest(eigenvals,K,fb);
            else
                D = phased.internal.mdltest(eigenvals,K,fb);
            end
        else
            D = obj.NumSignals;
        end
    end

    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractULASpectralDOA(obj);
        if strcmp(obj.NumSignalsSource,'Property')
            % Verify that requested number of signals is not greater than
            % the number allowed based on the array size and smoothing.
            maxNumSig = getMaxNumSignal(obj);
            cond = obj.NumSignals > maxNumSig;
            if cond
                coder.internal.errorIf(cond,'phased:phased:internal:AbstractULASpectralDOA:InvalidNumSources',maxNumSig);
            end
        end
    end

    function maxNumSig = getMaxNumSignal(obj)
        % NumSignals should be no greater than possible signal subspace dim - 1
        maxNumSig = getEffectiveChannel(obj) - 1;
    end

    function [eigenvals, eigenvects] = privEig(obj,Sx) %#ok<INUSL> % #ok<MANU>
        [eigenvects, eigenvalsC] = eig(Sx);
        eigenvals = real(eigenvalsC);
        [eigenvals,indx] = sort(diag(eigenvals),'descend');
        eigenvects= eigenvects(:,indx);
        eigenvals(eigenvals<0) = 0;
    end

    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        SourceAuto = strcmp(obj.NumSignalsSource,'Auto');
        if SourceAuto && strcmp(prop,'NumSignals')
            flag = true;
        end
        if ~SourceAuto && strcmp(prop,'NumSignalsMethod')
            flag = true;
        end
    end
end

methods (Static, Hidden)
    function titlestr = getSpectrumPlotTitleExtrinsic()
        titlestr = sprintf('%s %s',...
            getString(message('phased:doa:MUSIC')),...
            getString(message('phased:doa:SpatialSpectrum')));
    end
end

methods (Access = protected) %for Simulink
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('ULA MUSIC\nSpectrum');
    end  
end

methods (Static,Hidden,Access=protected)
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:MUSICEstimatorTitle')),...
            'Text',getString(message('phased:library:block:MUSICEstimatorDesc')));
    end

    function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractULASpectralDOA;
        pNumSignalsSource = matlab.system.display.internal.Property('NumSignalsSource', ...
            'IsGraphical', false, ...
            'UseClassDefault', false,'Default','Property');
        props = {...
            pNumSignalsSource ,...
            'NumSignalsMethod'};
        groups(1).PropertyList = [groups(1).PropertyList props];
    end
end

end


