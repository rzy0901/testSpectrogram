classdef (Sealed,StrictDefaults) MVDREstimator2D < phased.internal.AbstractSpectralDOA
%MVDREstimator2D 2-D MVDR spatial spectrum estimator
%   H = phased.MVDREstimator2D creates a 2-D MVDR spatial spectrum
%   estimator System object, H. The object estimates the signal's spatial
%   spectrum using a narrowband MVDR beamformer.
%
%   H = phased.MVDREstimator2D(Name,Value) creates a 2-D MVDR spatial
%   spectrum estimator object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) estimates the spatial spectrum from X using the estimator
%   H. X is a matrix whose columns correspond to channels. Y is a matrix
%   representing the magnitude of the estimated 2-D spatial spectrum. The
%   dimension of Y is [NumElevation NumAzimuth].
%
%   [Y,ANG] = step(H,X) returns additional output ANG as the signal's
%   direction of arrival (DOA) when the DOAOutputPort property is true. ANG
%   is a 2-row matrix where the first row represents estimated azimuth and
%   the 2nd row represents estimated elevation (in degrees).
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   MVDREstimator2D methods:
%
%   step         - Perform spatial spectrum estimation (see above)
%   release      - Allow property value and input characteristics changes
%   clone        - Create a 2-D MVDR spatial spectrum estimator object with
%                  same property values
%   isLocked     - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>        - Reset states of 2D MVDR spatial spectrum estimator 
%                  object
%   plotSpectrum - Plot spatial spectrum
%
%   MVDREstimator2D properties:
%
%   SensorArray              - Sensor array
%   PropagationSpeed         - Signal propagation speed
%   OperatingFrequency       - Operating frequency
%   NumPhaseShifterBits      - Number of bits in phase shifters
%   ForwardBackwardAveraging - Forward-backward averaging
%   AzimuthScanAngles        - Azimuth scan angles
%   ElevationScanAngles      - Elevation scan angles
%   DOAOutputPort            - Enable DOA output
%   NumSignals               - Number of signals
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
%   %   Estimate the DOAs of two signals received by a 50-element URA. The 
%   %   antenna operating frequency is 150 MHz. The actual direction of the
%   %   first signal is -37 degrees in azimuth and 0 degrees in elevation. 
%   %   The direction of the second signal is 17 degrees in azimuth and 20 
%   %   degrees in elevation.
%
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   array = phased.URA('Size',[5 10],'ElementSpacing',[1 0.6]);
%   array.Element.FrequencyRange = [100e6 300e6];
%   fc = 150e6;
%   x = collectPlaneWave(array,[x1 x2],[-37 0;17 20]',fc);
%   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
%   estimator = phased.MVDREstimator2D('SensorArray',array,...
%                  'OperatingFrequency',fc,...
%                  'DOAOutputPort',true,'NumSignals',2,...
%                  'AzimuthScanAngles',-50:50,...
%                  'ElevationScanAngles',-30:30);
%   [y,doas] = estimator(x+noise);
%   % Plot the spatial spectrum
%   surf(estimator.AzimuthScanAngles,estimator.ElevationScanAngles,y);
%   xlabel('Azimuth Angles (degrees)'); 
%   ylabel('Elevation Angles (degrees)');
%   zlabel('Magnitude Spectrum');
%
%   % Example 2:
%   %   Repeat the previous example and use plotSpectrum to plot the 
%   %   spatial spectrum.
%
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   array = phased.URA('Size',[5 10],'ElementSpacing',[1 0.6]);
%   array.Element.FrequencyRange = [100e6 300e6];
%   fc = 150e6;
%   x = collectPlaneWave(array,[x1 x2],[-37 0;17 20]',fc);
%   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
%   estimator = phased.MVDREstimator2D('SensorArray',array,...
%                  'OperatingFrequency',fc,...
%                  'AzimuthScanAngles',-50:50,...
%                  'ElevationScanAngles',-30:30);
%   y = estimator(x+noise);
%   plotSpectrum(estimator);
%
%   See also phased, phased.MVDREstimator.

%   Copyright 2009-2018 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable, Dependent = true)
    %NumPhaseShifterBits    Number of bits in phase shifters
    %   Specify the number of bits used in the phase shifter as a
    %   non-negative integer. The default value of this property is 0,
    %   indicating there is no quantization effect in the phase shifter.
    NumPhaseShifterBits = 0
end

properties (Nontunable, PositiveInteger) 
    %NumSignals Number of signals
    %   Specify the number of signals for DOA estimation as a positive
    %   scalar integer. The default value is 1. This property applies when
    %   you set the DOAOutputPort property to true.
    NumSignals = 1;
end

methods

    function obj = MVDREstimator2D(varargin)
        obj = obj@phased.internal.AbstractSpectralDOA(varargin{:});
    end
    
end

methods
    function set.NumPhaseShifterBits(obj,value)
        validateattributes( value, { 'double','single' }, ...
          { 'scalar', 'integer', 'nonnegative', 'finite' }, '', 'NumPhaseShifterBits');
        obj.pNumPhaseShifterBits = value;
    end
    function val = get.NumPhaseShifterBits(obj)
        val = obj.pNumPhaseShifterBits;
    end
end

methods (Access = protected)
        
    function privDOASpectrum(obj,Cx)

        if obj.pOneIterFlag
            sv = obj.pSteeringVectors;
            obj.pPattern = 1./real(sum(sv'.*(Cx\sv).',2));    % Eq. 9.13 in [1]
        else
            scanAngleIndex = obj.pScanAngleBlockIndex;
            hsv = obj.cSteeringVector;
            fc = obj.OperatingFrequency;
            for m = 0:obj.pNumIter-1
                tempidx = scanAngleIndex+(obj.pScanAngleBlockSize*m);
                tempangle = obj.pScanAngles(:,tempidx);
                sv = step(hsv,fc,tempangle);
                obj.pPattern(tempidx) = 1./real(sum(sv'.*(Cx\sv).',2));  % Eq. 9.13 in [1]
            end
        end
        
    end
    function titlestr = getSpectrumPlotTitle(obj) 
        coder.extrinsic('phased.MVDREstimator2D.getSpectrumPlotTitleExtrinsic');
        titlestr = obj.getSpectrumPlotTitleExtrinsic();
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('MVDR\nSpectrum');
    end    
end
methods (Static, Hidden) 
    function titlestr = getSpectrumPlotTitleExtrinsic()
        titlestr = sprintf('%s %s',...
            getString(message('phased:doa:MVDR2D')),...
            getString(message('phased:doa:SpatialSpectrum')));
    end
end
methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractSpectralDOA;
        props = {...
          'NumPhaseShifterBits'};
        iNew = find(strcmp(groups(1).PropertyList,'ForwardBackwardAveraging'));
        groups(1).PropertyList = [groups(1).PropertyList(1:iNew-1) props groups(1).PropertyList(iNew:end)];
        groups(1).DependOnPrivatePropertyList =  [groups(1).DependOnPrivatePropertyList {'NumPhaseShifterBits'}];
    end
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:MVDREstimator2DTitle')),...
            'Text',getString(message('phased:library:block:MVDREstimator2DDesc')));
    end    
end
end
