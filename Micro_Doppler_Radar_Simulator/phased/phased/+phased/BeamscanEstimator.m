classdef (Sealed,StrictDefaults) BeamscanEstimator < phased.internal.AbstractULASpectralDOA
%BeamscanEstimator Beamscan spatial spectrum estimator for ULA
%   H = phased.BeamscanEstimator creates a beamscan spatial spectrum
%   estimator System object, H. The object estimates the incoming signal's
%   spatial spectrum using a narrowband conventional beamformer for a
%   uniform linear array (ULA).
%
%   H = phased.BeamscanEstimator(Name,Value) creates a beamscan spatial
%   spectrum estimator object, H, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax:
%
%   Y = step(H,X) estimates the spatial spectrum from X using the estimator
%   H. X is a matrix whose columns correspond to channels. Y is a column
%   vector representing the magnitude of the estimated spatial spectrum.
%
%   [Y,ANG] = step(H,X) returns additional output ANG as the signal's
%   direction of arrival (DOA) when the DOAOutputPort property is true. ANG
%   is a row vector of the estimated broadside angles (in degrees).
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   BeamscanEstimator methods:
%
%   step         - Perform spatial spectrum estimation (see above)
%   release      - Allow property value and input characteristics changes
%   clone        - Create a beamscan spatial spectrum estimator object with
%                  same property values
%   isLocked     - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>        - Reset states of beamscan spatial spectrum estimator 
%                  object
%   plotSpectrum - Plot spatial spectrum
%
%   BeamscanEstimator properties:
%
%   SensorArray              - Sensor array
%   PropagationSpeed         - Signal propagation speed
%   OperatingFrequency       - Operating frequency
%   NumPhaseShifterBits      - Number of bits in phase shifters
%   ForwardBackwardAveraging - Forward-backward averaging
%   SpatialSmoothing         - Spatial smoothing
%   ScanAngles               - Scan angles
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
%   %   Estimate the DOAs of two signals received by a standard 10-element 
%   %   ULA with element spacing 1 meter. The antenna operating frequency 
%   %   is 150 MHz. The actual direction of the first signal is 10 degrees
%   %   in azimuth and 20 degrees in elevation. The direction of the second
%   %   signal is 60 degrees in azimuth and -5 degrees in elevation.
%
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   ha = phased.ULA('NumElements',10,'ElementSpacing',1);
%   ha.Element.FrequencyRange = [100e6 300e6];
%   fc = 150e6;
%   x = collectPlaneWave(ha,[x1 x2],[10 20;60 -5]',fc);
%   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
%   beamscan = phased.BeamscanEstimator('SensorArray',ha,...
%                  'OperatingFrequency',fc,...
%                  'DOAOutputPort',true,'NumSignals',2);
%   [y,doas] = beamscan(x+noise);
%   doas = broadside2az(sort(doas),[20 -5])
%   % Plot the spatial spectrum
%   plot(beamscan.ScanAngles,y);
%   xlabel('Broadside Angles (degrees)'); ylabel('Magnitude Spectrum');
%
%   % Example 2:
%   %   Repeat the previous example and use plotSpectrum to plot the 
%   %   spatial spectrum.
%
%   fs = 8000; t = (0:1/fs:1).';
%   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
%   ha = phased.ULA('NumElements',10,'ElementSpacing',1);
%   ha.Element.FrequencyRange = [100e6 300e6];
%   fc = 150e6;
%   x = collectPlaneWave(ha,[x1 x2],[10 20;60 -5]',fc);
%   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
%   beamscan = phased.BeamscanEstimator('SensorArray',ha,...
%                  'OperatingFrequency',fc);
%   y = beamscan(x+noise);
%   plotSpectrum(beamscan)
%
%
%   See also phased, phased.BeamscanEstimator2D, broadside2az.

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

    function obj = BeamscanEstimator(varargin)
        obj = obj@phased.internal.AbstractULASpectralDOA(varargin{:});
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
        sv = obj.pSteeringVector;
        obj.pPattern = real(sum((sv'*Cx).*sv.',2));    % Eq. 9.4 in [1]
                
    end
    
    function titlestr = getSpectrumPlotTitle(obj)
        coder.extrinsic('phased.BeamscanEstimator.getSpectrumPlotTitleExtrinsic');
        titlestr = obj.getSpectrumPlotTitleExtrinsic();
    end
    
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('ULA Beamscan\nSpectrum ');
    end    
end
methods (Static, Hidden) 
    function titlestr = getSpectrumPlotTitleExtrinsic()
        titlestr = sprintf('%s %s',...
            getString(message('phased:doa:Beamscan')),...
            getString(message('phased:doa:SpatialSpectrum')));
    end
end
methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractULASpectralDOA;
        props = {...
          'NumPhaseShifterBits'};
        iNew = find(strcmp(groups(1).PropertyList,'ForwardBackwardAveraging'));
        groups(1).PropertyList = [groups(1).PropertyList(1:iNew-1) props groups(1).PropertyList(iNew:end)];
        groups(1).DependOnPrivatePropertyList =  [groups(1).DependOnPrivatePropertyList {'NumPhaseShifterBits'}];

    end
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:BeamscanEstimatorTitle')),...
            'Text',getString(message('phased:library:block:BeamscanEstimatorDesc')));
    end    
end
end
