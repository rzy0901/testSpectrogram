classdef (Hidden) AbstractULASpectralDOA < phased.internal.AbstractULADOA  & ...
     matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%ABSTRACTULASPECTRALDOA Abstract class for estimating DOA from peaks of
%spatial spectrum using a ULA

%   Copyright 2009-2018 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Dependent, Nontunable)
    %SpatialSmoothing Spatial smoothing
    %   Specify the effective reduction in the size of the sensor array due
    %   to spatial smoothing as a nonnegative integer. If the array
    %   consists of M elements and the value of this property is L,
    %   maximally overlapped subarrays of M-L elements are formed.
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

properties (Nontunable)
    %ScanAngles Scan angles (deg)
    %   Specify the scan angles (in degrees) as a real row vector. The
    %   default value is -90:90. The angles are broadside angles and must
    %   be within [-90 90]. You must specify the angles in an ascending
    %   order.
    ScanAngles = -90:90;
end

properties (Abstract = true) 
    NumSignals 
end

properties (Nontunable, Logical) 
    %ForwardBackwardAveraging Forward-backward averaging
    %   Set this property to true to use forward-backward averaging to
    %   estimate the covariance matrix for sensor arrays with conjugate
    %   symmetric array manifold. The default value is false.
    ForwardBackwardAveraging = false;
    %DOAOutputPort Enable DOA output
    %   Set this property to true to output the signal's direction of
    %   arrival (DOA). The default value is false.
    DOAOutputPort = false;
end


properties (Access = protected, Nontunable)
    % the number of phase shifter bits
    pNumPhaseShifterBits=0;
end

properties (Access = protected)
    % the pattern of the spectrum
    pPattern;
    % precalculated steering vector matrix
    pSteeringVector;
end

methods (Abstract, Access = protected)
    pattern = privDOASpectrum(obj,Cx);
    titlestr = getSpectrumPlotTitle(obj)
end

methods
    
    function set.SpatialSmoothing(obj,value)
        validateattributes( value, { 'double','single' }, { 'scalar', 'integer', 'nonnegative', 'finite' }, '', 'SpatialSmoothing');
        obj.pSpatialSmoothing = value;
    end
    
    function value = get.SpatialSmoothing(obj)
        value = obj.pSpatialSmoothing;
    end
    
    function set.ScanAngles(obj,value)
        validateattributes( value, { 'double','single' }, { 'row', 'real', 'finite', '>=', -90, '<=', 90 }, '', 'ScanAngles');
        % findpeaks requires at least 3 samples.
        cond = length(value) < 3;
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractULASpectralDOA:NotEnoughSamples');
        end        
        % ensure ascending
        cond = any(diff(value)<=0);
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractULASpectralDOA:NotAscending');
        end
        obj.ScanAngles = value;
    end

end

methods (Access = protected)

    function obj = AbstractULASpectralDOA(varargin)
        obj = obj@phased.internal.AbstractULADOA(varargin{:});
    end
    
    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractULADOA(obj);
        cond = obj.DOAOutputPort && obj.NumSignals > length(obj.ScanAngles);
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractULASpectralDOA:InvalidNumSignals');
        end
    end
    
    function setupImpl(obj,X)
        setupImpl@phased.internal.AbstractULADOA(obj,X);
        obj.cCovEstimator.ForwardBackwardAveraging = ...
            obj.ForwardBackwardAveraging;
        ang = obj.ScanAngles;
        % calculate the matrix of steering vectors
        hSteeringVector = phased.SteeringVector(...
            'SensorArray',obj.SensorArray,...
            'PropagationSpeed',obj.PropagationSpeed,...
            'NumPhaseShifterBits',obj.pNumPhaseShifterBits);
        sv = cast(step(hSteeringVector,cast(obj.OperatingFrequency,'double'),...
            cast([ang; zeros(1,length(ang))],'double')),class(X));
        obj.pSteeringVector = sv(1:obj.pEffChannels,:);
        obj.pPattern = zeros(length(obj.ScanAngles),1,class(X));
    end
    
    function flag = isOutputComplexityLockedImpl(obj,index)
        flag = false;
        if index == 1
            flag = false;
        end
        if obj.DOAOutputPort && (index == 2)
            flag = true;
        end
            
    end
    
    function releaseImpl(obj)
        releaseImpl@phased.internal.AbstractULADOA(obj);
    end
    
    function num = getNumOutputsImpl(obj)
        if obj.DOAOutputPort
            num = 2;
        else
            num = 1;
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractULADOA(obj);
        if isLocked(obj)
            %s.SpatialSmoothing = obj.SpatialSmoothing;
            s.pPattern = obj.pPattern;
            s.pSteeringVector = obj.pSteeringVector;
            s.pNumPhaseShifterBits = obj.pNumPhaseShifterBits;
        end
    end

    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        s = loadSubObjects(obj,s);
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function [scanpattern,doasOut] = stepImpl(obj,X)
        
        classtouse=class(X);
        % estimate the covariance matrix
        Cx = cast(step(obj.cCovEstimator,X),classtouse);%g1812952
        
        % generate spatial magnitude spectrum
        privDOASpectrum(obj,Cx); 
        scanpattern = sqrt(abs(obj.pPattern));  
        
        if nargout > 1
            [~,locs] = findpeaks(scanpattern,'SortStr','descend');
            numSignals = obj.NumSignals;
            D = min(numSignals,length(locs));
            assert(D <= numSignals);            

            if D>0
                doasOut = NaN(1,numSignals,classtouse);%invalid angles
                doasOut(1:D) = obj.ScanAngles(locs(1:D).');
                if D < numSignals && isempty(coder.target)
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
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractULADOA(obj, prop);
        if ~obj.DOAOutputPort && strcmp(prop,'NumSignals')
            flag = true;
        end
    end
    
end

methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl
    groups = getPropertyGroupsImpl@phased.internal.AbstractULADOA;
    props = {...
      'ForwardBackwardAveraging',...
      'SpatialSmoothing',...
      'ScanAngles',...
      'DOAOutputPort',...
      'NumSignals'};
    groups(1).PropertyList = [groups(1).PropertyList props];
    % SpatialSmoothing as dependent on private-only
    groups(1).DependOnPrivatePropertyList =  {'SpatialSmoothing'};
  end
end

methods
    function hp = plotSpectrum(obj,varargin)
    %plotSpectrum Plot the spatial spectrum
    %   plotSpectrum(H) plots the spatial spectrum resulting from the last
    %   call of the step method. 
    %   
    %   plotSpectrum(...,'Unit',UNIT) plots the spatial spectrum using the
    %   unit specified in UNIT. UNIT can be one of 'mag'|'power'|'db'. 
    %   The default is 'db'.
    %
    %   plotSpectrum(...,'NormalizeResponse',NFLAG) plots the normalized
    %   spectrum if NFLAG is true. NFLAG can be true or false. The default
    %   is false.
    %
    %   plotSpectrum(...,'Title',TITLE) uses TITLE as the title of the
    %   resulting figure.
    %
    %   h = plotSpectrum(...) returns the line handle in the figure.
    %
    %   % Example:
    %   % Estimate the DOAs of two signals received by a standard  
    %   % 10-element ULA with element spacing 1 meter. The antenna 
    %   % operating frequency is 150 MHz. The actual direction of the first
    %   % signal is 10 degrees in azimuth and 20 degrees in elevation. The 
    %   % direction of the second signal is 60 degrees in azimuth and -5 
    %   % degrees in elevation.
    %
    %   fs = 8000; t = (0:1/fs:1).';
    %   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
    %   ha = phased.ULA('NumElements',10,'ElementSpacing',1);
    %   ha.Element.FrequencyRange = [100e6 300e6];
    %   fc = 150e6;
    %   x = collectPlaneWave(ha,[x1 x2],[10 20;60 -5]',fc);
    %   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
    %   hdoa = phased.BeamscanEstimator('SensorArray',ha,...
    %                  'OperatingFrequency',fc,...
    %                  'DOAOutputPort',true,'NumSignals',2);
    %   [y,doas] = step(hdoa,x+noise);
    %   doas = broadside2az(sort(doas),[20 -5])
    %   % Plot the spatial spectrum
    %   plotSpectrum(hdoa);
    
        coder.extrinsic('phased.internal.AbstractULASpectralDOA.plotSpectrumExtrinsic');
                
        cond = ~isLocked(obj);
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractSpectralDOA:NotLocked');
        end
        titlestr = getSpectrumPlotTitle(obj);
        if  nargout > 0
            if ~isempty(coder.target) 
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','plotSpectrum');
            end
             hp = obj.plotSpectrumExtrinsic(obj.ScanAngles,obj.pPattern,titlestr,varargin{:});
        else
            obj.plotSpectrumExtrinsic(obj.ScanAngles,obj.pPattern,titlestr,varargin{:});
        end
        
                
    end
end

methods(Static, Hidden)
    
    function  hp = plotSpectrumExtrinsic(scanAngles,pattern,titlestr,varargin)

        for m = 1:numel(varargin)
            if strcmp(varargin{m},'Unit')
                varargin{m} = 'Units';
            end
            if strcmp(varargin{m},'NormalizeResponse')
                varargin{m} = 'NormalizeResp';
            end
        end
        
        Hresp = phased.internal.RespPattern2D;
        Hresp.Angle = scanAngles(:);
        Hresp.Pattern = sqrt(pattern);  % magnitude spectrum
        h = plot(Hresp,'Title',titlestr,varargin{:});
        xlabel('Broadside Angle (degrees)');
        axis tight;
        
        if nargout > 0
            hp = h;
        end
        
    end
end
methods (Access = protected) %for Simulink
    function varargout = getOutputNamesImpl(~)
        varargout = {'Y','Ang'};
    end
    function varargout = getInputNamesImpl(~)
        varargout = {''};
    end
    function flag = isInputSizeLockedImpl(~,~)
        flag = false;
    end
    function varargout = getOutputSizeImpl(obj)
        varargout{1} = [size(obj.ScanAngles,2),1];
        if obj.DOAOutputPort
            varargout{2} = [1,obj.NumSignals];
        end
    end
    function varargout = isOutputFixedSizeImpl(~)
        varargout = {true, true};
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout = {propagatedInputDataType(obj,1),propagatedInputDataType(obj,1)};
    end
    function varargout = isOutputComplexImpl(~)
        varargout = {false, false};
    end
end
end

