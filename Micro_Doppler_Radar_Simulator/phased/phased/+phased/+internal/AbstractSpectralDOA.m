classdef (Hidden) AbstractSpectralDOA < phased.internal.AbstractDOA & ...
     matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%ABSTRACTSPECTRALDOA Abstract class for estimating DOA from peaks of
%spatial spectrum

%   Copyright 2009-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %AzimuthScanAngles Azimuth scan angles (deg)
    %   Specify the azimuth scan angles (in degrees) as a real row vector.
    %   The default value is -90:90. The angles must be within [-180 180].
    %   You must specify the angles in an ascending order.
    AzimuthScanAngles = -90:90;
    %ElevationScanAngles Elevation scan angles (deg)
    %   Specify the elevation scan angles (in degrees) as a real row vector
    %   or scalar. The default value is 0. The angles must be within [-90
    %   90]. You must specify the angles in an ascending order.
    ElevationScanAngles = 0;
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

properties(Access = protected, Nontunable)
    % the number of phase shifter bits
    pNumPhaseShifterBits=0;
end

properties(Access = protected)
    % the pattern of the spectrum
    pPattern;
    % steering vector when there is only one iteration
    pSteeringVectors;
    % the scanning angles
    pScanAngles;
end

properties(Access = protected, Nontunable)
    % handle to the steering vector object
    cSteeringVector;
    % the size of the resulting pattern
    pPatternSize;
    % the number of iterations to accommodate large computation
    pNumIter;
    % the number of scan angles for each iteration
    pScanAngleBlockSize;
    % the index for each scan angle block
    pScanAngleBlockIndex;
end

properties (Access = protected, Nontunable, Logical)
    % flag to identify whether there is only one iteration
    pOneIterFlag;
    % flag to identify whether we have a single cut
    pIsOneCut;
end

methods (Abstract, Access = protected)
    pattern = privDOASpectrum(obj,Cx,azang,elang);
    titlestr = getSpectrumPlotTitle(obj)
end

methods
    function set.AzimuthScanAngles(obj,value)
        validateattributes( value, { 'double','single' }, { 'row', 'real', 'finite', '>=', -180, '<=', 180 }, '', 'AzimuthScanAngles');
        %ensure ascending
        cond = any(diff(value) <=0);
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractSpectralDOA:NotAscending', 'AzimuthScanAngles');
        end
        obj.AzimuthScanAngles = value;
    end
    function set.ElevationScanAngles(obj,value)
        validateattributes( value, { 'double','single' }, { 'row', 'real', 'finite', '>=', -90, '<=', 90 }, '', 'ElevationScanAngles');
        %ensure ascending
        cond = any(diff(value) <=0);
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractSpectralDOA:NotAscending', 'ElevationScanAngles');
        end
        obj.ElevationScanAngles = value;
    end 
end

methods (Access = protected)
       
    function obj = AbstractSpectralDOA(varargin)
        obj = obj@phased.internal.AbstractDOA(varargin{:});
    end
    
    % function flag = isSubarraySupported(obj) %#ok<MANU>
    %     flag = true;
    % end

    function validatePropertiesImpl(obj)
        az = obj.AzimuthScanAngles;
        el = obj.ElevationScanAngles;
        len_az = numel(az);
        len_el = numel(el);
        cond = len_az == 1 && ...
            len_el < 3;
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractSpectralDOA:NotEnoughElevationSamples');
        end
        cond = len_el == 1 && ...
            len_az < 3;
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractSpectralDOA:NotEnoughAzimuthSamples');
        end
        cond = isa(obj.SensorArray,'phased.ULA') && ...
            obj.DOAOutputPort == true && ...
            len_el ~= 1;
        if cond
            coder.internal.errorIf(cond,'phased:phased:internal:AbstractSpectralDOA:ArrayAmbiguity');
        end
    end
    
    function setupImpl(obj,X)
        setupImpl@phased.internal.AbstractDOA(obj,X);
        obj.cCovEstimator = phased.internal.SpatialCovEstimator(...
            'NumSubarrays',1,... % No subarray smoothing for 2D.
            'ForwardBackwardAveraging',obj.ForwardBackwardAveraging);
        obj.cSteeringVector = phased.SteeringVector(...
            'SensorArray',obj.SensorArray,...
            'PropagationSpeed',obj.PropagationSpeed,...
            'NumPhaseShifterBits',obj.pNumPhaseShifterBits);
        az = obj.AzimuthScanAngles;
        el = obj.ElevationScanAngles;
        len_az = numel(az);
        len_el = numel(el);
        if (len_el ==  1) || (len_az == 1)
            obj.pIsOneCut = true;
        else
            obj.pIsOneCut = false;
        end
        NumEl = numel(el);
        NumAz = numel(az);
        obj.pPatternSize = [NumEl NumAz];
        [ScanAz, ScanEl] = meshgrid(az,el);
        obj.pScanAngles = [ScanAz(:) ScanEl(:)].';
        if NumEl*NumAz < 400
            obj.pNumIter = 1;
            obj.pScanAngleBlockSize = NumEl*NumAz;
            obj.pOneIterFlag = true;
        else
            obj.pNumIter = min(NumEl, NumAz);
            obj.pScanAngleBlockSize = max(NumEl, NumAz);
            obj.pOneIterFlag = false;
        end
        obj.pScanAngleBlockIndex = 1:obj.pScanAngleBlockSize;
        %obj.pScanAngleBlockIndex = [(0:obj.pNumIter-1).'*obj.pScanAngleBlockSize+1 ...
        %                           (1:obj.pNumIter).'*obj.pScanAngleBlockSize];
        %obj.pScanAngleBlockIndex(:,1) = obj.pScanAngleBlockIndex(:,1) + 1;
        obj.pPattern = zeros(prod(obj.pPatternSize),1);
        if obj.pOneIterFlag
            obj.pSteeringVectors = step(obj.cSteeringVector,...
                cast(obj.OperatingFrequency,'double'),cast(obj.pScanAngles,'double'));
        end
        
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
        release(obj.cSteeringVector);
        release(obj.cCovEstimator);
    end
    
    function resetImpl(obj)
        reset(obj.cSteeringVector);
        reset(obj.cCovEstimator);
    end
    
    function num = getNumOutputsImpl(obj)
        if obj.DOAOutputPort
            num = 2;
        else
            num = 1;
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractDOA(obj);
        if isLocked(obj)
            s.cSteeringVector = saveobj(obj.cSteeringVector);
            s.pPattern = obj.pPattern;
            s.pSteeringVectors = obj.pSteeringVectors;
            s.pPatternSize = obj.pPatternSize;
            s.pScanAngles = obj.pScanAngles;
            s.pNumIter = obj.pNumIter;
            s.pScanAngleBlockSize = obj.pScanAngleBlockSize;
            s.pScanAngleBlockIndex = obj.pScanAngleBlockIndex;
            s.pOneIterFlag = obj.pOneIterFlag;
            s.pIsOneCut = obj.pIsOneCut;
            s.pNumPhaseShifterBits = obj.pNumPhaseShifterBits;
        end
    end

    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractDOA(obj,s);
        if isfield(s,'isLocked')
            if s.isLocked
                obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                s = rmfield(s,'cSteeringVector');
            end
            s = rmfield(s,'isLocked');
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
        Cx = step(obj.cCovEstimator,X);
        
        % generate spatial spectrum
        privDOASpectrum(obj,Cx);
        scanpattern = cast(reshape(sqrt(abs(obj.pPattern)),obj.pPatternSize),classtouse);  
        
        azang = obj.AzimuthScanAngles;
        elang = obj.ElevationScanAngles;
        if nargout > 1
            numSignals = obj.NumSignals;
            if obj.pIsOneCut % 1D peak finding
                [~,locs] = findpeaks(scanpattern,'SortStr','descend');
                D = min(numSignals,length(locs));
                assert(D <= numSignals);
                if D>0
                    if length(elang) == 1 % peak in azimuth direction
                        doas = [azang(locs(1:D)); elang*ones(1,D)];
                    else % peak in elevation direction
                        doas = [azang*ones(1,D); elang(locs(1:D))];
                    end
                    doasOut = NaN(2,numSignals,classtouse);%invalid angles
                    doasOut(:,1:D) = doas;
                else % D == 0
                    if isempty(coder.target)
                        warning(message('phased:phased:doa:ZeroSourceNumber'));
                    end
                    doasOut = NaN(2,numSignals,classtouse);%invalid angles
                end
            else                % Find peaks in both azimuth and elevation
                % Note the dimension of scan pattern is (NumElev x NumAzim)
                [rowIndex,colIndex] = ...
                    obj.privFindPatternPeaks2D(scanpattern,numSignals);
                D = length(rowIndex);
                if D>0
                    doas = [azang(colIndex); elang(rowIndex)];
                    doasOut = NaN(2,numSignals,classtouse);%invalid angles
                    doasOut(:,1:D) = doas;
                else % D == 0
                    if isempty(coder.target)
                        warning(message('phased:phased:doa:ZeroSourceNumber'));
                    end
                    doasOut = NaN(2,numSignals,classtouse);%invalid angles
                end
            end
            
            if D < numSignals && D~=0 && isempty(coder.target)
                warning(message('phased:phased:internal:AbstractSpectralDOA:LessNumSources', D));
            end
        end
                
    end
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractDOA(obj, prop);
        if ~obj.DOAOutputPort && strcmp(prop, 'NumSignals')
            flag = true;
        end
    end    
        
end

methods (Static,Hidden,Access=protected)
  function groups = getPropertyGroupsImpl
    groups = getPropertyGroupsImpl@phased.internal.AbstractDOA('array');
    props = {...
      'ForwardBackwardAveraging',...
      'AzimuthScanAngles',...
      'ElevationScanAngles',...
      'DOAOutputPort',...
      'NumSignals'};
    groups(1).PropertyList = [groups(1).PropertyList props];
  end
end

methods
    function hp = plotSpectrum(obj,varargin)
    %plotSpectrum Plot the spatial spectrum
    %   plotSpectrum(H) plots the spatial spectrum resulting from the last
    %   call of the step method. 
    %   
    %   plotSpectrum(...,'Unit',UNIT) plots the spatial spectrum using the
    %   unit specified in UNIT. UNIT can be one of the following:
    %   'mag'|'power'|'db'. The default is 'db'.
    %
    %   plotSpectrum(...,'NormalizeResponse',NFLAG) plots the normalized
    %   spectrum if NFLAG is true. NFLAG can be true or false. The default
    %   is false.
    %
    %   plotSpectrum(...,'Title',TITLE) uses TITLE as the title of the
    %   resulting figure.
    %
    %   h = plotSpectrum(...) returns the line handle in a 1D figure, or
    %   the surface handle in a 2D figure.
    %
    %   % Example:
    %   % Estimate the DOAs of two signals received by a 50-element URA. 
    %   % The antenna operating frequency is 150 MHz. The actual direction 
    %   % of the first signal is -37 degrees in azimuth and 0 degrees in 
    %   % elevation. The direction of the second signal is 17 degrees in 
    %   % azimuth and 20 degrees in elevation.
    %
    %   fs = 8000; t = (0:1/fs:1).';
    %   x1 = cos(2*pi*t*300); x2 = cos(2*pi*t*400);
    %   ha = phased.URA('Size',[5 10],'ElementSpacing',[1 0.6]);
    %   ha.Element.FrequencyRange = [100e6 300e6];
    %   fc = 150e6;
    %   x = collectPlaneWave(ha,[x1 x2],[-37 0;17 20]',fc);
    %   noise = 0.1*(randn(size(x))+1i*randn(size(x)));
    %   hdoa = phased.BeamscanEstimator2D('SensorArray',ha,...
    %                  'OperatingFrequency',fc,...
    %                  'DOAOutputPort',true,'NumSignals',2,...
    %                  'AzimuthScanAngles',-50:50,...
    %                  'ElevationScanAngles',-30:30);
    %   [~,doas] = step(hdoa,x+noise);
    %   % Plot the spatial spectrum
    %   plotSpectrum(hdoa);

        coder.extrinsic('phased.internal.AbstractSpectralDOA.plotSpectrumExtrinsic');
        
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
             hp = obj.plotSpectrumExtrinsic(obj.AzimuthScanAngles,obj.ElevationScanAngles, ...
                                            obj.pPattern,obj.pPatternSize,titlestr,varargin{:});
        else
            obj.plotSpectrumExtrinsic(obj.AzimuthScanAngles,obj.ElevationScanAngles, ...
                                      obj.pPattern,obj.pPatternSize,titlestr,varargin{:});
        end        
    end
end
methods(Static, Hidden)
    
    function  hp = plotSpectrumExtrinsic(azimuthScanAngles,elevationScanAngles, ...
                                         pattern,patternSize,titlestr,varargin)

        for m = 1:numel(varargin)
            if strcmp(varargin{m},'Unit')
                varargin{m} = 'Units';
            end
            if strcmp(varargin{m},'NormalizeResponse')
                varargin{m} = 'NormalizeResp';
            end
        end
        
        
        azang = azimuthScanAngles(:);
        elang = elevationScanAngles(:);
        
        pattern = reshape(sqrt(pattern),patternSize);
        if isvector(pattern)
            Hresp = phased.internal.RespPattern2D;
            if length(azang) == 1
                Hresp.SliceDir = 'El'; % This is a weird thing of RespPattern3D
                Hresp.Angle = elang;
                ttstring = sprintf('%s %s',titlestr,...
                    getString(message('phased:phased:atAzimuth',sprintf('%2.2f',azang))));
            else
                Hresp.SliceDir = 'Az';
                Hresp.Angle = azang;
                ttstring = sprintf('%s %s',titlestr,...
                    getString(message('phased:phased:atElevation',sprintf('%2.2f',elang))));
            end
            Hresp.Pattern = pattern(:);  % magnitude spectrum
            h = plot(Hresp,'Title',ttstring,varargin{:});
            axis tight;
        else
            Hresp = phased.internal.RespPattern3D;
            Hresp.AzAngle = azang;
            Hresp.ElAngle = elang;
            Hresp.Pattern = pattern;  % magnitude spectrum
            h = plot(Hresp,'Title',titlestr,varargin{:});
        end
        
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
        varargout{1} = [size(obj.ElevationScanAngles,2), ...
                        size(obj.AzimuthScanAngles,2)];
        if obj.DOAOutputPort
            varargout{2} = [2,obj.NumSignals];
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

methods (Static, Access = protected)
    function [RowIndex, ColIndex] = privFindPatternPeaks2D(X,DArg)
    % find at most D dominant peaks from 2D matrix X and return their
    % corresponding row and column indexes. The returned peaks are in
    % descending order in peak heights.

        %find all regional peaks.
        xbk = phased.internal.findpeaks2D(X);
        [rows, cols] = find(xbk);

        %identify at most D dominant peaks.
        pkvalues = X(sub2ind(size(X),rows,cols));
        [~, locs] = sort(pkvalues,'descend');
        D = min(DArg, length(locs));
        assert(D <= DArg);
        locs = locs(1:D);

        % extract their indexes
        RowIndex = rows(locs);
        ColIndex = cols(locs);

    end
end
end


