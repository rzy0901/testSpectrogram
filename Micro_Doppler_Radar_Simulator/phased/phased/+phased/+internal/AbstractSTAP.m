classdef (Hidden) AbstractSTAP < phased.internal.AbstractNarrowbandArrayProcessing & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%AbstractSTAP   Define the AbstractSTAP class.

%   Copyright 2008-2017 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %PRFSource    Source of PRF
        %   Specify how to determine the PRF for the STAP processor as one
        %   of 'Property' | 'Input port', where the default is 'Property'.
        %   When you set this property to 'Property', the PRF is determined
        %   by the value of the PRF property. When you set this property to
        %   'Input port', the PRF is determined by the input argument.
        PRFSource = 'Property';
        %PRF    Pulse repetition frequency (Hz)
        %   Specify the pulse repetition frequency (PRF) (in Hz) of the
        %   received signal as a scalar. The default value of this property
        %   is 1.
        PRF = 1;
        %DirectionSource    Source of direction
        %   Specify how to determine the direction for the STAP processor
        %   as one of 'Property' | 'Input port', where the default is
        %   'Property'. When you set this property to 'Property', the
        %   direction is determined by the value of the Direction property.
        %   When you set this property to 'Input port', the direction is
        %   determined by the input argument.
        DirectionSource = 'Property';
        %DopplerSource    Source of targeting Doppler
        %   Specify how to determine the targeting Doppler for the STAP
        %   processor as one of 'Property' | 'Input port', where the
        %   default is 'Property'. When you set this property to
        %   'Property', the Doppler is determined by the value of the
        %   Doppler property. When you set this property to 'Input port',
        %   the Doppler is determined by the input argument.
        DopplerSource = 'Property';
        %Doppler    Targeting Doppler (Hz)
        %   Specify the targeting Doppler of the STAP processor as a
        %   scalar. This property applies when you set the DopplerSource
        %   property to 'Property'. The default value of this property is
        %   0.
        Doppler = 0;
    end

    properties (Nontunable, Logical) 
        %WeightsOutputPort    Enable weights output
        %   Set this property to true to output the weights used in the
        %   STAP processor. Set this property to false to not output the
        %   weights. The default value of this property is false.
        WeightsOutputPort = false;
    end
    
    properties(Nontunable)
        %NumPhaseShifterBits    Number of bits in phase shifters
        %   Specify the number of bits used in the phase shifter as a
        %   non-negative integer. The default value of this property is 0,
        %   indicating there is no quantization effect in the phase
        %   shifter.
        NumPhaseShifterBits = 0
    end
    
    properties(Abstract, Nontunable)
        %abstract property for direction. The interpretations for direction
        %are different between SMI and DPCA based algorithms
        Direction
    end
    
    properties(Constant, Hidden)
        PRFSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        DirectionSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        DopplerSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    properties(Access = protected)
        pMaximumCellIndex;
    end
    
    properties(Access = protected)
        pCubeDim = [-1 -1 -1];
    end
    
    methods (Access = protected)
        % Constructor
        function obj = AbstractSTAP(varargin)
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});
        end
        
        function privValidateSensorArray(obj,val) 
            if isSubarraySupported(obj)
                validateattributes( val, { 'phased.internal.AbstractArray',...
                    'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
            else
                validateattributes( val, { 'phased.internal.AbstractArray'},...
                    { 'scalar' }, '', 'SensorArray');
            end
        end
        
        function flag = isSubarraySupported(obj) %#ok<MANU>
            flag = false;
        end
    end
    
    methods
        function set.PRF(obj,val)
            sigdatatypes.validateFrequency(val,'','PRF',...
                {'scalar'});
            obj.PRF = val;
        end
        
        function set.Doppler(obj,val)
            validateattributes(val,{'double'},{'finite','nonempty','scalar','real'},...
                '','Doppler');
            obj.Doppler = val;
        end
    end
    
    methods (Static, Hidden)
        function [datavec, datacubeDim] = organizeSpaceTimeSnapshots(x)
        %organizeSpaceTimeSnapshots Organize data in space time snapshots
        %   DATA = organizeSpaceTimeSnapshots(X) transform the input X from
        %   a cube to a matrix DATA. The cube X is in the dimensions of
        %   [range, channels, pulses] while DATA is in the dimensions of
        %   [space time snapshots, range]. Each column of DATA is a space
        %   time snapshot, which concatenates channel measurements retained
        %   in different pulses.
        %
        %   [DATA, XDim] = organizeSpaceTimeSnapshots(X) also returns the
        %   dimensions of the input X.
        %
        %   % Example:
        %   %   Convert the received datacube x into space time snapshots
        %   %   for different ranges.
        %
        %   x = zeros(100,3,2);
        %   y = stap.AbstractSTAP.organizeSpaceTimeSnapshots(x);
        
            % validation of x is done in calling function
            datacubeDim = size(x);
            % Transform datacube in snapshots to vectors
            xShift = shiftdim(x,1);   % change to channel x doppler x range
            datavec = reshape(xShift,datacubeDim(2)*datacubeDim(3),...
                datacubeDim(1)); % (channel x doppler) x range
        end
        
    end
    
    methods (Access = protected)
        function validatePropertiesImpl(obj)
            if strncmpi(obj.DopplerSource,'p',1)
                if (obj.PRFSource(1) == 'P')
                    prf = obj.PRF;
                    cond = abs(obj.Doppler/prf) > 0.5;
                    if cond
                        coder.internal.errorIf(cond,'phased:stap:step:OutOfBoundDoppler', 'Doppler', sprintf( '%5.2f', -prf/2 ), sprintf( '%5.2f', prf/2 ));
                    end
                end
            end
        end
        
        function setupImpl(obj,x)
            sz_x = size(x);
            obj.pCubeDim = sz_x;
            obj.pMaximumCellIndex = sz_x(1);
        end
        
        function processInputSizeChangeImpl(obj,x,~,~,~,~)
            sz_x = size(x);
            obj.pCubeDim = sz_x;
            obj.pMaximumCellIndex = sz_x(1);    
        end
        
        function flag = isInputSizeLockedImpl(~,index)
            if index == 1
                flag = false;
            else
                flag = true;
            end
        end
        
        function flag = isInputComplexityLockedImpl(~,index) 
            flag = false;  % index == 1
            if (index == 2)
                flag  = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(~,~) 
            flag = false;  % (index == 1) || (index == 2)
        end
                
        function validateInputsImpl(obj,x,cutidx)
            %   validates X to be a 3-dimensional data cube. The dimensions
            %   of X is (range, channels, pulses). The second dimension of
            %   X (number of channels) must match the number of elements in
            %   the sensor array as specified in STAP processor H.          
            cond = ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','double');
            end
            sz_x = size(x);
            cond = ndims(x) ~= 3;
            if cond
                coder.internal.errorIf(cond,'phased:stap:InvalidData');
            end
            N = getDOF(obj.SensorArray);
            cond = sz_x(2) ~= N;
            if cond
                coder.internal.errorIf(cond,'phased:stap:InvalidDataDimension', N);
            end 
            
            if  obj.pCubeDim(3) ~= -1
                cond = obj.pCubeDim(3) ~= sz_x(3);
                if cond
                    coder.internal.errorIf(cond,'phased:step:NumInputPageNotConstant');
                end
            end
            
            cond = ~isa(cutidx,'double');
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','Idx','double');
            end
            cond = ~isscalar(cutidx);
            if cond
                coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeScalar','Idx');
            end
        end
        
        function num = getNumOutputsImpl(obj) 
            if obj.WeightsOutputPort
                num = 2;
            else
                num = 1;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            if isLocked(obj)
                s.pMaximumCellIndex = obj.pMaximumCellIndex; 
                s.pCubeDim = obj.pCubeDim;
            end
        end
        
   end
    
    methods(Access=protected)
        function validateInputPRFSpec(~,prfArg)
            cond = ~isa(prfArg,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','PRF','double');
            end
            cond = ~isscalar(prfArg);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeScalar','PRF');
            end
        end
        function validateInputAngleSpec(~,angle)
            sz_angle = size(angle);
            cond = ~isa(angle,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Ang','double');
            end
            cond = ~iscolumn(angle) || isempty(angle);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeColVector','Ang');
            end
            cond = sz_angle(1) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:stap:NeedTwoRows','Ang');
            end
            cond = ~isreal(angle);
            if cond
                coder.internal.errorIf(cond,'phased:stap:NeedReal', 'Ang');
            end
        end
        
        function validateInputDopplerSpec(~,dop)
            cond = ~isa(dop,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Dop','double');
            end
            cond = ~isscalar(dop);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeScalar','Dop');
            end
            cond = ~isreal(dop);
            if cond
                coder.internal.errorIf(cond,'phased:stap:NeedReal', 'Dop');
            end
        end
        
        function prf = validateInputPRF(~,prfArg)
            cond = ~isfinite(prfArg);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:step:expectedFinite','PRF');
            end
            cond = ~isreal(prfArg);
            if cond
                coder.internal.errorIf(cond,'phased:stap:NeedReal', 'PRF');
            end
            cond = prfArg < 0;
            if cond
                coder.internal.errorIf(cond,...
                    'phased:step:expectedPositive','PRF');
            end 
            prf = prfArg;
        end
        
        function angle = validateInputAngle(~,angleArg)
            if isscalar(angleArg)
                angle = [angleArg; 0];
            else
                angle = angleArg;
            end
            
            cond = angle(1) < -180 || angle(1) > 180;
            if cond
                coder.internal.errorIf(cond,'phased:stap:step:OutOfBoundAzimuth');
            end
            cond = angle(2) < -90 || angle(2) > 90;
            if cond
                coder.internal.errorIf(cond,'phased:stap:step:OutOfBoundElevation');
            end
        end
        
        function dop = validateInputDoppler(~,dop,prf)
            % Max Doppler within half PRF
            cond = abs(dop/prf) > 0.5;
            if cond
                coder.internal.errorIf(cond,'phased:stap:step:OutOfBoundDoppler', 'DOP', sprintf( '%5.2f', -prf/2 ), sprintf( '%5.2f', prf/2 ));
            end
        end
            
        function cutidx = validateInputCUTIdx(obj,cutidx)
            sigdatatypes.validateIndex(cutidx,'',...
                'Idx',{'<=',obj.pMaximumCellIndex});
        end
    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl(sensorType)
        groups = getPropertyGroupsImpl@phased.internal.AbstractNarrowbandArrayProcessing(sensorType);
        pPRFSource = matlab.system.display.internal.Property('PRFSource', ...
            'Description', 'Specify PRF as');
        pDirectionSource = matlab.system.display.internal.Property('DirectionSource', ...
                'Description', 'Specify direction as');
        pDopplerSource = matlab.system.display.internal.Property('DopplerSource', ...
                'Description', 'Specify targeting Doppler as');

        
        props = {pPRFSource,...
                 'PRF',...
                 pDirectionSource,...
                 'Direction',...
                 'NumPhaseShifterBits',...
                 pDopplerSource,...
                 'Doppler',...
                 'WeightsOutputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
    end
    methods (Access = protected) %For Simulink propagation and mask
        function varargout = getOutputNamesImpl(~)
            varargout = {'Y', 'W'};
        end
        
        function varargout = getInputNamesImpl(obj)
            varargout = {'X','Idx'};
            if (obj.PRFSource(1) == 'I') %Input Port
                varargout{end+1} = 'PRF';
            end
            if (obj.DirectionSource(1) == 'I') %Input Port
                varargout{end+1} = 'Ang'; 
            end
            if (obj.DopplerSource(1) == 'I') %Input Port
                varargout{end+1} = 'Dop'; 
            end
        end
        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            varargout{1} = [szX(1) 1];
            if obj.WeightsOutputPort
                varargout{2} = [szX(2)*szX(3) 1];
            end
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj, 1);
            varargout{2} = true;
        end
        function varargout = getOutputDataTypeImpl(obj)  %#ok<MANU>
            varargout = {'double','double'};
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout = {true,true};
        end
    end

end




