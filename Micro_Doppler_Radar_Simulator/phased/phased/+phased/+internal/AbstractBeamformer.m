classdef (Hidden) AbstractBeamformer < phased.internal.AbstractArrayOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%AbstractBeamformer Abstract class for beamformers

%   Copyright 2009-2016 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %DirectionSource    Source of beamforming direction
        %   Specify how to determine the beamforming direction for the
        %   beamformer as one of 'Property' | 'Input port', where the
        %   default is 'Property'. When you set this property to
        %   'Property', the direction is determined by the value of the
        %   Direction property. When you set this property to 'Input port',
        %   the direction is determined by the input argument.
        DirectionSource = 'Property';
        %Direction    Beamforming direction (deg)
        %   Specify the beamforming direction of the beamformer as a 2-row
        %   matrix. Each column of Direction is specified in the format of
        %   [AzimuthAngle; ElevationAngle] (in degrees). Azimuth angle
        %   should be between -180 and 180. Elevation angle should be
        %   between -90 and 90. This property applies when you set the
        %   DirectionSource property to 'Property'. The default value of
        %   this property is [0; 0].
        Direction = [0; 0];
    end

    properties(Constant, Hidden)
        DirectionSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    properties(Access = protected, Nontunable)
        pDOF
    end
    
    methods
        function set.Direction(obj,val)
            if isMultipleInputAnglesAllowed(obj)
                sigdatatypes.validateAzElAngle(val,'Beamformer',...
                    'Direction',{'double','single'},{'real'});
            else
                sigdatatypes.validateAzElAngle(val,'Beamformer',...
                    'Direction',{'double','single'},{'size',[2 1]});
            end
            obj.Direction = val;
        end
    end
       
    methods (Access = protected)
        function privValidateSensorArray(obj,val) 
            if isSubarraySupported(obj)
                validateattributes( val, { 'phased.internal.AbstractArray',...
                    'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
            else
                validateattributes( val, { 'phased.internal.AbstractArray'},...
                    { 'scalar' }, '', 'SensorArray');
            end
        end
        
        function flag = isMultipleInputAnglesAllowed(obj) %#ok<MANU>
            flag = false;
        end
        
        function flag = isSubarraySupported(obj) %#ok<MANU>
            flag = false;
        end
    end
    
    methods (Access = protected)
        function obj = AbstractBeamformer(varargin)
            obj@phased.internal.AbstractArrayOperation(varargin{:});
        end

    end

    methods (Access = protected)
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if strncmpi(obj.DirectionSource,'Input port',1) && ...
                    strcmp(prop,'Direction')
                flag = true;
            end
        end
        
        function setupImpl(obj,x)
            obj.pDOF = getDOF(obj.SensorArray);
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)  %#ok<INUSL>
            flag = false; 
            if index == 1
                flag = false;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function validateInputsImpl(obj,x)
            validateInputSignal(obj,x,'X');
        end

        function num = getNumInputsImpl(obj)
            if strncmpi(obj.DirectionSource, 'p',1)
                num = 1;
            else
                num = 2;
            end
        end
        
        function num = getNumOutputsImpl(obj) %#ok<MANU>
            num = 1;
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractArrayOperation(obj);
            if isLocked(obj)
                s.pDOF = obj.pDOF;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractArrayOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    if isfield(s,'pSensorArrayNumElements')
                        s = rmfield(s,'pSensorArrayNumElements');
                    end
                end
            end
        end
        
        function validateInputSignal(obj,x,str_x)
            cond = ~(isa(x,'float'));
            if cond 
                coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType',str_x,'float');
            end
            cond = ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                          'MATLAB:system:inputMustBeMatrix',str_x);
            end
            sz_x = size(x);
            cond = sz_x(2) ~= getDOF(obj.SensorArray);
            if cond
                coder.internal.errorIf(cond,'phased:beamformer:InvalidDataDimension', str_x, getDOF( obj.SensorArray ));
            end
            validateNumChannels(obj,x)
        end
            
        function validateInputAngle(obj,ang)
          coder.extrinsic('mat2str')
          cond = ~(isa(ang,'float'));
          if cond
              coder.internal.errorIf(cond, ...
                  'MATLAB:system:invalidInputDataType','Ang','float');
          end
          sz_ang = size(ang);
          if isMultipleInputAnglesAllowed(obj)
              cond = sz_ang(1) ~= 2;
              if cond
                  coder.internal.errorIf(cond,'phased:phased:invalidRowNumbers','Ang',2);
              end
          else
              cond = ~isequal(sz_ang,[2 1]);
              if cond
                  coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDimensions','Ang',...
                    '[2 1]', coder.internal.const(mat2str(sz_ang)));
              end
          end
          cond = ~isreal(ang);
          if cond
              coder.internal.errorIf(cond,'phased:beamformer:NeedReal','Ang');
          end
        end

        function validateAngleRange(obj,ang)  %#ok<INUSL>
            cond = ang(1) < -180 || ang(1) > 180;
            if cond
                coder.internal.errorIf(cond,'phased:beamformer:step:OutOfBoundAzimuth');
            end
            cond = ang(2) < -90 || ang(2) > 90;
            if cond
                coder.internal.errorIf(cond,'phased:beamformer:step:OutOfBoundElevation');
            end
        end
                
    end
    methods (Access = protected) %for Simulink
        function varargout = getOutputNamesImpl(~)
            varargout = {'Y','W'};
        end
        function varargout = getInputNamesImpl(~)
            varargout = {'X','Ang'};
        end
        function varargout = getOutputSizeImpl(obj)
            szX = propagatedInputSize(obj,1);
            if strncmpi(obj.DirectionSource, 'p',1)
                szAng = size(obj.Direction);
            else
                szAng = propagatedInputSize(obj,getAngInputIdx(obj));
            end
            varargout{1} = [szX(1) szAng(2)];
            % Output weights, if specified
            varargout{2} = [szX(2) szAng(2)];
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj, 1);
            if strncmpi(obj.DirectionSource, 'p',1)
                varargout{2} = true;
            else
                varargout{2} = ...
                  propagatedInputFixedSize(obj, getAngInputIdx(obj));
            end
        end
        function varargout = getOutputDataTypeImpl(obj)
            dt = propagatedInputDataType(obj,1);
            varargout = {dt dt};
        end
        function varargout = isOutputComplexImpl(obj) %#ok<MANU>
            varargout = {true, true};
        end
        function aIdx = getAngInputIdx(obj)  %#ok<MANU>
            aIdx = 2;
        end
    end   
end

