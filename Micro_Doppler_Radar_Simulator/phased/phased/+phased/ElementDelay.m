classdef (Sealed,StrictDefaults) ElementDelay < phased.internal.AbstractArrayOperation & ...
     matlab.system.mixin.Propagates
%ElementDelay   Sensor array element delay estimator
%   H = phased.ElementDelay creates an element delay estimator System
%   object, H. This object calculates the signal delay for elements in an
%   array when the signal arrives the array from specified directions. By
%   default a 2-element uniform linear array (ULA) is used.
%
%   H = phased.ElementDelay(Name,Value) creates an element delay estimator
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   TAU = step(H,ANG) returns the delay (in seconds), TAU, of each element
%   relative to the array's phase center for given signal incident
%   directions specified in ANG (in degrees). ANG can be either a length M
%   row vector or a 2xM matrix. TAU is an NxM matrix where N is the number
%   of elements in the array. Each column of TAU contains the delays of
%   array elements for the corresponding directions specified in ANG.
%
%   When ANG is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in the [azimuth; elevation] form. The azimuth
%   angle should be between [-180 180] degrees and the elevation angle
%   should be between [-90 90] degrees. If ANG is a length M row vector,
%   each element specifies a direction's azimuth angle and the
%   corresponding elevation angle is assumed to be 0.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ElementDelay methods:
%
%   step     - Calculate the delay for the elements (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create an element delay object with same property values
%   isLocked - Locked status (logical)
%
%   ElementDelay properties:
%
%   SensorArray      - Sensor array
%   PropagationSpeed - Signal propagation speed 
%
%   % Example:
%   %   Calculate the element delay for a 4-element uniform linear array 
%   %   when the input is impinging the array from 30 degrees azimuth and 
%   %   20 degrees elevation.
%
%   ha = phased.ULA(4);
%   elementdelay = phased.ElementDelay('SensorArray',ha);
%   tau = elementdelay([30; 20])
%
%   See also phased, phased.SteeringVector, phased.ArrayResponse,
%   phased.ArrayGain.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Access = private, Nontunable)
        % private property to hold array positions. Calculated at compile
        % time since this is a a constant.
        pArrayPos
        cSensorArray
    end

    methods

        function obj = ElementDelay(varargin)
            %ElementDelay   Construct the ElementDelay class.

            obj@phased.internal.AbstractArrayOperation(varargin{:});
            
        end

    end
    
    methods (Access = protected)
        
      function validateInputsImpl(~,angle)        
        cond = ~isa(angle,'double');
        if cond
            coder.internal.errorIf(cond, ...
            'MATLAB:system:invalidInputDataType','Ang','double');
        end
        cond = ~ismatrix(angle) || isempty(angle);
        if cond
            coder.internal.errorIf(cond, ...
            'MATLAB:system:inputMustBeMatrix','Ang');
        end
        sz_angle = size(angle);
        cond = sz_angle(1) > 2;
        if cond
            coder.internal.errorIf(cond,'phased:system:ElementDelay:NeedTwoRows','Ang');
        end
        cond = ~isreal(angle);
        if cond
            coder.internal.errorIf(cond,'phased:step:NeedReal', 'Ang');
        end
      end
        
        function setupImpl(obj,ang) %#ok<INUSD>

            if isempty(coder.target)
                obj.cSensorArray = cloneSensor(obj.SensorArray);
            else
                if isElementFromAntenna(obj.SensorArray)
                    coder.internal.errorIf(true, ...
                        'phased:system:element:AntennaToolboxCodegenNotSupported','em.Antenna','phased.CustomAntennaElement');
                else
                    obj.cSensorArray = clonecg(obj.SensorArray);
                end
            end
            
            obj.pArrayPos = getElementPosition(obj.cSensorArray);
        
        end
        
        function flag = isInputComplexityLockedImpl(obj,~) %#ok<MANU>
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~) %#ok<MANU>
            flag = true;
        end
        
        function tau = stepImpl(obj,angArg)
            if size(angArg,1) == 1;
                ang = [angArg; zeros(size(angArg))];
            else
                ang = angArg;
            end
            azang = ang(1,:);
            elang = ang(2,:);
            cond = any(azang>180) || any(elang>90);
            if cond
                coder.internal.errorIf(cond,'phased:step:AzElNotLessEqual');
            end
            cond = any(azang<-180) || any(elang<-90);
            if cond
                coder.internal.errorIf(cond,'phased:step:AzElNotGreaterEqual');
            end
            % angles defined in local coordinate system
            incidentdir = [-cosd(elang).*cosd(azang);...
                -cosd(elang).*sind(azang);...
                -sind(elang)];
            tau = obj.pArrayPos.'*incidentdir/obj.PropagationSpeed;
            
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractArrayOperation(obj);
            if isLocked(obj)
                s.cSensorArray = saveobj(obj.cSensorArray);
                s.pArrayPos = obj.pArrayPos;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractArrayOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSensorArray = phased.internal.AbstractArray.loadobj(s.cSensorArray);
                    s = rmfield(s,'cSensorArray');
                end
                s = rmfield(s,'isLocked');
            end
        end        
        
        function loadObjectImpl(obj,s,wasLocked)
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end       
    end
    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractArrayOperation('array');
      end
    end

    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)  %#ok<MANU>
            varargout = {''};
        end
        
        function varargout = getOutputNamesImpl(obj) %#ok<MANU>
            varargout = {''};
        end
        
        function varargout = getOutputSizeImpl(obj)
            szAng = propagatedInputSize(obj,1);
            numDOF = getDOF(obj.SensorArray);
            varargout{1} = [numDOF szAng(2)];
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj,1);
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout{1} = 'double';
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = false;
        end
    end        
end


