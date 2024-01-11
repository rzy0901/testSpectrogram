classdef (Sealed, StrictDefaults) Radiator < phased.internal.AbstractSensorOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%Radiator Narrowband signal radiator
%   H = phased.Radiator creates a narrowband signal radiator System object,
%   H. The object returns radiated narrowband signals for given directions
%   using a sensor array or a single element.
%
%   H = phased.Radiator(Name,Value) creates a radiator object, H, with the
%   specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax when the CombineRadiatedSignals property is true:
%
%   Y = step(H,X,ANG) radiates signal X to multiple directions
%   specified in ANG, and combines the radiated signal for each direction
%   in Y. The combination is performed using the phase approximation of
%   time delays across radiating elements at the far field.
%
%   X can be either a vector or a matrix. If X is a vector, it is
%   radiated through all elements. If X is a matrix, its column number
%   must equal the number of subarrays if Sensor contains subarrays, or the
%   number of elements otherwise. Each column of X is radiated by
%   corresponding element/subarray.
%
%   ANG is a 2-row matrix. Each column of ANG specifies a radiating
%   direction in the form of an [AzimuthAngle; ElevationAngle] pair (in
%   degrees). Y is a matrix whose number of columns equals the number of
%   radiating directions. Each column of Y contains the combined far field
%   output from all elements/subarrays to the corresponding direction.
%
%   Y = step(H,X,ANG,LAXES) specifies the radiator's local coordinate
%   system in LAXES when you set the EnablePolarization property to true.
%   LAXES is a 3x3 matrix whose columns specify the local coordinate
%   system's orthonormal x, y, and z axes, respectively. Each axis is
%   specified in [x;y;z] form measured in the global coordinate system. Y
%   is a 1xM struct array where M is the number of entries in ANG. Each
%   struct in the array has three fields: X, Y, and Z. Each field contains
%   the X, Y, and Z component, also measured in the global coordinate
%   system, of the polarized far field output, respectively. Within each
%   field is a column vector containing the combined far field output from
%   all elements/subarrays to the corresponding direction.
%
%   Y = step(H,X,ANG,W) uses W as the weight vector when
%   the WeightsInputPort property is true. W must be a column vector
%   whose length is the same as the number of radiating elements or
%   subarrays.
%
%   Y = step(H,X,ANG,STEER) uses STEER as the subarray
%   steering angle (in degrees). STEER must be a length-2 column
%   vector in the form of [AzimuthAngle; ElevationAngle]. This syntax is
%   only applicable when you use subarrays in the Sensor property and set
%   the SubarraySteering property in the Sensor to either 'Phase' or
%   'Time'.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANG,LAXES,W,STEER)
%
%   Step method syntax when the CombineRadiatedSignals property is false:
%
%   Y = step(H,X,ANG) radiates the input signal to different angles for
%   different radiating elements. The number of columns of ANG is the
%   same as the number of radiating elements and each column of ANG
%   specifies the radiating direction for the corresponding input signal
%   through the corresponding radiating element. Note that this syntax does
%   not allow the Sensor to contain subarrays. Each column of Y represents
%   the radiated signal from the corresponding radiating element.
%
%   Y = step(H,X,ANG,LAXES) specifies the radiator's local coordinate
%   system in LAXES when you set the EnablePolarization property to true.
%   LAXES is a 3x3 matrix whose columns specify the local coordinate
%   system's orthonormal x, y, and z axes, respectively. Each axis is
%   specified in [x;y;z] form measured in the global coordinate system.  Y
%   is a 1xN struct array where N is the number of elements/subarrays in
%   the sensor array. Each struct in the array has three fields: X, Y, and
%   Z. Each field contains the X, Y, and Z component of the polarized far
%   field output, also measured in the global coordinate system,
%   respectively. Within each field is a column vector representing the
%   radiated signal from the corresponding radiating element.
%  
%   Y = step(H,X,ANG,W) uses W as the weight vector when
%   the WeightsInputPort property is true. W must be a column vector
%   whose length is the same as the number of radiating elements.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   Y = step(H,X,ANG,LAXES,W)
%
%   Radiator methods:
%
%   step     - Radiate signals (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a radiator object with same property values
%   isLocked - Locked status (logical)
%
%   Radiator properties:
%
%   Sensor                 - Handle of the sensor
%   PropagationSpeed       - Propagation speed
%   OperatingFrequency     - Operating frequency
%   SignalType             - Signal type
%   CombineRadiatedSignals - Combine radiated signals
%   EnablePolarization     - Enable polarization
%   WeightsInputPort       - Enable weights input
%
%   % Examples:
%
%   % Example 1: 
%   %   Radiate signal with a single antenna.
%
%   ha = phased.IsotropicAntennaElement;
%   hr = phased.Radiator('Sensor',ha,'OperatingFrequency',300e6);
%   x = [1;1];
%   radiatingAngle = [30 10]';
%   y = step(hr,x,radiatingAngle);
%
%   % Example 2: 
%   %   Radiate a far field signal with a 5-element array.
%
%   ha = phased.ULA('NumElements',5);
%   hr = phased.Radiator('Sensor',ha,'OperatingFrequency',300e6);
%   x = [1;1];
%   radiatingAngle = [30 10; 20 0]'; % two directions
%   y = step(hr,x,radiatingAngle);
%
%   % Example 3: 
%   %   Radiate signal with a 3-element antenna array. Each antenna 
%   %   radiates a separate signal to a separate direction.
%
%   ha = phased.ULA('NumElements',3);
%   hr = phased.Radiator('Sensor',ha,'OperatingFrequency',1e9,...
%                        'CombineRadiatedSignals',false);
%   x = [1 2 3;1 2 3];
%   radiatingAngle = [10 0; 20 5; 45 2]'; % One angle for one antenna
%   y = step(hr,x,radiatingAngle);
%
%   See also phased, phased.Collector.

%   Copyright 2009-2013 The MathWorks, Inc.
%     

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
properties (Nontunable)
    %OperatingFrequency     Operating frequency (Hz)
    %   Specify the operating frequency (in Hz) as a positive scalar. The
    %   default value is 3e8.
    OperatingFrequency = 3e8;    
    %SignalType     Signal type
    %   Specify whether the signal should be treated as one of 'Field' |
    %   'Power', where 'Field' is the default. When you set SignalType to
    %   'Field', the radiated signal is scaled by the field at the
    %   corresponding direction. When you set SignalType to 'Power', the
    %   radiated signal is scaled by the directivity at the corresponding
    %   direction.
    SignalType = 'Field';
end

properties (Nontunable, Logical) 
    %CombineRadiatedSignals Combine radiated signals
    %   Set this property to true to combine radiated signals from all
    %   radiating elements. Set this property to false to obtain the
    %   radiated signal for each radiating element. The default value is
    %   true. CombineRadiatedSignals must be true if you set SignalType to
    %   'Power'.
    CombineRadiatedSignals = true;
    %EnablePolarization  Enable polarization
    %   Set this property to true to enable polarization. Set this property
    %   to false to ignore polarization. The default value of this property
    %   is false. This property applies when the sensor specified in the
    %   Sensor property is capable of simulating polarization.
    EnablePolarization = false;
end

properties(Constant, Hidden)
    SignalTypeSet = matlab.system.StringSet({'Field','Power'});
end

properties (Access = private, Nontunable)
    cSteeringVector;
    cIntegratedPattern;
    pDOF;
    pCommonSignalForChannels;
    pPowerDistributionFactor;
end

properties (Access = private, Logical, Nontunable)
    pIsPowerSignal
end

methods
    function set.OperatingFrequency(obj,value)
        validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'OperatingFrequency');
        obj.OperatingFrequency = value;
    end    
end

methods
    function obj = Radiator(varargin)
        obj@phased.internal.AbstractSensorOperation(varargin{:});
    end
end
    
methods (Access = protected)
      
    function num = getNumInputsImpl(obj)
        num = getNumInputsImpl@phased.internal.AbstractSensorOperation(obj);
        if obj.EnablePolarization
            num = num+1;
        end
    end
    
    function validatePropertiesImpl(obj)
        cond = isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~obj.CombineRadiatedSignals;
        if cond
            coder.internal.errorIf(cond,'phased:phased:collector:invalidSettingsForSubarray','CombineRadiatedSignals','true','Sensor');
        end
        
        cond = strcmp(obj.SignalType,'Power') && ~obj.CombineRadiatedSignals;
        if cond
            coder.internal.errorIf(cond,'phased:system:array:IrrelevantSetting','Power','CombineRadiatedSignals','false');
        end

        cond = obj.EnablePolarization && ~isPolarizationCapable(obj.Sensor);
        if cond
            coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
        end
    end
    
    function validateInputsImpl(obj,x,ang,lclaxes,wArg,stangArg)
        coder.extrinsic('mat2str');
        % input signal
        cond = ~isa(x,'double');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','X','double');
        end
        cond = ~ismatrix(x) || isempty(x);
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:inputMustBeMatrix','X');
        end
            
        xsize = size(x);
        UseArrayFlag = isa(obj.Sensor,'phased.internal.AbstractArray') || ....
            isa(obj.Sensor,'phased.internal.AbstractSubarray');
        if UseArrayFlag
            M = getDOF(obj.Sensor);
        else
            M = 1;
        end
       
        cond = xsize(2)>1 && xsize(2)~=M;
        if cond
            coder.internal.errorIf(cond,'phased:phased:Radiator:DimMismatch');
        end
        
        % For near field, the angle must match the channels
        angsize = size(ang);
        cond = UseArrayFlag && ~obj.CombineRadiatedSignals && angsize(2)~=M;
        if cond
            coder.internal.errorIf(cond,'phased:phased:Radiator:InvalidAngle','CombineRadiatedSignals');
        end
        cond = ~isa(ang,'double');
        if cond
            coder.internal.errorIf(cond, ...
                'MATLAB:system:invalidInputDataType','Ang','double');
        end
        cond = ~isreal(ang);
        if cond
            coder.internal.errorIf(cond,'phased:step:NeedReal', 'Ang');
        end
        
        if obj.EnablePolarization
            % check lclaxes to be 3x3  
            sigdatatypes.validate3DCartCoord(lclaxes,'','LAxes',...
                {'size',[3 3]});
        end
        
        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmp(obj.Sensor.SubarraySteering,'None',1)
            if ~obj.WeightsInputPort
                if obj.EnablePolarization
                    stang = wArg;
                else
                    stang = lclaxes;
                end
            else
                if ~obj.EnablePolarization
                    stang = wArg;
                else
                    stang = stangArg;
                end
            end
            sz_stang = size(stang);
            cond = ~ismatrix(stang) || isempty(stang);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','Steer');
            end
            cond = sz_stang(1) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:NeedTwoRows','Steer');
            end
            cond = sz_stang(2) > 1;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:NeedOneColumn','Steer');
            end
            cond = ~isreal(stang);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidAngle','Steer');
            end
            cond = ~isa(stang,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Steer','double');
            end

        end
        % weights
        if obj.WeightsInputPort
            if ~obj.EnablePolarization
                w = lclaxes;
            else
                w = wArg;
            end
            wsize = size(w);
            cond = ~isequal(wsize, [M 1]);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDimensions','W',...
                    coder.const(mat2str([M 1])), ...
                    coder.const(mat2str(wsize)));
            end
            cond = ~isa(w,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','W','double');
            end
        end
            
            
    end
    
    function setupImpl(obj,x,~,~,~,~)
        setupImpl@phased.internal.AbstractSensorOperation(obj);
        if obj.pUseArray
            M = getDOF(obj.cSensor);
        else
            M = 1;            
        end
        obj.pDOF = M;
        xsize = size(x);
        if obj.pUseArray && (xsize(2) == 1)
            obj.pCommonSignalForChannels = true;
        else
            obj.pCommonSignalForChannels = false;
        end
        obj.pIsPowerSignal = strcmp(obj.SignalType,'Power');
        if obj.pUseArray 
            if obj.CombineRadiatedSignals
                obj.cSteeringVector = phased.SteeringVector(...
                    'SensorArray',obj.cSensor,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'IncludeElementResponse',true,...
                    'EnablePolarization',obj.EnablePolarization);
                if obj.pIsPowerSignal
                    sensorElem = getElementHandle(obj.cSensor);
                    if isElementFromAntenna(obj.cSensor) || ...
                            isa(sensorElem,'phased.internal.AntennaAdapter')
                        obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                            'Sensor',obj.cSensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',obj.WeightsInputPort,...
                            'EnablePolarization',obj.EnablePolarization);
                    else
                        obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                            'Sensor',obj.cSensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',obj.WeightsInputPort,...
                            'EnablePolarization',obj.EnablePolarization);
                    end
                end
            end
        elseif obj.pIsPowerSignal
            if isElementFromAntenna(obj.cSensor) || ...
                    isa(obj.cSensor,'phased.internal.AntennaAdapter')
                obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                    'Sensor',obj.cSensor,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'WeightsInputPort',false,...  % single element weights is just a scalar, separable
                    'EnablePolarization',obj.EnablePolarization);
            else
                obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                    'Sensor',obj.cSensor,...
                    'PropagationSpeed',obj.PropagationSpeed,...
                    'WeightsInputPort',false,...  % single element weights is just a scalar, separable
                    'EnablePolarization',obj.EnablePolarization);
            end
        end
            
        if isa(obj.Sensor,'phased.internal.AbstractSubarray')
            obj.pPowerDistributionFactor = getNumElements(obj.cSensor,...
                1:getNumSubarrays(obj.cSensor));
        else
            obj.pPowerDistributionFactor = 1;
        end
    end
    
    function flag = isInputComplexityLockedImpl(obj,index)
        flag = false;  % index == 1
        if index == 2
            flag = true;
        end
        if obj.WeightsInputPort && (index == 3)
            flag = false;
        end
    end

    function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
        flag = false;
    end

    function resetImpl(obj)
        if obj.pUseArray
            if obj.CombineRadiatedSignals
                reset(obj.cSteeringVector);
                if obj.pIsPowerSignal
                    reset(obj.cIntegratedPattern);
                end
            else
                reset(obj.cSensor);
            end
        else
            reset(obj.cSensor);
            if obj.pIsPowerSignal
                reset(obj.cIntegratedPattern);
            end
        end
    end
    
    function releaseImpl(obj)
        if obj.pUseArray
            if obj.CombineRadiatedSignals
                release(obj.cSteeringVector);
                if obj.pIsPowerSignal
                    release(obj.cIntegratedPattern);
                end
            else
                release(obj.cSensor)
            end
        else
            release(obj.cSensor)
            if obj.pIsPowerSignal
                release(obj.cIntegratedPattern);
            end
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractSensorOperation(obj);
        if isLocked(obj)
            s.pDOF = obj.pDOF;
            s.cSteeringVector = saveobj(obj.cSteeringVector);
            s.cIntegratedPattern = saveobj(obj.cIntegratedPattern);
            s.pCommonSignalForChannels = obj.pCommonSignalForChannels;
            s.pPowerDistributionFactor = obj.pPowerDistributionFactor;
            s.pIsPowerSignal = obj.pIsPowerSignal;
        end
    end

    function s = loadSubObjects(obj,s)
        s = loadSubObjects@phased.internal.AbstractSensorOperation(obj,s);
        if isfield(s,'isLocked')
            if s.isLocked
                obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                s = rmfield(s,'cSteeringVector');
                if isfield(s,'cIntegratedPattern') && ...
                            isfield(s.cIntegratedPattern,'ClassNameForLoadTimeEval')
                    obj.cIntegratedPattern = eval(...
                        sprintf('%s.loadobj(s.cIntegratedPattern)',s.cIntegratedPattern.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cIntegratedPattern');
                end
            end
            s = rmfield(s,'isLocked');
        end
    end
    
    function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
        s = loadSubObjects(obj,s);
        if isfield(s,'pNumSensorArrayElements')
            s = rmfield(s,'pNumSensorArrayElements');
        end
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end

    function y = stepImpl(obj,x,ang,lclaxes,wArg,stangArg)
        
        fc = obj.OperatingFrequency;
        if obj.EnablePolarization
            if obj.pUseArray
                if obj.CombineRadiatedSignals
                    if obj.pNeedSteeringAngle
                        if ~obj.WeightsInputPort
                            stang = wArg;
                            if obj.pIsPowerSignal
                                intresp = step(obj.cIntegratedPattern,fc,stang);
                            end
                        else
                            stang = stangArg;
                            w = wArg;
                            if obj.pIsPowerSignal
                                intresp = step(obj.cIntegratedPattern,fc,w,stang);
                            end
                        end
                        sv = step(obj.cSteeringVector,fc,ang,stang);
                    else
                        sv = step(obj.cSteeringVector,fc,ang);
                        if obj.WeightsInputPort
                            w = wArg;
                            if obj.pIsPowerSignal
                                intresp = step(obj.cIntegratedPattern,fc,w);
                            end
                        elseif obj.pIsPowerSignal
                            intresp = step(obj.cIntegratedPattern,fc);
                        end
                    end
                    pfactor = 1./sqrt(obj.pPowerDistributionFactor(:));
                    sv_h = bsxfun(@times,pfactor,sv.H);
                    sv_v = bsxfun(@times,pfactor,sv.V);
                    if obj.WeightsInputPort
                        for m = 1:obj.pDOF
                            sv_h(m,:) = w(m)*sv_h(m,:);
                            sv_v(m,:) = w(m)*sv_v(m,:);
                        end
                    end
                    if obj.pCommonSignalForChannels 
                        svc_h = sum(sv_h,1);
                        svc_v = sum(sv_v,1);
                    else
                        svc_h =  sv_h;
                        svc_v =  sv_v;
                    end
                    if obj.pIsPowerSignal
                        svc_h = phased.internal.normalizeIntegratedPower(svc_h,intresp.H,false);
                        svc_v = phased.internal.normalizeIntegratedPower(svc_v,intresp.V,false);
                    end
                    y_h = x*svc_h;
                    y_v = x*svc_v;
                else
                    g = step(obj.cSensor,fc,ang);
                    g_h = diag(g.H);
                    g_v = diag(g.V);
                    if obj.WeightsInputPort
                        w = wArg;
                        g_h = w.*g_h;
                        g_v = w.*g_v;
                    end
                    if obj.pCommonSignalForChannels
                        y_h = x*g_h.';
                        y_v = x*g_v.';
                    else
                        y_h = complex(x);
                        y_v = complex(x);
                        for m = 1:obj.pDOF
                            y_h(:,m) = y_h(:,m)*g_h(m);
                            y_v(:,m) = y_v(:,m)*g_v(m);
                        end
                    end
                end
            else
                % Single sensor
                g = step(obj.cSensor,fc,ang);
                g_h = g.H;
                g_v = g.V;
                if obj.WeightsInputPort
                    w = wArg;
                    g_h = w*g_h;
                    g_v = w*g_v;
                end
                if obj.pIsPowerSignal
                    intresp = step(obj.cIntegratedPattern,fc);
                    if obj.WeightsInputPort
                        intresp.H = intresp.H*abs(w)^2;
                        intresp.H = intresp.V*abs(w)^2;
                    end
                    g_h = phased.internal.normalizeIntegratedPower(g_h,intresp.H);
                    g_v = phased.internal.normalizeIntegratedPower(g_v,intresp.V);
                end
                y_h = x*g_h.';
                y_v = x*g_v.';
            end
            initYfields = complex(zeros(size(y_h,1),1));
            y = repmat(struct('X',initYfields, 'Y', initYfields,'Z', initYfields), ...
                1, size(ang,2));
                
            for m = size(ang,2):-1:1
                y_temp = sph2cartvec(...
                    [y_h(:,m).';y_v(:,m).';zeros(1,size(x,1))],...
                    ang(1,m),ang(2,m));
                y_temp = phased.internal.local2globalvec(...
                    y_temp,lclaxes);
                y(m).X = y_temp(1,:).';
                y(m).Y = y_temp(2,:).';
                y(m).Z = y_temp(3,:).';
            end
            
        else  % no polarization
            if obj.pUseArray
                if obj.CombineRadiatedSignals
                    if obj.pNeedSteeringAngle
                        if ~obj.WeightsInputPort
                            stang = lclaxes;
                            if obj.pIsPowerSignal
                                intresp = step(obj.cIntegratedPattern,fc,stang);
                            end
                        else
                            stang = wArg;
                            w = lclaxes;
                            if obj.pIsPowerSignal
                                intresp = step(obj.cIntegratedPattern,fc,w,stang);
                            end
                        end
                        sv = step(obj.cSteeringVector,fc,ang,stang);
                    else
                        sv = step(obj.cSteeringVector,fc,ang);
                        if obj.WeightsInputPort
                            w = lclaxes;
                            if obj.pIsPowerSignal
                                intresp = step(obj.cIntegratedPattern,fc,w);
                            end
                        else
                            if obj.pIsPowerSignal
                                intresp = step(obj.cIntegratedPattern,fc);
                            end
                        end
                    end
                    sv = bsxfun(@times,1./sqrt(obj.pPowerDistributionFactor(:)),sv);
                    if obj.WeightsInputPort
                        for m = 1:obj.pDOF
                            sv(m,:) = w(m)*sv(m,:);
                        end
                    end
                    if obj.pIsPowerSignal
                        sv = phased.internal.normalizeIntegratedPower(sv,intresp,false);
                    end
                    if obj.pCommonSignalForChannels 
                       y = x*sum(sv,1);
                    else
                       y = x*sv; 
                    end
                    
                else
                    g_temp = step(obj.cSensor,fc,ang);
                    if isPolarizationEnabled(obj.cSensor)
                        g = diag(hypot(g_temp.H,g_temp.V));
                    else
                        g = diag(g_temp);
                    end
                    if obj.WeightsInputPort
                        w = lclaxes;
                        g = w.*g;
                    end
                    if obj.pCommonSignalForChannels
                        y = x*g.';
                    else
                        y = complex(zeros(size(x)));
                        for m = 1:obj.pDOF
                            y(:,m) = x(:,m)*g(m);
                        end
                    end
                end
            else
                % Single sensor
                g_temp = step(obj.cSensor,fc,ang);
                if isPolarizationEnabled(obj.cSensor)
                    g = hypot(g_temp.H,g_temp.V);
                else
                    g = g_temp;
                end
                if obj.WeightsInputPort
                    w = lclaxes;
                    g = w*g;
                end
                if obj.pIsPowerSignal
                    intresp = step(obj.cIntegratedPattern,fc);
                    if obj.WeightsInputPort
                        intresp = intresp*abs(w)^2;
                    end
                    g = phased.internal.normalizeIntegratedPower(g,intresp,false);
                end
                y = x*g.';
            end
        end
        
    end
           
end

methods (Static,Hidden,Access=protected)
    function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractSensorOperation;
        
        pCombineRadiatedSignals = ...
            matlab.system.display.internal.Property('CombineRadiatedSignals', ...
            'IsGraphical', false);
        dEnablePolarization = ...
            matlab.system.display.internal.Property('EnablePolarization', ...
            'IsGraphical', false);
        
        props = {'OperatingFrequency',...
            'SignalType',...
            pCombineRadiatedSignals,...
            dEnablePolarization,...
            'WeightsInputPort'};
        groups(1).PropertyList = [groups(1).PropertyList props];
        action = matlab.system.display.Action(@phased.Radiator.onActionCalled, ...
            'Label', 'Analyze');
        matlab.system.display.internal.setCallbacks(action, ...
            'DialogAppliedFcn', @phased.Radiator.onDialogApplied, ...
            'SystemDeletedFcn', @phased.Radiator.onSystemDeleted);
        groups(2).Actions = action;
    end
    function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'Title',getString(message('phased:library:block:RadiatorTitle')),...
            'Text',getString(message('phased:library:block:RadiatorDesc')));
    end
end
methods (Static,Hidden)
    function onActionCalled(actionData, sysObj)
        if ~phased.apps.internal.SensorArrayViewer.SettingsDialog.verifySensorArray(sysObj.Sensor)
            return;
        end
        f = actionData.UserData;
        if isempty(f) || isStale(f);
            actionData.UserData = ...
                phased.apps.internal.SensorArrayViewer.SensorArrayViewer(sysObj);
        end
    end
    function onDialogApplied(actionData, sysObj)
 
        f = actionData.UserData;
        if ~isempty(f)
            if isStale(f)
                actionData.UserData = [];
            else
                if ~phased.apps.internal.SensorArrayViewer.SettingsDialog.verifySensorArray(sysObj.Sensor)
                    return;
                end
                updateSettingsWithSystemObject(f,sysObj);
            end
        end
    end
    function onSystemDeleted(actionData)
      f = actionData.UserData;
      if ~isempty(f) && ~isStale(f)
         close(f);
         actionData.UserData = [];
      end
   end
end 

methods (Access = protected) %for Simulink
    function varargout = getOutputNamesImpl(~)
        varargout = {''};
    end
    function varargout = getInputNamesImpl(obj)
        varargout = {'X','Ang','LAxes','W','Steer'};
        lastIdx = 3;
        if obj.EnablePolarization
            varargout(lastIdx) = {'Axes'};
            lastIdx = lastIdx+1;
        end
        if obj.WeightsInputPort
            varargout(lastIdx) = {'W'};
            lastIdx = lastIdx+1;
        end
        if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                ~strncmpi(obj.Sensor.SubarraySteering,'None',1);
            varargout(lastIdx) = {'Steer'};
        end
    end
    function str = getIconImpl(obj) %#ok<MANU>
        str = sprintf('Narrowband\n Tx Array');
    end        
    function varargout = getOutputSizeImpl(obj)
        szX = propagatedInputSize(obj,1);
        szAng = propagatedInputSize(obj,2);
        %Output is number of snapshots * number of angles
        varargout{1} = [szX(1) szAng(2)];
    end
    function varargout = isOutputFixedSizeImpl(obj)
        varargout{1} = propagatedInputFixedSize(obj, 1);
    end
    function varargout = getOutputDataTypeImpl(obj)
        varargout{1} = propagatedInputDataType(obj,1);
    end
    function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
        varargout{1} = true;
    end
end
end

