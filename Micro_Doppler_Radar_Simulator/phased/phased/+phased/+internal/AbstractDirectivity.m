classdef (Hidden) AbstractDirectivity < phased.internal.AbstractSensorOperation
%This class is for internal use only. It may be removed in the future.

%AbstractDirectivity   Base computation needed for directivity 

%   Copyright 2013-2017 The MathWorks, Inc.

%   Reference
%   [1] Constantine A. Balanis, Antenna Theory, 3rd ed. Wiley-Interscience,
%       2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable, Logical) 
        %EnablePolarization  Enable polarization simulation
        %   Set this property to true to enable polarization. Set this
        %   property to false to ignore polarization. The default value of
        %   this property is false. This property applies when the array
        %   specified in the Sensor property is capable of simulating
        %   polarization.
        EnablePolarization = false
    end
    
    properties (Access = private, Nontunable)
        cArrayResponse
        cIntegratedPattern
        pDOF
    end
    
    methods (Access = protected, Abstract)
        resp = extractResponse(obj,respIn) 
        % extract response from computed responseIn depending on whether
        % polarization is enabled, may need to combine polarization
        % components
        g = extractOutput(obj,resp,intresp)
        % compute either directivity or directivity response based on the
        % response and integrated response
    end

    methods (Access = protected)

        function obj = AbstractDirectivity(varargin)
            obj@phased.internal.AbstractSensorOperation(varargin{:});
        end

    end
    
    methods (Access = protected)
        
        function validatePropertiesImpl(obj)
            cond = obj.EnablePolarization && ~isPolarizationCapable(obj.Sensor);
            if cond
                coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
            end
        end
        
        function setupImpl(obj,~,~, ~, ~) 
            setupImpl@phased.internal.AbstractSensorOperation(obj)
            
            if obj.pUseArray
                obj.pDOF = getDOF(obj.Sensor);
                if obj.WeightsInputPort
                    if isPolarizationCapable(obj.Sensor)
                        if obj.EnablePolarization
                            obj.cArrayResponse = phased.ArrayResponse(...
                                'SensorArray',obj.Sensor,...
                                'PropagationSpeed',obj.PropagationSpeed,...
                                'WeightsInputPort',true,...
                                'EnablePolarization',true);
                        else
                            obj.cArrayResponse = phased.ArrayResponse(...
                                'SensorArray',obj.Sensor,...
                                'PropagationSpeed',obj.PropagationSpeed,...
                                'WeightsInputPort',true,...
                                'EnablePolarization',false);
                        end
                    else
                        obj.cArrayResponse = phased.ArrayResponse(...
                            'SensorArray',obj.Sensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',true,...
                            'EnablePolarization',false);
                    end
                else
                    if isPolarizationCapable(obj.Sensor)
                        if obj.EnablePolarization
                            obj.cArrayResponse = phased.ArrayResponse(...
                                'SensorArray',obj.Sensor,...
                                'PropagationSpeed',obj.PropagationSpeed,...
                                'WeightsInputPort',false,...
                                'EnablePolarization',true);
                        else
                            obj.cArrayResponse = phased.ArrayResponse(...
                                'SensorArray',obj.Sensor,...
                                'PropagationSpeed',obj.PropagationSpeed,...
                                'WeightsInputPort',false,...
                                'EnablePolarization',false);
                        end
                    else
                        obj.cArrayResponse = phased.ArrayResponse(...
                            'SensorArray',obj.Sensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',false,...
                            'EnablePolarization',false);
                    end
                end
                
                sensorElem = getElementHandle(obj.Sensor);
                if isElementFromAntenna(obj.Sensor) || ...
                        isa(sensorElem,'phased.internal.AntennaAdapter')
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                        'Sensor',obj.Sensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',obj.WeightsInputPort,...
                        'EnablePolarization',obj.EnablePolarization);
                else
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                        'Sensor',obj.Sensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',obj.WeightsInputPort,...
                        'EnablePolarization',obj.EnablePolarization);
                end
            else
                if isempty(coder.target)
                    obj.cArrayResponse = clone(obj.cSensor);
                else
                    obj.cArrayResponse = clonecg(obj.cSensor);
                end
                obj.pDOF = 1;
                if isElementFromAntenna(obj.Sensor) || ...
                        isa(obj.Sensor,'phased.internal.AntennaAdapter') 
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPattern(...
                        'Sensor',obj.Sensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',false,...
                        'EnablePolarization',obj.EnablePolarization);
                else
                    obj.cIntegratedPattern = phased.internal.IntegratedPowerPatternReference(...
                        'Sensor',obj.Sensor,...
                        'PropagationSpeed',obj.PropagationSpeed,...
                        'WeightsInputPort',false,...
                        'EnablePolarization',obj.EnablePolarization);
                end
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index) 
            flag = true;
            if obj.WeightsInputPort && (index == 3)
                flag = false;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(~,~) 
            flag = true;
        end
        
        function releaseImpl(obj)
            release(obj.cArrayResponse);
            if ~isa(obj.cIntegratedPattern,'double')
                release(obj.cIntegratedPattern);
            end
        end
        
        function resetImpl(obj)
            reset(obj.cArrayResponse);
            if ~isa(obj.cIntegratedPattern,'double')
                reset(obj.cIntegratedPattern);
            end
        end

        function validateInputsImpl(obj, freq, angle, w, stang)
            cond = ~isrow(freq) || isempty(freq);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeRowVector','Freq');
            end
            cond = ~isa(freq,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Freq','double');
            end
            cond = ~isreal(freq);
            if cond
                coder.internal.errorIf(cond,'phased:measure:NeedReal', 'Freq');
            end
            
            cond = ~ismatrix(angle) || isempty(angle);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','Ang');
            end
            sz_angle = size(angle);
            cond = sz_angle(1) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:system:measure:ArrayGain:NeedTwoRows','Ang');
            end
            cond = ~isa(angle,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Ang','double');
            end
            cond = ~isreal(angle);
            if cond
                coder.internal.errorIf(cond,'phased:measure:NeedReal', 'Ang');
            end
            
            if obj.WeightsInputPort
                cond = ~ismatrix(w) || isempty(w);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeMatrix','W');
                end
                sz_w = size(w);
                cond = ~isa(w,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','W','double');
                end
                cond = (size(freq,2)~=sz_w(2)) && (sz_w(2)~=1);
                if cond
                    coder.internal.errorIf(cond,'phased:measure:MatrixVectorDimensionMismatch','W','Freq');
                end
                if ~isa(obj.Sensor,'phased.internal.AbstractElement')
                    N = getDOF(obj.Sensor);
                else
                    N = 1;
                end
                cond = sz_w(1)~=N;
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidRowNumbers','W',N);
                end
            end
            
            if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    ~strncmp(obj.Sensor.SubarraySteering,'None',1)
                if ~obj.WeightsInputPort
                    stang = w;
                end
                if strncmp(obj.Sensor.SubarraySteering,'Phase',1) || ...
                        strncmp(obj.Sensor.SubarraySteering,'Time',1)
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
                else
                    % does not support multiple weights for multiple
                    % frequency yet because this is still an analog
                    % behavior so at any moment, there is only one set of
                    % weights can be applied.
                    
                    ws = stang;  % weights
                    Ns = getNumSubarrays(obj.Sensor);
                    cond = (~iscell(ws) && ~ismatrix(ws)) || isempty(ws);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:phased:expectedCellOrMatrix','WS');
                    end
                    Nse = zeros(1,Ns);
                    for m = 1:Ns
                        Nse(m) = getNumElements(obj.Sensor,m);
                    end
                    if iscell(ws)
                        cond = ~isrow(ws) || (numel(ws)~= Ns);
                        if cond
                            coder.internal.errorIf(cond, ...
                                'phased:phased:expectedMatrixSize','WS',1,Ns);
                        end
                        for m = 1:Ns
                            cond = ~iscolumn(ws{m}) || (numel(ws{m})~=Nse(m));
                            if cond
                                coder.internal.errorIf(cond, ...
                                    'phased:system:array:SubarrayElementWeightsSizeMismatch',...
                                    m,'WS',Nse(m));
                            end
                            cond = ~isa(ws{m},'double');
                            if cond
                                coder.internal.errorIf(cond, ...
                                    'phased:system:array:SubarrayElementWeightsInvalidDataType',...
                                    m,'WS','double');
                            end
                        end
                    else
                        sz_ws = size(ws);
                        Nsemax = max(Nse);
                        cond = ~isequal(sz_ws,[Nsemax Ns]);
                        if cond
                            coder.internal.errorIf(cond, ...
                                'phased:phased:expectedMatrixSize','WS',Nsemax,Ns);
                        end
                        cond = ~isa(ws,'double');
                        if cond
                            coder.internal.errorIf(cond, ...
                                'MATLAB:system:invalidInputDataType','WS','double');
                        end
                    end
                end
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSensorOperation(obj);
            if isLocked(obj)
                s.cArrayResponse = saveobj(obj.cArrayResponse);
                s.cIntegratedPattern = saveobj(obj.cIntegratedPattern);
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSensorOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cArrayResponse = phased.ArrayResponse.loadobj(s.cArrayResponse);
                    s = rmfield(s,'cArrayResponse');
                    obj.cIntegratedPattern = eval(...
                        sprintf('%s.loadobj(s.cIntegratedPattern)',s.cIntegratedPattern.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cIntegratedPattern');
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
        
        function g = stepImpl(obj,freq,ang, weights, stang)

            if ~obj.WeightsInputPort
                if obj.pNeedSteeringAngle
                    stang = weights;
                end
                weights = ones(obj.pDOF,1);
            end
            
            if obj.pUseArray
                if obj.pNeedSteeringAngle
                    if ~obj.WeightsInputPort
                        resp_temp = step(obj.cArrayResponse,freq,ang,stang);
                        intresp = step(obj.cIntegratedPattern,freq,stang);
                    else
                        resp_temp = step(obj.cArrayResponse,freq,ang,weights,stang);
                        intresp = step(obj.cIntegratedPattern,freq,weights,stang);
                    end
                else
                    if ~obj.WeightsInputPort
                        resp_temp = step(obj.cArrayResponse,freq,ang);
                        intresp = step(obj.cIntegratedPattern,freq);
                    else
                        resp_temp = step(obj.cArrayResponse,freq,ang,weights);
                        intresp = step(obj.cIntegratedPattern,freq,weights);
                    end
                end
                
                resp = extractResponse(obj,resp_temp);
                
                g = extractOutput(obj,resp,intresp);
            elseif ~isElementFromAntenna(obj.Sensor)
                resp_temp = step(obj.cArrayResponse,freq,ang);  % no weights because it cancels
                resp = extractResponse(obj,resp_temp);
                intresp = step(obj.cIntegratedPattern,freq);
                             
                g = extractOutput(obj,resp,intresp);
            else
                g = zeros(size(ang,2),numel(freq));
                if size(ang,1) == 1
                    angM = [ang;zeros(size(ang))];
                else
                    angM = ang;
                end
                for m = 1:numel(freq)
                    d_temp = pattern(obj.Sensor,freq(m),angM(1,:),angM(2,:));
                    g(:,m) = diag(d_temp);
                end
            end
            
        end
        
    end
    
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            props = {'Sensor','PropagationSpeed',...
                     'WeightsInputPort','EnablePolarization'};
            groups = matlab.system.display.Section('Title', 'Parameters', ...
                                                   'PropertyList', props);
        end
    end

end

