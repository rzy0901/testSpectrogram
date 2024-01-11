classdef (Hidden, StrictDefaults) IntegratedPowerPatternReference < phased.internal.AbstractSensorOperation
%This class is for internal use only. It may be removed in the future.

%AbstractDirectivity   Base computation needed for directivity 

%   Copyright 2014-2017 The MathWorks, Inc.

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
        cIntegratedPattern
        pDOF
        pNumFreqs
        pNumAzSteps = 360
        pNumElSteps = 180
    end
    
    properties (Access = private)
        pIntegratedResponse
        pComputedFreq
        pComputedWeights
        pComputedSteerAngle
        pSteerVecCov
        pLastRunNumFreqs = 0
    end
    
    properties (Access = private, Nontunable, Logical)
        pIsSingleWeights
    end
    
    methods 

        function obj = IntegratedPowerPatternReference(varargin)
            obj@phased.internal.AbstractSensorOperation(varargin{:});
        end

    end
    
    methods (Access = protected)
        
        function num = getNumInputsImpl(obj)
            num = 1;
            if obj.WeightsInputPort
                num = num+1;
            end
            if isa(obj.Sensor,'phased.internal.AbstractSubarray') && ...
                    ~strncmpi(obj.Sensor.SubarraySteering,'None',1)
                num = num+1;
            end
        end
    
        function validatePropertiesImpl(obj)
            cond = obj.EnablePolarization && ~isPolarizationCapable(obj.Sensor);
            if cond
                coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
            end
        end
        
        function setupImpl(obj,freq,weights, stang) 
            setupImpl@phased.internal.AbstractSensorOperation(obj)
            Nfreq = numel(freq);
            if obj.pUseArray
                obj.pDOF = getDOF(obj.Sensor);
                if obj.WeightsInputPort
                    if isPolarizationCapable(obj.Sensor)
                        obj.cIntegratedPattern = phased.SteeringVector(...
                            'SensorArray',obj.Sensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'IncludeElementResponse',true,...
                            'EnablePolarization',obj.EnablePolarization);
                    else
                        obj.cIntegratedPattern = phased.SteeringVector(...
                            'SensorArray',obj.Sensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'IncludeElementResponse',true,...
                            'EnablePolarization',false);
                    end
                else
                    if isPolarizationCapable(obj.Sensor)
                        obj.cIntegratedPattern = phased.ArrayResponse(...
                            'SensorArray',obj.Sensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',false,...
                            'EnablePolarization',obj.EnablePolarization);
                    else
                        obj.cIntegratedPattern = phased.ArrayResponse(...
                            'SensorArray',obj.Sensor,...
                            'PropagationSpeed',obj.PropagationSpeed,...
                            'WeightsInputPort',false,...
                            'EnablePolarization',false);
                    end
                end
                
            else
                obj.pDOF = 1;
            end
            obj.pNumFreqs = Nfreq;

            if ~obj.WeightsInputPort
                obj.pIsSingleWeights = true;
                obj.pComputedWeights = ones(obj.pDOF,1);
            else
                obj.pIsSingleWeights = (size(weights,2)==1);
                obj.pComputedWeights = weights;
            end
            
            if obj.pUseArray 
                % use radians to save computation
                Naz = obj.pNumAzSteps;
                Nel = obj.pNumElSteps;
                azstep = 2*pi/Naz;
                elstep = pi/Nel;
                azrad = -pi:azstep:pi-azstep;
                elrad = -pi/2:elstep:pi/2-elstep;
                
                obj.pComputedFreq = freq;
                
                if obj.WeightsInputPort
                    if obj.pNeedSteeringAngle
                        obj.pSteerVecCov = calcSteeringVectorCovariance(...
                            obj,freq,azrad,elrad,stang);
                    else
                        obj.pSteerVecCov = calcSteeringVectorCovariance(...
                            obj,freq,azrad,elrad);
                    end
                
                    if isPolarizationCapable(obj.Sensor) && ...
                            obj.EnablePolarization
                        obj.pIntegratedResponse = struct(...
                            'H',zeros(1,Nfreq),...
                            'V',zeros(1,Nfreq));
                    else
                        obj.pIntegratedResponse = zeros(1,Nfreq);
                    end
                    
                    for m = Nfreq:-1:1
                        if obj.pIsSingleWeights
                            w = weights;
                        else
                            w = weights(:,m);
                        end
                        if isPolarizationCapable(obj.Sensor)
                            if obj.EnablePolarization
                                obj.pIntegratedResponse.H(m) = ...
                                    pi/Nel*2*pi/Naz*...
                                    real(w'*obj.pSteerVecCov.H(:,:,m)*w);
                                obj.pIntegratedResponse.V(m) = ...
                                    pi/Nel*2*pi/Naz*...
                                    real(w'*obj.pSteerVecCov.V(:,:,m)*w);
                            else
                                obj.pIntegratedResponse(m) = ...
                                    pi/Nel*2*pi/Naz*...
                                    real(w'*obj.pSteerVecCov(:,:,m)*w);
                            end
                        else
                            obj.pIntegratedResponse(m) = ...
                                pi/Nel*2*pi/Naz*...
                                real(w'*obj.pSteerVecCov(:,:,m)*w);
                        end
                    end
                        
                else
                    
                    if obj.pNeedSteeringAngle
                        stang = weights;
                        obj.pSteerVecCov = calcRespPattern(...
                            obj,freq,azrad,elrad,stang);
                    else
                        obj.pSteerVecCov = calcRespPattern(...
                            obj,freq,azrad,elrad);
                    end
                
                    if isPolarizationCapable(obj.Sensor)
                        if obj.EnablePolarization
                            obj.pIntegratedResponse.H = ...
                                phased.internal.integratePattern(...
                                obj.pSteerVecCov.H,...
                                elrad(:),azstep,elstep);
                            obj.pIntegratedResponse.V = ...
                                phased.internal.integratePattern(...
                                obj.pSteerVecCov.V,...
                                elrad(:),azstep,elstep);
                        else
                            obj.pIntegratedResponse = ...
                                phased.internal.integratePattern(...
                                obj.pSteerVecCov,...
                                elrad(:),azstep,elstep);
                        end
                    else
                        obj.pIntegratedResponse = ...
                            phased.internal.integratePattern(...
                            obj.pSteerVecCov,elrad(:),azstep,elstep);
                    end
                    
                end
                
            else
                angstep = [1 1];
                azstep = angstep(1);
                elstep = angstep(2);
                azang = -180:azstep:180-azstep;
                elang = -90:elstep:90-elstep;
                
                [resp,obj.pComputedFreq] = ...
                    getPowerPattern(obj.Sensor,azang,elang); 
                if isPolarizationCapable(obj.Sensor)
                    if obj.EnablePolarization
                        obj.pIntegratedResponse.H = ...
                            phased.internal.integratePattern(...
                            resp.H,phased.internal.deg2rad(elang(:)),...
                            phased.internal.deg2rad(azstep),...
                            phased.internal.deg2rad(elstep));
                        obj.pIntegratedResponse.V = ...
                            phased.internal.integratePattern(...
                            resp.V,phased.internal.deg2rad(elang(:)),...
                            phased.internal.deg2rad(azstep),...
                            phased.internal.deg2rad(elstep));
                    else
                        obj.pIntegratedResponse = ...
                            phased.internal.integratePattern(...
                            resp.H+resp.V,...
                            phased.internal.deg2rad(elang(:)),...
                            phased.internal.deg2rad(azstep),...
                            phased.internal.deg2rad(elstep));
                    end
                else
                    obj.pIntegratedResponse = ...
                        phased.internal.integratePattern(...
                        resp,phased.internal.deg2rad(elang(:)),...
                        phased.internal.deg2rad(azstep),...
                        phased.internal.deg2rad(elstep));
                end
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index) 
            flag = true;
            if obj.WeightsInputPort && (index == 2)
                flag = false;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(~,~) 
            flag = true;
        end
        
        function releaseImpl(obj)
            if obj.pUseArray && ~isa(obj.cIntegratedPattern,'double')
                release(obj.cIntegratedPattern);
            end
            if coder.target('MATLAB')
                obj.pIntegratedResponse = [];
            end
        end
        
        function resetImpl(obj)
            if obj.pUseArray && ~isa(obj.cIntegratedPattern,'double')
                reset(obj.cIntegratedPattern);
            end
        end

        function validateInputsImpl(obj, freq, w, stang)
            cond = ~isrow(freq) || isempty(freq);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeRowVector','FREQ');
            end
            cond = ~isa(freq,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','FREQ','double');
            end
            cond = ~isreal(freq);
            if cond
                coder.internal.errorIf(cond,'phased:measure:NeedReal', 'FREQ');
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
                    coder.internal.errorIf(cond,'phased:measure:MatrixVectorDimensionMismatch','W','FREQ');
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
                            'MATLAB:system:inputMustBeMatrix','STEERANGLE');
                    end
                    cond = sz_stang(1) > 2;
                    if cond
                        coder.internal.errorIf(cond,'phased:system:array:NeedTwoRows','STEERANGLE');
                    end
                    cond = sz_stang(2) > 1;
                    if cond
                        coder.internal.errorIf(cond,'phased:system:array:NeedOneColumn','STEERANGLE');
                    end
                    cond = ~isreal(stang);
                    if cond
                        coder.internal.errorIf(cond,'phased:system:array:InvalidAngle','STEERANGLE');
                    end
                    cond = ~isa(stang,'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'MATLAB:system:invalidInputDataType','STEERANGLE','double');
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
                s.cIntegratedPattern = saveobj(obj.cIntegratedPattern);
                s.pDOF = obj.pDOF;
                s.pNumFreqs = obj.pNumFreqs;
                s.pIsSingleWeights = obj.pIsSingleWeights;
                s.pNumAzSteps = obj.pNumAzSteps;
                s.pNumElSteps = obj.pNumElSteps;
                s.pIntegratedResponse = obj.pIntegratedResponse;
                s.pComputedFreq = obj.pComputedFreq;
                s.pComputedWeights = obj.pComputedWeights;
                s.pComputedSteerAngle = obj.pComputedSteerAngle;
                s.pSteerVecCov = obj.pSteerVecCov;
                s.pLastRunNumFreqs = obj.pLastRunNumFreqs;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSensorOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    if s.WeightsInputPort
                        obj.cIntegratedPattern = phased.SteeringVector.loadobj(s.cIntegratedPattern);
                    else
                        obj.cIntegratedPattern = phased.ArrayResponse.loadobj(s.cIntegratedPattern);
                    end
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
        
        function g = stepImpl(obj,freq, weights, stang)

            if ~obj.WeightsInputPort
                if obj.pNeedSteeringAngle
                    stang = weights;
                end
                weights = ones(obj.pDOF,1);
            end
            
            if obj.pUseArray
                if ~obj.WeightsInputPort
                    
                    if obj.pNeedSteeringAngle
                        intPatChangeIdx = (freq~=obj.pComputedFreq) | ...
                            ~isequal(stang,obj.pComputedSteerAngle);
                        obj.pComputedFreq = freq;
                        obj.pComputedSteerAngle = stang;
                    else
                        intPatChangeIdx = (freq~=obj.pComputedFreq);
                        obj.pComputedFreq = freq;
                    end
                    
                    if any(intPatChangeIdx)
                        Naz = obj.pNumAzSteps;
                        Nel = obj.pNumElSteps;
                        azstep = 2*pi/Naz;
                        elstep = pi/Nel;
                        azrad = -pi:azstep:pi-azstep;
                        elrad = -pi/2:elstep:pi/2-elstep;
                        
                        if isPolarizationCapable(obj.Sensor)
                            if coder.target('MATLAB')
                                if obj.pNeedSteeringAngle
                                    temp = calcRespPattern(...
                                        obj,freq(intPatChangeIdx),azrad,elrad,stang);
                                else
                                    temp = calcRespPattern(...
                                        obj,freq(intPatChangeIdx),azrad,elrad);
                                end
                            else
                                if obj.pNeedSteeringAngle
                                    tempcg = calcRespPattern(...
                                        obj,freq,azrad,elrad,stang);
                                    if obj.EnablePolarization
                                        temp.H = tempcg.H(:,:,intPatChangeIdx);
                                        temp.V = tempcg.V(:,:,intPatChangeIdx);
                                    else
                                        temp = tempcg(:,:,intPatChangeIdx);
                                    end
                                else
                                    tempcg = calcRespPattern(...
                                        obj,freq(intPatChangeIdx),azrad,elrad);
                                    if obj.EnablePolarization
                                        temp.H = tempcg.H(:,:,intPatChangeIdx);
                                        temp.V = tempcg.V(:,:,intPatChangeIdx);
                                    else
                                        temp = tempcg(:,:,intPatChangeIdx);
                                    end
                                end
                            end
                            
                            if obj.EnablePolarization
                                obj.pSteerVecCov.H(:,:,intPatChangeIdx) = temp.H;
                                obj.pSteerVecCov.V(:,:,intPatChangeIdx) = temp.V;
                                obj.pIntegratedResponse.H(intPatChangeIdx) = ...
                                    phased.internal.integratePattern(...
                                    temp.H,elrad(:),azstep,elstep);
                                obj.pIntegratedResponse.V(intPatChangeIdx) = ...
                                    phased.internal.integratePattern(...
                                    temp.V,elrad(:),azstep,elstep);
                            else
                                obj.pSteerVecCov(:,:,intPatChangeIdx) = temp;
                                obj.pIntegratedResponse(intPatChangeIdx) = ...
                                    phased.internal.integratePattern(...
                                    temp,elrad(:),azstep,elstep);
                            end
                        else
                            if coder.target('MATLAB')
                                if obj.pNeedSteeringAngle
                                    temp = calcRespPattern(...
                                        obj,freq(intPatChangeIdx),azrad,elrad,stang);
                                else
                                    temp = calcRespPattern(...
                                        obj,freq(intPatChangeIdx),azrad,elrad);
                                end
                            else
                                if obj.pNeedSteeringAngle
                                    tempcg = calcRespPattern(...
                                        obj,freq,azrad,elrad,stang);
                                    temp = tempcg(:,:,intPatChangeIdx);
                                else
                                    tempcg = calcRespPattern(...
                                        obj,freq,azrad,elrad);
                                    temp = tempcg(:,:,intPatChangeIdx);
                                end
                            end
                            obj.pSteerVecCov(:,:,intPatChangeIdx) = temp;
                            obj.pIntegratedResponse(intPatChangeIdx) = ...
                                phased.internal.integratePattern(...
                                temp,elrad(:),azstep,elstep);
                        end
                    end
                else
                    if obj.pNeedSteeringAngle
                        intPatChangeIdx = (freq~=obj.pComputedFreq) | ...
                            ~isequal(stang,obj.pComputedSteerAngle); 
                        obj.pComputedFreq = freq;
                        obj.pComputedSteerAngle = stang;
                    else
                        intPatChangeIdx = (freq~=obj.pComputedFreq);
                        obj.pComputedFreq = freq;
                    end
                    
                    weightChangeIdx = any(weights~=obj.pComputedWeights) & ...
                        ~intPatChangeIdx;
                    obj.pComputedWeights = weights;
                                                               
                    if any(intPatChangeIdx)
                        Naz = obj.pNumAzSteps;
                        Nel = obj.pNumElSteps;
                        azstep = 2*pi/Naz;
                        elstep = pi/Nel;
                        azrad = -pi:azstep:pi-azstep;
                        elrad = -pi/2:elstep:pi/2-elstep;
                        
                        if obj.pNeedSteeringAngle
                            if isPolarizationCapable(obj.Sensor)
                                if coder.target('MATLAB')
                                    temp = calcSteeringVectorCovariance(...
                                        obj,freq(intPatChangeIdx),azrad,elrad,stang);
                                else
                                    tempcg = calcSteeringVectorCovariance(...
                                        obj,freq,azrad,elrad,stang);
                                    if obj.EnablePolarization
                                        temp.H = tempcg.H(:,:,intPatChangeIdx);
                                        temp.V = tempcg.V(:,:,intPatChangeIdx);
                                    else
                                        temp = tempcg(:,:,intPatChangeIdx);
                                    end
                                end
                                if obj.EnablePolarization
                                    obj.pSteerVecCov.H(:,:,intPatChangeIdx) = ...
                                        temp.H;
                                    obj.pSteerVecCov.V(:,:,intPatChangeIdx) = ...
                                        temp.V;
                                else
                                    obj.pSteerVecCov(:,:,intPatChangeIdx) = ...
                                        temp;
                                end
                            else
                                if coder.target('MATLAB')
                                    obj.pSteerVecCov(:,:,intPatChangeIdx) = ...
                                        calcSteeringVectorCovariance(...
                                        obj,freq(intPatChangeIdx),azrad,elrad,stang);
                                else
                                    tempcg = ...
                                        calcSteeringVectorCovariance(...
                                        obj,freq,azrad,elrad,stang);
                                    obj.pSteerVecCov(:,:,intPatChangeIdx) = tempcg(:,:,intPatChangeIdx);
                                end
                            end
                        else
                            if isPolarizationCapable(obj.Sensor)
                                if coder.target('MATLAB')
                                    temp = calcSteeringVectorCovariance(...
                                        obj,freq(intPatChangeIdx),azrad,elrad);
                                else
                                    tempcg = calcSteeringVectorCovariance(...
                                        obj,freq,azrad,elrad);
                                    if obj.EnablePolarization
                                        temp.H = tempcg.H(:,:,intPatChangeIdx);
                                        temp.V = tempcg.V(:,:,intPatChangeIdx);
                                    else
                                        temp = tempcg(:,:,intPatChangeIdx);
                                    end
                                end
                                if obj.EnablePolarization
                                    obj.pSteerVecCov.H(:,:,intPatChangeIdx) = ...
                                        temp.H;
                                    obj.pSteerVecCov.V(:,:,intPatChangeIdx) = ...
                                        temp.V;
                                else
                                    obj.pSteerVecCov(:,:,intPatChangeIdx) = ...
                                        temp;
                                end
                            else
                                if coder.target('MATLAB')
                                    obj.pSteerVecCov(:,:,intPatChangeIdx) = ...
                                        calcSteeringVectorCovariance(...
                                        obj,freq(intPatChangeIdx),azrad,elrad);
                                else
                                    tempcg = ...
                                        calcSteeringVectorCovariance(...
                                        obj,freq,azrad,elrad);
                                    obj.pSteerVecCov(:,:,intPatChangeIdx) = tempcg(:,:,intPatChangeIdx);
                                end
                            end
                        end
                
                        if obj.pIsSingleWeights
                            if isPolarizationCapable(obj.Sensor)
                                if obj.EnablePolarization
                                    obj.pIntegratedResponse.H(intPatChangeIdx) = ...
                                        pi/Nel*2*pi/Naz*real(...
                                        sum(bsxfun(@times,sum(bsxfun(@times,conj(weights),...
                                        obj.pSteerVecCov.H(:,:,intPatChangeIdx))),...
                                        permute(weights,[3 1 2])),2));
                                    obj.pIntegratedResponse.V(intPatChangeIdx) = ...
                                        pi/Nel*2*pi/Naz*real(...
                                        sum(bsxfun(@times,sum(bsxfun(@times,conj(weights),...
                                        obj.pSteerVecCov.V(:,:,intPatChangeIdx))),...
                                        permute(weights,[3 1 2])),2));
                                else
                                    obj.pIntegratedResponse(intPatChangeIdx) = ...
                                        pi/Nel*2*pi/Naz*real(...
                                        sum(bsxfun(@times,sum(bsxfun(@times,conj(weights),...
                                        obj.pSteerVecCov(:,:,intPatChangeIdx))),...
                                        permute(weights,[3 1 2])),2));
                                end
                            else
                                obj.pIntegratedResponse(intPatChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,conj(weights),...
                                    obj.pSteerVecCov(:,:,intPatChangeIdx))),...
                                    permute(weights,[3 1 2])),2));
                            end
                        else
                            if isPolarizationCapable(obj.Sensor)
                                if obj.EnablePolarization
                                    obj.pIntegratedResponse.H(intPatChangeIdx) = ...
                                        pi/Nel*2*pi/Naz*real(...
                                        sum(bsxfun(@times,sum(bsxfun(@times,...
                                        permute(conj(weights(:,intPatChangeIdx)),[1 3 2]),...
                                        obj.pSteerVecCov.H(:,:,intPatChangeIdx))),...
                                        permute(weights(:,intPatChangeIdx),[3 1 2])),2));
                                    obj.pIntegratedResponse.V(intPatChangeIdx) = ...
                                        pi/Nel*2*pi/Naz*real(...
                                        sum(bsxfun(@times,sum(bsxfun(@times,...
                                        permute(conj(weights(:,intPatChangeIdx)),[1 3 2]),...
                                        obj.pSteerVecCov.V(:,:,intPatChangeIdx))),...
                                        permute(weights(:,intPatChangeIdx),[3 1 2])),2));
                                else
                                    obj.pIntegratedResponse(intPatChangeIdx) = ...
                                        pi/Nel*2*pi/Naz*real(...
                                        sum(bsxfun(@times,sum(bsxfun(@times,...
                                        permute(conj(weights(:,intPatChangeIdx)),[1 3 2]),...
                                        obj.pSteerVecCov(:,:,intPatChangeIdx))),...
                                        permute(weights(:,intPatChangeIdx),[3 1 2])),2));
                                end
                            else
                                obj.pIntegratedResponse(intPatChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,...
                                    permute(conj(weights(:,intPatChangeIdx)),[1 3 2]),...
                                    obj.pSteerVecCov(:,:,intPatChangeIdx))),...
                                    permute(weights(:,intPatChangeIdx),[3 1 2])),2));
                            end
                        end
                    end
                    
                    if any(weightChangeIdx)
                        
                        Naz = obj.pNumAzSteps;
                        Nel = obj.pNumElSteps;

                        if obj.pIsSingleWeights
                            if obj.EnablePolarization
                                obj.pIntegratedResponse.H(weightChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,conj(weights),...
                                    obj.pSteerVecCov.H(:,:,weightChangeIdx))),...
                                    permute(weights,[3 1 2])),2));
                                obj.pIntegratedResponse.V(weightChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,conj(weights),...
                                    obj.pSteerVecCov.V(:,:,weightChangeIdx))),...
                                    permute(weights,[3 1 2])),2));
                            else
                                obj.pIntegratedResponse(weightChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,conj(weights),...
                                    obj.pSteerVecCov(:,:,weightChangeIdx))),...
                                    permute(weights,[3 1 2])),2));
                            end
                        else
                            if obj.EnablePolarization
                                obj.pIntegratedResponse.H(weightChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,...
                                    permute(conj(weights(:,weightChangeIdx)),[1 3 2]),...
                                    obj.pSteerVecCov.H(:,:,weightChangeIdx))),...
                                    permute(weights(:,weightChangeIdx),[3 1 2])),2));
                                obj.pIntegratedResponse.V(weightChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,...
                                    permute(conj(weights(:,weightChangeIdx)),[1 3 2]),...
                                    obj.pSteerVecCov.V(:,:,weightChangeIdx))),...
                                    permute(weights(:,weightChangeIdx),[3 1 2])),2));
                            else
                                obj.pIntegratedResponse(weightChangeIdx) = ...
                                    pi/Nel*2*pi/Naz*real(...
                                    sum(bsxfun(@times,sum(bsxfun(@times,...
                                    permute(conj(weights(:,weightChangeIdx)),[1 3 2]),...
                                    obj.pSteerVecCov(:,:,weightChangeIdx))),...
                                    permute(weights(:,weightChangeIdx),[3 1 2])),2));
                            end
                        end
                        
                    
                    end
                    
                end
                
                g = obj.pIntegratedResponse;
            else

                if obj.EnablePolarization
                    if isscalar(obj.pComputedFreq)
                        intresp.H = repmat(obj.pIntegratedResponse.H,...
                            1,numel(freq));
                        intresp.V = repmat(obj.pIntegratedResponse.V,...
                            1,numel(freq));
                    else
                        intresp.H = interp1(obj.pComputedFreq,...
                            obj.pIntegratedResponse.H,freq,'nearest','extrap');
                        intresp.V = interp1(obj.pComputedFreq,...
                            obj.pIntegratedResponse.V,freq,'nearest','extrap');
                    end
                else
                    if isscalar(obj.pComputedFreq)
                        intresp = repmat(obj.pIntegratedResponse,...
                            1,numel(freq));
                    else
                        intresp = interp1(obj.pComputedFreq,...
                            obj.pIntegratedResponse,freq,'nearest','extrap');
                    end
                end
                
                g = intresp;
            end
            
        end
        
    end
    
    methods (Access = private)
        function stvcov = calcSteeringVectorCovariance(obj,freq,azrad,elrad,varargin)
            Nfreq = numel(freq);
            az = phased.internal.rad2deg(azrad);
            el = phased.internal.rad2deg(elrad);
            Naz = numel(az);
            Nel = numel(el);
            
            myArrayStv = obj.cIntegratedPattern;
            if isPolarizationCapable(obj.Sensor) 
                if obj.EnablePolarization
                    stvcov = struct(...
                        'H',complex(zeros(obj.pDOF,obj.pDOF,Nfreq)),...
                        'V',complex(zeros(obj.pDOF,obj.pDOF,Nfreq)));
                else
                    stvcov = complex(zeros(obj.pDOF,obj.pDOF,Nfreq));
                end
            else
                stvcov = complex(zeros(obj.pDOF,obj.pDOF,Nfreq));
            end
            if obj.pNeedSteeringAngle
                stang = varargin{1};
                obj.pComputedSteerAngle = stang;
                for n = 1:Nfreq
                    for m = 1:Nel
                        tempV = step(myArrayStv,freq(n),...
                            [az;el(m)*ones(1,Naz)],...
                            obj.pComputedSteerAngle);
                        if isPolarizationCapable(obj.Sensor) 
                            if obj.EnablePolarization
                                stvcov.H(:,:,n) = stvcov.H(:,:,n) + ...
                                    tempV.H*tempV.H'*cos(elrad(m));
                                stvcov.V(:,:,n) = stvcov.V(:,:,n) + ...
                                    tempV.V*tempV.V'*cos(elrad(m));
                            else
                                stvcov(:,:,n) = stvcov(:,:,n) + ...
                                    tempV*tempV'*cos(elrad(m));
                            end
                        else
                            stvcov(:,:,n) = stvcov(:,:,n) + ...
                                tempV*tempV'*cos(elrad(m));
                        end
                    end
                end

            else
                for n = 1:Nfreq
                    for m = 1:Nel
                        tempV = step(myArrayStv,freq(n),...
                            [az;el(m)*ones(1,Naz)]);
                        if isPolarizationCapable(obj.Sensor) 
                            if obj.EnablePolarization
                                stvcov.H(:,:,n) = stvcov.H(:,:,n) + ...
                                    tempV.H*tempV.H'*cos(elrad(m));
                                stvcov.V(:,:,n) = stvcov.V(:,:,n) + ...
                                    tempV.V*tempV.V'*cos(elrad(m));
                            else
                                stvcov(:,:,n) = stvcov(:,:,n) + ...
                                    tempV*tempV'*cos(elrad(m));
                            end
                        else
                            stvcov(:,:,n) = stvcov(:,:,n) + ...
                                tempV*tempV'*cos(elrad(m));
                        end
                    end
                end

            end
            

        end
        
        function resp = calcRespPattern(obj,freq,azrad,elrad,varargin)
            Nfreq = numel(freq);
            az = phased.internal.rad2deg(azrad);
            el = phased.internal.rad2deg(elrad);
            Naz = numel(az);
            Nel = numel(el);
            
            myArrayResp = obj.cIntegratedPattern;
            if Nfreq ~= obj.pLastRunNumFreqs
                release(myArrayResp);
                obj.pLastRunNumFreqs = Nfreq;
            end
            if isPolarizationCapable(obj.Sensor)
                if obj.EnablePolarization
                    resp = struct(...
                        'H',complex(zeros(Nel,Naz,Nfreq)),...
                        'V',complex(zeros(Nel,Naz,Nfreq)));
                else
                    resp = complex(zeros(Nel,Naz,Nfreq));
                end
            else
                resp = complex(zeros(Nel,Naz,Nfreq));
            end

            if obj.pNeedSteeringAngle
                stang = varargin{1};
                obj.pComputedSteerAngle = stang;
                for m = 1:Nel
                    tempV = step(myArrayResp,freq,...
                        [az;el(m)*ones(1,Naz)],...
                        obj.pComputedSteerAngle);
                    if isPolarizationCapable(obj.Sensor)
                        if obj.EnablePolarization
                            resp.H(m,:,:) = resp.H(m,:,:) + ...
                                permute(abs(tempV.H).^2,[3 1 2]);
                            resp.V(m,:,:) = resp.V(m,:,:) + ...
                                permute(abs(tempV.V).^2,[3 1 2]);
                        else
                            resp(m,:,:) = resp(m,:,:) + ...
                                permute(abs(tempV).^2,[3 1 2]);
                        end
                    else
                        resp(m,:,:) = resp(m,:,:) + ...
                            permute(abs(tempV).^2,[3 1 2]);
                    end
                end

            else
                for m = 1:Nel
                    tempV = step(myArrayResp,freq,...
                        [az;el(m)*ones(1,Naz)]);
                    if isPolarizationCapable(obj.Sensor)
                        if obj.EnablePolarization
                            resp.H(m,:,:) = resp.H(m,:,:) + ...
                                permute(abs(tempV.H).^2,[3 1 2]);
                            resp.V(m,:,:) = resp.V(m,:,:) + ...
                                permute(abs(tempV.V).^2,[3 1 2]);
                        else
                            resp(m,:,:) = resp(m,:,:) + ...
                                permute(abs(tempV).^2,[3 1 2]);
                        end
                    else
                        resp(m,:,:) = resp(m,:,:) + ...
                            permute(abs(tempV).^2,[3 1 2]);
                    end
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

