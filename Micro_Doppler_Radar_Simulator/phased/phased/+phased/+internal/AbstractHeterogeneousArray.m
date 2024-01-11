classdef (Hidden) AbstractHeterogeneousArray < phased.internal.AbstractArray
%This class is for internal use only. It may be removed in the future.

%AbstractHeterogeneousArray   Define the AbstractHeterogeneousArray class.

%   Copyright 2012-2016 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %ElementSet  Set of elements used in the array
        %   Specify the set of different elements used in the sensor array
        %   as a row MATLAB cell array. Each member of the cell array
        %   contains an element object in the phased package. The default
        %   value of this property is a cell containing one isotropic
        %   antenna element.
        %
        %   Elements specified in the ElementSet property must be either
        %   all antennas or all microphones. In addition, all specified
        %   antenna elements should have same polarization capability.
        ElementSet
    end
    
    properties (Abstract, Nontunable)
        %ElementIndices Mapping of elements in the array
        %   Need to document and process according to the array.
        ElementIndices
    end

    properties (Access = protected, Nontunable)
        pMap
        pisUsedType        
    end
    
    methods
        function set.ElementSet(obj,val)
            cond = ~(iscell(val) && isrow(val) && isscalar(val{1}) &&  ...
                 (isa(val{1},'phased.internal.AbstractElement') || isa(val{1},'em.Antenna')));
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidElementSet','ElementSet');
            end

            if isa(val{1}, 'phased.internal.AbstractMicrophoneElement')
                firstSensorType = 'Microphone';
            else
                firstSensorType = 'Antenna';
            end
            for idx = 2:length(val)
                cond = (strcmp(firstSensorType,'Microphone') && ...
                    ~isa(val{idx},'phased.internal.AbstractMicrophoneElement')) || ...
                    (strcmp(firstSensorType,'Antenna') && ...
                    isa(val{idx},'phased.internal.AbstractMicrophoneElement'));
                if cond
                    coder.internal.errorIf(cond,'phased:system:array:InvalidRemainingElementTypes','ElementSet');
                end
            end

            if strcmp(firstSensorType,'Antenna')
                if isElementFromAntenna(val{1})
                    flag = true;
                else
                    flag = isPolarizationCapable(val{1});
                end
                for idx = 2:numel(val)
                    if isa(val{idx},'phased.internal.AbstractAntennaElement')
                        polflag = isPolarizationCapable(val{idx});
                    else
                        polflag = true;
                    end
                    cond = xor(flag,polflag);
                    if cond
                        coder.internal.errorIf(cond,'phased:system:array:InvalidAntennaPolarizationCapability',...
                            'isPolarizationCapable');
                    end
                end
            end
            obj.ElementSet = val;
        end
        
    end
    
    methods

        function obj = AbstractHeterogeneousArray(varargin)

            obj@phased.internal.AbstractArray(varargin{:});
            if isempty(coder.target)
                if isempty(obj.ElementSet)
                    obj.ElementSet = {phased.IsotropicAntennaElement};
                end
            else
                if ~coder.internal.is_defined(obj.ElementSet)
                    obj.ElementSet = {phased.IsotropicAntennaElement};
                end
            end

        end

    end
    
    methods (Access = protected)
        function validatePropertiesImpl(obj)
            numEl = getNumElements(obj);
            numElTypes = numel(obj.ElementSet);

            cond = any(obj.ElementIndices(:) > numElTypes) || ...
                numel(obj.ElementIndices) ~= numEl;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:InvalidElementMap', ...
                    'ElementIndices',numEl,numElTypes);
            end
            
        end
        
        function setupImpl(obj,freq,ang) 
            sz_freq = size(freq);
            obj.pNumFreqs = sz_freq(2);
            sz_angle = size (ang);
            obj.pNumAngles = sz_angle(2);
            num_angles = sz_angle(2);
            num_elements = getNumElements(obj);
            
            num_patterns = numel(obj.ElementSet); 
            tempElement = cell(1,num_patterns);
            for idx = 1:num_patterns
                if isElementFromAntenna(obj.ElementSet{idx})
                    tempElement{idx} = phased.internal.AntennaAdapter(...
                        'Antenna',obj.ElementSet{idx});
                else
                    if isempty(coder.target)
                        tempElement{idx} =  cloneSensor(obj.ElementSet{idx});
                    else
                        tempElement{idx} =  clonecg(obj.ElementSet{idx});
                    end
                end
                release(tempElement{idx});
            end
            obj.cElement = tempElement;
            lElementMap = obj.ElementIndices(:).';
            
            [lMap,isUsedType] = getUsageMap(lElementMap,num_patterns,num_elements,num_angles);
            obj.pMap = lMap;
            obj.pisUsedType = isUsedType;
                               
            num_freqs = sz_freq(2);
            obj.pNumFreqs = num_freqs;
            obj.pNumAngles = num_angles;
            obj.pAzimuthOnly = (sz_angle(1) == 1);
            
            lTaper = getTaper(obj);
            if isscalar(lTaper)
                obj.pTaper = lTaper;
            else
                obj.pTaper =repmat(lTaper,num_angles,num_freqs);
            end
   
        end
        
        function resetImpl(obj)
            for idx = 1:length(obj.cElement)
                reset(obj.cElement{idx});
            end
        end
        
        function releaseImpl(obj)
            for idx = 1:length(obj.cElement)
                release(obj.cElement{idx});
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractArray(obj);
            for idx=1:length(obj.ElementSet)
                s.ElementSet{idx} = saveobj(obj.ElementSet{idx});
            end
            if isLocked(obj)
                for idx=1:length(obj.cElement)
                    s.cElement{idx} = saveobj(obj.cElement{idx});
                end
                s.pMap = obj.pMap;
                s.pisUsedType = obj.pisUsedType;
            end
        end
        
        function s = loadSubObjects(obj,s)
            elemset = cell(size(s.ElementSet));
            for idx=1:length(s.ElementSet)
                if isa(s.ElementSet{idx},'em.Antenna')
                    elemset{idx} = em.EmStructures.loadobj(s.ElementSet{idx});
                else
                    elemset{idx} = phased.internal.AbstractElement.loadobj(s.ElementSet{idx});
                end
            end
            obj.ElementSet = elemset;
            s = rmfield(s,'ElementSet');
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cElement = cell(size(s.cElement));
                    for idx=1:length(s.cElement)
                        obj.cElement{idx} = ...
                            phased.internal.AbstractElement.loadobj(s.cElement{idx});
                    end
                    s = rmfield(s,'cElement');
                end
                s = rmfield(s,'isLocked');
            end
        end
    end
        
    methods (Access = protected)
        
        function resp = getPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles)
            num_elements = getNumElements(obj);
            inc_angle = convertIncidentToAzEl(obj,ang);
            
            lMap = obj.pMap;
            isUsedType = obj.pisUsedType;
            resp_h = zeros(num_elements*num_angles,obj.pNumFreqs,'like',1+1i);
            resp_v = zeros(num_elements*num_angles,obj.pNumFreqs,'like',1+1i);
            for idx = 1:numel(obj.cElement)
                if isUsedType(idx)
                    currentIdxMap = lMap(idx,:);
                    temp = ...
                        step(obj.cElement{idx},freq,...
                        inc_angle(:,currentIdxMap));
                    resp_h(currentIdxMap,:) = temp.H;
                    resp_v(currentIdxMap,:) = temp.V;
                end
            end
            
            resp_h = resp_h.*obj.pTaper;
            resp.H = reshape(resp_h,num_elements,num_angles,[]);
            resp_v = resp_v.*obj.pTaper;
            resp.V = reshape(resp_v,num_elements,num_angles,[]);
        end
        
        function resp = getNonPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles)
            num_elements = getNumElements(obj);
            inc_angle = convertIncidentToAzEl(obj,ang);
            
            lMap = obj.pMap;
            isUsedType = obj.pisUsedType;
            resp = zeros(num_elements*num_angles,obj.pNumFreqs,'like',1+1i);
            for idx = 1:numel(obj.cElement)
                if isUsedType(idx)
                    currentIdxMap = lMap(idx,:);
                    resp(currentIdxMap,:) = ...
                        step(obj.cElement{idx},freq,...
                             inc_angle(:,currentIdxMap));
                end
            end
            resp = resp.*obj.pTaper;
            resp = reshape(resp,num_elements,num_angles,[]);
        end
    end
    
    methods(Static, Hidden, Access=protected)
        function groups = getPropertyGroupsImpl
            props = {...
              'ElementSet'};
            groups = matlab.system.display.Section('Title', 'Parameters', ...
                                                   'PropertyList', props);
        end
    end
    
    methods
        function flag = isPolarizationCapable(obj)
        %isPolarizationCapable Indicate if the array is capable of 
        %simulating polarization
        %   F = isPolarizationCapable(H) returns the flag, F, which
        %   indicates whether the array H is capable of simulating
        %   polarization.
        %
        %   % Example:
        %   %   Determine whether an array of IsotropicAntennaElement 
        %   %   is capable of simulating polarization.
        %   
        %   h = phased.HeterogeneousULA;
        %   f = isPolarizationCapable(h)
        %
        %   See also phased.

            flag = true;
            for idx = 1:numel(obj.ElementSet)
                if ~isElementFromAntenna(obj.ElementSet{idx})
                    flag = flag && ...
                        isPolarizationCapable(obj.ElementSet{idx});
                end
            end
        end
        
    end
    
    methods (Hidden)
        function flag = isPolarizationEnabled(obj)
        %isPolarizationEnabled Indicate if the array is enabled to 
        %simulate polarization
        %   F = isPolarizationEnabled(H) returns the flag, F, which
        %   indicates whether the array H is enabled to simulate
        %   polarization.
        %
        %   % Example:
        %   %   Determine whether an array of IsotropicAntennaElement 
        %   %   is enabled to simulate polarization.
        %   
        %   h = phased.HeterogeneousULA;
        %   f = isPolarizationEnabled(h)
        %
        %   See also phased.

            flag = true;
            for idx = 1:numel(obj.ElementSet)
                if isElementFromAntenna(obj.ElementSet{idx}) 
                    if isLocked(obj)
                        flag = flag && ...
                            isPolarizationEnabled(obj.cElement{idx});
                    end
                else
                    flag = flag && ...
                        isPolarizationEnabled(obj.ElementSet{idx});
                end
            end
        end
        
    end
    
    methods (Hidden)
        function h = getElementHandle(obj)
            h = obj.ElementSet;
        end
        function newObj = cloneSensor(obj)
            if isElementFromAntenna(obj)
                newObj = clone(obj); 
                release(newObj); 
                for idx = 1:numel(obj.ElementSet)
                    if isElementFromAntenna(obj.ElementSet{idx})
                        newObj.ElementSet{idx} = phased.internal.AntennaAdapter(...
                            'Antenna',obj.ElementSet{idx});
                    end
                end
            else
                newObj = clone(obj);
                release(newObj);
            end
        end
        function flag = isElementFromAntenna(obj) 
            flag = false;
            for idx = 1:numel(obj.ElementSet)
                if isElementFromAntenna(obj.ElementSet{idx})
                    flag = true;
                    break;
                end
            end
        end
    end
    
    methods (Hidden, Access = {?phased.internal.AbstractArray, ?phased.gpu.internal.AbstractClutterSimulator})
        %Methods used by the GPU ConstantGammaClutter model.
        function xscaled = scaleByElemResponse(obj, azin, elin, freq, x, idx)
            %scale x - the steering vector, by the element pattern for the
            %elements selected by the indices in idx
            xscaled = gpuArray.zeros(size(x));

            %For heterogeneous arrays scale x by the element pattern
            %for each different element. use partitionElements to get
            %the indices for identical elements in the ElementIndices. 
            part = partitionElements(obj, obj.ElementIndices(idx));
            for ii=1:numel(part)
                if ~isempty(part{ii})
                    hidx = part{ii}; %homogeneous index
                    if ~isvector(azin) %slice out only the indices needed.
                        aztmp = azin(:,hidx,:);
                        eltmp = elin(:,hidx, :);
                    else
                        aztmp = azin;
                        eltmp = elin;
                    end
                    if isElementFromAntenna(obj.ElementSet{ii})
                        hele = cloneSensor(obj.ElementSet{ii});
                    else
                        hele = obj.ElementSet{ii};
                    end
                    elPat = getgpuElemResponse(hele, aztmp, eltmp, freq);
                    if isscalar(obj.Taper)
                        elPatwTaper = elPat.*obj.Taper;
                    else
                        elPatwTaper = bsxfun(@times, elPat, reshape(obj.Taper(hidx), 1, [], 1));
                    end
                    xscaled(:, hidx, :) = bsxfun(@times, elPatwTaper, x(:,hidx,:));
                end
            end
        end
    end
    
    methods (Access = private)
        function idxOut = partitionElements(obj, idx)
            %Returns a cell array of indices. Each element of the cell array is an
            %index vector to extract elements of the same type from the element map,
            %emap.
            
            %Calculate how many indices to return == number of antenna
            %types
            emap = obj.ElementIndices(idx);
            numtypes = max(emap);
            
            %Count the number of each type and build start and end index
            %using scan (cumsum)
            cnt = histc(emap, (1:numtypes));
            endIdx = cumsum(cnt);
            startIdx = endIdx - cnt + 1;
            
            %Use the indices to pick apart the index vector that results
            %from sorting the element map. These new indices will extract
            %identical elements from the ElementIndices
            [~,ei] = sort(emap);
            idxOut = cell(1,numtypes);
            for ii=1:numel(idxOut)
                idxOut{ii} = ei(startIdx(ii):endIdx(ii));
            end
        end
    end
        
end

function [lMap,isUsedType] = getUsageMap(lElementMap,num_patterns,num_elements,num_angles)
            lMap = false(num_patterns,num_elements*num_angles);
            isUsedType = false(1,num_patterns);
            for idx = 1:num_patterns
                currentMap = (lElementMap == idx);
                if any(currentMap)
                    isUsedType(idx) = true;
                    lMap(idx,:) = repmat(currentMap,1,num_angles);
                end
            end
end


% [EOF]
