classdef (Sealed,StrictDefaults) PartitionedArray < phased.internal.AbstractSubarray
%PartitionedArray   Phased array partitioned into subarrays
%   H = phased.PartitionedArray creates a partitioned array System object,
%   H. This object models an array that is partitioned into subarrays. The
%   default array is a 4-element ULA partitioned into two 2-element ULAs.
%
%   H = phased.PartitionedArray(Name,Value) creates a partitioned array
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   RESP = step(H,FREQ,ANGLE,V) returns the subarrays' responses, RESP,
%   given the array's operating frequency FREQ (in Hz) and the directions
%   specified in ANGLE (in degrees). FREQ is a row vector of length L and
%   ANGLE can be either a row vector of length M or a 2xM matrix. V is a
%   scalar representing the propagation speed (in m/s). 
%
%   If the subarray is not capable of simulating polarization, RESP has a
%   dimension of NxMxL where N is the number of subarrays in the phased
%   array. Each column of RESP contains the responses of the array elements
%   for the corresponding direction specified in ANGLE. Each page of RESP
%   contains the responses of the array for the given frequency specified
%   in FREQ. If the subarray is capable of simulating polarization, RESP is
%   a structure containing two fields, H and V. H represents the array's
%   response in horizontal polarization and V represents the array's
%   response in vertical polarization. Each field is an NxMxL array whose
%   columns contains the responses of the array elements for the
%   corresponding direction and frequency.
%
%   When ANGLE is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in the [azimuth; elevation] form. The azimuth
%   angle should be between [-180 180] degrees and the elevation angle
%   should be between [-90 90] degrees. If ANGLE is a length-M row vector,
%   then each element specifies a direction's azimuth angle, and the
%   corresponding elevation angle is assumed to be 0.
%   
%   Note that the elements within each subarray are connected to the
%   subarray phase center using an equal-path feed.
%   
%   RESP = step(H,FREQ,ANGLE,V,STEERANGLE) uses STEERANGLE as the
%   subarray's steering direction. STEERANGLE can be either a scalar or a
%   2-element column vector. This syntax is only applicable if you set the
%   SubarraySteering property to either 'Phase' or 'Time'.
%
%   If STEERANGLE is a column vector, it specifies the steering direction
%   in the [azimuth; elevation] form. The azimuth angle should be between
%   [-180 180] degrees and the elevation angle should be between [-90 90]
%   degrees. If STEERANGLE is a scalar, it specifies the direction's
%   azimuth angle, and the corresponding elevation angle is assumed to be
%   0.
%   
%   RESP = step(H,FREQ,ANGLE,V,WS) uses WS as the weights applied to each
%   element in the subarray. WS can be either a matrix or a cell array.
%   This syntax is only applicable when you set the SubarraySteering
%   property to 'Custom'.
%   
%   If each individual subarray has same number of elements, WS must be an
%   NSExN matrix where NSE is the number of elements in each individual
%   subarray and N is the number of subarrays. Each column in WS specifies
%   the weights for the elements in the corresponding subarray.
%
%   If subarrays can have different number of elements, WS can be either an
%   NSExN matrix, where NSE indicates the number of elements in the largest
%   subarray and N is the number of subarrays, or a 1xN cell array, where N
%   is the number of subarrays and each cell contains a column vector whose
%   length is the same as the number of elements of the corresponding
%   subarray.  If WS is a matrix, the first K entries in each column, where
%   K is the number of elements in the corresponding subarray, specifies
%   the weights for the elements in the corresponding subarray. If WS is a
%   cell array, each cell in the array is a column vector specifying the
%   weights for the elements in the corresponding subarray.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   RESP = step(H,FREQ,ANGLE,V,STEERANGLE)
%
%   or
%
%   RESP = step(H,FREQ,ANGLE,V,WS)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   PartitionedArray methods:
%
%   step                  - Output responses of the subarrays (see above)
%   release               - Allow property value and input characteristics
%                           changes
%   clone                 - Create a replicated subarray object with same 
%                           property values
%   collectPlaneWave      - Simulate received plane waves
%   isLocked              - Locked status (logical)
%   isPolarizationCapable - Indicate if the subarray is capable of 
%                           simulating polarization
%   getNumElements        - Number of elements in array
%   getNumSubarrays       - Number of subarrays in array
%   getElementPosition    - Element positions in array
%   getSubarrayPosition   - Subarray positions in array
%   directivity           - Compute array directivity
%   pattern               - Plot subarray response pattern
%   patternAzimuth        - Plot azimuth pattern
%   patternElevation      - Plot elevation pattern
%   viewArray             - View array geometry
%
%   PartitionedArray properties:
%
%   Array                 - Array aperture
%   SubarraySelection     - Subarray definition matrix
%   SubarraySteering      - Subarray steering method
%   PhaseShifterFrequency - Phase shifter frequency
%   NumPhaseShifterBits   - Number of bits in phase shifters
%
%   % Examples:
%
%   % Example 1:
%   %   Construct a 4-element ULA and then partition it into two 2-element 
%   %   ULAs. Plot its azimuth response. Assume the operating frequency is 
%   %   1 GHz and the wave propagation speed is 3e8 m/s.
%
%   array = phased.ULA(4,0.5);
%   subarray = phased.PartitionedArray('Array',array,...
%                   'SubarraySelection',[1 1 0 0;0 0 1 1]);
%   fc = 1e9; c = 3e8;
%   pattern(subarray,fc,-180:180,0,'CoordinateSystem','polar',...
%       'PropagationSpeed',c);
%
%   % Example 2:
%   %   Find the response of each subarray in the above array at the 
%   %   boresight.
%
%   array = phased.ULA(4,0.5);
%   subarray = phased.PartitionedArray('Array',array,...
%                   'SubarraySelection',[1 1 0 0;0 0 1 1]);
%   fc = 1e9; c = 3e8; ang = [0;0];
%   resp = subarray(fc,ang,c)
%
%   See also phased.ReplicatedSubarray, phased.ULA, phased.URA,
%   phased.ConformalArray.

%   Copyright 2011-2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %Array  Array aperture
        %   Specify the array aperture. Array must be of type phased.ULA,
        %   phased.URA, or phased.ConformalArray. The default value of this
        %   property is a 4-element ULA.
        Array
        %SubarraySelection  Subarray definition matrix
        %   Specify the subarray selection as an MxN matrix, where N is the
        %   number of elements in the array and M is the number of
        %   subarrays. Each row of the matrix indicates which elements
        %   belong to the corresponding subarray. Each entry in the matrix
        %   represents whether the corresponding element in the subarray.
        %   If the entry is 0, the element is not present in the array. If
        %   the entry is not 0, it indicates the element is in the subarray
        %   and the value of the entry represents the weights applied to
        %   this element. The default value of this property is [1 1 0 0;0
        %   0 1 1].
        %
        %   Note that the phase center of each subarray is at its geometric
        %   center.
        SubarraySelection = [1 1 0 0;0 0 1 1];
        
    end
    
    methods
        function set.SubarraySelection(obj,value)
            validateattributes(value,{'double'},{'2d','nonnan','finite'},'','SubarraySelection');
            cond = any(sum(value,2)==0);
            if cond
                coder.internal.errorIf(cond,'phased:system:array:expectedNonZeroRows','SubarraySelection');
            end
            obj.SubarraySelection = value;
        end
        
        function set.Array(obj,value)
            validateattributes(value,{'phased.internal.AbstractArray'},...
                {'scalar'},'','Array');
            cond = getNumElements(value) < 2;
            if cond
                coder.internal.errorIf(cond,'phased:system:array:expectMoreElements',...
                    'Array',2);
            end
            obj.Array = value;
        end
    end

    methods

        function obj = PartitionedArray(varargin)
            obj@phased.internal.AbstractSubarray(varargin{:});
            if isempty(coder.target)
                if isempty(obj.Array)
                    obj.Array = phased.ULA('NumElements',4,...
                                           'ElementSpacing',0.5);
                end
            else
                if ~coder.internal.is_defined(obj.Array)
                    obj.Array = phased.ULA('NumElements',4,...
                                           'ElementSpacing',0.5);
                end
            end
        end                
    end
    methods(Hidden)
            function cl = clonecg(obj)
            if strcmp(obj.SubarraySteering,'Phase')
                cl = phased.PartitionedArray(...
                    'Array',clonecg(obj.Array), ...
                    'SubarraySelection',obj.SubarraySelection, ...                
                    'SubarraySteering',obj.SubarraySteering, ...
                    'PhaseShifterFrequency',obj.PhaseShifterFrequency, ...
                    'NumPhaseShifterBits',obj.NumPhaseShifterBits);
            else
                cl = phased.PartitionedArray(...
                    'Array',clonecg(obj.Array), ...
                    'SubarraySelection',obj.SubarraySelection, ...                
                    'SubarraySteering',obj.SubarraySteering);
            end
        end

    end
    methods (Access = protected)
        function num = calcNumSubarrays(obj)
            num = size(obj.SubarraySelection,1);
        end
        
        function num = calcNumElements(obj,subarrayidx)
            if nargin > 1
                sigdatatypes.validateIndex(subarrayidx,'getNumElements',...
                    'SubarrayIndex',{'row','<=',getNumSubarrays(obj)});
                num = (sum(logical(obj.SubarraySelection(subarrayidx,:)),2)).';
            else
                num = getNumElements(obj.Array);
            end
        end
        
        function pos = calcSubarrayPosition(obj)
            N = coder.internal.const(getNumSubarrays(obj));
            elem_pos = coder.internal.const(getElementPosition(obj));
            sSelection = coder.internal.const(obj.SubarraySelection);
            pos = coder.internal.const(...
              calcSubarrayPositionCG(N,elem_pos,sSelection));
        end
         
        function pos = calcElementPosition(obj)
            pos = getElementPosition(obj.Array);
        end
    end
    
    methods (Access = protected)

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSubarray(obj);
            s.Array = saveobj(obj.Array);
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSubarray(obj,s);
            obj.Array = eval(...
                sprintf('%s.loadobj(s.Array)',s.Array.ClassNameForLoadTimeEval));
            s = rmfield(s,'Array');
            if isfield(s,'isLocked')
                s = rmfield(s,'isLocked');
            end
            if isfield(s,'pPropagationSpeed')  % obsolete property
                s = rmfield(s,'pPropagationSpeed');
            end
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function validatePropertiesImpl(obj)
            N = getNumElements(obj.Array);
            cond = size(obj.SubarraySelection,2) ~= N;
            if cond
                coder.internal.errorIf(cond,'phased:phased:invalidColumnNumbers',...
                    'SubarraySelection',N);
            end
        end

        function setupImpl(obj,freq,ang,varargin)

            if isempty(coder.target)
                obj.cArray = clone(obj.Array);
                release(obj.cArray);
            else
                obj.cArray = clonecg(obj.Array);
            end

            setupImpl@phased.internal.AbstractSubarray(obj,freq,ang,varargin{:});
            if isPolarizationEnabled(obj)
                obj.cModuleResponse = phased.SteeringVector(...
                    'SensorArray',obj.cArray,...
                    'IncludeElementResponse',true,...
                    'EnablePolarization',true);
                obj.cModuleResponse.PropagationSpeedSource = 'Input port';
            else
                obj.cModuleResponse = phased.SteeringVector(...
                    'SensorArray',obj.cArray,...
                    'IncludeElementResponse',true);
                obj.cModuleResponse.PropagationSpeedSource = 'Input port';
            end
            obj.cModuleSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.cArray,...
                'IncludeElementResponse',false);
            if obj.pPhaseSteeringFlag
                obj.cModuleSteeringVector.NumPhaseShifterBits = ...
                    obj.NumPhaseShifterBits;
            end
            obj.cModuleSteeringVector.PropagationSpeedSource = 'Input port';
       end

       function resp = stepImpl(obj, freq, angArg, c, stang)

            sigdatatypes.validateSpeed(c,'','V',{'positive'});
                        
            if size(angArg,1) == 1
                ang = [angArg;zeros(size(angArg))];
            else
                ang = angArg;
            end
            
            if obj.pTimeSteeringFlag || obj.pPhaseSteeringFlag
                [f,theta] = getSteeringConfiguration(obj,freq,stang);
            end
            
            [Na,N,M,L] = calcRespDimensions(obj);
            
            % For partitioning, if no steering is used, because the
            % element is already at its final location, no compensation
            % needed for subarray phase center. If there is steering,
            % because steering vector cancels out the phase center, it
            % needs to compensate for the subarray phase center in the
            % steering vector because they should be within each
            % subarray.

            if obj.pNoSteeringFlag
                w = ones(Na,1);
                wstv = ones(N,1);
                subarraystv = phased.internal.steeringvec(...
                    obj.pSubarrayPosition,freq,c,ang);
            elseif obj.pPhaseSteeringFlag || obj.pTimeSteeringFlag
                w = step(obj.cModuleSteeringVector,f,theta,c);
                wstv = phased.internal.steeringvec(...
                    obj.pSubarrayPosition,f,c,theta);
                subarraystv = phased.internal.steeringvec(...
                    obj.pSubarrayPosition,freq,c,ang);
            else % custom weights
                w = stang;  % this is like a taper, needs conj when passed in
                if iscell(w)
                    for m = 1:N
                        cond = ~all(isfinite(w{m}));
                        if cond        
                            coder.internal.errorIf(cond,'phased:step:expectedFinite','W');
                        end
                    end
                else
                    cond = ~all(isfinite(w(:)));
                    if cond        
                        coder.internal.errorIf(cond,'phased:step:expectedFinite','W');
                    end
                end
                wstv = ones(N,1);
                subarraystv = phased.internal.steeringvec(...
                    obj.pSubarrayPosition,freq,c,ang);
            end
            if isPolarizationEnabled(obj)
                elemresp = ...
                    step(obj.cModuleResponse,freq,ang,c);
                resp_h = complex(zeros(N,M,L));
                resp_v = complex(zeros(N,M,L));
            else
                elemresp = step(obj.cModuleResponse,freq,ang,c);
                resp = complex(zeros(N,M,L));
            end
            selmatrix = logical(obj.SubarraySelection);
            selmatweights = obj.SubarraySelection.';
            for m = N:-1:1
                idx = selmatrix(m,:);
                if obj.pIsSingleWeights
                    if ~(obj.pNoSteeringFlag || obj.pPhaseSteeringFlag || ...
                            obj.pTimeSteeringFlag)
                        w_temp = conj(getElementWeightsForSubarray(obj,w,m))*conj(wstv(m)).*selmatweights(idx,m);
                    else
                        w_temp = w(idx)*conj(wstv(m)).*selmatweights(idx,m);
                    end
                    if isPolarizationEnabled(obj)
                        elemresp_h_temp = bsxfun(@times,elemresp.H(idx,:),conj(subarraystv(m,:)));
                        elemresp_v_temp = bsxfun(@times,elemresp.V(idx,:),conj(subarraystv(m,:)));
                        if obj.pIsSingleFrequency
                            resp_h(m,:) = w_temp'*elemresp_h_temp;
                            resp_v(m,:) = w_temp'*elemresp_v_temp;
                        else
                            resp_h(m,:,:) = reshape(...
                                w_temp'*elemresp_h_temp,[],M,L);
                            resp_v(m,:,:) = reshape(...
                                w_temp'*elemresp_v_temp,[],M,L);
                        end
                    else
                        elemresp_temp = bsxfun(@times,elemresp(idx,:),conj(subarraystv(m,:)));
                        if obj.pIsSingleFrequency
                            resp(m,:) = w_temp'*elemresp_temp;
                        else
                            resp(m,:,:) = reshape(...
                                w_temp'*elemresp_temp,[],M,L);
                        end
                    end
                else % multiple weights
                    % Note Custom steering does not apply here
                    if isPolarizationEnabled(obj)
                        for l = L:-1:1
                            w_temp = w(idx,l)*conj(wstv(m,l)).*selmatweights(idx,m);
                            elemresp_h_temp = bsxfun(@times,elemresp.H(idx,:,l),conj(subarraystv(m,:,l)));
                            elemresp_v_temp = bsxfun(@times,elemresp.V(idx,:,l),conj(subarraystv(m,:,l)));
                            resp_h(m,:,l) = w_temp'*elemresp_h_temp;
                            resp_v(m,:,l) = w_temp'*elemresp_v_temp;
                        end
                    else
                        for l = L:-1:1
                            w_temp = w(idx,l)*conj(wstv(m,l)).*selmatweights(idx,m);
                            elemresp_temp = bsxfun(@times,elemresp(idx,:,l),conj(subarraystv(m,:,l)));
                            resp(m,:,l) = w_temp'*elemresp_temp;
                        end
                    end
                end
            end

            if isPolarizationEnabled(obj)
                resp.H = reshape(resp_h,N,M,L);
                resp.V = reshape(resp_v,N,M,L);
            end
        end

        function validateElementWeights(obj,ws)
            % for ParitionedArray, it is possible to have unequal number of
            % elements in each subarray. That requires cell array support.
            cond = (~iscell(ws) && ~ismatrix(ws)) || isempty(ws);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:phased:expectedCellOrMatrix','WS');
            end
            Ns = getNumSubarrays(obj);
            Nse = sum(logical(obj.SubarraySelection),2);
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
                            'phased:system:array:SubarrayElementWeightsSizeMismatch',m,'WS',Nse(m));
                    end
                    cond = ~isa(ws{m},'double');
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:system:array:SubarrayElementWeightsInvalidDataType',m,'WS','double');
                    end
                end
            else % matrix
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
        
        function wele = getElementWeightsForSubarray(obj,w,subidx)

            % subidx is the index of the subarray
            if iscell(w)
                wele = w{subidx};
            else % ismatrix
                wele = w(1:getNumElements(obj,subidx),subidx);
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractSubarray(obj, prop);
        end
        
    end
    
    methods (Hidden, Access = {?phased.gpu.internal.AbstractClutterSimulator, ?phased.internal.AbstractSubarray})
        function [pos, az, el] = getSubarrayPosAzEl(obj, subArrIdx, azin, elin)
            %Get the element positions for this subarray. az and el remain
            %the same.
            posAll = getElementPosition(obj);
            selmatrix = logical(obj.SubarraySelection(subArrIdx, :));
            pos = posAll(:, selmatrix);
            az = azin;
            el = elin;
        end   
        function sv = getSubArraySteeringVec(obj, subArrIdx, freq, c, stang )
            %Do this on the CPU, it's small\
            if strncmp(obj.SubarraySteering,'Phase',1) || ...
                    strncmp(obj.SubarraySteering,'Time',1)
                %Get the steering vector for this subarray at steering angle
                %stang.
                if strncmp(obj.SubarraySteering,'Phase',1)
                    f = obj.PhaseShifterFrequency;
                    theta = stang;
                else % if strncmp(obj.SubarraySteering,'Time',1)
                    f = freq;
                    theta = stang;
                end
            
                pos = getElementPosition(obj);
                subArrElemPos = pos(:,logical(obj.SubarraySelection(subArrIdx,:)));
                w = phased.internal.steeringvec(subArrElemPos, f,c,theta);
                
                posSub = getSubarrayPosition(obj);
                wstv = phased.internal.steeringvec(...
                    posSub(:, subArrIdx),f,c,theta);
                
                %second conjugation to match Hermitian transpose in resp
                %computation above.
                sv = conj(w*conj(wstv));
                
            else % Custom
                sv = getElementWeightsForSubarray(obj,stang,subArrIdx);
            end
            
        end
            
    end
    
    methods (Access = {?phased.ArrayGain, ?phased.gpu.internal.AbstractClutterSimulator})
        function Rn = getNoiseGainMatrix(obj,freq,c,stangArg)
            if strcmp(obj.SubarraySteering,'None')
                w = ones(getNumElements(obj.Array),1);
                ws = ones(getDOF(obj.Array),1);
            elseif strcmp(obj.SubarraySteering,'Phase')
                if size(stangArg,1) == 1
                    stang = [stangArg; zeros(size(stangArg))];
                else
                    stang = stangArg;
                end
                hstv = phased.SteeringVector(...
                    'SensorArray',obj.Array,...
                    'PropagationSpeed',c,...
                    'NumPhaseShifterBits',obj.NumPhaseShifterBits);
                w = step(hstv,obj.PhaseShifterFrequency,stang);
                ws = phased.internal.steeringvec(...
                    getSubarrayPosition(obj),...
                    obj.PhaseShifterFrequency,c,stang);
            elseif strcmp(obj.SubarraySteering,'Time')
                if size(stangArg,1) == 1
                    stang = [stangArg; zeros(size(stangArg))];
                else
                    stang = stangArg;
                end
                hstv = phased.SteeringVector(...
                    'SensorArray',obj.Array,...
                    'PropagationSpeed',c);
                w = step(hstv,freq,stang);
                ws = phased.internal.steeringvec(...
                    getSubarrayPosition(obj),freq,c,stang);
            else % custom weights
                w = stangArg;
                ws = ones(getDOF(obj.Array),1);
            end
            M = getDOF(obj);
            Rn = complex(zeros(M));
            SelMat = obj.SubarraySelection;
            for m = 1:M
                for n = 1:M
                    idx1 = find(SelMat(m,:));
                    idx2 = find(SelMat(n,:));
                    RnTemp = bsxfun(@eq,idx1(:),idx2);
                    if strcmp(obj.SubarraySteering,'Custom')
                        w1 = getElementWeightsForSubarray(obj,w,m);
                        w2 = getElementWeightsForSubarray(obj,w,n);
                    else
                        w1 = w(idx1);
                        w2 = w(idx2);
                    end
                    Rn(m,n) = ws(m)*w1'*RnTemp*w2*conj(ws(n));
                end
            end
        end
    end
    
    methods 
        function flag = isPolarizationCapable(obj)
        %isPolarizationCapable Indicate if the subarray is capable of 
        %simulating polarization
        %   F = isPolarizationCapable(H) returns the flag, F, which
        %   indicates whether the subarray H is capable of simulating
        %   polarization.
        %
        %   % Example: 
        %   %   Determine whether a subarray is capable of simulating 
        %   %   polarization.
        %   
        %   h = phased.PartitionedArray;
        %   f = isPolarizationCapable(h)
        %
        %   See also phased.
            flag = isPolarizationCapable(obj.Array);
        end
        
    end
    
    methods (Hidden)
        function flag = isPolarizationEnabled(obj)
        %isPolarizationEnabled Indicate if the subarray is enabled to 
        %simulate polarization
        %   F = isPolarizationEnabled(H) returns the flag, F, which
        %   indicates whether the subarray H is enabled to simulate
        %   polarization.
        %
        %   % Example: 
        %   %   Determine whether a subarray is enabled to simulate 
        %   %   polarization.
        %   
        %   h = phased.PartitionedArray;
        %   f = isPolarizationEnabled(h)
        %
        %   See also phased.
            flag = isPolarizationEnabled(obj.Array);
        end
        
    end
    
    methods (Hidden)
        function h = getElementHandle(obj)
            if isa(obj.Array,'phased.internal.AbstractHeterogeneousArray')
                h = obj.Array.ElementSet;
            else
                h = obj.Array.Element;
            end
        end
        function newObj = cloneSensor(obj)
            if isElementFromAntenna(obj)
                newObj = clone(obj); 
                release(newObj); 
                newObj.Array = cloneSensor(obj.Array);
            else
                newObj = clone(obj);
                release(newObj);
            end
        end
        function flag = isElementFromAntenna(obj) 
            flag = isElementFromAntenna(obj.Array);
        end
    end

    methods(Static, Hidden, Access=protected)
        function groups = getPropertyGroupsImpl
            sensorClassSet = matlab.system.display.internal.ClassStringSet(...
                {
                    'phased.ULA',...
                    'phased.URA',...
                    'phased.UCA',...                    
                    'phased.ConformalArray'...
                }, ...
                'ConstructorExpressions', ...
                {
                    'phased.ULA(''NumElements'',4)',...
                    'phased.URA',...
                    'phased.UCA',...                    
                    'phased.ConformalArray(''ElementPosition'',[0 0 0 0;-0.25 -0.25 0.25 0.25;0.25 -0.25 0.25 -0.25])' ...
                }, ...                
                'PropertiesTitle', '', ...
                'NestDisplay', false, ...
                'Labels', { 
                    'ULA',...
                    'URA',...
                    'UCA',...                    
                    'Conformal Array'...
                          });
            sensorProp = matlab.system.display.internal.Property('Array', ...
                'Description', 'Geometry', ...
                'ClassStringSet', sensorClassSet);  
            props = {sensorProp,...
                'SubarraySelection',...
                'SubarraySteering',...
                'PhaseShifterFrequency',...
                'NumPhaseShifterBits'};
            groups = matlab.system.display.Section('Title', 'Parameters', ...
                'PropertyList', props);
        end
    end

end

function pos = calcSubarrayPositionCG(N,elem_pos,subarraySelection)
pos = zeros(3,N);
%for m = 1:N
for m = coder.unroll(1:N)
    subSelect = subarraySelection(m,:);
    lsubSelect = logical(subSelect);
    selemPos = elem_pos(:,lsubSelect);
    pos(:,m) = mean(selemPos,2);
end
end


% [EOF]
