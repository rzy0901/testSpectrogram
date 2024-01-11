classdef (Sealed,StrictDefaults) ArrayResponse < phased.internal.AbstractArrayOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%ArrayResponse  Sensor array response
%   H = phased.ArrayResponse creates an array response System object, H.
%   This object calculates the response of a sensor array for the specified
%   directions. By default a 2-element uniform linear array (ULA) is used.
%
%   H = phased.ArrayResponse(Name,Value) creates an array response object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   RESP = step(H,FREQ,ANG) returns the array response, RESP, for the
%   directions specified in ANG (in degrees) at the given operating
%   frequency FREQ (in Hz). FREQ is a row vector of length L and ANG can
%   be either a row vector of length M or a 2xM matrix. 
%
%   When ANG is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in the [azimuth; elevation] form. The azimuth
%   angle should be between [-180 180] degrees and the elevation angle
%   should be between [-90 90] degrees. If ANG is a length M row vector,
%   each element specifies a direction's azimuth angle and the
%   corresponding elevation angle is assumed to be 0.
%
%   If you set the EnablePolarization property to false, RESP is an MxL
%   matrix whose elements contain the responses of the sensor array at
%   angles specified in ANG and frequencies specified in FREQ. If you set
%   the EnablePolarization property to true, RESP is a structure containing
%   two fields, H and V. H represents the array's response in horizontal
%   polarization and V represents the array's response in vertical
%   polarization. Note that when you set the EnablePolarization to false,
%   the polarization information is discarded and the combined pattern from
%   both H and V polarizations is used at each element to compute the array
%   response.
%
%   RESP = step(H,FREQ,ANG,W) uses W as the weights applied on the sensor
%   array when you set the WeightsInputPort property to true. W can be a
%   length-N column vector or an NxL matrix where N is the number of
%   subarrays if SensorArray contains subarrays, or the number of elements
%   otherwise. L is the number of frequencies specified in FREQ. If W is a
%   vector, the weights are applied at all frequencies. If W is a matrix,
%   each column of W represents the weights used at the corresponding
%   frequency specified in FREQ.
%
%   RESP = step(H,FREQ,ANG,STEER) uses STEER as the subarray steering angle
%   (in degrees). STEER can be a scalar or a length-2 column vector. If
%   STEER is a vector, it is in the form of [AzimuthAngle; ElevationAngle].
%   If STEER is a scalar, it represents the azimuth angle and the elevation
%   angle is assumed to be 0. This syntax is only applicable when you use
%   subarrays in the SensorArray property and set the SubarraySteering
%   property in the SensorArray to either 'Phase' or 'Time'.
%
%   RESP = step(H,FREQ,ANG,WS) uses WS as the weights applied to each
%   element in the subarray. WS can be either a matrix or a cell array.
%   This syntax is only applicable when you use subarrays in the
%   SensorArray property, set the SubarraySteering property in the
%   SensorArray to 'Custom' and set the IncludeElementResponse property to
%   true.
%   
%   If the Sensor property is a phased.ReplicatedSubarray, WS must be an
%   NSExN matrix where NSE is the number of elements in each individual
%   subarray and N is the number of subarrays. Each column in WS specifies
%   the weights for the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its individual
%   subarrays have same number of elements, WS must be an NSExN matrix
%   where NSE is the number of elements in each individual subarray and N
%   is the number of subarrays. Each column in WS specifies the weights for
%   the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its subarrays
%   can have different number of elements, WS can be either an NSExN
%   matrix, where NSE indicates the number of elements in the largest
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
%   RESP = step(H,FREQ,ANG,W,STEER)
%
%   or
%
%   RESP = step(H,FREQ,ANG,W,WS)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ArrayResponse methods:
%
%   step     - Calculate the array response of the sensor array (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create an array response object with same property values
%   isLocked - Locked status (logical)
%
%   ArrayResponse properties:
%
%   SensorArray        - Sensor array 
%   PropagationSpeed   - Signal propagation speed 
%   WeightsInputPort   - Add input to specify the weights
%   EnablePolarization - Enable polarization simulation
%
%   % Examples:
%
%   % Example 1:
%   %   Calculate and plot the azimuth response between -90 and 90 degrees
%   %   for a 4-element ULA. Assume the array's operating frequency is 300
%   %   MHz.
%
%   ha = phased.ULA(4);
%   response = phased.ArrayResponse('SensorArray',ha);
%   ang = -90:90; fc = 300e6; 
%   resp = response(fc,ang);
%   plot(ang,abs(resp)); xlabel('Angle (degrees)'); ylabel('Magnitude');
%   title('Azimuth Response');
%
%   % Example 2:
%   %   Calculate the array response for a 4-element uniform linear array 
%   %   at the direction of 30 degrees azimuth and 20 degrees elevation.
%   %   Assume the array's operating frequency is 300 MHz.
%
%   ha = phased.ULA(4);
%   response = phased.ArrayResponse('SensorArray',ha);
%   fc = 3e8; ang = [30;20];
%   resp = response(fc,ang)
%
%   See also phased, phased.ElementDelay, phased.SteeringVector,
%   phased.ArrayGain, phased.ULA/pattern, phased.URA/pattern,
%   phased.ConformalArray/pattern.

%   Copyright 2010-2017 The MathWorks, Inc.

%   Reference
%   [1] Harry Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable, Logical) 
        %WeightsInputPort    Add input to specify the weights
        %   Set this property to true to add input to specify the weights.
        %   Set this property to false to not add input to specify the
        %   weights. The default value of this property is false.
        WeightsInputPort = false;
        %EnablePolarization  Enable polarization simulation
        %   Set this property to true to enable polarization. Set this
        %   property to false to ignore polarization. The default value of
        %   this property is false. This property applies when the array
        %   specified in the SensorArray property is capable of simulating
        %   polarization.
        EnablePolarization = false
    end

    properties (Hidden, Nontunable, Access={?phased.ArrayResponse,...
            ?phased.internal.AbstractSubarray})
        %PropagationSpeedSource     Source of propagation speed
        %   Specify the source of propagation speed as one of 'Property' |
        %   'Input port', where the default is 'Property'. 
        PropagationSpeedSource = 'Property'
    end
    
    properties (Constant,Hidden)
        PropagationSpeedSourceSet = dsp.CommonSets.getSet(...
            'PropertyOrInputPort');
    end
    
    properties (Access = private, Nontunable)
        cSteeringVector;
        pDOF;
        pNumInputAngles;
        pNumFreqs;
    end
    
    properties (Access = private, Nontunable, Logical)
        pIsFreqScalar;
        pIsSingleWeights;
        pNeedSteeringAngle;
        pIsPropagationSpeedViaProperty
    end

    methods

        function obj = ArrayResponse(varargin)
            %ArrayResponse   Construct the ArrayResponse class.
            obj@phased.internal.AbstractArrayOperation(varargin{:});

        end

    end
    
    methods (Access = protected)
        function privValidateSensorArray(obj,val)  %#ok<INUSL>
            validateattributes( val, { 'phased.internal.AbstractArray',...
                'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
        end
    end
    
    methods (Access = protected)
        
        function validatePropertiesImpl(obj)
            cond = obj.EnablePolarization && ...
                    ~isPolarizationCapable(obj.SensorArray);
            if cond
                coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
            end
        end
        
        function setupImpl(obj,freq,ang, weights,~) 
            obj.pIsPropagationSpeedViaProperty = strcmp(...
                obj.PropagationSpeedSource,'Property');
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'IncludeElementResponse',true,...
                'EnablePolarization',obj.EnablePolarization);
            if obj.pIsPropagationSpeedViaProperty
                obj.cSteeringVector.PropagationSpeed = ...
                    obj.PropagationSpeed;
            else
                obj.cSteeringVector.PropagationSpeedSource = ...
                    obj.PropagationSpeedSource;
            end
            obj.pDOF = getDOF(obj.SensorArray);
            sz_freq = size(freq);
            sz_angle = size(ang);
            obj.pNumInputAngles = sz_angle(2);
            obj.pNumFreqs = sz_freq(2);
            obj.pIsFreqScalar = (sz_freq(2)==1);
            if ~obj.WeightsInputPort
                obj.pIsSingleWeights = true;
            else
                obj.pIsSingleWeights = (size(weights,2)==1);
            end
            obj.pNeedSteeringAngle = ...
                isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                ~strncmp(obj.SensorArray.SubarraySteering,'None',1);
        end
        
        function flag = isInputComplexityLockedImpl(obj,index) 
            flag = true;
            if obj.WeightsInputPort && (index == 3)
                flag = false;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function resetImpl(obj)
            reset(obj.cSteeringVector);
        end
        
        function releaseImpl(obj)
            release(obj.cSteeringVector);
        end

        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if strcmp(obj.PropagationSpeedSource,'Input port') && ...
                    strcmp(prop,'PropagationSpeed')
                flag = true;
            end
        end
        
        function validateInputsImpl(obj,freq,angle,weights,stang)
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
                coder.internal.errorIf(cond,'phased:system:measure:ArrayResponse:NeedTwoRows','Ang');
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
                cond = ~ismatrix(weights) || isempty(weights);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeMatrix','W');
                end
                sz_w = size(weights);
                cond = ~isa(weights,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','W','double');
                end
                cond = (size(freq,2)~=sz_w(2)) && (sz_w(2)~=1);
                if cond
                    coder.internal.errorIf(cond,'phased:measure:MatrixVectorDimensionMismatch','W','Freq');
                end
                N = getDOF(obj.SensorArray);
                cond = sz_w(1)~=N;
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidRowNumbers','W',N);
                end
            end
            
            if isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                    ~strncmp(obj.SensorArray.SubarraySteering,'None',1)
                if ~obj.WeightsInputPort
                    stang = weights;
                end
                if strncmp(obj.SensorArray.SubarraySteering,'Phase',1) || ...
                        strncmp(obj.SensorArray.SubarraySteering,'Time',1)
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
                    Ns = getNumSubarrays(obj.SensorArray);
                    cond = (~iscell(ws) && ~ismatrix(ws)) || isempty(ws);
                    if cond
                        coder.internal.errorIf(cond, ...
                            'phased:phased:expectedCellOrMatrix','WS');
                    end
                    Nse = zeros(1,Ns);
                    for m = 1:Ns
                        Nse(m) = getNumElements(obj.SensorArray,m);
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
            s = saveObjectImpl@phased.internal.AbstractArrayOperation(obj);
            if isLocked(obj)
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.pDOF = obj.pDOF;
                s.pNumInputAngles = obj.pNumInputAngles;
                s.pIsFreqScalar = obj.pIsFreqScalar;
                s.pIsSingleWeights = obj.pIsSingleWeights;
                s.pNumFreqs = obj.pNumFreqs;
                s.pNeedSteeringAngle = obj.pNeedSteeringAngle;
                s.pIsPropagationSpeedViaProperty = obj.pIsPropagationSpeedViaProperty;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractArrayOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
                    if ~isfield(s,'pIsPropagationSpeedViaProperty')
                        obj.pIsPropagationSpeedViaProperty = true;
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
        
        function resp = stepImpl(obj,freq,ang, weightsArg, stangArg, c)
            
            if ~obj.WeightsInputPort
                if obj.pNeedSteeringAngle
                    stang = weightsArg;
                    if ~obj.pIsPropagationSpeedViaProperty
                        pspeed = stangArg;
                    end
                else
                    if ~obj.pIsPropagationSpeedViaProperty
                        pspeed = weightsArg;
                    end
                end
                weights = ones(obj.pDOF,1);
            else
                weights = weightsArg;
                if obj.pNeedSteeringAngle
                    stang = stangArg;
                    if ~obj.pIsPropagationSpeedViaProperty
                        pspeed = c;
                    end
                else
                    if ~obj.pIsPropagationSpeedViaProperty
                        pspeed = stangArg;
                    end
                end
            end
            
            if obj.EnablePolarization
                if obj.pNeedSteeringAngle
                    if obj.pIsPropagationSpeedViaProperty
                        sv = step(obj.cSteeringVector,freq,ang,stang);
                    else
                        sv = step(obj.cSteeringVector,freq,ang,stang,pspeed);
                    end
                else
                    if obj.pIsPropagationSpeedViaProperty
                        sv = step(obj.cSteeringVector,freq,ang);
                    else
                        sv = step(obj.cSteeringVector,freq,ang,pspeed);
                    end
                end
                
                sv_h = sv.H;
                sv_v = sv.V;

                if obj.pIsFreqScalar
                    resp_h = (weights'*sv_h).';
                    resp_v = (weights'*sv_v).';
                else
                    if obj.pIsSingleWeights
                        resp_h = reshape(weights'*...
                            reshape(sv_h,obj.pDOF,[]),...
                            obj.pNumInputAngles,[]);
                        resp_v = reshape(weights'*...
                            reshape(sv_v,obj.pDOF,[]),...
                            obj.pNumInputAngles,[]);
                    else
                        resp_h = ...
                            coder.nullcopy(complex(zeros(obj.pNumInputAngles,obj.pNumFreqs)));
                        resp_v = coder.nullcopy(resp_h);
                        
                        for m = obj.pNumFreqs:-1:1
                            resp_h(:,m) = weights(:,m)'*sv_h(:,:,m);
                            resp_v(:,m) = weights(:,m)'*sv_v(:,:,m);
                        end
                    end

                end
                resp.H = resp_h;
                resp.V = resp_v;
            else
                if obj.pNeedSteeringAngle
                    if obj.pIsPropagationSpeedViaProperty
                        sv = step(obj.cSteeringVector,freq,ang,stang);
                    else
                        sv = step(obj.cSteeringVector,freq,ang,stang,pspeed);
                    end
                else
                    if obj.pIsPropagationSpeedViaProperty
                        sv = step(obj.cSteeringVector,freq,ang);
                    else
                        sv = step(obj.cSteeringVector,freq,ang,pspeed);
                    end
                end

                if obj.pIsFreqScalar
                    resp = (weights'*sv).';
                else
                    if obj.pIsSingleWeights
                        resp = reshape(weights'*...
                            reshape(sv,obj.pDOF,[]),...
                            obj.pNumInputAngles,[]);
                    else
                        resp = ...
                            coder.nullcopy(complex(zeros(obj.pNumInputAngles,obj.pNumFreqs)));
                        for m = obj.pNumFreqs:-1:1
                            resp(:,m) = weights(:,m)'*sv(:,:,m);
                        end
                    end

                end
            end
            
        end
        
        function num = getNumInputsImpl(obj) 
            num = 2;
            if obj.WeightsInputPort
                num = num+1;
            end
            if isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                    ~strncmp(obj.SensorArray.SubarraySteering,'None',1)
                num = num+1;
            end
            if strcmp(obj.PropagationSpeedSource,'Input port')
                num = num+1;
            end
        end
        
    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractArrayOperation('subarray');
        dEnablePolarization = ...
                matlab.system.display.internal.Property('EnablePolarization', ...
                                                        'IsGraphical', false);
        props = {...
          'WeightsInputPort',...
          dEnablePolarization};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:ArrayResponseTitle')),...
              'Text',getString(message('phased:library:block:ArrayResponseDesc')));
      end
    end
    methods (Access = protected)
        function varargout = getInputNamesImpl(obj) 
            varargout = {'Freq','Ang','W','Steer'};
            if ~obj.WeightsInputPort
                varargout(3) = {'Steer'};
            end
        end
        
        function varargout = getOutputNamesImpl(obj) %#ok<MANU>
            varargout = {''};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Array\nResponse');
        end
        
        function varargout = getOutputSizeImpl(obj)
            szFreq = propagatedInputSize(obj,1);
            szAng = propagatedInputSize(obj,2);
            varargout{1} = [szAng(2) szFreq(2)];
        end
        function varargout = isOutputFixedSizeImpl(~)
            varargout{1} = true;
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout{1} = 'double';
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = true;
        end
    end    
    
end

