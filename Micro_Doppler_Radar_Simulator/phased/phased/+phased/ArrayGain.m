classdef (Sealed,StrictDefaults) ArrayGain < phased.internal.AbstractArrayOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%ArrayGain   Sensor array gain
%   H = phased.ArrayGain creates an array gain System object, H. This
%   object calculates the array gain of a sensor array for the specified
%   directions. By default a 2-element uniform linear array (ULA) is used.
%
%   H = phased.ArrayGain(Name,Value) creates an array gain object, H, with
%   the specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN). 
%
%   Step method syntax:
%
%   G = step(H,FREQ,ANG) returns the gain (in dB), G, of the array for the
%   directions specified in ANG (in degrees), at the given operating
%   frequency FREQ (in Hz). FREQ is a row vector of length L and ANG can be
%   either a row vector of length M or a 2xM matrix. G is an MxL matrix
%   whose elements contain the gains of the sensor array at angles
%   specified in ANG and frequencies specified in FREQ.
%
%   When ANG is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in the [azimuth; elevation] form. The azimuth
%   angle should be between [-180 180] degrees and the elevation angle
%   should be between [-90 90] degrees. If ANG is a length M row vector,
%   each element specifies a direction's azimuth angle and the
%   corresponding elevation angle is assumed to be 0.
%
%   The array gain is defined as the signal to noise ratio (SNR)
%   improvement between the array output and the individual channel input,
%   assuming the noise is spatially white. It is related to the array
%   response but is not the same. If a rectangular taper is used in the
%   array, then the array gain is the square of the array response
%   normalized by the number of elements in the array.
%
%   G = step(H,FREQ,ANG,W) uses W as the weights applied on the sensor
%   array when you set the WeightsInputPort property to true. W can be a
%   length-N column vector or an NxL matrix where N is the number of
%   subarrays if SensorArray contains subarrays, or the number of elements
%   otherwise. L is the number of frequencies specified in FREQ. If W is a
%   vector, the weights are applied at all frequencies. If W is a matrix,
%   each column of W represents the weights used at the corresponding
%   frequency specified in FREQ.
%
%   G = step(H,FREQ,ANG,STEER) uses STEER as the subarray steering angle
%   (in degrees). STEER can be a scalar or a length-2 column vector. If
%   STEER is a vector, it is in the form of [AzimuthAngle; ElevationAngle].
%   If STEER is a scalar, it represents the azimuth angle and the elevation
%   angle is assumed to be 0. This syntax is only applicable when you use
%   subarrays in the SensorArray property and set the SubarraySteering
%   property in the SensorArray to either 'Phase' or 'Time'.
%
%   G = step(H,FREQ,ANG,WS) uses WS as the weights applied to each element
%   in the subarray. WS can be either a matrix or a cell array. This syntax
%   is only applicable when you use subarrays in the SensorArray property,
%   set the SubarraySteering property in the SensorArray to 'Custom' and
%   set the IncludeElementResponse property to true.
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
%   G = step(H,FREQ,ANG,W,STEER)
%
%   or
%
%   G = step(H,FREQ,ANG,W,WS)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   ArrayGain methods:
%
%   step     - Calculate the array gain of the sensor array (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create an array gain object with same property values
%   isLocked - Locked status (logical)
%
%   ArrayGain properties:
%
%   SensorArray      - Sensor array 
%   PropagationSpeed - Signal propagation speed 
%   WeightsInputPort - Add input to specify the weights
%
%   % Examples:
%
%   % Example 1:
%   %   Calculate and plot the array gain for a 4-element ULA between -90 
%   %   and 90 degrees in azimuth.  
%   %   Assume the array's operating frequency is 300 MHz.
%
%   ha = phased.ULA(4);
%   gain = phased.ArrayGain('SensorArray',ha);
%   fc = 300e6; ang = -90:90;
%   g = gain(fc,ang);
%   plot(ang,g); title('Array Gain'); 
%   xlabel('Angle (degrees)'); ylabel('Array Gain (dB)');
%
%   % Example 2:
%   %   Calculate the array gain for the above array at the direction of
%   %   30 degrees azimuth and 20 degrees elevation. Assuming a Hamming 
%   %   taper is used.
%
%   ha = phased.ULA(4);
%   gain = phased.ArrayGain('SensorArray',ha,'WeightsInputPort',true);
%   fc = 300e6; ang = [30;20]; w = hamming(4);
%   g = gain(fc,ang,w)
%
%   See also phased, phased.ElementDelay, phased.SteeringVector,
%   phased.ArrayResponse.

%   Copyright 2010-2017 The MathWorks, Inc.

%   Reference
%   [1] J. R. Guerci, Space-Time Adaptive Processing for Radar, Artech
%       House, 2003
%   [2] Harry Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable, Logical) 

        %WeightsInputPort    Add input to specify the weights
        %   Set this property to true to add input to specify the weights.
        %   Set this property to false to not add input to specify the
        %   weights. The default value of this property is false.
        WeightsInputPort = false;
    end
    
    properties (Access = private, Nontunable)
        cArrayResponse
        pDOF
        pNumFreqs
    end
    
    properties (Access = private, Nontunable, Logical)
        pNeedSteeringAngle
        pIsSingleWeights
    end

    methods

        function obj = ArrayGain(varargin)
            %ArrayGain   Construct the ArrayGain class.
            obj@phased.internal.AbstractArrayOperation(varargin{:});

        end

    end
    
    methods (Access = protected)
        function privValidateSensorArray(~,val) 
            validateattributes( val, { 'phased.internal.AbstractArray',...
                'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
        end
    end
    
    methods (Access = protected)
        
        function setupImpl(obj,freq,~, weights, ~) 
            obj.cArrayResponse = phased.ArrayResponse(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'WeightsInputPort',true,...
                'EnablePolarization',false);
            obj.pNeedSteeringAngle = ...
                isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                ~strncmp(obj.SensorArray.SubarraySteering,'None',1);
            obj.pNumFreqs = size(freq,2);
            obj.pDOF = getDOF(obj.SensorArray);

            if ~obj.WeightsInputPort
                obj.pIsSingleWeights = true;
            else
                obj.pIsSingleWeights = (size(weights,2)==1);
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
        end
        
        function resetImpl(obj)
            reset(obj.cArrayResponse);
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
                N = getDOF(obj.SensorArray);
                cond = sz_w(1)~=N;
                if cond
                    coder.internal.errorIf(cond,'phased:phased:invalidRowNumbers','W',N);
                end
            end
            
            if isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                    ~strncmp(obj.SensorArray.SubarraySteering,'None',1)
                if ~obj.WeightsInputPort
                    stang = w;
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
                s.cArrayResponse = saveobj(obj.cArrayResponse);
                s.pNeedSteeringAngle = obj.pNeedSteeringAngle;
                s.pDOF = obj.pDOF;
                s.pIsSingleWeights = obj.pIsSingleWeights;
                s.pNumFreqs = obj.pNumFreqs;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractArrayOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cArrayResponse = phased.ArrayResponse.loadobj(s.cArrayResponse);
                    s = rmfield(s,'cArrayResponse');
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
        
        function g = stepImpl(obj,freq,ang, weightsArg, stang)

            if ~obj.WeightsInputPort
                if obj.pNeedSteeringAngle
                    stang = weightsArg;
                end
                weights = ones(obj.pDOF,1);
            else
                weights = weightsArg;
            end
               
            if obj.pNeedSteeringAngle
                resp = step(obj.cArrayResponse,freq,ang,weights,stang);
            else
                resp = step(obj.cArrayResponse,freq,ang,weights);
            end
            g = zeros(size(resp));
            if obj.pIsSingleWeights
                for m = obj.pNumFreqs:-1:1
                    if obj.pNeedSteeringAngle
                        Rn = getNoiseGainMatrix(...
                            obj.SensorArray,freq(m),obj.PropagationSpeed,stang);
                    else
                        Rn = getNoiseGainMatrix(...
                            obj.SensorArray,freq(m),obj.PropagationSpeed);
                    end
                    weightnorm = real(weights'*Rn*weights);

                    if weightnorm
                        g(:,m) = abs(resp(:,m)).^2/weightnorm;  % ensure real
                    end
                end
            else
                for m = obj.pNumFreqs:-1:1
                    temp_w = weights(:,m);
                    if obj.pNeedSteeringAngle
                        Rn = getNoiseGainMatrix(...
                            obj.SensorArray,freq(m),obj.PropagationSpeed,stang);
                    else
                        Rn = getNoiseGainMatrix(...
                            obj.SensorArray,freq(m),obj.PropagationSpeed);
                    end
                    weightnorm = real(temp_w'*Rn*temp_w);

                    if weightnorm
                        g(:,m) = abs(resp(:,m)).^2/weightnorm;
                    end
                end
            end
            g = pow2db(g); % return dB
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
        end
        
    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractArrayOperation('subarray');
        props = 'WeightsInputPort';
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:ArrayGainTitle')),...
              'Text',getString(message('phased:library:block:ArrayGainDesc')));
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
            str = sprintf('Array Gain');
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
            varargout{1} = false;
        end
    end    

end
