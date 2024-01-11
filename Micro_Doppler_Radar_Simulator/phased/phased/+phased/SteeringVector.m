classdef (Sealed,StrictDefaults) SteeringVector < phased.internal.AbstractArrayOperation & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%SteeringVector   Sensor array steering vector
%   H = phased.SteeringVector creates a steering vector System object, H.
%   This object calculates the steering vector of a sensor array for the
%   specified directions. By default a 2-element uniform linear array (ULA)
%   is used.
%
%   H = phased.SteeringVector(Name,Value) creates a steering vector object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   SV = step(H,FREQ,ANG) returns the steering vector, SV, of the array for
%   the directions specified in ANG (in degrees), at the given operating
%   frequency FREQ (in Hz). FREQ is a row vector of length L and ANG can be
%   either a row vector of length M or a 2xM matrix. When ANG is a 2xM
%   matrix, each column of the matrix specifies the direction in space in
%   the [azimuth; elevation] form. The azimuth angle should be between
%   [-180 180] degrees and the elevation angle should be between [-90 90]
%   degrees. If ANG is a length M row vector, each element specifies a
%   direction's azimuth angle and the corresponding elevation angle is
%   assumed to be 0.
%
%   The dimensions of SV are NxMxL where N is the number of subarrays if
%   SensorArray contains subarrays, or the number of elements otherwise.
%   Each column of SV contains the steering vector of the array for the
%   corresponding directions specified in ANG. Each page of SV contains
%   the steering vectors of the array for the given frequency specified in
%   FREQ.
%
%   If you set the IncludeElementResponse property to true, the resulting
%   steering vector SV includes the individual element responses. If you
%   set the IncludeElementResponse property to false, the elements are
%   assumed to be isotropic so the steering vector SV does not include the
%   individual element responses.
%
%   When you set the EnablePolarization property to true, the resulting SV
%   is a structure containing two fields, H and V. H represents the array's
%   response in horizontal polarization and V represents the array's
%   response in vertical polarization. Each field is an NxMxL array whose
%   columns contain the steering vectors of the sensor array for the
%   corresponding directions and frequencies. If you set the
%   EnablePolarization to false, then the polarization information is
%   discarded and the combined pattern from both H and V polarizations is
%   used at each element to compute the steering vector. This syntax is
%   only applicable when the sensor array is capable of simulating
%   polarization and when you set the IncludeElementResponse property to
%   true.
%
%   SV = step(H,FREQ,ANG,STEER) uses STEER as the subarray steering angle
%   (in degrees). STEER can be a scalar or a length-2 column vector. If
%   STEER is a vector, it is in the form of [AzimuthAngle; ElevationAngle].
%   If STEER is a scalar, it represents the azimuth angle and the elevation
%   angle is assumed to be 0. This syntax is only applicable when you use
%   subarrays in the SensorArray property, set the SubarraySteering
%   property in the SensorArray to either 'Phase' or 'Time' and set the
%   IncludeElementResponse property to true.
%
%   SV = step(H,FREQ,ANG,WS) uses WS as the weights applied to each element
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
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   SteeringVector methods:
%
%   step     - Calculate the steering vector (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a steering vector object with same property values
%   isLocked - Locked status (logical)
%
%   SteeringVector properties:
%
%   SensorArray            - Sensor array 
%   PropagationSpeed       - Signal propagation speed 
%   IncludeElementResponse - Include element response in steering vector
%   NumPhaseShifterBits    - Number of bits in phase shifters
%   EnablePolarization     - Enable polarization simulation
%
%   % Example:
%   %   Calculate the steering vector for a 4-element uniform linear array 
%   %   at the direction of 30 degrees azimuth and 20 degrees elevation. 
%   %   Assume the array's operating frequency is 300 MHz. Compare the beam
%   %   pattern before and after the steering.
%
%   array = phased.ULA(4);
%   steervector = phased.SteeringVector('SensorArray',array);
%   sv = steervector(3e8,[30; 20])
%   c = steervector.PropagationSpeed;
%   subplot(211)
%   pattern(array,3e8,-180:180,0,'PropagationSpeed',c,...
%       'CoordinateSystem','rectangular'); 
%   title('Before steering');
%   subplot(212)
%   pattern(array,3e8,-180:180,0,'PropagationSpeed',c,'Weights',sv,...
%       'CoordinateSystem','rectangular'); 
%   title('After steering');
%
%   See also phased, phased.ElementDelay, phased.ArrayResponse,
%   phased.ArrayGain.

%   Copyright 2009-2017 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen    

    properties (Logical, Nontunable) 
        %IncludeElementResponse Include element response in steering vector
        %   Set this property to true to include individual element
        %   responses in the steering vector. Setting this property to
        %   false assumes that the elements are isotropic so the steering
        %   vector does not include the individual element responses. The
        %   default value of this property is false. 
        %
        %   When IncludeElementResponse is false, the resulting steering
        %   vector is the array factor among the subarrays if the
        %   SensorArray contains subarrays or the array factor among the
        %   elements otherwise.
        IncludeElementResponse = false
    end
    
    properties (Nontunable)
        %NumPhaseShifterBits    Number of bits in phase shifters
        %   Specify the number of bits used in the phase shifter as a
        %   non-negative integer. The default value of this property is 0,
        %   indicating there is no quantization effect in the phase
        %   shifter.
        NumPhaseShifterBits = 0
    end
    
    properties (Logical, Nontunable)    
        %EnablePolarization  Enable polarization simulation
        %   Set this property to true to enable polarization. Set this
        %   property to false to ignore polarization. The default value of
        %   this property is false. This property applies when the array
        %   specified in the SensorArray property is capable of simulating
        %   polarization and you set the IncludeElementResponse property to
        %   true.
        EnablePolarization = false
    end
    
    properties (Hidden, Nontunable, Access={?phased.ArrayResponse,...
            ?phased.internal.AbstractSubarray,?phased.SteeringVector})
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
        cSensorArray;
        pPortPosition;
        pDOF;
        pNumInputAngles;
    end
    
    properties (Access = private, Nontunable, Logical)
        pIsFreqScalar
        pUseSubarray
        pNeedSteeringAngle
        pNeedCustomSteering
        pEnablePolarization
        pIsPropagationSpeedViaProperty
    end
    
    methods

        function obj = SteeringVector(varargin)
            %SteeringVector   Construct the SteeringVector class.
            
            obj@phased.internal.AbstractArrayOperation(varargin{:});

        end
        
    end
    
    methods
        function set.NumPhaseShifterBits(obj,val)
            validateattributes(val,{'double'},...
                {'scalar','integer','nonnegative'},...
                '','NumPhaseShifterBits');
            obj.NumPhaseShifterBits = val;
        end
    end
    
    methods (Access = protected)
        function privValidateSensorArray(obj,val)  %#ok<INUSL>
            validateattributes( val, { 'phased.internal.AbstractArray',...
                'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
        end
    end
    
    methods (Access = protected)
        
        function setupImpl(obj,freq,ang,stang) %#ok<INUSD>

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
            
            obj.pIsPropagationSpeedViaProperty = strcmp(...
                obj.PropagationSpeedSource,'Property');
            %pEnablePolarization can only be assigned once for codegen
            if obj.IncludeElementResponse 
                if obj.EnablePolarization
                    if isPolarizationCapable(obj.cSensorArray)
                        obj.pEnablePolarization = true;
                    else
                        obj.pEnablePolarization = false;
                    end
                else
                    obj.pEnablePolarization = false;
                    hele = getElementHandle(obj.cSensorArray);
                    if iscell(hele)
                        if isa(obj.cSensorArray,'phased.ReplicatedSubarray')
                            for m = 1:numel(hele)
                                % disablePolarization(hele{m});
                                disablePolarization(obj.cSensorArray.Subarray.ElementSet{m});
                            end
                        elseif isa(obj.cSensorArray,'phased.PartitionedArray')
                            for m = 1:numel(hele)
                                % disablePolarization(hele{m});
                                disablePolarization(obj.cSensorArray.Array.ElementSet{m});
                            end
                        else
                            for m = 1:numel(hele)
                                % disablePolarization(hele{m});
                                disablePolarization(obj.cSensorArray.ElementSet{m});
                            end
                        end
                    else
%                         disablePolarization(hele)
                        if isa(obj.cSensorArray,'phased.ReplicatedSubarray')
                            disablePolarization(obj.cSensorArray.Subarray.Element);
                        elseif isa(obj.cSensorArray,'phased.PartitionedArray')
                            disablePolarization(obj.cSensorArray.Array.Element);
                        else
                            disablePolarization(obj.cSensorArray.Element);
                        end
                    end
                end
            else
                obj.pEnablePolarization = false;
            end
            obj.pDOF = getDOF(obj.cSensorArray);
            if isa(obj.cSensorArray,'phased.internal.AbstractArray')
                obj.pUseSubarray = false;
                a = coder.internal.const(getElementPosition(obj.cSensorArray));
                obj.pPortPosition =a;
            else  % Subarray
                obj.pUseSubarray = true;
                obj.pPortPosition = getSubarrayPosition(obj.cSensorArray);
            end
            sz_freq = size(freq);
            sz_angle = size (ang);
            obj.pNumInputAngles = sz_angle(2);
            obj.pIsFreqScalar = (sz_freq(2)==1);
            obj.pNeedSteeringAngle = obj.pUseSubarray && ...
                obj.IncludeElementResponse && ...
                (strncmp(obj.cSensorArray.SubarraySteering,'Phase',1) || ...
                strncmp(obj.cSensorArray.SubarraySteering,'Time',1));
            obj.pNeedCustomSteering = obj.pUseSubarray && ...
                obj.IncludeElementResponse && ...
                strncmp(obj.cSensorArray.SubarraySteering,'Custom',1);
        end
        
        function flag = isInputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = false;
        end
        
        function releaseImpl(obj)
            %release Steering Vector
            if obj.IncludeElementResponse 
                release(obj.cSensorArray);
            end
        end
        
        function resetImpl(obj)
            %reset Steering Vector
            if obj.IncludeElementResponse
                 reset(obj.cSensorArray);
            end    
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if ~obj.IncludeElementResponse && ...
                    strcmp(prop,'EnablePolarization')
                flag = true;
            end
            if strcmp(obj.PropagationSpeedSource,'Input port') && ...
                    strcmp(prop,'PropagationSpeed')
                flag = true;
            end
        end
        
        function validatePropertiesImpl(obj)
            cond = obj.IncludeElementResponse && obj.EnablePolarization && ...
                    ~isPolarizationCapable(obj.SensorArray);
            if cond
                coder.internal.errorIf(cond,'phased:polarization:invalidElementPolarizationSetting');
            end
        end
        
        function validateInputsImpl(obj,freq,ang,stang)
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
            cond = ~ismatrix(ang) || isempty(ang);
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:inputMustBeMatrix','Ang');
            end
            sz_ang = size(ang);
            cond = sz_ang(1) > 2;
            if cond
                coder.internal.errorIf(cond,'phased:system:measure:SteeringVector:NeedTwoRows','Ang');
            end
            cond = ~isa(ang,'double');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','Ang','double');
            end
            cond = ~isreal(ang);
            if cond
                coder.internal.errorIf(cond,'phased:measure:NeedReal', 'Ang');
            end
            if isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                    obj.IncludeElementResponse && ...
                    ~strncmp(obj.SensorArray.SubarraySteering,'None',1)
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
                s.cSensorArray = saveobj(obj.cSensorArray);
                s.pDOF = obj.pDOF;
                s.pPortPosition = obj.pPortPosition;
                s.pNumInputAngles = obj.pNumInputAngles;
                s.pIsFreqScalar = obj.pIsFreqScalar;
                s.pUseSubarray = obj.pUseSubarray;
                s.pNeedSteeringAngle = obj.pNeedSteeringAngle;
                s.pNeedCustomSteering = obj.pNeedCustomSteering;
                s.pEnablePolarization = obj.pEnablePolarization;
                s.pIsPropagationSpeedViaProperty = obj.pIsPropagationSpeedViaProperty;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractArrayOperation(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSensorArray = eval(...
                        sprintf('%s.loadobj(s.cSensorArray)',s.cSensorArray.ClassNameForLoadTimeEval));
                    s = rmfield(s,'cSensorArray');
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
            if isfield(s,'cElementDelay')
                s = rmfield(s,'cElementDelay');
            end
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function sv = stepImpl(obj,freq,angArg,stangArg,c)
            cond = any(freq < 0);
            if cond
                coder.internal.errorIf(cond,'phased:step:expectedNonnegative', 'Freq');
            end
            if size(angArg,1)==1
                ang = [angArg;zeros(size(angArg))];
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
            
            if obj.pNeedSteeringAngle 
                if isscalar(stangArg)
                    stang = [stangArg; 0];
                else
                    stang = stangArg;
                end
                azstang = stang(1,:);
                elstang = stang(2,:);
                cond = any(azstang>180) || any(elstang>90);
                if cond
                    coder.internal.errorIf(cond,'phased:step:AzElNotLessEqual');
                end
                cond = any(azstang<-180) || any(elstang<-90);
                if cond
                    coder.internal.errorIf(cond,'phased:step:AzElNotGreaterEqual');
                end
                if obj.pIsPropagationSpeedViaProperty
                    pspeed = obj.PropagationSpeed;
                else
                    pspeed = c;
                end
            elseif obj.pNeedCustomSteering
                w = stangArg; % weights
                if iscell(w)
                    for m = 1:numel(w)
                        cond = ~all(isfinite(w{m}));
                        if cond        
                            coder.internal.errorIf(cond,'phased:step:expectedFinite','WS');
                        end
                    end
                else
                    cond = ~all(isfinite(w(:)));
                    if cond        
                        coder.internal.errorIf(cond,'phased:step:expectedFinite','WS');
                    end
                end
                if obj.pIsPropagationSpeedViaProperty
                    pspeed = obj.PropagationSpeed;
                else
                    pspeed = c;
                end                
            else
                if obj.pIsPropagationSpeedViaProperty
                    pspeed = obj.PropagationSpeed;
                else
                    pspeed = stangArg;
                end
            end

            sv_temp = phased.internal.steeringvec(obj.pPortPosition,...
                freq,pspeed,ang,obj.NumPhaseShifterBits);
            if obj.pEnablePolarization
                if obj.pUseSubarray
                    if obj.pNeedSteeringAngle
                        elempattern = step(obj.cSensorArray,freq,ang,...
                            pspeed,stang);
                    elseif obj.pNeedCustomSteering
                        elempattern = step(obj.cSensorArray,freq,ang,...
                            pspeed,w);
                    else
                        elempattern = step(obj.cSensorArray,freq,ang,...
                            pspeed);
                    end
                else
                    elempattern = ...
                        step(obj.cSensorArray,freq,ang);
                end
                sv.H = sv_temp.*elempattern.H;
                sv.V = sv_temp.*elempattern.V;
            elseif obj.IncludeElementResponse
                if obj.pUseSubarray
                    if obj.pNeedSteeringAngle
                        elempattern = step(obj.cSensorArray,freq,ang,...
                            pspeed,stang);
                    elseif obj.pNeedCustomSteering
                        elempattern = step(obj.cSensorArray,freq,ang,...
                            pspeed,w);
                    else
                        elempattern = step(obj.cSensorArray,freq,ang,...
                            pspeed);
                    end
                else
                    elempattern = ...
                        step(obj.cSensorArray,freq,ang);
                end
                sv = sv_temp.*elempattern;
            else
                sv = sv_temp;
            end
        end
        
        function num = getNumInputsImpl(obj) 
            if isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                    obj.IncludeElementResponse && ...
                    ~strncmp(obj.SensorArray.SubarraySteering,'None',1)
                num = 3;
            
            else
                num = 2;
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
          'IncludeElementResponse',...
          'NumPhaseShifterBits',...
          dEnablePolarization};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:SteeringVectorTitle')),...
              'Text',getString(message('phased:library:block:SteeringVectorDesc')));
      end
    end
    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)  %#ok<MANU>
            varargout = {'Freq','Ang','Steer'};
        end
        
        function varargout = getOutputNamesImpl(obj) %#ok<MANU>
            varargout = {''};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Steering\nVector');
        end
        
        function varargout = getOutputSizeImpl(obj)
            szFreq = propagatedInputSize(obj,1);
            szAng = propagatedInputSize(obj,2);
            numDOF = getDOF(obj.SensorArray);
            varargout{1} = [numDOF szAng(2) szFreq(2)];
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj,1) && ...
                propagatedInputFixedSize(obj,2);
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout{1} = 'double';
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = true;
        end
    end    
end

