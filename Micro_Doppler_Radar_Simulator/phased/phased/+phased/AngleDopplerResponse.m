classdef (Sealed,StrictDefaults) AngleDopplerResponse < phased.internal.AbstractNarrowbandArrayProcessing & ...
        matlab.system.mixin.CustomIcon & ...
        matlab.system.mixin.Propagates
%AngleDopplerResponse   Angle-Doppler response
%   H = phased.AngleDopplerResponse creates an angle-Doppler response
%   System object, H. This object calculates the angle-Doppler response of
%   the input data.
%
%   H = phased.AngleDopplerResponse(Name,Value) creates an angle-Doppler
%   response object, H, with the specified property Name set to the
%   specified Value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1,...,NameN,ValueN).
%
%   AngleDopplerResponse object generates the response using a conventional
%   beamformer and an FFT-based Doppler filter.
%
%   Step method syntax:
%   
%   [RESP,ANG,DOP] = step(H,X) calculates the angle-Doppler
%   response of the input data X.
%
%   X must be a matrix or a column vector. If X is a matrix, its number of
%   rows must equal the number of subarrays if SensorArray contains
%   subarrays, or the number of elements otherwise. If X is a vector, the
%   number of rows must be an integer multiple of the number of
%   elements/subarrays of the array specified in the SensorArray property.
%   In addition, the multiplier must be at least 2.
%
%   RESP is a PxQ matrix containing the complex angle-Doppler response of
%   the input X, where P is determined by the NumDopplerSamples property
%   and Q is determined by the NumAngleSamples property. ANG is a
%   length Q column vector containing the angle samples at which the
%   angle-Doppler response is evaluated and DOP is a length P column
%   vector containing the Doppler samples at which the angle-Doppler
%   response is evaluated.
%
%   [RESP,ANG,DOP] = step(H,X,PRF) uses input PRF as the pulse repetition
%   frequency (in Hz) to calculate the angle-Doppler response when the
%   PRFSource property is set to 'Input port'.
%
%   [RESP,ANG,DOP] = step(H,X,EL) uses input EL as the
%   elevation angle (in degrees) to calculate the angle-Doppler response
%   when the ElevationAngleSource property is set to 'Input port'.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   AngleDopplerResponse methods:
%
%   step         - Calculate angle-Doppler response
%   release      - Allow property value and input characteristics changes
%   clone        - Create angle-Doppler response object with same property
%                  values
%   isLocked     - Locked status (logical)
%   plotResponse - Plot angle-Doppler response
%
%   AngleDopplerResponse properties:
%
%   SensorArray          - Sensor array
%   PropagationSpeed     - Signal propagation speed
%   OperatingFrequency   - Operating frequency
%   PRFSource            - Source of PRF
%   PRF                  - Pulse repetition frequency
%   ElevationAngleSource - Source of elevation angle
%   ElevationAngle       - Elevation angle
%   NumAngleSamples      - Number angular bins
%   NumDopplerSamples    - Number Doppler bins
%
%   % Example:
%   %   Calculate the angle Doppler response of the 190th cell of a
%   %   collected data cube.
%
%   load STAPExampleData;
%   x = shiftdim(STAPEx_ReceivePulse(190,:,:));
%   response = phased.AngleDopplerResponse(...
%               'SensorArray',STAPEx_HArray,...
%               'OperatingFrequency',STAPEx_OperatingFrequency,...
%               'PropagationSpeed',STAPEx_PropagationSpeed,...
%               'PRF',STAPEx_PRF);
%   [resp,ang_grid,dop_grid] = response(x);
%   contour(ang_grid,dop_grid,abs(resp)),xlabel('Angle'),ylabel('Doppler');
%
%   See also phased, phased.STAPSMIBeamformer, phased.DPCACanceller,
%   phased.ADPCACanceller.

%   Copyright 2010-2018 The MathWorks, Inc.

%   Reference
%   [1] J. R. Guerci, Space-Time Adaptive Processing for Radar, Artech
%       House, 2003


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        
        %PRFSource    Source of PRF
        %   Specify how to determine the PRF for the STAP processor as one
        %   of 'Property' | 'Input port', where the default is 'Property'.
        %   When you set this property to 'Property', the PRF is determined
        %   by the value of the PRF property. When you set this property to
        %   'Input port', the PRF is determined by the input argument.
        PRFSource = 'Property';
        %PRF     Pulse repetition frequency (Hz)
        %   Specify the pulse repetition frequency (PRF) (in Hz) of the
        %   input signal as a positive scalar. The default value of this
        %   property is 1.
        PRF = 1;
        %ElevationAngleSource   Source of elevation angle
        %   Specify how to determine the elevation angle when calculating
        %   the angle-Doppler response as one of 'Property' | 'Input port',
        %   where the default is 'Property'. When you set this property to
        %   'Property', the elevation angle is determined by the
        %   ElevationAngle property. When you set this property to 'Input
        %   port', the elevation angle is determined by the input argument.
        ElevationAngleSource = 'Property';
        %ElevationAngle     Elevation angle (deg)
        %   Specify the elevation angle (in degrees) used to calculate the
        %   angle-Doppler response as a scalar. The angle must be between
        %   -90 and 90. This property applies when you set the
        %   ElevationAngleSource property to 'Property'. The default value
        %   of this property is 0.
        ElevationAngle = 0;
        %NumAngleSamples    Number of angle bins
        %   Specify the number of samples in angular domain used to
        %   calculate the angle-Doppler response as a positive integer.
        %   This value must be greater than 2. The default value of this
        %   property is 256.
        NumAngleSamples = 256;
        %NumDopplerSamples    Number of Doppler bins
        %   Specify the number of samples in Doppler domain used to
        %   calculate the angle-Doppler response as a positive integer.
        %   This value must be greater than 2. The default value of this
        %   property is 256.
        NumDopplerSamples = 256;
    end
    
    properties(Constant, Hidden)
        PRFSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
        ElevationAngleSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    properties(Access = private, Nontunable)
        pNumPulses
        cSteeringVector
    end
    
    methods
        function set.PRF(obj,val)
            sigdatatypes.validateFrequency(val,'phased.AngleDopplerResponse',...
                'PRF',{'scalar'});
            obj.PRF = val;
        end
        function set.ElevationAngle(obj,val)
            sigdatatypes.validateAngle(val,'phased.AngleDopplerResponse',...
                'ElevationAngle',{'scalar','>=',-90,'<=',90});
            obj.ElevationAngle = val;
        end
        function set.NumAngleSamples(obj,val)
            sigdatatypes.validateIndex(val,'phased.AngleDopplerResponse',...
                'NumAngleSamples',{'scalar','>=',2});
            obj.NumAngleSamples = val;
        end
        function set.NumDopplerSamples(obj,val)
            sigdatatypes.validateIndex(val,'phased.AngleDopplerResponse',...
                'NumDopplerSamples',{'scalar','>=',2});
            obj.NumDopplerSamples = val;
        end
    end
    
    methods

        function obj = AngleDopplerResponse(varargin)
            
            obj@phased.internal.AbstractNarrowbandArrayProcessing(varargin{:});

        end

        function varargout = plotResponse(obj,x,varargin)
        %plotResponse   Plot angle-Doppler response
        %   plotResponse(H,X) plots the angle-Doppler response in dB scale.
        %
        %   plotResponse(H,X,PRF) plots the angle-Doppler response
        %   calculated using the specified pulse repetition frequency (in
        %   Hz) PRF. This option is available when you set the PRFSource
        %   property to 'Input port'
        %
        %   plotResponse(H,X,EL) plots the angle-Doppler response
        %   calculated using the specified elevation angle (in degrees)
        %   EL. This option is available when you set the
        %   ElevationAngleSource property to 'Input port'.
        %
        %   plotResponse(...,Name,Value) plots the angle-Doppler response
        %   with the specified parameter Name set to the specified value.
        %   The parameter Names are
        %                   Unit: The unit of the plot, using one of 
        %                         | 'db' | 'mag' | 'pow' |. The default 
        %                         value is 'db'.
        %       NormalizeDoppler: Set this to true to normalize the Doppler
        %                         frequency. Set this to false to not
        %                         normalize the Doppler frequency. The
        %                         default value is false.
        %
        %   % Example:
        %   %   Plot the angle-Doppler response of 190th Cell of a
        %   %   collected data cube.
        %
        %   load STAPExampleData;
        %   x = shiftdim(STAPEx_ReceivePulse(190,:,:));
        %   hadresp = phased.AngleDopplerResponse(...
        %               'SensorArray',STAPEx_HArray,...
        %               'OperatingFrequency',STAPEx_OperatingFrequency,...
        %               'PropagationSpeed',STAPEx_PropagationSpeed,...
        %               'PRF',STAPEx_PRF);
        %   plotResponse(hadresp,x,'NormalizeDoppler',true);
        %
        %   See also phased.AngleDopplerResponse, phased.ULA/pattern, 
        %   phased.URA/pattern, phased.ConformalArray/pattern.
        
            if ~isempty(coder.target)
                coder.internal.assert(false, ...
                                      'phased:Waveform:CodegenNotSupported','plotResponse');
            end

            narginchk(2,inf);

            hadresp = clone(obj);
            release(hadresp);
            if (obj.PRFSource(1) == 'P')
                prf = obj.PRF;
                if (obj.ElevationAngleSource(1) == 'P')
                    [resp,ang_grid,dop_grid] = step(hadresp,x);
                else
                    if isempty(varargin)
                        error(message('phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                            'El','ElevationAngleSource','Input port'));
                    end
                    elang = varargin{1};
                    sigdatatypes.validateAngle(elang,'plotResponse','El',...
                        {'scalar','>=',-90,'<=',90});
                    varargin(1) = [];
                    [resp,ang_grid,dop_grid] = step(hadresp,x,elang);
                    
                end
                
            else
                prf = varargin{1};
                sigdatatypes.validateFrequency(prf,'','PRF',...
                    {'scalar'});
                if (obj.ElevationAngleSource(1) == 'P')
                    varargin(1)=[];
                    [resp,ang_grid,dop_grid] = step(hadresp,x,prf);
                else
                    if isempty(varargin)
                        error(message('phased:AngleDopplerResponse:plotResponse:MissingParameter',...
                            'El','ElevationAngleSource','Input port'));
                    end
                    elang = varargin{2};
                    sigdatatypes.validateAngle(elang,'plotResponse','El',...
                        {'scalar','>=',-90,'<=',90});
                    varargin(2) = [];
                    varargin(1)=[];
                    [resp,ang_grid,dop_grid] = step(hadresp,x,prf,elang);
                    
                end
            end
            
            unit = 'db';
            normalizedoppler = false;
            sigutils.pvparse(varargin{:});
            unit = validatestring(unit,{'db','mag','pow'},...
                'plotResponse','Unit');
            validateattributes(normalizedoppler,{'logical'},{'scalar'},...
                'plotResponse','NormalizeDoppler');
            
            hresp = phased.internal.AngleDopplerPattern(...
                'Doppler',dop_grid(:)',...
                'Angle',ang_grid(:)',...
                'Pattern',abs(resp),...
                'PRF',prf);
            
            if nargout
                varargout{1} = ...
                    plot(hresp,'Units',unit,'NormalizeDoppler',normalizedoppler);
            else
                plot(hresp,'Units',unit,'NormalizeDoppler',normalizedoppler);
            end
        
        end
    end
    
    methods(Access = protected)
        function privValidateSensorArray(~,val) 
            validateattributes( val, { 'phased.internal.AbstractArray',...
                'phased.internal.AbstractSubarray'}, { 'scalar' }, '', 'SensorArray');
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if (obj.ElevationAngleSource(1) == 'I') && ...
                    strcmp(prop, 'ElevationAngle')
                flag = true;
            end
            if(obj.PRFSource(1) == 'I'&& strcmp(prop, 'PRF'))
                flag = true;
            end
        end
        
        function num = getNumInputsImpl(obj)
            if (obj.PRFSource(1) == 'P')
                if (obj.ElevationAngleSource(1) == 'P')
                    num = 1;
                else
                    num = 2;
                end
            else
                if (obj.ElevationAngleSource(1) == 'P')
                    num = 2;
                else
                    num = 3;
                end
            end
        end
        
        function num = getNumOutputsImpl(obj) %#ok<MANU>
            num = 3;
        end
        
        function validateInputsImpl(obj,x,elang,stang) %#ok<INUSD>
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
            sz_x = size(x);
            N = getDOF(obj.SensorArray);
            if sz_x(2) == 1
                temp = sz_x(1)/N;
                cond = rem(temp,1) || temp < 2;
                if cond
                    coder.internal.errorIf(cond,'phased:AngleDopplerResponse:InvalidWeightDimension', N, 2*N);
                end
            else
                cond = sz_x(1) ~= N;
                if cond
                    coder.internal.errorIf(cond,'phased:AngleDopplerResponse:InvalidDataDimension', N);
                end
            end
            if (obj.ElevationAngleSource(1) == 'I')
                cond = ~isa(elang,'double');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','El','double');
                end
                cond = ~isscalar(elang);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeScalar','El');
                end
                cond = ~isreal(elang);
                if cond
                    coder.internal.errorIf(cond,'phased:AngleDopplerResponse:NeedReal','El');
                end
            end
        end
        
        function setupImpl(obj,x,varargin) 
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed);
            sz_x = size(x);
            N = getDOF(obj.SensorArray);
            if sz_x(2) == 1 
                obj.pNumPulses = sz_x(1)/N;
            else
                obj.pNumPulses = sz_x(2);
            end
        end
        
        function flag = isInputComplexityLockedImpl(~,index)
            flag = true;  % index == 2 || index == 3
            if index == 1  
                flag = false;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(~,index) 
            if index == 1
                flag = false;
            else % (index == 2) || (index == 3)
                flag = true;
            end
        end
        
        function releaseImpl(obj)
            release(obj.cSteeringVector);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractNarrowbandArrayProcessing(obj);
            if isLocked(obj)
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.pNumPulses = obj.pNumPulses;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractNarrowbandArrayProcessing(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
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

        function [respC,ang_gridC,dop_grid] = stepImpl(obj,xArg,varargin)
            if (obj.PRFSource(1) == 'P')
                prf = obj.PRF;
                if (obj.ElevationAngleSource(1) == 'P')
                    elang = obj.ElevationAngle;
                else
                    elang = varargin{1};
                end
            else
                prf = varargin{1};
                sigdatatypes.validateFrequency(prf,'','PRF',...
                    {'scalar'});
                if (obj.ElevationAngleSource(1) == 'P')
                    elang = obj.ElevationAngle;
                else
                    elang = varargin{2};
                end
            end
            cond = (elang > 90) || (elang < -90);
            if cond
                coder.internal.errorIf(cond,'phased:AngleDopplerResponse:Step:OutOfBoundElAng','El');
            end
            x = reshape(xArg,[],obj.pNumPulses);
            
            ang_grid = linspace(-90,90,obj.NumAngleSamples);
            ang_stv = step(obj.cSteeringVector,obj.OperatingFrequency,...
                [ang_grid; elang*ones(size(ang_grid))]);

            ang_gridC = ang_grid.';  % output column vector for angle grid
            
            maxdop = prf/2;
            dop_grid = linspace(-maxdop,maxdop,obj.NumDopplerSamples).';
            dop_stv = dopsteeringvec(dop_grid.',obj.pNumPulses,prf);
            
            resp = ang_stv'*x*conj(dop_stv);
            respC = resp.';  % make angle in x axis and Doppler in y axis
            
        end
        
    end
    
    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractNarrowbandArrayProcessing('subarray');
        pPRFSource = matlab.system.display.internal.Property('PRFSource', ...
            'Description', 'Specify PRF as');
        props = {pPRFSource,...
          'PRF',...
          'ElevationAngleSource',...
          'ElevationAngle',...
          'NumAngleSamples',...
          'NumDopplerSamples'};
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:AngleDopplerResponseTitle')),...
              'Text',getString(message('phased:library:block:AngleDopplerResponseDesc')));
      end      
    end
    methods (Access = protected) %For Simulink propagation and mask
        function varargout = getOutputNamesImpl(~)
            varargout = {'Resp','Ang','Dop'};
        end
        
        function varargout = getInputNamesImpl(obj)
            varargout = {'X'};
            if (obj.PRFSource(1) == 'I') %Input Port
                varargout{end+1} = 'PRF';
            end
            if (obj.ElevationAngleSource(1) == 'I') %Input Port
                varargout{end+1} = 'El';
            end
            
        end
        
        function varargout = getOutputSizeImpl(obj)
            P = obj.NumDopplerSamples;
            Q = obj.NumAngleSamples;
            varargout = {[P Q],[Q 1],[P 1]};
        end
        function varargout = isOutputFixedSizeImpl(obj) %#ok<MANU>
            varargout = {true, true, true};
        end
        function varargout = getOutputDataTypeImpl(obj)  %#ok<MANU>
            varargout = {'double','double','double'};
        end
        function varargout = isOutputComplexImpl(obj)  %#ok<MANU>
            varargout = {true,false,false};
        end    
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Angle Doppler\nResponse');
        end                               
    end
    
end
