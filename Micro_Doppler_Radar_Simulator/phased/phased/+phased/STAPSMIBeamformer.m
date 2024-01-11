classdef (Sealed,StrictDefaults) STAPSMIBeamformer < phased.internal.AbstractSTAP
%STAPSMIBeamformer    Sample matrix inversion (SMI) STAP beamformer
%   H = phased.STAPSMIBeamformer creates a sample matrix inversion (SMI)
%   STAP beamformer System object, H. This object performs SMI space-time
%   adaptive processing (STAP) on the input data.
%
%   H = phased.STAPSMIBeamformer(Name,Value) creates an SMI object, H, with
%   the specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,IDX) applies SMI processing to the input data X. X must
%   be a 3-dimensional MxNxP numeric array whose dimensions are (range,
%   channels, pulses). The number of channels is the number of subarrays if
%   SensorArray contains subarrays, or the number of elements otherwise.
%   The processing weights are calculated according to the range cell
%   specified by IDX. The targeting direction and the targeting Doppler
%   are specified by the Direction and Doppler properties, respectively. Y
%   is a length M column vector. This option is available when
%   DirectionSource property is 'Property' and DopplerSource property is
%   'Property'.
%
%   [Y,W] = step(H,X,IDX) returns additional output W as the processing
%   weights, when you set the WeightsOutputPort property to true. W is a
%   length (N*P) column vector.
%
%   Y = step(H,X,IDX,PRF) uses PRF as the pulse repetition frequency (in
%   Hz), when you set the PRFSource property to 'Input port'. PRF must be a
%   scalar.
%
%   [Y,W] = step(H,X,IDX,PRF) uses PRF as the pulse repetition frequency
%   (in Hz), when you set the PRFSource property to 'Input port'.
%
%   Y = step(H,X,IDX,ANG) uses ANG as the targeting direction, when
%   you set the DirectionSource property to 'Input port'. ANG must be a
%   2x1 vector in the form of [AzimuthAngle; ElevationAngle] (in degrees).
%   Azimuth angle should be between [-180 180]. Elevation angle should be
%   between [-90 90].
%
%   [Y,W] = step(H,X,IDX,ANG) uses ANG as the targeting direction,
%   when you set the DirectionSource property to 'Input port'.
%
%   Y = step(H,X,IDX,DOP) uses DOP as the targeting Doppler frequency
%   (in Hz), when you set the DopplerSource property to 'Input port'. DOP
%   must be a scalar.
%
%   [Y,W] = step(H,X,IDX,DOP) uses DOP as the targeting Doppler
%   frequency (in Hz), when you set the DopplerSource property to 'Input
%   port'.
%
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   [Y,W] = step(H,X,IDX,PRF,ANG,DOP)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   STAPSMIBeamformer methods:
%
%   step     - Perform SMI STAP processing on input data (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a STAP SMI beamformer object with same property 
%              values
%   isLocked - Locked status (logical)
%
%   STAPSMIBeamformer properties:
%
%   SensorArray         - Sensor array
%   PropagationSpeed    - Signal propagation speed
%   OperatingFrequency  - Operating frequency
%   PRFSource           - Source of PRF
%   PRF                 - Pulse repetition frequency
%   DirectionSource     - Source of direction
%   Direction           - Targeting direction
%   NumPhaseShifterBits - Number of bits in phase shifters
%   DopplerSource       - Source of targeting Doppler
%   Doppler             - Targeting Doppler
%   NumGuardCells       - Number of guard cells
%   NumTrainingCells    - Number of training cells
%   WeightsOutputPort   - Enable weights output
%   
%   % Example:
%   %   Process the data cube using an SMI processor. The weights are 
%   %   calculated for the 71st cell of a collected data cube pointing 
%   %   to the direction of [45; -35] degrees and the Doppler of 12980 Hz.
%
%   load STAPExampleData;   % load data 
%   stap = phased.STAPSMIBeamformer('SensorArray',STAPEx_HArray,...
%           'PRF',STAPEx_PRF,'PropagationSpeed',STAPEx_PropagationSpeed,...
%           'OperatingFrequency',STAPEx_OperatingFrequency,...
%           'NumTrainingCells',100,'WeightsOutputPort',true,...
%           'DirectionSource','Input port','DopplerSource','Input port');
%   [y,w] = stap(STAPEx_ReceivePulse,71,[45; -35],12980);
%   adresp = phased.AngleDopplerResponse('SensorArray',stap.SensorArray,...
%           'OperatingFrequency',stap.OperatingFrequency,'PRF',stap.PRF,...
%           'PropagationSpeed',stap.PropagationSpeed);
%   plotResponse(adresp,w);
%
%   See also phased, phased.DPCACanceller, phased.ADPCACanceller,
%   phased.PhaseShiftBeamformer, phased.AngleDopplerResponse.

%   Copyright 2008-2017 The MathWorks, Inc.

% References:
% [1] James Ward, "Space-Time Adaptive Processing for Airborne
%     Radar". MIT Lincoln Lab tech report 1015, 1994.
% [2] J. R. Guerci, "Space-Time Adaptive Processing for Radar".
%     Artech House, 2003.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties (Nontunable)
        %Direction    Targeting direction (deg)
        %   Specify the targeting direction of the SMI processor as a
        %   length 2 column vector. The direction is specified in the
        %   format of [AzimuthAngle; ElevationAngle] (in degrees). Azimuth
        %   angle should be between -180 and 180. Elevation angle should be
        %   between -90 and 90. This property applies when you set the
        %   DirectionSource property to 'Property'. The default value of
        %   this property is [0; 0].
        Direction = [0; 0];
        %NumGuardCells  Number of guard cells
        %   Specify the number of guard cells used in the training as an
        %   even integer. It specifies the total number of cells on both
        %   sides of the cell under test. The default value of this
        %   property is 2 indicating there is one guard cell at both the
        %   front and back of the cell under test.
        NumGuardCells = 2;
        %NumTrainingCells   Number of training cells
        %   Specify the number of training cells used in the training as an
        %   even integer. Whenever possible, the training cells are equally
        %   divided before and after the cell under test. The default value
        %   of this property is 2 indicating there is one training cell at
        %   both the front and back of the cell under test.
        NumTrainingCells = 2;
    end

    properties (Access = private, Nontunable)
        cTraining
        cSteeringVector
    end
    
    methods

        function set.Direction(obj,val)
            sigdatatypes.validateAzElAngle(val,'phased.STAPSMIBeamformer',...
                'Direction',{'size',[2 1]});
            obj.Direction = val;
        end
        
        function set.NumGuardCells(obj,val)
            validateattributes(val,{'double'},...
                {'finite','nonnan','nonnegative','even','scalar'},...
                'phased.SMI','NumGuardCells');
            obj.NumGuardCells = val;
        end
        
        function set.NumTrainingCells(obj,val)
            validateattributes(val,{'double'},...
                {'finite','nonnan','positive','even','scalar'},...
                'phased.SMI','NumTrainingCells');
            obj.NumTrainingCells = val;
        end
    end
    
    methods
        
        % Constructor
        function obj = STAPSMIBeamformer(varargin)
        %SMI    Construct SMI object
            obj@phased.internal.AbstractSTAP(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        
        function flag = isSubarraySupported(obj) %#ok<MANU>
            flag = true;
        end
        
        function num = getNumInputsImpl(obj)
            if (obj.DirectionSource(1) == 'P')
                if (obj.DopplerSource(1) == 'P')
                    num = 2;
                else
                    num = 3;
                end
            else
                if (obj.DopplerSource(1) == 'P')
                    num = 3;
                else
                    num = 4;
                end
            end
            if (obj.PRFSource(1) == 'I')
                num = num+1 ; %5
            end
        end
                
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if (obj.DirectionSource(1) == 'P')
                if (obj.DopplerSource(1) ~= 'P') && ...
                        strcmp(prop, 'Doppler')
                    flag = true;
                end
            else
                if (obj.DopplerSource(1) == 'P')
                    if strcmp(prop, 'Direction')
                        flag = true;
                    end
                else
                    if strcmp(prop,'Direction') || ...
                            strcmp(prop,'Doppler')
                        flag = true;
                    end
                end
            end
            if(obj.PRFSource(1) == 'I'&& strcmp(prop, 'PRF'))
                flag = true;
            end  
        end
        
        function validateInputsImpl(obj,x,cutidx,varargin)
            validateInputsImpl@phased.internal.AbstractSTAP(obj,x,cutidx);
            
            if(obj.PRFSource(1)=='P')
                prf = [];
                if (obj.DirectionSource(1) == 'P')
                    adir = [];
                    if (obj.DopplerSource(1) == 'I')
                        dop = varargin{1};
                    else
                        dop = [];
                    end
                else
                    adir = varargin{1};
                    if (obj.DopplerSource(1) == 'I')
                        dop = varargin{2};
                    else
                        dop = [];
                    end
                end
            else
                prf = varargin{1};
                if (obj.DirectionSource(1) == 'P')
                    adir = [];
                    if (obj.DopplerSource(1) == 'I')
                        dop = varargin{2};
                    else
                        dop = [];
                    end
                else
                    adir = varargin{2};
                    if (obj.DopplerSource(1) == 'I')
                        dop = varargin{3};
                    else
                        dop = [];
                    end
                end
            end
            
            if ~isempty(prf)
                validateInputPRFSpec(obj,prf);
            end
            if ~isempty(adir)
                validateInputAngleSpec(obj,adir);
            end
            if ~isempty(dop)
                validateInputDopplerSpec(obj,dop);
            end
        end
        
        function setupImpl(obj, x, CUTIdx, varargin) %#ok<INUSD>
            setupImpl@phased.internal.AbstractSTAP(obj,x);
            obj.cTraining = phased.internal.STAPTraining(...
                obj.NumGuardCells,obj.NumTrainingCells);
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'NumPhaseShifterBits',obj.NumPhaseShifterBits);
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = isInputComplexityLockedImpl@phased.internal.AbstractSTAP(obj,index);
            if(obj.PRFSource(1) == 'P')
                if (obj.DirectionSource(1) == 'P')
                    if (obj.DopplerSource(1) == 'I') && (index == 3)
                        flag = true;
                    end
                else
                    if index == 3
                        flag = true;
                    end
                    if (obj.DopplerSource(1) == 'I') && (index == 4)
                        flag = true;
                    end
                end
                
            else
                if index == 3
                    flag = true;
                end
                if (obj.DirectionSource(1) == 'P')
                    if (obj.DopplerSource(1) == 'I') && (index == 4)
                        flag = true;
                    end
                else
                    if index == 4
                        flag = true;
                    end
                    if (obj.DopplerSource(1) == 'I') && (index == 5)
                        flag = true;
                    end
                end
            end
        end
        
        function releaseImpl(obj)
            release(obj.cTraining);
            release(obj.cSteeringVector);
        end
        
        function resetImpl(obj)
            reset(obj.cTraining);
            reset(obj.cSteeringVector);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSTAP(obj);
            if isLocked(obj)
                s.cTraining = saveobj(obj.cTraining);
                s.cSteeringVector = saveobj(obj.cSteeringVector);
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSTAP(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cTraining = phased.internal.STAPTraining.loadobj(s.cTraining);
                    s = rmfield(s,'cTraining');
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
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

        function [y, w] = stepImpl(obj, x, CUTIdx, varargin)
            
            if(obj.PRFSource(1) == 'P')
                prf = obj.PRF;
                if (obj.DirectionSource(1) == 'P')
                    angle = obj.Direction;
                    if (obj.DopplerSource(1) == 'P')
                        dop = obj.Doppler;
                    else
                        dop = varargin{1};
                    end
                else
                    angle = varargin{1};
                    if (obj.DopplerSource(1) == 'P')
                        dop = obj.Doppler;
                    else
                        dop = varargin{2};
                    end
                end
                
            else
                prf = varargin{1};
                if (obj.DirectionSource(1) == 'P')
                    angle = obj.Direction;
                    if (obj.DopplerSource(1) == 'P')
                        dop = obj.Doppler;
                    else
                        dop = varargin{2};
                    end
                else
                    angle = varargin{2};
                    if (obj.DopplerSource(1) == 'P')
                        dop = obj.Doppler;
                    else
                        dop = varargin{3};
                    end
                end
            end
            
            prf = validateInputPRF(obj,prf);
            angle = validateInputAngle(obj,angle);
            dop = validateInputDoppler(obj,dop,prf);
            CUTIdx = validateInputCUTIdx(obj,CUTIdx);
            
            [w, datavec] = privCalcWeights(obj,x,CUTIdx,prf,angle,dop);
                       
            y = w'*datavec;
            y = y(:);
            
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('SMI\nBeamformer');
        end                               
    end
    
    methods (Access = private)
        function [w, datavec] = privCalcWeights(obj,x,CUTIdx,prf,angle,dp)
            % Normalize Doppler
            dp = dp/prf;
 
            [datavec, datadim]  = obj.organizeSpaceTimeSnapshots(x);
            dopplerSteeringVec = dopsteeringvec(dp,datadim(3));
            angleSteeringVec = step(obj.cSteeringVector,obj.OperatingFrequency,angle);
            stSteeringVec = kron(dopplerSteeringVec, angleSteeringVec);

            TrainingData = step(obj.cTraining,datavec,CUTIdx);
            
            w = phased.internal.lcmvweights(TrainingData.',...
                stSteeringVec,1,0);

        end
    end
    
    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractSTAP('subarray');
        props = {...
          'NumGuardCells',...
          'NumTrainingCells'};
        groups(1).PropertyList = [groups(1).PropertyList(1:end-1) props ...
                                  groups(1).PropertyList(end)];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:STAPSMIBeamformerTitle')),...
              'Text',getString(message('phased:library:block:STAPSMIBeamformerDesc')));
      end      
    end
    
end



