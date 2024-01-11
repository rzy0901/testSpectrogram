classdef  (Sealed,StrictDefaults) ADPCACanceller < phased.internal.AbstractDPCA 
%ADPCACanceller  Adaptive DPCA (ADPCA) pulse canceller
%   H = phased.ADPCACanceller creates an adaptive displaced phase center
%   array (ADPCA) canceller System object, H. This object performs 2-pulse
%   ADPCA processing on the input data for a ULA.
%
%   H = phased.ADPCACanceller(Name,Value) creates an ADPCA object, H,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,IDX) applies the ADPCA processing to the input data X.
%   X must be a 3-dimensional MxNxP numeric array whose dimensions are
%   (range, channels, pulses). The processing weights are calculated
%   according to the range cell specified by IDX. The receiving mainlobe
%   direction and the targeting Doppler are specified by Direction and
%   Doppler properties, respectively. 
%
%   If PreDopplerOutput property is true, Y is an Mx(P-1) matrix that
%   contains the pre-Doppler data. Each column in Y represents the result
%   obtained by cancelling the two successive pulses. If PreDopplerOutput
%   property is false, an FFT based Doppler processing is applied to the
%   pre-Doppler data and Y is a length M column vector. Note that the
%   targeting Doppler only applies when PreDopplerOutput property is false.
%   This option is available when DirectionSource property is 'Property'
%   and DopplerSource property is 'Property'.
%
%   [Y,W] = step(H,X,IDX) returns the additional output W as the
%   processing weights, when you set the WeightsOutputPort property to
%   true. When you set the PreDopplerOutput property to true, W is a
%   (N*2)x(P-1) matrix containing the weights used to obtain the
%   pre-Doppler data. Each column in W represents the weights corresponding
%   to successive pulses in X. When you set the PreDopplerOutput property
%   to false, W is a length (N*P) column vector.
%
%   Y = step(H,X,IDX,PRF) uses PRF as the pulse repetition frequency (in
%   Hz), when you set the PRFSource property to 'Input port'. PRF must be a
%   scalar.
%
%   [Y,W] = step(H,X,IDX,PRF) uses PRF as the pulse repetition frequency
%   (in Hz), when you set the PRFSource property to 'Input port'.
%
%   Y = step(H,X,IDX,ANG) uses ANG as the receiving mainlobe
%   direction, when you set the DirectionSource property to 'Input port'.
%   ANG must be a 2x1 vector in the form of [AzimuthAngle;
%   ElevationAngle] (in degrees). Azimuth angle should be between [-180
%   180]. Elevation angle should be between [-90 90].
%
%   [Y,W] = step(H,X,IDX,ANG) uses ANG as the receiving mainlobe
%   direction, when you set the DirectionSource property to 'Input port'.
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
%   ADPCACanceller methods:
%
%   step     - Perform ADPCA processing on input data (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create an ADPCA object with same property values
%   isLocked - Locked status (logical)
%
%   ADPCACanceller properties:
%
%   SensorArray         - Sensor array
%   PropagationSpeed    - Signal propagation speed
%   OperatingFrequency  - Operating frequency
%   PRFSource           - Source of pulse repetition frequency
%   PRF                 - Pulse repetition frequency
%   DirectionSource     - Source of direction
%   Direction           - Receiving mainlobe direction
%   NumPhaseShifterBits - Number of bits in phase shifters
%   DopplerSource       - Source of targeting Doppler
%   Doppler             - Targeting Doppler
%   NumGuardCells       - Number of guard cells
%   NumTrainingCells    - Number of training cells
%   WeightsOutputPort   - Enable weights output
%   PreDopplerOutput    - Output pre-Doppler result
%   
%   % Example:
%   %    Process the data cube using an ADPCA processor.  The weights are 
%   %    calculated for the 71st cell of a collected data cube. The look 
%   %    direction is [0; 0] degrees and the Doppler is 12980 Hz.
%
%   load STAPExampleData;   % load data 
%   canceller = phased.ADPCACanceller('SensorArray',STAPEx_HArray,...
%           'PRF',STAPEx_PRF,'PropagationSpeed',STAPEx_PropagationSpeed,...
%           'OperatingFrequency',STAPEx_OperatingFrequency,...
%           'NumTrainingCells',100,'WeightsOutputPort',true,...
%           'DirectionSource','Input port','DopplerSource','Input port');
%   [y,w] = canceller(STAPEx_ReceivePulse,71,[0; 0],12980);
%   Hresp = phased.AngleDopplerResponse('SensorArray',...
%           canceller.SensorArray,'OperatingFrequency',...
%           canceller.OperatingFrequency,'PRF',canceller.PRF,...
%           'PropagationSpeed',canceller.PropagationSpeed);
%   plotResponse(Hresp,w);
%
%   See also phased, phased.DPCACanceller, phased.STAPSMIBeamformer,
%   phased.AngleDopplerResponse.

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
    
    properties(Access = private, Nontunable)
        cTraining
    end
     
    methods

        function set.NumGuardCells(obj,val)
            validateattributes(val,{'double'},...
                {'finite','nonnan','nonnegative','even','scalar'},...
                'phased.ADPCA','NumGuardCells');
            obj.NumGuardCells = val;
        end
        
        function set.NumTrainingCells(obj,val)
            validateattributes(val,{'double'},...
                {'finite','nonnan','positive','even','scalar'},...
                'phased.ADPCA','NumTrainingCells');
            obj.NumTrainingCells = val;
        end
    end
    

    methods

        % Constructor
        function obj = ADPCACanceller(varargin)
            %ADPCA  Construct ADPCA object
            obj@phased.internal.AbstractDPCA(varargin{:});
        end

    end

    methods (Access = protected)
        
        function setupImpl(obj,varargin)
            setupImpl@phased.internal.AbstractDPCA(obj,varargin{:});
            obj.cTraining = phased.internal.STAPTraining(...
                obj.NumGuardCells,obj.NumTrainingCells);
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractDPCA(obj);
            release(obj.cTraining);
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractDPCA(obj);
            reset(obj.cTraining);
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractDPCA(obj);
            if isLocked(obj)
                s.cTraining = saveobj(obj.cTraining);
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractDPCA(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cTraining = phased.internal.STAPTraining.loadobj(s.cTraining);
                    s = rmfield(s,'cTraining');
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

        function w = calcPreDopplerWeights(obj,angle,tempdata,CUTIdx)
        %calcPreDopplerWeights Calculate pre-Doppler weights
        %   W = calcPreDopplerWeights(Hs,ANG,X,IDX) calculates the
        %   pre-Doppler weights W of Hs for the IDX-th range cell of the
        %   input data cube X. X is a 3-dimensional numeric array whose
        %   dimensions are (range, channels, pulses). ANG specifies the
        %   receiving mainlobe direction.  ANG must be a 2x1 vector with
        %   the format of [AzimuthAngle;ElevationAngle] (in degrees).
        %
        %   Because ADPCA is a reduced dimension method, a pre-Doppler
        %   weight vector is generated for every two consecutive pulses.
        %   Therefore, W is a matrix whose columns are the pre-Doppler
        %   weights for each 2-consecutive-pulse space time snapshot.

            dpcaSteeringVec = formDPCAWeights(obj,angle);
            
            % return w as a matrix, dimPulseSnapshot * dimReducedDoppler
            dimPulseSnapshot = obj.NumDPCAPulses * obj.pArrayNumElements;
            dimReducedDoppler = obj.pReducedDimension;
            tempdata = reshape(tempdata,...
                dimPulseSnapshot*dimReducedDoppler,[]);
            trnData = step(obj.cTraining,tempdata,CUTIdx);
            % in the form of[pulse1;pulse2;pulse2;pulse3;pulse3;pulse4,...]
            
            w = complex(zeros(dimPulseSnapshot,dimReducedDoppler));
            for m = 1:dimReducedDoppler
                tempdata = trnData((m-1)*dimPulseSnapshot+1:m*dimPulseSnapshot,:);
                w(:,m) = phased.internal.lcmvweights(tempdata.',...
                    dpcaSteeringVec,1,0);
            end
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('ADPCA\nCanceller');
        end                       
    end

    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractDPCA;
        props = {...
          'NumGuardCells',...
          'NumTrainingCells'};
        groups(1).PropertyList = [groups(1).PropertyList(1:end-2) ...
                                  props ...
                                  groups(1).PropertyList(end-1:end)];
      end
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:ADPCACancellerTitle')),...
              'Text',getString(message('phased:library:block:ADPCACancellerDesc')));
      end      
    end
end





