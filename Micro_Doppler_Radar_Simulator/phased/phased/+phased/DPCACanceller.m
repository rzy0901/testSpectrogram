classdef (Sealed,StrictDefaults) DPCACanceller < phased.internal.AbstractDPCA
%DPCACanceller  Displaced phase center array (DPCA) pulse canceller
%   H = phased.DPCACanceller creates a displaced phase center array (DPCA)
%   canceller System object, H. This object performs 2-pulse DPCA
%   processing on the input data for a ULA.
%
%   H = phased.DPCACanceller(Name,Value) creates a DPCA object, H, with
%   the specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X,IDX) applies the DPCA processing to the input data X. X
%   must be a 3-dimensional MxNxP numeric array whose dimensions are
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
%   [Y,W] = step(H,X,IDX) returns additional output W as the processing
%   weights, when you set the WeightsOutputPort property to true. When you
%   set the PreDopplerOutput property to true, W is a (N*2)x(P-1) matrix
%   containing the weights used to obtain the pre-Doppler data. Each column
%   in W represents the weights corresponding to successive pulses in X.
%   When you set the PreDopplerOutput property to false, W is a length
%   (N*P) column vector.
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
%   DPCACanceller methods:
%
%   step     - Perform DPCA processing on input data (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a DPCA object with same property values
%   isLocked - Locked status (logical)
%
%   DPCACanceller properties:
%
%   SensorArray         - Sensor array
%   PropagationSpeed    - Signal propagation speed
%   OperatingFrequency  - Operating frequency
%   PRFSource           - Source of PRF
%   PRF                 - Pulse repetition frequency
%   DirectionSource     - Source of direction
%   Direction           - Receiving mainlobe direction
%   NumPhaseShifterBits - Number of bits in phase shifters
%   DopplerSource       - Source of targeting Doppler
%   Doppler             - Targeting Doppler
%   WeightsOutputPort   - Enable weights output
%   PreDopplerOutput    - Output pre-Doppler result
%   
%   % Example:
%   %    Process the data cube using a DPCA processor.  The weights are 
%   %    calculated for the 71st cell of a collected data cube. The look 
%   %    direction is [0; 0] degrees and the Doppler is 12980 Hz.
%
%   load STAPExampleData;   % load data 
%   canceller = phased.DPCACanceller('SensorArray',STAPEx_HArray,...
%           'PRF',STAPEx_PRF,'PropagationSpeed',STAPEx_PropagationSpeed,...
%           'OperatingFrequency',STAPEx_OperatingFrequency,...
%           'WeightsOutputPort',true,...
%           'DirectionSource','Input port','DopplerSource','Input port');
%   [y,w] = canceller(STAPEx_ReceivePulse,71,[0; 0],12980);
%   resp = phased.AngleDopplerResponse(...
%           'SensorArray',canceller.SensorArray,...
%           'OperatingFrequency',canceller.OperatingFrequency,...
%           'PRF',canceller.PRF,...
%           'PropagationSpeed',canceller.PropagationSpeed);
%   plotResponse(resp,w);
%
%   See also phased, phased.STAPSMIBeamformer, phased.ADPCACanceller,
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
    methods

        function obj = DPCACanceller(varargin)
        %DPCA   Construct the DPCA object
            obj@phased.internal.AbstractDPCA(varargin{:});
        end
    end
    
    methods(Access = protected)
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractDPCA(obj,s);
            if isfield(s,'isLocked')
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

    end
    
    methods(Access = protected)
        
        function wOut = calcPreDopplerWeights(obj,angle,~,~) 
        %calcPreDopplerWeights Calculate pre-Doppler weights
        %   W = calcPreDopplerWeights(H,ANG) calculates the pre-Doppler
        %   weights W of H. ANG specifies the receiving mainlobe
        %   direction.  ANG must be a 2x1 vector with the format of
        %   [AzimuthAngle; ElevationAngle] (in degrees).
        %
        %   Because DPCA is a reduced dimension method, a pre-Doppler
        %   weight vector is generated for every two consecutive pulses.
        %   Therefore, W is a matrix whose columns are the pre-Doppler
        %   weights for each 2-consecutive-pulse space time snapshot.
            
            % return w as a matrix, dimPulseSnapshot * dimReducedDoppler
            w = formDPCAWeights(obj,angle);
            wOut = repmat(w,1,obj.pReducedDimension);
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('DPCA\nCanceller');
        end                
    end   
    methods (Static,Hidden,Access=protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:DPCACancellerTitle')),...
                'Text',getString(message('phased:library:block:DPCACancellerDesc')));
        end
    end
end




