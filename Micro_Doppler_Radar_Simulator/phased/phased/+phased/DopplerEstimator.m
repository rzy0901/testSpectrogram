classdef (Sealed,StrictDefaults) DopplerEstimator < ...
        phased.internal.AbstractParameterEstimator
%DopplerEstimator    Doppler estimation for radar detections
%   H = phased.DopplerEstimator creates a Doppler estimator System object,
%   H. This object estimates Doppler measurements from detections found in
%   the radar response data by fitting a quadratic curve to the return
%   signal received from the target along the Doppler dimension of the
%   response data.
%
%   H = phased.DopplerEstimator(Name,Value) creates a Doppler estimator
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   DOPEST = step(H,RESP,DOPGRID,DETIDX) estimates Doppler measurements
%   from detections in the Doppler response data of RESP. Doppler estimates
%   are computed for each detection position reported in DETIDX.
%
%   RESP is the complex, Doppler processed, radar response data. RESP is
%   either a Px1 column vector, an MxP matrix, an NxP matrix, or an MxNxP
%   array containing the radar response, where M indicates the number of
%   range samples, N indicates the number of beams, and P is the number of
%   Doppler bins in the radar response data.
% 
%   DOPGRID is a Px1 column vector defining the Doppler frequencies or the
%   radial speed values assigned to the Doppler axis for the sampled radar
%   data in RESP. When DOPGRID is a vector of Doppler frequencies,
%   DopplerEstimator will return the estimated Doppler frequency for each
%   detection. When DOPGRID is a vector of radial speeds, DopplerEstimator
%   will return the estimated radial speed for each detection.
% 
%   DETIDX is a DxL matrix where L is the number of detections and D
%   represents the number of dimensions in the radar response data. Each
%   column of DETIDX contains the indices for each of the D-dimensions of
%   RESP where that detection was found.
% 
%   DOPEST is an Lx1 column vector of Doppler estimates for each of the
%   detection locations in DETIDX.
% 
%   [DOPEST,DOPVAR] = STEP(H,RESP,DOPGRID,DETIDX,NOISEPOWER) uses the
%   values in NOISEPOWER to estimate the variance of the Doppler estimates.
%   This syntax applies when you set VarianceOutputPort to true and
%   NoisePowerSource to 'Input port'.
% 
%   NOISEPOWER is either a scalar or a 1xL row vector. When NOISEPOWER is a
%   scalar, the same noise power value is applied to all detections. When
%   NOISEPOWER is a vector, each value represents the noise power estimate
%   at the corresponding detection in DETIDX. NOISEPOWER units are linear
%   (not decibels) and must use the same scale factor as the response data.
% 
%   DOPVAR is an Lx1 column vector of variances for each of the Doppler
%   estimates in DOPEST. The estimator variance is computed using the
%   Ziv-Zakai bound.
% 
%   DOPEST = step(H,RESP,DOPGRID,DETIDX,CLUSTERIDS) uses the cluster ID
%   values to map multiple detections in DETIDX to a single hypothetical
%   target in the radar response data. This syntax applies when you set
%   ClusterInputPort to true.
% 
%   CLUSTERIDS is a 1xL row vector of positive integers. Each value in
%   CLUSTERIDS is a unique identifier indicating a hypothetical target to
%   which the corresponding detection crossing in DETIDX is mapped.
% 
%   In this syntax, the length L of DOPEST and DOPVAR corresponds to the
%   number of unique cluster IDs defined in CLUSTERIDS.
%
%   You can combine optional input and output arguments when their enabling
%   properties are set. Optional inputs and outputs must be listed in the
%   same order as the order of the enabling properties. For example,
%
%   [DOPEST,DOPVAR] = STEP(H,RESP,DOPGRID,DETIDX,NOISEPOWER,CLUSTERIDS)
% 
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   DopplerEstimator methods:
%
%   step     - Calculate the Doppler estimates
%   release  - Allow property value and input characteristics changes
%   clone    - Create Doppler estimator object with same property values
%   isLocked - Locked status (logical)
%
%   DopplerEstimator properties:
%
%   NumEstimatesSource - Source of number of Doppler estimates to report
%   NumEstimates       - Maximum number of Doppler estimates to report
%   VarianceOutputPort - Output variance for Doppler estimates
%   NumPulses          - Number of pulses processed in Doppler response
%   NoisePowerSource   - Source of noise power used to compute variance of
%                        Doppler estimates
%   NoisePower         - Value used for all detections to compute variance
%                        of Doppler estimates
%   ClusterInputPort   - Enable cluster ID input
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Estimate the radial speed from the detections found in the range-
%   %   Doppler response of a pulsed radar signal. The signal includes
%   %   three target returns. The targets have radial speeds of 60, -20,
%   %   and -40 meters per second.
%
%   load('RangeDopplerEstimatorData', ...
%       'resp','spdgrid','detidx','noise','numPulse');
% 
%   dopestimator = phased.DopplerEstimator(...
%       'VarianceOutputPort',true,...
%       'NoisePowerSource','Input port',...
%       'NumPulses',numPulse);
% 
%   [spdest,spdvar] = dopestimator(resp,spdgrid,detidx,noise)
%
%   See also phased, phased.RangeEstimator, phased.RangeDopplerResponse,
%   phased.CFARDetector, phased.CFARDetector2D.

%   Copyright 2016 The MathWorks, Inc.

%   References
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, 2nd ed.,
%       McGraw-Hill Professional Engineering, 2014
%   [2] Mark Richards, James Scheer, William Holm, Principles of Modern
%       Radar: Basic Principles, SciTech Publishing, 2010

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable, PositiveInteger)
        %NumPulses  Number of pulses in Doppler processed waveform
        %   Specify a scalar value representing the number of pulses in
        %   the Doppler processed waveform. The default value of this
        %   property is 2. This property only applies when you set the
        %   VarianceOutputPort property to true.
        NumPulses = 2
    end
    
    methods 
        function obj = DopplerEstimator(varargin)
            % DopplerEstimator Constructor for phased.DopplerEstimator class
            obj@phased.internal.AbstractParameterEstimator(varargin{:});
        end
    end
        
    methods (Access = protected)
        
        function setParameterDimension(obj)
            % Doppler is always the last dimension
            if obj.pCubeSize(3) ~= -1
                obj.pDimension = 3;
            elseif obj.pCubeSize(2) > 1
                obj.pDimension = 2;
            else
                obj.pDimension = 1;
            end
        end
        
        function var = estimateVariance(obj,snr,dopgrid)
            
            numDopbin = length(dopgrid);
            prf = (max(dopgrid)-min(dopgrid))*numDopbin/(numDopbin-1);
            
            apb = prf^2/12; % a priori bound
            
            numPulses = obj.NumPulses;
            crlb = 6*prf^2./((2*pi)^2*(numPulses^2-1)*snr);
            
            var = apb*erfc(sqrt(snr/2))+crlb.*gammainc(3/2,snr/(2*sqrt(2)),'upper');
            var = real(var);
        end
        
        function validateInputsImpl(obj,resp,dopgrid,detpos,varargin)
            validateInputsImpl@phased.internal.AbstractParameterEstimator(obj,resp,dopgrid,detpos,varargin{:});

            % dopgrid must have same length as last dim in resp
            szResp = size(resp);
            if numel(szResp)==3
                exLen = szResp(3);
            else
                if szResp(2)>1
                    exLen = szResp(2);
                else
                    exLen = szResp(1);
                end
            end
            
            validateattributes(dopgrid,{'double','single'},{'real','column','numel',exLen},'','DOPGRID');
        end
        
        function num = getNumInputsImpl(obj)
            if obj.VarianceOutputPort
                num = 5;
                
                if strcmp(obj.NoisePowerSource,'Property')
                    num = num-1;
                end
                
                if ~obj.ClusterInputPort
                    num = num-1;
                end
            else
                num = 4;
                
                if ~obj.ClusterInputPort
                    num = num-1;
                end
            end
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@phased.internal.AbstractParameterEstimator(obj, prop);
            switch prop
                case 'NumPulses'
                    if ~obj.VarianceOutputPort
                        flag = true;
                    end
            end
        end
        
        function varargout = getInputNamesImpl(obj)
            varargout = {'Resp','Doppler','DetIdx','NoisePower','Clusters'};
            
            if ~obj.VarianceOutputPort
                varargout = varargout(~ismember(varargout,{'NoisePower'}));
            end
            
            if strcmp(obj.NoisePowerSource,'Property')
                varargout = varargout(~ismember(varargout,{'NoisePower'}));
            end
            
            if ~obj.ClusterInputPort
                varargout = varargout(~ismember(varargout,{'Clusters'}));
            end
        end
        
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Doppler Estimator');
        end
    end
    
    methods (Access = protected, Static, Hidden)
        function groups = getPropertyGroupsImpl
            pNumEstimatesSource = matlab.system.display.internal.Property(...
                'NumEstimatesSource', ...
                'IsGraphical', false, ...
                'UseClassDefault', false,'Default','Property');
            
            props = {...
                pNumEstimatesSource ,...
                'NumEstimates'};
    
            groups = matlab.system.display.Section(...
                'Title','Parameters',...
                'PropertyList',props);
            
            props = {...
                'VarianceOutputPort',...
                'NumPulses',...
                'NoisePowerSource',...
                'NoisePower',...
                'ClusterInputPort'
                };
            
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:DopplerEstimatorTitle')),...
              'Text',getString(message('phased:library:block:DopplerEstimatorDesc')));
      end
    end
end


