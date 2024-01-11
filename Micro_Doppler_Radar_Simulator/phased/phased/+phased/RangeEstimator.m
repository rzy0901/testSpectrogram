classdef (Sealed,StrictDefaults) RangeEstimator < ...
        phased.internal.AbstractParameterEstimator
%RangeEstimator    Range estimation for radar detections
%   H = phased.RangeEstimator creates a range estimator System object, H.
%   This object estimates range measurements from detections found in the
%   radar response data by fitting a quadratic curve to the return signal
%   received from the target along the range dimension of the response
%   data.
%
%   H = phased.RangeEstimator(Name,Value) creates a range estimator object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   RNGEST = step(H,RESP,RANGEGRID,DETIDX) estimates range measurements
%   from detections in the range response data of RESP. Range estimates are
%   computed for each detection position reported in DETIDX.
%
%   RESP is the complex, range processed, radar response data. RESP is
%   either an Mx1 column vector, an MxN matrix, an MxP matrix, or an MxNxP
%   array containing the radar response, where M indicates the number of
%   range samples, N indicates the number of beams, and P is the number of
%   Doppler bins or pulses in the radar response data.
% 
%   RANGEGRID is an Mx1 column vector defining the range values assigned to
%   the range axis for the sampled radar data in RESP.
% 
%   DETIDX is a DxL matrix where L is the number of detections and D
%   represents the number of dimensions in the radar response data. Each
%   column of DETIDX contains the indices for each of the D-dimensions of
%   RESP where that detection was found.
% 
%   RNGEST is an Lx1 column vector of range estimates for each of the
%   detection locations in DETIDX.
% 
%   [RNGEST,RNGVAR] = step(H,RESP,RANGEGRID,DETIDX,NOISEPOWER) uses the
%   values in NOISEPOWER to estimate the variance of the range estimates.
%   This syntax applies when you set VarianceOutputPort to true and
%   NoisePowerSource to 'Input port'.
% 
%   NOISEPOWER is either a scalar or a 1xL row vector. When NOISEPOWER is a
%   scalar, the same noise power value is applied to all detections. When
%   NOISEPOWER is a vector, each value represents the noise power estimate
%   at the corresponding detection in DETIDX. NOISEPOWER units are linear
%   (not decibels) and must use the same scale factor as the response data.
% 
%   RNGVAR is an Lx1 column vector of variances for each of the range
%   estimates in RNGEST. The estimator variance is computed using the
%   Ziv-Zakai bound.
% 
%   RNGEST = step(H,RESP,RANGEGRID,DETIDX,CLUSTERIDS) uses the cluster ID
%   values to map multiple detections in DETIDX to a single hypothetical
%   target in the radar response data. This syntax applies when you set
%   ClusterInputPort to true.
% 
%   CLUSTERIDS is a 1xL row vector of positive integers. Each value in
%   CLUSTERIDS is a unique identifier indicating a hypothetical target to
%   which the corresponding detection crossing in DETIDX is mapped.
% 
%   In this syntax, the length L of RNGEST and RNGVAR corresponds to the
%   number of unique cluster IDs defined in CLUSTERIDS.
%
%   You can combine optional input and output arguments when their enabling
%   properties are set. Optional inputs and outputs must be listed in the
%   same order as the order of the enabling properties. For example,
%
%   [RNGEST,RNGVAR] = step(H,RESP,RANGEGRID,DETIDX,NOISEPOWER,CLUSTERIDS)
% 
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   RangeEstimator methods:
%
%   step     - Calculate the range estimates (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create range estimator object with same property values
%   isLocked - Locked status (logical)
%
%   RangeEstimator properties:
%
%   NumEstimatesSource - Source of number of range estimates to report
%   NumEstimates       - Maximum number of range estimates to report
%   VarianceOutputPort - Output variance for range estimates
%   RMSResolution      - Root-mean-square (RMS) range resolution of the
%                        waveform after range processing has been applied
%   NoisePowerSource   - Source of noise power used to compute variance of
%                        range estimates
%   NoisePower         - Value used for all detections to compute variance
%                        of range estimates
%   ClusterInputPort   - Enable cluster ID input
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data RESP is single precision,
%   the output data is single precision. If the input data RESP is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Estimate range measurements from the detections found in the
%   %   range-Doppler response of a pulsed radar signal. The signal
%   %   includes three target returns. The targets are located 500, 520,
%   %   and 750 meters away from the radar.
%
%   load('RangeDopplerEstimatorData', ...
%       'resp','rnggrid','detidx','noise','rngrms');
% 
%   rngestimator = phased.RangeEstimator(...
%       'VarianceOutputPort',true,...
%       'NoisePowerSource','Input port',...
%       'RMSResolution',rngrms);
% 
%   [rngest,rngvar] = rngestimator(resp,rnggrid,detidx,noise)
%
%   See also phased, bw2range, phased.DopplerEstimator,
%   phased.RangeResponse, phased.RangeDopplerResponse, phased.CFARDetector,
%   phased.CFARDetector2D.

%   Copyright 2016 The MathWorks, Inc.

%   References
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, 2nd ed.,
%       McGraw-Hill Professional Engineering, 2014
%   [2] Mark Richards, James Scheer, William Holm, Principles of Modern
%       Radar: Basic Principles, SciTech Publishing, 2010

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %RMSResolution  Root-mean-square range resolution
        %   Specify a scalar value representing the root-mean-square range
        %   resolution of the processed waveform. The default value of this
        %   property is 1. RMSResolution should be specified in the same
        %   units as were used for the RANGEGRID input argument. This
        %   property only applies when you set the VarianceOutputPort
        %   property to true.
        RMSResolution = 1
    end
    
    methods 
        function obj = RangeEstimator(varargin)
            % RangeEstimator Constructor for phased.RangeEstimator class
            obj@phased.internal.AbstractParameterEstimator(varargin{:});
        end
    end
    
    methods
        function set.RMSResolution(obj,val)
            validateattributes( val, { 'double','single' }, { 'real', 'scalar', 'positive', 'finite' }, '', 'RMSResolution');
            obj.RMSResolution = val;
        end
    end
    
    methods (Access = protected)
        function setParameterDimension(obj)
            % Range is always the first dimension
            obj.pDimension = 1;
        end
        
        function var = estimateVariance(obj,snr,rangegrid)
            rtot = max(rangegrid)-min(rangegrid);
            apb = rtot^2/12; % a priori bound
            
            rrms = obj.RMSResolution;
            crlb = rrms^2./(8*pi^2*snr);
            
            var = apb*erfc(sqrt(snr/2))+crlb.*gammainc(3/2,snr/(2*sqrt(2)),'upper');
            var = real(var);
        end

        function validateInputsImpl(obj,resp,rangegrid,detpos,varargin)
            validateInputsImpl@phased.internal.AbstractParameterEstimator(obj,resp,rangegrid,detpos,varargin{:});
            
            validateattributes(rangegrid,{'double','single'},{'real','column','numel',size(resp,1)},'','RANGEGRID');
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
                case 'RMSResolution'
                    if ~obj.VarianceOutputPort
                        flag = true;
                    end
            end
        end
        
        function varargout = getInputNamesImpl(obj)
            varargout = {'Resp','Range','DetIdx','NoisePower','Clusters'};
            
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
            str = sprintf('Range Estimator');
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
                'RMSResolution',...
                'NoisePowerSource',...
                'NoisePower',...
                'ClusterInputPort'
                };
            
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:RangeEstimatorTitle')),...
                'Text',getString(message('phased:library:block:RangeEstimatorDesc')));
        end
    end
end


