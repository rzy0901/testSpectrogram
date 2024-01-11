classdef (Hidden) AbstractDetector <  phased.internal.AbstractVarSizeEngine & ...
        matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates
%This class is for internal use only. It may be removed in the future.

%   Copyright 2009-2015 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Method     CFAR algorithm
        %   Specify the algorithm of CFAR detector using one of 'CA' |
        %   'SOCA' | 'GOCA' | 'OS', where the default is 'CA'. When you set
        %   the Method property to 'CA', the CFAR detector uses the
        %   cell-averaging algorithm. When you set the Method property to
        %   'SOCA', the CFAR detector uses the smallest-of cell-averaging
        %   algorithm. When you set the Method property to 'GOCA', the CFAR
        %   detector uses the greatest-of cell-averaging algorithm. When
        %   you set the Method property to 'OS', the CFAR detector uses the
        %   order statistic algorithm.
        Method = 'CA'
        
        %Rank   Rank of order statistic
        %   Specify the rank of the order statistic used in the order
        %   statistic CFAR algorithm as a positive integer. The value of
        %   the Rank property must be between 1 and N, where N is the total
        %   number of training cells. The default value of this property is
        %   1. This property only applies when you set the Method property
        %   to 'OS'.
        Rank = 1
        
        %ThresholdFactor    Threshold factor method
        %   Specify the method of obtaining the threshold factor using one
        %   of 'Auto' | 'Input port' | 'Custom', where the default is
        %   'Auto'. When you set the ThresholdFactor property to 'Auto',
        %   the threshold factor is calculated based on the desired
        %   probability of false alarm specified in the
        %   ProbabilityFalseAlarm property. The calculation assumes that
        %   each independent signal in the input is a single pulse coming
        %   out of a square law detector with no pulse integration. In
        %   addition, the noise is assumed to be white Gaussian. When you
        %   set the ThresholdFactor property to 'Input port', the threshold
        %   factor is specified through an input argument. When you set the
        %   ThresholdFactor property to 'Custom', the threshold factor is
        %   the value of the CustomThresholdFactor property.
        ThresholdFactor = 'Auto'
        
        %ProbabilityFalseAlarm  Probability of false alarm
        %   Specify the desired probability of false alarm as a scalar
        %   between 0 and 1 (not inclusive). This property only applies
        %   when you set the ThresholdFactor property to 'Auto'. The
        %   default value of this property is 0.1.
        ProbabilityFalseAlarm = 0.1
    end
    
    properties
        %CustomThresholdFactor  Custom threshold factor
        %   Specify the custom threshold factor as a positive scalar. This
        %   property only applies when you set the ThresholdFactor property
        %   to 'Custom'. This property is tunable. The default value of
        %   this property is 1.
        CustomThresholdFactor = 1
    end
    
    properties (Nontunable)
        %OutputFormat  Output format
        %   Specify the format used to report detections as one of 'CUT
        %   result' | 'Detection index', where the default is 'CUT result'.
        %   When set to 'CUT result', logical detection results will be
        %   reported for each of the cells under test. When set to
        %   'Detection index', a DxQ matrix is returned with the indices of
        %   the location of the detection for each of the input data's
        %   D-dimensions specified along each column. Here, Q is the
        %   product of the number of cells under test with the number of
        %   independent signals in the input data. Columns without valid
        %   detections are set to NaN.
        OutputFormat = 'CUT result'
    end
    
    properties (Nontunable, Logical)
        %ThresholdOutputPort    Output detection threshold
        %   Set this property to true to output the detection threshold.
        %   Set this property to false to not output the detection
        %   threshold. The default value of this property is false.
        ThresholdOutputPort = false;
        
        %NoisePowerOutputPort    Output estimated noise power
        %   Set this property to true to output the estimated noise power.
        %   Set this property to false to not output the estimated noise
        %   power. The default value of this property is false.
        NoisePowerOutputPort = false;
    end
    
    properties (Nontunable)
        %NumDetectionsSource  Source of the number of detections
        % Specify the source of the number of detections as one of 'Auto' |
        % 'Property'. The default is 'Auto'. If you set this property to
        % 'Auto', the number of detection indices reported is the total
        % number of cells under test that have a detection. If you set this
        % property to 'Property', the number of reported detections is
        % determined by the value of the NumDetections property. This
        % property only applies when you set the OutputFormat property to
        % 'Detection index'.
        NumDetectionsSource = 'Auto'
    end
    
    properties (Nontunable, PositiveInteger)
        %NumDetections  Maximum number of detections
        %   Specify the maximum number of detection indices to report as a
        %   positive integer scalar. The default value is 1. This property
        %   applies when you set the OutputFormat property to 'Detection
        %   index' and the NumDetectionsSource property to 'Property'.
        NumDetections = 1
    end
    
    properties (Access=protected, Nontunable)
        pNumChannels;
    end
    
    properties (Access=protected)
        pFactor;
    end
    
       
    properties(Constant, Hidden)
        ThresholdFactorSet = matlab.system.StringSet(...
            {'Auto','Input port','Custom'});
        MethodSet = matlab.system.StringSet(...
            {'CA','GOCA','SOCA','OS'});
        OutputFormatSet = matlab.system.StringSet(...
            {'CUT result','Detection index'});
        NumDetectionsSourceSet = dsp.CommonSets.getSet('AutoOrProperty');
    end
    
    methods
        function obj = AbstractDetector(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods
        function set.ProbabilityFalseAlarm(obj,val)
            % Pfa cannot be either 0 or 1 which makes Pd also 0 and 1.
            sigdatatypes.validateProbability(val,'phased.CFARDetector',...
                'ProbabilityFalseAlarm',{'double','single'},{'scalar','>',0,'<',1});
            obj.ProbabilityFalseAlarm = val;
        end
        
        function set.CustomThresholdFactor(obj,val)
            validateattributes( val, { 'double','single' }, ...
              { 'nonempty', 'finite', 'positive', 'scalar' }, '', 'CustomThresholdFactor');
            obj.CustomThresholdFactor = val;
        end
        
        function set.Rank(obj,val)
            validateattributes(val, {'double','single'}, ...
              {'nonempty','finite','positive','scalar','integer'},'','Rank');
            obj.Rank = val;
        end
        
    end
    
    methods (Access = protected)
        
        function validatePropertiesImpl(obj)
            NumTrainingCells = getNumTrainingCells(obj); 
      
            if obj.Method(1) == 'O' % OS
                sigdatatypes.validateIndex(obj.Rank,...
                    '','Rank',{'double','single'},{'<=',NumTrainingCells});
            end
        end
        
        function num = getNumInputsImpl(obj)
            if obj.ThresholdFactor(1) == 'I' %Input port
                num = 3;
            else
                num = 2;
            end
        end
        
        function num = getNumOutputsImpl(obj)
            num = 3;
            
            if ~obj.ThresholdOutputPort
                num = num-1;
            end
            
            if ~obj.NoisePowerOutputPort
                num = num-1;
            end
        end
       
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if  (obj.ThresholdFactor(1) ~= 'A') && ... %Auto
                    strcmp(prop, 'ProbabilityFalseAlarm')
                flag = true;
            end
            if (obj.ThresholdFactor(1) ~= 'C') && ... %Custom
                    strcmp(prop, 'CustomThresholdFactor')
                flag = true;
            end
            if (obj.Method(1) ~='O') && ...  %OS
                    strcmp(prop, 'Rank')
                flag = true;
            end
            if strcmp(obj.OutputFormat,'CUT result') && ...
                    strcmp(prop, 'NumDetectionsSource')
                flag = true;
            end
            if ( strcmp(obj.OutputFormat,'CUT result') || ...
                    strcmp(obj.NumDetectionsSource,'Auto') ) && ...
                    strcmp(prop, 'NumDetections')
                flag = true;
            end
        end
        
        function setupImpl(obj) 
            setupImpl@phased.internal.AbstractVarSizeEngine(obj);
            if obj.ThresholdFactor(1) == 'A' %Auto
                obj.pFactor = calcWGNThresholdFactor(obj);
            elseif obj.ThresholdFactor(1) == 'C' %Custom
                processTunedPropertiesImpl(obj);
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = true;
            if (obj.ThresholdFactor(1) == 'I') && (index == 3)
                flag = true;
            end
        end
        
        function flag = isOutputComplexityLockedImpl(obj,~)  %#ok<INUSD>
            flag = true;
        end
        
        function validateInputsImpl(obj,x,idx,varargin)
            cond = ~isa(x,'float');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','X','float');
            end
            
            cond = ~isreal(x);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:CFARDetector:ComplexInput', 'X');
            end
            
            cond = ~isa(idx,'float');
            if cond
                coder.internal.errorIf(cond, ...
                    'MATLAB:system:invalidInputDataType','IDX','float');
            end
            
            cond = ~isreal(idx);
            if cond
                coder.internal.errorIf(cond, ...
                    'phased:CFARDetector:ComplexInput', 'IDX');
            end
            
            if obj.ThresholdFactor(1) == 'I'
                thfac = varargin{1};
                cond = ~isa(thfac,'float');
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:invalidInputDataType','K','float');
                end
                
                cond = ~isscalar(thfac);
                if cond
                    coder.internal.errorIf(cond, ...
                        'MATLAB:system:inputMustBeScalar','K');
                end
                
                cond = ~isreal(thfac);
                if cond
                    coder.internal.errorIf(cond, ...
                        'phased:CFARDetector:ComplexInput', 'K');
                end
            end
        end
        
        function processTunedPropertiesImpl(obj)
            if obj.ThresholdFactor(1) == 'C'
                obj.pFactor = obj.CustomThresholdFactor;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
            if isLocked(obj)
                s.pFactor = obj.pFactor;
                s.pNumChannels = obj.pNumChannels;
            end
        end
        
        function flag = isInputSizeLockedImpl(~,idx)
            if idx == 1 || idx == 2
                flag = false;
            else
                flag = true;
            end
        end
     end
     
    methods (Access = protected)
        function varargout = getInputNamesImpl(~)
            varargout = {'X','Idx','K'};
        end
        
        function varargout = getOutputNamesImpl(obj)
            varargout = {'Y','Th','N'};
            if ~obj.ThresholdOutputPort
                varargout = varargout(~ismember(varargout,{'Th'}));
            end
            if ~obj.NoisePowerOutputPort
                varargout = varargout(~ismember(varargout,{'N'}));
            end
        end
        
        function varargout = isOutputFixedSizeImpl(obj)
            %Fixed if both X and IDX are fixed.
            varargout{1} = propagatedInputFixedSize(obj, 1) && ...
                propagatedInputFixedSize(obj, 2);
            varargout{2} = varargout{1};
            varargout{3} = varargout{1};
        end
        function varargout = getOutputDataTypeImpl(obj)
            if strcmp(obj.OutputFormat,'CUT result')
                varargout{1} = 'logical';
            else
                varargout{1} = propagatedInputDataType(obj,1);
            end
            varargout{2} = propagatedInputDataType(obj,1);
            varargout{3} = varargout{2};
        end
        function varargout = isOutputComplexImpl(~)
            varargout = {false, false, false};
        end
    end

    methods (Access = protected, Abstract)
        NumTrainingCells = getNumTrainingCells(obj) 
    end
    
    methods (Access = private)
        function alpha = calcWGNThresholdFactor(obj)
        %calcWGNThresholdFactor calculate threshold factor for white
        %                       Gaussian noise
        
            % currently we can only calculate threshold when the input is
            % single pulses, with no pulse integration performed.
            % single pulse threshold factor, see [1]
            Nc = getNumTrainingCells(obj); 
            Pfa = obj.ProbabilityFalseAlarm;
            CAThreshold = Nc*(Pfa^(-1/Nc)-1); 
            if obj.Method(1) == 'C'
                alpha = CAThreshold; 
            elseif obj.Method(1) == 'S'
                %alpha = fzero(@(x) SOCAWGNThresholdFactor(x,Nc,Pfa), CAThreshold);
                SOCAWGNThresholdFactor(Nc,Pfa);
                alpha = fzero(@SOCAWGNThresholdFactor,cast(CAThreshold,'double'));
            elseif obj.Method(1) == 'G'
                %alpha = fzero(@(x) GOCAWGNThresholdFactor(x,Nc,Pfa), CAThreshold);
                GOCAWGNThresholdFactor(Nc,Pfa);
                alpha = fzero(@GOCAWGNThresholdFactor,cast(CAThreshold,'double'));
            elseif obj.Method(1) == 'O'
                OSWGNThresholdFactor(Nc,obj.Rank,Pfa);
                
                if ~isempty(coder.target) %Codegen
                    alpha = fzero(@OSWGNThresholdFactor,cast(CAThreshold,'double'));
                else % In MATLAB
                
                   % Try to find a threshold with initial guess. If fzero does 
                   % not converge, fzero throws an error. We configure fzero to 
                   % avoid throwing the error and increase the guess by 10 times 
                   % and try another call to fzero until max iteration reached.
                   % If fzero can still not find a solution, we throw a meaningful 
                   % error.
                
                   flag = -3;
                   initguess = 0.1*CAThreshold;
                   foption = optimset('Display','off');
                   maxiter = 50;
                   iter = 1;
                   while flag == -3 && iter <= maxiter
                       initguess = 10*initguess;
                       [alpha,~,flag] = fzero(@OSWGNThresholdFactor,...
                                              cast(initguess,'double'),foption);
                       iter = iter+1;
                   end
                   if flag == -3
                       error(message('phased:CFARDetector:MaxIterForOS','Rank'));
                   end
                end
            end
        end
    end
end

function c = GOCASOCAThresholdCore(x,N)

temp = 0;
for k = 0:N/2-1
    tempval = gammaln(N/2+k)-gammaln(k+1)-gammaln(N/2);
    temp = temp+exp(tempval)*(2+x/(N/2))^(-k);
end
c = temp*(2+x/(N/2))^(-N/2);

end

function y = SOCAWGNThresholdFactor(varargin)

persistent N pfa;
if nargin > 1
    N = varargin{1};
    pfa = varargin{2};
else
    % Lines in isempty conditionals are dead code
    % The isempty conditionals were added since codegen
    % was complaining that those variables were undefined
    % in certain code paths. Compiler could not detect
    % that the persistent variables were always set first
    % by calling the function with two arguments. The latter
    % workaround was added since codegen did not support
    % anonymous functions.
    
    if isempty(N)
        N = 0;
    end
    if isempty(pfa)
        pfa = 0;
    end
    
    % [1] (7.35)
    y = GOCASOCAThresholdCore(varargin{1},N)-pfa/2;
end
end

function y = GOCAWGNThresholdFactor(varargin)

persistent N pfa;
if nargin > 1
    N = varargin{1};
    pfa = varargin{2};
else
    if isempty(N)
        N = 0;
    end
    if isempty(pfa)
        pfa = 0;
    end
    
    x = varargin{1};
    % [1] (7.37)
    y = (1+x/(N/2))^(-N/2)-GOCASOCAThresholdCore(x,N)-pfa/2;
end
end

function y = OSWGNThresholdFactor(varargin)
persistent N k pfa;
if nargin > 1
    N = varargin{1};
    k = varargin{2};
    pfa = varargin{3};
else
    if isempty(N)
        N = 0;
    end
    if isempty(pfa)
        pfa = 0;
    end
    if isempty(k)
        k = 0;
    end
    x = varargin{1};
    % [1] (7.49)
    temp1 = gammaln(N+1)-gammaln(k)-gammaln(N-k+1);
    c = x+N-k+1;
    if c >= 0
        temp2 = betaln(c,k);
        y = exp(temp1+temp2)-pfa;
    else % betaln used to return inf, now error out. Restore the behavior
        y = nan;
    end
end
end

