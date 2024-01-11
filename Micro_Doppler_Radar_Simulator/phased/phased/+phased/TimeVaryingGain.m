classdef (Sealed,StrictDefaults) TimeVaryingGain < phased.internal.AbstractVarSizeEngine & ...
        matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates
%TimeVaryingGain    Time varying gain control
%   H = phased.TimeVaryingGain creates a time varying gain control System
%   object, H. The object applies a time varying gain to the input signal
%   to compensate for the signal power loss due to the range.
%
%   H = phased.TimeVaryingGain(Name,Value) creates a time varying gain
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(H,X) applies time varying gains to the input signal X. The
%   process equalizes power levels across all samples to match a given
%   reference range. The compensated signal is returned in Y. X can be
%   either a column vector, a matrix or a cube. The gain is applied to each
%   column in X independently. The number of rows in X must be less than or
%   equal to the length of the loss vector specified in the RangeLoss
%   property. Y has the same dimensionality as X.
%
%   Y = step(H,X,L) specifies the range loss, L (in dB), corresponding to
%   each sample as a column vector. The number of rows in L must match the
%   number of rows in X. Each entry in L is the range loss for each sample
%   in X along the first dimension.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   TimeVaryingGain methods:
%
%   step     - Apply time varying gains to the input signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create time varying gain object with same property values
%   isLocked - Locked status (logical)
%
%   TimeVaryingGain properties:
%
%   RangeLossSource - Source of range losses
%   RangeLoss       - Range losses
%   ReferenceLoss   - Reference range loss
%
%   This System object supports single and double precision for input data,
%   properties and arguments. If the input data X is single precision, the
%   output data is single precision. If the input data X is double
%   precision, the output data is double precision. The precision of the
%   output is independent of the precision of the properties and other
%   arguments.
%
%   % Example:
%   %   Apply time varying gain to a signal to compensate for signal power
%   %   loss due to range. 
%
%   rngloss = 10:22;    refloss = 16;   % in dB
%   t = (1:length(rngloss))';  x = 1./db2mag(rngloss(:));
%   tvg = phased.TimeVaryingGain('RangeLoss',rngloss,...
%          'ReferenceLoss',refloss);
%   y = tvg(x);
%
%   % Plot signals
%   tref = find(rngloss==refloss);
%   stem([t t],[abs(x) abs(y)]);
%   hold on; stem(tref,x(tref),'filled','r'); hold off;
%   xlabel('Time (s)'); ylabel('Magnitude (V)'); grid on;
%   legend('Before time varying gain',...
%          'After time varying gain',...
%          'Reference range');
%
%   See also phased, phased.MatchedFilter, pulsint.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001
%   [2] Byron Edde, Radar: Principles, Technology, Applications, Prentice
%       Hall, 1993


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %RangeLossSource  Source of range losses
        %   Specify the source of range losses as one of 'Property' |
        %   'Input port', where the default is 'Property'. When you specify
        %   the RangeLossSource as 'Property', the range loss for each
        %   sample is given in the RangeLoss property. When you specify the
        %   RangeLossSource as 'Input port', the range loss for each
        %   sample is specified as an input argument.
        RangeLossSource = 'Property'
        %RangeLoss  Range losses (dB)
        %   Specify the loss (in dB) due to the range for each sample in
        %   the input signal as a vector. The default value of this
        %   property is 0.
        RangeLoss = 0;
        %ReferenceLoss  Reference range loss (dB)
        %   Specify the loss (in dB) at a given reference range as a
        %   scalar. The default value of this property is 0.
        ReferenceLoss = 0;
    end
    
    properties(Constant, Hidden)
        RangeLossSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    properties (Access = private)
        pCubeDim = [-1 -1 -1];
        pNumMtxChannels;
        pNumSamplesPerChannel;
    end

    properties (Access = private, Nontunable)
        pInput3DFlag 
        pNormalizeVector;
    end
    
    methods
        function set.RangeLoss(obj,val)
            validateattributes( val, { 'double','single' }, { 'real', 'vector' }, '', 'RangeLoss');
            obj.RangeLoss = val;
        end
        function set.ReferenceLoss(obj,val)
            validateattributes( val, { 'double','single' }, { 'finite', 'real', 'scalar' }, '', 'ReferenceLoss');
            obj.ReferenceLoss = val;
        end
    end

    methods 
        function obj = TimeVaryingGain(varargin)
        %TimeVaryingGain Constructor for phased.TimeVaryingGain class
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods (Access = protected)
        
        function num = getNumInputsImpl(obj)
            if strcmp(obj.RangeLossSource,'Property')
                num = 1;
            else
                num = 2;
            end
        end
        
        function validateInputsImpl(obj,x,L)
            size_x = size(x);
            cond =  ~isa(x,'float');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','X','float');
            end
            cond =  numel(size_x) > 3;
            if cond
                coder.internal.errorIf(cond, ...
                     'phased:TimeVaryingGain:TooManyDimensions');            
            end
            
            if strcmp(obj.RangeLossSource,'Property')
                Nrow = numel(obj.RangeLoss);
                cond =  (size_x(1) > Nrow);
                if cond
                    coder.internal.errorIf(cond, ...
                         'phased:phased:expectedLessRows', 'X', Nrow);
                end
            else
                validateattributes(L,{'double','single'},{'real','nonempty','nonnan',...
                    'column','numel',size_x(1)},'','L');
            end
            
            validateNumChannels(obj,x);
            if ndims(x) == 3 && obj.pCubeDim(3) ~= -1
                validateNumPages(obj,x,obj.pCubeDim(3));
            end
            
        end
        
        function setupImpl(obj,x,~)
            if strcmp(obj.RangeLossSource,'Property')
                rngloss = cast(obj.RangeLoss(:)-obj.ReferenceLoss,'like',x);
                obj.pNormalizeVector = db2mag(rngloss);
            end
            
            sz_x = size(x);
            obj.pNumInputChannels = sz_x(2);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        
            if numel(sz_x) == 3
                obj.pInput3DFlag = true;
                obj.pCubeDim = sz_x;
                obj.pNumMtxChannels = sz_x(2)*sz_x(3);
            else
                obj.pInput3DFlag = false;
                obj.pNumMtxChannels = sz_x(2);
            end
            obj.pNumSamplesPerChannel = sz_x(1);
        end
        
        function processInputSizeChangeImpl(obj,x,~)
            sz_x = size(x);
            if obj.pInput3DFlag
                obj.pCubeDim = sz_x;
                obj.pNumMtxChannels = sz_x(2)*sz_x(3);
            else
                obj.pNumMtxChannels = sz_x(2);
            end
        end

        function flag = isInputComplexityLockedImpl(~,~) 
            flag = false;
        end

        function flag = isOutputComplexityLockedImpl(~,~) 
            flag = false;
        end

        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractVarSizeEngine(obj);
        end
        
        function flag = isInputSizeLockedImpl(obj,~) %#ok<INUSD>
            flag = false;
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if strcmp(prop,'RangeLoss') && ...
                    strcmp(obj.RangeLossSource,'Input port')
                flag = true;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
            if isLocked(obj)
                s.pNormalizeVector = obj.pNormalizeVector;
                s.pInput3DFlag = obj.pInput3DFlag;
                s.pCubeDim = obj.pCubeDim;
                s.pNumMtxChannels = obj.pNumMtxChannels;
                s.pNumSamplesPerChannel = obj.pNumSamplesPerChannel;
            end
        end

        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function yout = stepImpl(obj,xin,L)
        %   compensates the radar received signal X so that all the echos
        %   from the target have the same power level as if the targets are
        %   located at a given reference range. The compensated signal is
        %   returned in Y. If X is a matrix or a cube, the operation is
        %   along the first dimension, i.e., the operation works on each
        %   column independently.

            classtouse = class(xin);
            if obj.pInput3DFlag
                sigsize = size(xin);
                % input is a cube
                % need to reshape, process and transform it back
                x = reshape(xin,sigsize(1),[]);
            else
                x = xin;
            end
            
            if obj.RangeLossSource(1) == 'P'
                normvec = obj.pNormalizeVector(1:size(x,1));
            else
                normvec = db2mag(cast((L-obj.ReferenceLoss),classtouse));
            end
            y = bsxfun(@times,normvec,x);

            if obj.pInput3DFlag
                % input is a cube, so we need to transform the result back
                % to its original shape.
                yout = reshape(y,sigsize);
            else
                yout = y;
            end
        end
    end    
  
    methods (Static,Hidden,Access=protected)  
      function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:TimeVaryingGainTitle')),...
              'Text',getString(message('phased:library:block:TimeVaryingGainDesc')));
      end
    end
    
    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)  
            if strcmp(obj.RangeLossSource,'Property')
                varargout = {'X'};
            else
                varargout = {'X','L'};
            end
        end
        
        function varargout = getOutputNamesImpl(obj) %#ok<MANU>
            varargout = {''};
        end
      
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('TVG');
        end
        function varargout = getOutputSizeImpl(obj)
           varargout{1} = propagatedInputSize(obj,1);
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj, 1);
        end
        function varargout = getOutputDataTypeImpl(obj)
            varargout{1} = propagatedInputDataType(obj,1);
        end
        function varargout = isOutputComplexImpl(obj)
            varargout{1} = propagatedInputComplexity(obj,1);
        end
    end
end


