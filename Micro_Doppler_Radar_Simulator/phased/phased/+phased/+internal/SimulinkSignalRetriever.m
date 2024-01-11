classdef (Sealed, StrictDefaults) SimulinkSignalRetriever < phased.internal.AbstractVarSizeEngine & ...
     matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.
 
%   Copyright 2015-2018 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods

        function obj = SimulinkSignalRetriever(varargin)
            
            setProperties(obj, nargin, varargin{:});
        end
        
    end
    
    methods (Access = protected)
        
        function loadObjectImpl(obj,s,~)
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function validateInputsImpl(obj,x,szref)
            validateattributes(x,{'double'},{'nonempty','nonnan','2d','finite'},...
                '','X');
            validateNumChannels(obj,x);
            validateattributes(szref,{'double'},{'nonempty','nonnan','2d','finite'},...
                '','SzRef');
            validateNumChannels(obj,szref);
        end
        
        function setupImpl(obj,x,~)
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
        end
        
        function num = getNumInputsImpl(obj) %#ok<MANU>
            num = 2;
        end
        
        function y = stepImpl(obj,x,sz_in) %#ok<INUSL>

            sz = sz_in(end);
            cond = (size(x,1) < sz);
            if cond
                coder.internal.errorIf(cond,'phased:phased:expectedMoreRows','X',sz);
            end
            y = x(1:sz,:);
           
        end        
    end

    methods (Static,Hidden,Access=protected)
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:SimulinkSignalRetrieverTitle')),...
              'Text',getString(message('phased:library:block:SimulinksignalRetrieverDesc')));
        end
        
    end
    
    methods (Access = protected) %For Simulink propagation and mask
        function varargout = getOutputNamesImpl(~)
            varargout = {''};
        end
        
        function varargout = getInputNamesImpl(obj) %#ok<MANU>
            varargout = {'X','Sz'};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Signal Retriever');
        end
        
        function flag = isInputSizeLockedImpl(obj,ind)  %#ok<INUSL>
            if ind == 1
                flag = true;
            else
                flag = false;
            end
        end
        
        function varargout = getOutputSizeImpl(obj)
            varargout{1} = propagatedInputSize(obj,1);
        end
        function varargout = isOutputFixedSizeImpl(obj) %#ok<MANU>
            varargout{1} = false;
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout{1} = 'double';
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = true;
        end
    end
end




