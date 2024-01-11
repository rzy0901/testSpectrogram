classdef (Sealed, StrictDefaults) SimulinkSignalPadder < phased.internal.AbstractVarSizeEngine & ...
     matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.
 
%   Copyright 2015-2018 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable, Logical)
        %SignalLengthOutputPort     Output input signal length
        %   Set this property to true to output input signal's length. The
        %   signal length is defined as number of rows of the input signal.
        %   Set this property to false to not output input signal's length.
        %   The default value of this property is false.
        SignalLengthOutputPort = false
    end
    
    properties (Access = private)
        pSignalLength
    end

    methods

        function obj = SimulinkSignalPadder(varargin)
            
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
        
        function validateInputsImpl(obj,x)
            validateattributes(x,{'double'},{'nonempty','nonnan','2d','finite'},...
                '','X');
            validateNumChannels(obj,x);
        end
        
        function setupImpl(obj,x)
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            obj.pSignalLength = getPropagatedNumInputSamples(obj,x);
            
        end
        
        function [y,sigsz] = stepImpl(obj,x)

            y = zeros(obj.pSignalLength,obj.pValidatedNumInputChannels,'like',x);
            sigsz = size(x,1);
            y(1:sigsz,:) = x;
           
        end 
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);
            if isLocked(obj)
                s.pSignalLength = obj.pSignalLength;
            end
        end
        
        function num = getNumOutputsImpl(obj) 
            if obj.SignalLengthOutputPort
                num = 2;
            else
                num = 1;
            end
        end
    end

    methods (Static,Hidden,Access=protected)
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:SimulinkSignalPadderTitle')),...
              'Text',getString(message('phased:library:block:SimulinkSignalPadderDesc')));
        end
        
    end
    methods (Access = protected) %For Simulink propagation and mask
        function varargout = getOutputNamesImpl(obj)
            if obj.SignalLengthOutputPort
                varargout = {'Xp','Sz'};
            else
                varargout = {'Xp'};
            end
        end
        
        function varargout = getInputNamesImpl(obj) %#ok<MANU>
            varargout = {'X'};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Signal Padder');
        end
        
        function flag = isInputSizeLockedImpl(obj,~) %#ok<INUSD>
            flag = false;
        end
        
        function varargout = getOutputSizeImpl(obj)
            if obj.SignalLengthOutputPort
                varargout = {propagatedInputSize(obj,1),...
                    [1 1]};
            else
                varargout = {propagatedInputSize(obj,1)};
            end
        end
        
        function varargout = isOutputFixedSizeImpl(obj) 
            %varargout{1} = propagatedInputFixedSize(obj, 1);
            if obj.SignalLengthOutputPort
                varargout = {true,true};
            else
                varargout = {true};
            end
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout = {'double','double'};
        end
        function varargout = isOutputComplexImpl(~)
            varargout = {true,false};
        end
    end
end




