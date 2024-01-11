classdef (StrictDefaults)SignalValidator < phased.internal.AbstractVarSizeEngine & matlab.system.mixin.Propagates ...
        & matlab.system.mixin.CustomIcon
    % SignalValidator  Validate input signal

    % Validate the following attributes of the input signal
    % Type: double
    % Dimension: row, column, or matrix

    %   Copyright 2017 The MathWorks, Inc.
    
    %#ok<*EMCLS>
    %#ok<*EMCA>
    %#codegen

    % Public, non-tunable properties
    properties(Nontunable)
        %Dimension     Expected dimensions
        Dimension   
    end

    properties(Constant, Hidden)
        DimensionSet = matlab.system.StringSet(...
            {'Row','Column','Matrix'});
    end
    
    methods
        % Constructor
        function obj = SignalValidator(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
        end
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj,x)
            % Perform one-time calculations, such as computing constants
            setupImpl@phased.internal.AbstractVarSizeEngine(obj);
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
        end

        function y = stepImpl(obj,u) %#ok<INUSL>
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            y = u;
        end

        function resetImpl(obj) %#ok<MANU>
            % Initialize / reset discrete-state properties
        end

        function validateInputsImpl(obj,x)
            % Validate inputs to the step method at initialization
            validateNumChannels(obj,x);
            
            cond =  ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','X','double');
            end
            
            switch obj.Dimension(1)
                case 'C'
                    cond =  ~iscolumn(x) || isempty(x);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeColVector','X');
                    end
                case 'R'
                    cond =  ~isrow(x) || isempty(x);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeRowVector','X');
                    end
                case 'M'
                    cond = ~ismatrix(x) || isempty(x);
                    if cond
                        coder.internal.errorIf(cond, ...
                             'MATLAB:system:inputMustBeMatrix','X');
                    end
            end
                    
        end

        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractVarSizeEngine(obj);

            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end

        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            % Set properties in object obj to values in structure s

            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        %% Simulink functions
        function flag = isInputSizeLockedImpl(obj,index) %#ok<INUSD>
            % Return true if input size is not allowed to change while
            % system is running
            flag = false;
        end

        function out = getOutputSizeImpl(obj)
            % Return size for each output port
            out = propagatedInputSize(obj,1);
        end

        function out = getOutputDataTypeImpl(obj)
            % Return data type for each output port
            out = propagatedInputDataType(obj,1);
        end

        function out = isOutputComplexImpl(obj)
            % Return true for each output port with complex data
            out = propagatedInputComplexity(obj,1);
        end

        function icon = getIconImpl(obj) %#ok<MANU>
            % Define icon for System block
            icon = 'Signal\nValidator'; 
        end

        function name = getInputNamesImpl(obj) %#ok<MANU>
            % Return input port names for System block
            name = 'x';
        end
    end

    methods(Static, Access = protected)
        %% Simulink customization functions
        function header = getHeaderImpl
            % Define header panel for System block dialog
            header = matlab.system.display.Header(mfilename('class'));
        end

        function group = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
            group = matlab.system.display.Section(mfilename('class'));
        end
    end
end
