classdef (Sealed,StrictDefaults) BarrageJammer < phased.internal.AbstractBarrageJammer
%BarrageJammer Barrage jammer
%   H = phased.BarrageJammer creates a barrage jammer System object, H.
%   This object generates a complex white Gaussian noise jamming signal.
%
%   H = phased.BarrageJammer(Name,Value) creates a barrage jammer object, 
%   H, with the specified property Name set to the specified Value. You 
%   can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   H = phased.BarrageJammer(E,Name,Value) creates a barrage jammer
%   object, H, with the ERP property set to E and other specified property
%   Names set to the specified Values. 
%
%   Step method syntax:
%
%   Y = step(H) returns a column vector, Y, that is a complex white
%   Gaussian noise jamming signal. The power of the jamming signal is
%   specified by the ERP property. The length of the jamming signal is
%   specified by the SamplesPerFrame property. This option is available
%   when the SamplesPerFrameSource property is 'Property'.
%
%   Y = step(H,N) returns the jamming signal with length N. This option is
%   available when the SamplesPerFrameSource property is 'Input port'.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H) and y = H() are
%   equivalent.
%
%   BarrageJammer methods:
%
%   step     - Generate noise jamming signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a barrage jammer object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset random number generator for noise generation
%
%   BarrageJammer properties:
%
%   ERP                   - Effective radiated power 
%   SamplesPerFrameSource - Source of number of samples per frame
%   SamplesPerFrame       - Number of samples per frame
%   SeedSource            - Source of seed for random number generator
%   Seed                  - Seed for random number generator
%   
%   % Example:
%   %   Create a barrage jammer with an effective radiated power of 1000 
%   %   Watts and plot the magnitude of its output.
%
%   jammer = phased.BarrageJammer(1000);
%   x = jammer();
%   plot(abs(x)); xlabel('Samples'); ylabel('Magnitude');
%
%   See also phased, phased.RadarTarget, phased.Platform.

%   Copyright 2009-2016 The MathWorks, Inc.

% References:
% [1] James Ward, "Space-Time Adaptive Processing for Airborne
%     Radar". MIT Lincoln Lab tech report 1015, 1994.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %SamplesPerFrameSource  Source of number of samples per frame
        %   Specify how to determine the number of samples of the output
        %   jamming signal as one of 'Property' | 'Input port', where the
        %   default is 'Property'. When you set this property to
        %   'Property', the number of samples of the jamming signal is
        %   determined by the value of the SamplesPerFrame property. When
        %   you set this property to 'Input port', the number of samples of
        %   the jamming signal is determined by the input argument.
        SamplesPerFrameSource = 'Property';
    end
    
    properties(Constant, Hidden)
        SamplesPerFrameSourceSet = dsp.CommonSets.getSet('PropertyOrInputPort');
    end
    
    methods

        function obj = BarrageJammer(varargin)
            
            obj@phased.internal.AbstractBarrageJammer(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        function validateInputsImpl(obj, N)
            if getNumInputs(obj) == 1
                cond =  ~isa(N,'double') ;
                if cond
                    coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType','N','double');
                end
                cond =  ~isscalar(N);
                if cond
                    coder.internal.errorIf(cond, ...
                                 'MATLAB:system:inputMustBeScalar','N');
                end
                cond =  ~isreal(N);
                if cond
                    coder.internal.errorIf(cond, ...
                                 'phased:step:NeedReal', 'N');
                end
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractBarrageJammer(obj);
        end

        function loadObjectImpl(obj,s,~)
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function y = stepImpl(obj,Narg)
            if obj.pSamplesPerFrameViaProp
                N = obj.SamplesPerFrame;
            else
                N = Narg;
                sigdatatypes.validateIndex(N,'step','N');           
            end
            
            y = step(obj.cNoiseSource,obj.ERP,[N 1]);
           
        end
        
    end

    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            dSPF = matlab.system.display.internal.Property(...
                'SamplesPerFrameSource','IsGraphical',false);
            dSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false,'UseClassDefault',false,...
                'Default','Property');
            dSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false,'UseClassDefault',false,...
                'Default','randi(65535,1)');
            props = {...
                'ERP',...
                dSPF,...
                'SamplesPerFrame',...
                dSeedSource,...
                dSeed };
            
            groups = matlab.system.display.Section('Title', 'Parameters', ...
                                                   'PropertyList', props);
        end
    end
    methods (Access = protected) 
        function varargout = getOutputNamesImpl(~)
            varargout = {''};
        end
        
        function varargout = getInputNamesImpl(obj)
            if ~strcmp(obj.SamplesPerFrameSource, 'Property')
                varargout = {'N'};
            else
                varargout = {};
            end      
        end
        function flag = isInputSizeLockedImpl(~,~)
            flag = true;
        end
        
        %SamplesPerFrameSource from Input Port is disabled for Simulink
        function varargout = getOutputSizeImpl(obj)
           varargout{1} = [obj.SamplesPerFrame,1];
        end
        function varargout = isOutputFixedSizeImpl(~)
            varargout{1} = true;
        end
    end
    
    methods (Static,Hidden)
        function a = getAlternateBlock
            a = 'phasedenvlib/Barrage Jammer';
        end
    end
end




