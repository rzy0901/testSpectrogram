classdef (Sealed, StrictDefaults) SimulinkBarrageJammer < phased.internal.AbstractBarrageJammer & ...
     matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.
 
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
%   BarrageJammer methods:
%
%   step     - Generate noise jamming signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a barrage jammer object with same property values
%   isLocked - Locked status (logical)
%   reset    - Reset random number generator for noise generation
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
%   Hjammer = phased.BarrageJammer(1000);
%   x = step(Hjammer);
%   plot(abs(x)); xlabel('Samples'); ylabel('Magnitude');
%
%   See also phased, phased.RadarTarget, phased.Platform.

%   Copyright 2009-2012 The MathWorks, Inc.

% References:
% [1] James Ward, "Space-Time Adaptive Processing for Airborne
%     Radar". MIT Lincoln Lab tech report 1015, 1994.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %SamplesPerFrameSource  Source of number of samples per frame
        %   Specify how to determine the number of samples of the output
        %   jamming signal as one of 'Property' | 'Derive from reference
        %   input port', where the default is 'Property'. When you set this
        %   property to 'Property', the number of samples of the jamming
        %   signal is determined by the value of the SamplesPerFrame
        %   property. When you set this property to 'Derive from reference
        %   input port', the number of samples of the jamming signal is
        %   determined by the number of rows of the reference input signal.
        SamplesPerFrameSource = 'Property';
    end
    
    properties(Constant, Hidden)
        SamplesPerFrameSourceSet = matlab.system.StringSet(...
            {'Property','Derive from reference input port'});
    end
    
    methods

        function obj = SimulinkBarrageJammer(varargin)
            
            obj@phased.internal.AbstractBarrageJammer(varargin{:});
        end
        
    end
    
    methods (Access = protected)
        
        function validateInputsImpl(obj, x) 
            if getNumInputs(obj) == 1
                cond =  ~isa(x,'double') ;
                if cond
                    coder.internal.errorIf(cond, ...
                                 'MATLAB:system:invalidInputDataType','Ref','double');
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
        
        
        function y = stepImpl(obj,x)

            if obj.pSamplesPerFrameViaProp
                N = obj.SamplesPerFrame;
            else
                N = size(x,1);
            end
            
            y = step(obj.cNoiseSource,obj.ERP,[N 1]);
           
        end        
    end

    methods (Static,Hidden,Access=protected)
        function header = getHeaderImpl
          header = matlab.system.display.Header(...
              'Title',getString(message('phased:library:block:BarrageJammerTitle')),...
              'Text',getString(message('phased:library:block:BarrageJammerDesc')));
        end
        
        function groups = getPropertyGroupsImpl
            dSeedSource = matlab.system.display.internal.Property(...
                'SeedSource','IsGraphical',false,'UseClassDefault',false,...
                'Default','Property');
            dSeed = matlab.system.display.internal.Property(...
                'Seed','IsGraphical',false,'UseClassDefault',false,...
                'Default','randi(65535,1)');
            props = {...
                'ERP',...
                'SamplesPerFrameSource',...
                'SamplesPerFrame',...
                dSeedSource,...
                dSeed };
            
            groups = matlab.system.display.Section('Title', 'Parameters', ...
                                                   'PropertyList', props);
        end
    end
    methods (Access = protected) %For Simulink propagation and mask
        function varargout = getOutputNamesImpl(~)
            varargout = {''};
        end
        
        function varargout = getInputNamesImpl(obj)
            if ~strcmp(obj.SamplesPerFrameSource, 'Property')
                varargout = {'Ref'};
            else
                varargout = {};
            end      
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Barrage');
        end
        
        function flag = isInputSizeLockedImpl(obj,~)
            if strcmp(obj.SamplesPerFrameSource, 'Property')
                flag = true;
            else
                flag = false;
            end
        end
        
        %SamplesPerFrameSource from Input Port is disabled for Simulink
        function varargout = getOutputSizeImpl(obj)
            if strcmp(obj.SamplesPerFrameSource, 'Property')
                varargout{1} = [obj.SamplesPerFrame,1];
            else
                varargout{1} = propagatedInputSize(obj,1);
            end
        end
        function varargout = isOutputFixedSizeImpl(obj)
            if strcmp(obj.SamplesPerFrameSource, 'Property')
                varargout{1} = true;
            else
                varargout{1} = propagatedInputFixedSize(obj, 1);
            end
        end
        function varargout = getOutputDataTypeImpl(~)
            varargout{1} = 'double';
        end
        function varargout = isOutputComplexImpl(~)
            varargout{1} = true;
        end
    end
end




