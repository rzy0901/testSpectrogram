classdef GridSpacingProperty < matlab.system.display.internal.Property
    %phased.internal.GridSpacingProperty   Reserved for MathWorks internal use only
    
    %   Copyright 2014 The MathWorks, Inc.
    
    methods
        function obj = GridSpacingProperty(varargin)
            obj@matlab.system.display.internal.Property(varargin{:});
        end
        
        function addDialogValue(obj, paramValue, builder)
            
            % GridSpacing acts like a string literal if value is 'Auto'
            % but otherwise it acts like a numerical property
            if strcmp(paramValue, 'Auto')
                builder.addStringParameterValue(obj.Name, 'Auto');
            else
                builder.addLiteralParameterValue(obj.Name, paramValue);
            end
        end
        
        function addParsedExpression(obj, expression, builder)
            if strcmp(expression, '''Auto''') 
                expression = 'Auto';
            end
            obj.addDialogValue(expression, builder);
        end
    end
end