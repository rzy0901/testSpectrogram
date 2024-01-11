classdef (StrictDefaults)SensorArrayDialog < matlab.system.display.internal.Presenter

    %   Copyright 2014 The MathWorks, Inc.

    properties (Nontunable)
        SensorSource = 'Array (no subarrays)';
        SensorElement;
        SensorArray;
        SensorPartitionedArray;
        SensorReplicatedSubArray;
        SensorReplicatedArray;
        SensorExpression = '<matlab expression>';
    end

    properties(Constant, Hidden)
        SensorSourceSet = matlab.system.StringSet({'Single element', 'Array (no subarrays)', 'Partitioned array', 'Replicated subarray', 'MATLAB expression'});
    end

    methods
        function obj = SensorArrayDialog(varargin)
            setProperties(obj, nargin, varargin{:});

            % Set default values, which are used when the user updates
            % SensorSource
            if isempty(obj.SensorArray)
                obj.SensorArray = phased.ULA;
            end
            if isempty(obj.SensorElement)
                obj.SensorElement = phased.IsotropicAntennaElement;
            end
            if isempty(obj.SensorPartitionedArray)
                obj.SensorPartitionedArray =  phased.PartitionedArray;
            end
            if isempty(obj.SensorReplicatedSubArray)
                obj.SensorReplicatedSubArray = phased.ULA;
            end
            if isempty(obj.SensorReplicatedArray)
                obj.SensorReplicatedArray = phased.ReplicatedSubarray;
            end
        end
    end

    methods(Access=protected)
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            switch (prop)
                case 'SensorElement'
                    flag = strcmp(obj.SensorSource, 'MATLAB expression');
                case 'SensorArray'
                    flag = ~strcmp(obj.SensorSource, 'Array (no subarrays)');
                case 'SensorPartitionedArray'
                    flag = ~strcmp(obj.SensorSource, 'Partitioned array');
                case 'SensorReplicatedSubArray'
                    flag = ~strcmp(obj.SensorSource, 'Replicated subarray');
                case 'SensorReplicatedArray'
                    flag = ~strcmp(obj.SensorSource, 'Replicated subarray');
                case 'SensorExpression'
                    flag = ~strcmp(obj.SensorSource, 'MATLAB expression');
            end
        end
    end

    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl(sensorType)
            
            %if called by disp there will be no argument
            %use a default one
            if nargin < 1
                sensorType = 'all';
            end

            [cSensorSource, cSensorElement, cSensorArray, ...
             cSensorPartitionedArray, cSensorReplicatedSubArray, cSensorReplicatedArray] =  ...
                phased.internal.SensorArrayDialog.getSensorClassSet(sensorType);

            pSensorSource = matlab.system.display.internal.Property('SensorSource', ...
                            'StringSetValues', cSensorSource, ...
                                     'Description', 'Specify sensor array as');

            pSensorExpression = matlab.system.display.internal.Property('SensorExpression', ...
                'Description', 'Expression');


            pSensorElement = matlab.system.display.internal.Property('SensorElement', ...
                'Description', 'Element type', ...
                'ClassStringSet', cSensorElement);
            sSensorElement = matlab.system.display.Section('Title', 'Element', ...
                'PropertyList', {pSensorElement});

            pSensorArray = matlab.system.display.internal.Property('SensorArray', ...
                'Description', 'Geometry', ...
                'ClassStringSet', cSensorArray);
             arrayLabels = cSensorArray.Labels;
            if length(arrayLabels) == 1
                arraySectionLabel = arrayLabels{1};
            else
                arraySectionLabel = 'Array';
            end
            sSensorArray = matlab.system.display.Section('Title', arraySectionLabel, ...
                'PropertyList', {pSensorArray});

            pSensorPartitionedArray = matlab.system.display.internal.Property('SensorPartitionedArray', ...
                'Description', 'Array', ...
                'ClassStringSet', cSensorPartitionedArray);
            sSensorPartitionedArray = matlab.system.display.Section('Title', 'Array', ...
                'PropertyList', {pSensorPartitionedArray});

            pSensorReplicatedSubArray = matlab.system.display.internal.Property('SensorReplicatedSubArray', ...
                'Description', 'Geometry', ...
                'ClassStringSet', cSensorReplicatedSubArray);
            sSensorReplicatedSubArray = matlab.system.display.Section('Title', 'Subarray', ...
                'PropertyList', {pSensorReplicatedSubArray});

            pSensorReplicatedArray = matlab.system.display.internal.Property('SensorReplicatedArray', ...
                'Description', 'Array', ...
                'ClassStringSet', cSensorReplicatedArray);
            sSensorReplicatedArray = matlab.system.display.Section('Title', 'Array', ...
                'PropertyList', {pSensorReplicatedArray});

            groups = matlab.system.display.SectionGroup('Title', 'Sensor Array', ...
                'PropertyList', {pSensorSource, pSensorExpression}, ...
                'Sections', [sSensorElement, sSensorArray, sSensorPartitionedArray, sSensorReplicatedSubArray, sSensorReplicatedArray]);
        end

        function [cSensorSource, cSensorElement, cSensorArray, cSensorPartitionedArray, ...
                  cSensorReplicatedSubArray, cSensorReplicatedArray] = getSensorClassSet(query)
            query = validatestring(query,{'all','antenna','array','subarray','ula','ura','ulauca'});

            if strcmp(query,'antenna')
                query = 'all';
                elementValues = { ...
                    'phased.IsotropicAntennaElement', ...
                    'phased.CosineAntennaElement', ...
                    'phased.CustomAntennaElement'};

                elementLabels = { ...
                    'Isotropic Antenna', ...
                    'Cosine Antenna', ...
                    'Custom Antenna'};
            
            else
                elementValues = { ...
                    'phased.IsotropicAntennaElement', ...
                    'phased.CosineAntennaElement', ...
                    'phased.CustomAntennaElement', ...
                    'phased.OmnidirectionalMicrophoneElement', ...
                    'phased.CustomMicrophoneElement'};

                elementLabels = { ...
                    'Isotropic Antenna', ...
                    'Cosine Antenna', ...
                    'Custom Antenna', ...
                    'Omni Microphone', ...
                    'Custom Microphone'};
            end
             
            if strcmp(query,'array')
                cSensorSource = {'Array (no subarrays)', 'MATLAB expression'};
                arrayValues = {'phased.ULA', 'phased.URA', 'phased.UCA','phased.ConformalArray'};
                arrayLabels = {'ULA','URA','UCA','Conformal Array'};
            elseif strcmp(query,'ula')
                cSensorSource = {'Array (no subarrays)', 'MATLAB expression'};
                arrayValues = {'phased.ULA'};
                arrayLabels = {'ULA'};
            elseif strcmp(query,'ulauca')
                cSensorSource = {'Array (no subarrays)', 'MATLAB expression'};
                arrayValues = {'phased.ULA','phased.UCA'};
                arrayLabels = {'ULA','UCA'};
            elseif strcmp(query,'ura')
                cSensorSource = {'Array (no subarrays)', 'MATLAB expression'};
                arrayValues = {'phased.URA'};
                arrayLabels = {'URA'};
            elseif strcmp(query,'subarray')
                cSensorSource = {'Array (no subarrays)', 'Partitioned array', ...
                                 'Replicated subarray', 'MATLAB expression'};
                arrayValues = {'phased.ULA', 'phased.URA', 'phased.UCA','phased.ConformalArray'};
                arrayLabels = {'ULA','URA','UCA','Conformal Array'};               
            else % all
                cSensorSource = {'Single element', 'Array (no subarrays)', 'Partitioned array', ...
                                 'Replicated subarray', 'MATLAB expression'};
                arrayValues = {'phased.ULA', 'phased.URA', 'phased.UCA','phased.ConformalArray'};
                arrayLabels = {'ULA','URA','UCA','Conformal Array'};               
            end

            cSensorElement = matlab.system.display.internal.ClassStringSet(...
                elementValues, ...
                'PropertiesTitle', 'Element', ...
                'NestDisplay', false, ...
                'Labels', elementLabels);
            cSensorArray = matlab.system.display.internal.ClassStringSet(...
                arrayValues, ...
                'PropertiesTitle', 'Array', ...
                'NestDisplay', false, ...
                'Labels', arrayLabels);

            cSensorPartitionedArray = matlab.system.display.internal.ClassStringSet(...
                {'phased.PartitionedArray'}, ...
                'PropertiesTitle', 'Array', ...
                'NestDisplay', false);

            cSensorReplicatedSubArray = matlab.system.display.internal.ClassStringSet(...
                arrayValues, ...
                'PropertiesTitle', 'SubArray', ...
                'Labels', arrayLabels,...
                'NestDisplay', false);

            cSensorReplicatedArray = matlab.system.display.internal.ClassStringSet(...
                {'phased.ReplicatedSubarray'}, ...
                'PropertiesTitle', 'Array', ...
                'NestDisplay', false);

        end

    end

    methods(Static)
        function dialogExpression = getDialogExpression(parameterExpression) % On dialog open

            dialogExpressionBuilder = matlab.system.ui.ConstructorBuilder('phased.internal.SensorArrayDialog');

            paramBuilder = matlab.system.ui.ConstructorBuilder.parse(parameterExpression);
            if isempty(paramBuilder) 
                dialogExpressionBuilder.addStringParameterValue('SensorSource','MATLAB expression');
                dialogExpressionBuilder.addStringParameterValue('SensorExpression', parameterExpression);
                dialogExpression = dialogExpressionBuilder.buildExpression;
                return;
            end

            switch (paramBuilder.ClassName)
                case {'phased.ULA', 'phased.URA', 'phased.UCA', 'phased.ConformalArray'}
                    dialogExpressionBuilder.addStringParameterValue('SensorSource','Array (no subarrays)');
                    dialogExpressionBuilder.addObjectParameterValue('SensorArray', paramBuilder);
                    if paramBuilder.isParameter('Element')
                        dialogExpressionBuilder.addObjectParameterValue('SensorElement', paramBuilder.getParameterBuilder('Element'));
                    end
                case 'phased.PartitionedArray'
                    dialogExpressionBuilder.addStringParameterValue('SensorSource','Partitioned array');
                    dialogExpressionBuilder.addObjectParameterValue('SensorPartitionedArray', paramBuilder);
                    if paramBuilder.isParameter('Array')
                        arrayBuilder = paramBuilder.getParameterBuilder('Array');
                        if arrayBuilder.isParameter('Element')
                            dialogExpressionBuilder.addObjectParameterValue('SensorElement', arrayBuilder.getParameterBuilder('Element'));
                        end
                    end
                case 'phased.ReplicatedSubarray'
                    dialogExpressionBuilder.addStringParameterValue('SensorSource','Replicated subarray');
                    dialogExpressionBuilder.addObjectParameterValue('SensorReplicatedArray', paramBuilder);
                    if paramBuilder.isParameter('Subarray')
                        arrayBuilder = paramBuilder.getParameterBuilder('Subarray');
                        dialogExpressionBuilder.addObjectParameterValue('SensorReplicatedSubArray', arrayBuilder);
                        if arrayBuilder.isParameter('Element')
                            dialogExpressionBuilder.addObjectParameterValue('SensorElement', arrayBuilder.getParameterBuilder('Element'));
                        end
                    end
                case {'phased.IsotropicAntennaElement', 'phased.CosineAntennaElement', 'phased.CustomAntennaElement', ...
                       'phased.OmnidirectionalMicrophoneElement', 'phased.CustomMicrophoneElement'}
                    dialogExpressionBuilder.addStringParameterValue('SensorSource','Single element');
                    dialogExpressionBuilder.addObjectParameterValue('SensorElement', paramBuilder);
                otherwise
                    dialogExpressionBuilder.addStringParameterValue('SensorSource','MATLAB expression');
                    dialogExpressionBuilder.addStringParameterValue('SensorExpression', parameterExpression);
                    dialogExpression = dialogExpressionBuilder.buildExpression;
                    return;
            end

            dialogExpression = dialogExpressionBuilder.buildExpression;
        end

        function parameterExpression = getParameterExpression(dialogExpression) % On dialog update

            dialogExpressionBuilder = matlab.system.ui.ConstructorBuilder.parse(dialogExpression);
            switch dialogExpressionBuilder.getLiteralParameterValue('SensorSource')
                case 'Array (no subarrays)'
                    if dialogExpressionBuilder.isParameter('SensorArray')
                        arrayBuilder = dialogExpressionBuilder.getParameterBuilder('SensorArray');
                        arrayBuilder.addObjectParameterValue('Element', dialogExpressionBuilder.getParameterValue('SensorElement'));
                        parameterExpression = arrayBuilder.buildExpression;
                    else
                        arrayBuilder = matlab.system.ui.ConstructorBuilder('phased.ULA');
                        if dialogExpressionBuilder.isParameter('SensorElement')
                            arrayBuilder.addObjectParameterValue('Element',dialogExpressionBuilder.getParameterValue('SensorElement'));
                        end
                        parameterExpression = arrayBuilder.buildExpression;
                    end
                case 'Partitioned array'
                    if dialogExpressionBuilder.isParameter('SensorPartitionedArray')
                        arrayBuilder = dialogExpressionBuilder.getParameterBuilder('SensorPartitionedArray');
                        elementBuilder = arrayBuilder.getParameterBuilder('Array');
                        elementBuilder.addObjectParameterValue('Element', dialogExpressionBuilder.getParameterValue('SensorElement'));
                        parameterExpression = arrayBuilder.buildExpression;
                    else
                        paramExpressionBuilder = matlab.system.ui.ConstructorBuilder('phased.PartitionedArray');
                        if dialogExpressionBuilder.isParameter('SensorElement')
                            arrayBuilder = matlab.system.ui.ConstructorBuilder('phased.ULA');
                            arrayBuilder.addObjectParameterValue('Element',dialogExpressionBuilder.getParameterValue('SensorElement'));
                            % When used in partitioned array, default NumElements for phased.ULA is 4 (instead of 2)
                            arrayBuilder.addLiteralParameterValue('NumElements','4');
                            paramExpressionBuilder.addObjectParameterValue('Array',arrayBuilder);
                        end
                        parameterExpression = paramExpressionBuilder.buildExpression;
                    end
                case 'Replicated subarray'
                    if dialogExpressionBuilder.isParameter('SensorReplicatedArray')
                        arrayBuilder = dialogExpressionBuilder.getParameterBuilder('SensorReplicatedArray');
                        arrayBuilder.addObjectParameterValue('Subarray', dialogExpressionBuilder.getParameterBuilder('SensorReplicatedSubArray'));
                        elementBuilder = arrayBuilder.getParameterBuilder('Subarray');
                        elementBuilder.addObjectParameterValue('Element', dialogExpressionBuilder.getParameterValue('SensorElement'));
                        parameterExpression = arrayBuilder.buildExpression;
                    else
                        paramExpressionBuilder = matlab.system.ui.ConstructorBuilder('phased.ReplicatedSubarray');
                        if dialogExpressionBuilder.isParameter('SensorElement')
                            subarrayBuilder = matlab.system.ui.ConstructorBuilder('phased.ULA');
                            subarrayBuilder.addObjectParameterValue('Element',dialogExpressionBuilder.getParameterValue('SensorElement'));
                            paramExpressionBuilder.addObjectParameterValue('Subarray',subarrayBuilder);
                        end
                        parameterExpression = paramExpressionBuilder.buildExpression;
                    end
                case 'Single element'
                    if dialogExpressionBuilder.isParameter('SensorElement')
                        arrayBuilder = dialogExpressionBuilder.getParameterBuilder('SensorElement');
                        parameterExpression = arrayBuilder.buildExpression;
                    else
                        parameterExpression = 'phased.IsotropicAntennaElement';
                    end
                case 'MATLAB expression'
                    if dialogExpressionBuilder.isParameter('SensorExpression')
                        parameterExpression = dialogExpressionBuilder.getLiteralParameterValue('SensorExpression');
                    else
                        parameterExpression = '<matlab expression>';
                    end
            end
        end
    end
end
