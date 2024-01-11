classdef (Hidden) AbstractHomogeneousArray < phased.internal.AbstractArray
%This class is for internal use only. It may be removed in the future.

%AbstractHomogeneousArray   Define the AbstractHomogeneousArray class.

%   Copyright 2012-2014 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Element  Element of the array
        %   Specify the element used in the sensor array as a scalar
        %   element objects in the phased package. The default value of
        %   this property is an isotropic antenna element.
        Element
    end
    
    methods
        function set.Element(obj,val)
            if ~isa(val,'phased.internal.AbstractElement') && ~isa(val,'em.Antenna')
                error(message('phased:system:array:InvalidElement','Element'));
            end
            obj.Element = val;
        end
        
    end
    
    methods (Access = protected)

        function obj = AbstractHomogeneousArray(varargin)
            obj@phased.internal.AbstractArray(varargin{:});
            if isempty(coder.target)
                if isempty(obj.Element)
                    obj.Element = phased.IsotropicAntennaElement;
                end
            else
                if ~coder.internal.is_defined(obj.Element)
                    obj.Element = phased.IsotropicAntennaElement;
                end
            end

        end

    end
    
    methods (Access = protected)
        
        function setupImpl(obj,freq,ang)
            sz_freq = size(freq);
            obj.pNumFreqs = sz_freq(2);
            sz_angle = size(ang);
            obj.pNumAngles = sz_angle(2);
            num_angles = sz_angle(2);
            if isempty(coder.target)
                if ~isElementFromAntenna(obj.Element)
                    obj.cElement = clone(obj.Element);
                else
                    obj.cElement = phased.internal.AntennaAdapter('Antenna',clone(obj.Element));
                end
                release(obj.cElement);
            else
                if isElementFromAntenna(obj.Element)
                    coder.internal.errorIf(true, ...
                        'phased:system:element:AntennaToolboxCodegenNotSupported','em.Antenna','phased.CustomAntennaElement');
                end
                obj.cElement = clonecg(obj.Element);
            end

            num_freqs =  sz_freq(2);
            obj.pNumFreqs = num_freqs;
            obj.pNumAngles = num_angles;
            obj.pAzimuthOnly = (sz_angle(1) == 1);
            
            lTaper = getTaper(obj);
            if isscalar(lTaper)
                obj.pTaper = lTaper;
            else
                obj.pTaper = ...
                    repmat(lTaper,num_angles,num_freqs);
            end
   
        end
        
        function resetImpl(obj)
            %Reset AbstractHomogeneous
            reset(obj.cElement);         
        end
        
        function releaseImpl(obj)
            %Release AbstractHomogeneous
            release(obj.cElement);   
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractArray(obj);
            s.Element = saveobj(obj.Element);
            if isLocked(obj)
                s.cElement = saveobj(obj.cElement);
            end
        end
        
        function s = loadSubObjects(obj,s)
            if isa(s.Element,'em.Antenna')
                obj.Element = em.EmStructures.loadobj(s.Element);
            else
                obj.Element = phased.internal.AbstractElement.loadobj(s.Element);
            end
            s = rmfield(s,'Element');
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cElement = phased.internal.AbstractElement.loadobj(s.cElement);
                    s = rmfield(s,'cElement');
                end
                s = rmfield(s,'isLocked');
            end
        end
    end
    
    methods (Access = protected)
        function respOut = getPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles)
            num_elements = getNumElements(obj);
            inc_angle = convertIncidentToAzEl(obj,ang);
            resp = step(obj.cElement,freq,inc_angle);
            
            resp.H = resp.H.*obj.pTaper;
            respOut.H = reshape(resp.H,num_elements,num_angles,[]);
            resp.V = resp.V.*obj.pTaper;
            respOut.V = reshape(resp.V,num_elements,num_angles,[]);
        end
        
        function resp = getNonPolarizationOutputInArrayAzEl(obj,freq,ang,num_angles)
            num_elements = getNumElements(obj);
            inc_angle = convertIncidentToAzEl(obj,ang);
            
            resp = step(obj.cElement,freq,inc_angle);
            resp = resp.*obj.pTaper;
            resp = reshape(resp,num_elements,num_angles,[]);
        end

    end
    
    methods
        function flag = isPolarizationCapable(obj)
        %isPolarizationCapable Indicate if the array is capable of 
        %simulating polarization
        %   F = isPolarizationCapable(H) returns the flag, F, which
        %   indicates whether the array H is capable of simulating
        %   polarization.
        %
        %   % Example:
        %   %   Determine whether an array of IsotropicAntennaElement 
        %   %   is capable of simulating polarization.
        %   
        %   h = phased.ULA;
        %   f = isPolarizationCapable(h)
        %
        %   See also phased.

            if isElementFromAntenna(obj.Element)
                flag = true;  % always true
            else
                flag = isPolarizationCapable(obj.Element);
            end
        end
        
    end
    
    methods (Hidden)
        function flag = isPolarizationEnabled(obj)
        %isPolarizationEnabled Indicate if the array is enabled to 
        %simulate polarization
        %   F = isPolarizationEnabled(H) returns the flag, F, which
        %   indicates whether the array H is enabled to simulate
        %   polarization.
        %
        %   % Example:
        %   %   Determine whether an array of IsotropicAntennaElement 
        %   %   is enabled to simulate polarization.
        %   
        %   h = phased.ULA;
        %   f = isPolarizationEnabled(h)
        %
        %   See also phased.

            if isElementFromAntenna(obj.Element)
                % only accessible via step
                if isLocked(obj)
                    flag = isPolarizationEnabled(obj.cElement);
                else
                    flag = true;
                end
            else
                flag = isPolarizationEnabled(obj.Element);
            end
        end
        
    end
    
    methods(Static, Hidden, Access=protected)
        function groups = getPropertyGroupsImpl
            classSet = matlab.system.display.internal.ClassStringSet(...
                {
                    'phased.IsotropicAntennaElement',...
                    'phased.CosineAntennaElement',...
                    'phased.ShortDipoleAntennaElement',...
                    'phased.CrossedDipoleAntennaElement',...
                    'phased.CustomAntennaElement',...
                    'phased.OmnidirectionalMicrophoneElement',...
                    'phased.CustomMicrophoneElement' ...
                }, ...
                'PropertiesTitle', '', ...
                'Labels', { ...
                    'Isotropic Antenna', ...
                    'Cosine Antenna', ...
                    'Short Dipole Antenna', ...
                    'Crossed Dipole Antenna', ...
                    'Custom Antenna', ...
                    'Omni Microphone', ...
                    'Custom Microphone' ...
                          });        
            elementProp = matlab.system.display.internal.Property('Element', ...
                'IsObjectDisplayOnly', true, ...
                'Description', 'Sensor element', ...
                'ClassStringSet', classSet);
            
            props = {elementProp};
            groups = matlab.system.display.Section('Title', 'Parameters', ...
                'PropertyList', props);
        end
    end
    
    methods (Hidden)
        function h = getElementHandle(obj)
            h = obj.Element;
        end
        function newObj = cloneSensor(obj)
            if isElementFromAntenna(obj)
                newObj = clone(obj); 
                release(newObj); 
                newObj.Element = phased.internal.AntennaAdapter(...
                    'Antenna',obj.Element);
            else
                newObj = clone(obj);
                release(newObj);
            end
        end
        function flag = isElementFromAntenna(obj) 
            flag = isElementFromAntenna(obj.Element);
        end
    end
    
    methods (Hidden, Access = {?phased.internal.AbstractArray, ?phased.gpu.internal.AbstractClutterSimulator})
        %Methods used by the GPU ConstantGammaClutter model.
        function xscaled = scaleByElemResponse(obj, azin, elin, freq, x, idx)
            %scale x - the steering vector, by the element pattern for the
            %elements selected by the indices in idx
   
            %For homogeneous elements scale by the element pattern.
            if isElementFromAntenna(obj.Element)
                hele = cloneSensor(obj.Element);
            else
                hele = obj.Element;
            end
            elPat = getgpuElemResponse(hele, azin, elin, freq);
            if isscalar(obj.Taper),
                elPatwTaper = elPat.*obj.Taper;
            else
                elPatwTaper = bsxfun(@times, elPat, reshape(obj.Taper(idx), 1, [], 1));
            end
            xscaled = bsxfun(@times, elPatwTaper, x);
        end
    end
    
end

% [EOF]
