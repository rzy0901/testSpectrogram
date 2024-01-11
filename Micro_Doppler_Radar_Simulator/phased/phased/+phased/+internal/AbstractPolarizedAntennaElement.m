classdef(Hidden) AbstractPolarizedAntennaElement < phased.internal.AbstractAntennaElement
%This class is for internal use only. It may be removed in the future.

%   Copyright 2012 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    methods

        function obj = AbstractPolarizedAntennaElement(varargin)
            %AbstractPolarizedAntennaElement   Construct the
            %AbstractPolarizedAntennaElement class.
            obj@phased.internal.AbstractAntennaElement(varargin{:});
        end
        
    end
    
    methods (Abstract, Access = protected)
        resp = getHResponse(obj,freq,angle);
        %getSpatialResponse Returns the spatial response of the element at
        %a given direction for a given freq that is inside the element's
        %frequency range.
        resp = getVResponse(obj,freq,angle);
        %getFrequencyResponse Returns the frequency response of the element
        %at a given frequency.
    end
    
    methods (Access = protected)

        function num = getNumOutputsImpl(obj)
            num = getNumOutputsImpl@phased.internal.AbstractAntennaElement(obj);
        end
        
        function num = getNumInputsImpl(obj)
            num = getNumInputsImpl@phased.internal.AbstractAntennaElement(obj);
        end
        
        function resp = stepImpl(obj,freq,angArg)
            [ang, num_ang] = validateInputValues(obj,freq,angArg);
            
            if isPolarizationEnabled(obj) 
                resp.H = getOutputResponse(obj,freq,ang,'horizontal',num_ang);
                resp.V = getOutputResponse(obj,freq,ang,'vertical',num_ang);
            else
                resp = getOutputResponse(obj,freq,ang,'combined',num_ang);
            end
                
        end
        
        function resp = getSpatialResponse(obj,freq,ang)
            resp = hypot(...
                getHResponse(obj,freq,ang),...
                getVResponse(obj,freq,ang));
        end
        
        function resp = getOutputSpatialResponse(obj,freq,ang,option)
            if strcmp(option,'horizontal')   % horizontal
                fh = @getHResponse;
            elseif strcmp(option,'vertical')   % vertical
                fh = @getVResponse;
            else 
                fh = @getSpatialResponse;
            end
            resp = fh(obj,freq,ang);
        end
        
    end
    
    methods 
        function flag = isPolarizationCapable(obj)   %#ok<MANU>
            flag = true;
        end
    end
    
    
   
end

% [EOF]
