classdef (Sealed,StrictDefaults) SimulinkSteeringVector < phased.internal.AbstractSteeringVector & ...
        matlab.system.mixin.CustomIcon 
%This class is for internal use only. It may be removed in the future.

%SteeringVector   Sensor array steering vector
%   H = phased.SteeringVector creates a steering vector System object, H.
%   This object calculates the steering vector of a sensor array for the
%   specified directions. By default a 2-element uniform linear array (ULA)
%   is used.
%
%   H = phased.SteeringVector(Name,Value) creates a steering vector object,
%   H, with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   SV = step(H,FREQ,ANG) returns the steering vector, SV, of the array for
%   the directions specified in ANG (in degrees), at the given operating
%   frequency FREQ (in Hz). FREQ is a row vector of length L and ANG can be
%   either a row vector of length M or a 2xM matrix. When ANG is a 2xM
%   matrix, each column of the matrix specifies the direction in space in
%   the [azimuth; elevation] form. The azimuth angle should be between
%   [-180 180] degrees and the elevation angle should be between [-90 90]
%   degrees. If ANG is a length M row vector, each element specifies a
%   direction's azimuth angle and the corresponding elevation angle is
%   assumed to be 0.
%
%   The dimensions of SV are NxMxL where N is the number of subarrays if
%   SensorArray contains subarrays, or the number of elements otherwise.
%   Each column of SV contains the steering vector of the array for the
%   corresponding directions specified in ANG. Each page of SV contains
%   the steering vectors of the array for the given frequency specified in
%   FREQ.
%
%   If you set the IncludeElementResponse property to true, the resulting
%   steering vector SV includes the individual element responses. If you
%   set the IncludeElementResponse property to false, the elements are
%   assumed to be isotropic so the steering vector SV does not include the
%   individual element responses.
%
%   When you set the EnablePolarization property to true, the resulting SV
%   is a structure containing two fields, H and V. H represents the array's
%   response in horizontal polarization and V represents the array's
%   response in vertical polarization. Each field is an NxMxL array whose
%   columns contain the steering vectors of the sensor array for the
%   corresponding directions and frequencies. If you set the
%   EnablePolarization to false, then the polarization information is
%   discarded and the combined pattern from both H and V polarizations is
%   used at each element to compute the steering vector. This syntax is
%   only applicable when the sensor array is capable of simulating
%   polarization and when you set the IncludeElementResponse property to
%   true.
%
%   SV = step(H,FREQ,ANG,STEER) uses STEER as the subarray steering angle
%   (in degrees). STEER can be a scalar or a length-2 column vector. If
%   STEER is a vector, it is in the form of [AzimuthAngle; ElevationAngle].
%   If STEER is a scalar, it represents the azimuth angle and the elevation
%   angle is assumed to be 0. This syntax is only applicable when you use
%   subarrays in the SensorArray property, set the SubarraySteering
%   property in the SensorArray to either 'Phase' or 'Time' and set the
%   IncludeElementResponse property to true.
%
%   SV = step(H,FREQ,ANG,WS) uses WS as the weights applied to each element
%   in the subarray. WS can be either a matrix or a cell array. This syntax
%   is only applicable when you use subarrays in the SensorArray property,
%   set the SubarraySteering property in the SensorArray to 'Custom' and
%   set the IncludeElementResponse property to true.
%   
%   If the Sensor property is a phased.ReplicatedSubarray, WS must be an
%   NSExN matrix where NSE is the number of elements in each individual
%   subarray and N is the number of subarrays. Each column in WS specifies
%   the weights for the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its individual
%   subarrays have same number of elements, WS must be an NSExN matrix
%   where NSE is the number of elements in each individual subarray and N
%   is the number of subarrays. Each column in WS specifies the weights for
%   the elements in the corresponding subarray.
%
%   If the Sensor property is a phased.PartitionedArray and its subarrays
%   can have different number of elements, WS can be either an NSExN
%   matrix, where NSE indicates the number of elements in the largest
%   subarray and N is the number of subarrays, or a 1xN cell array, where N
%   is the number of subarrays and each cell contains a column vector whose
%   length is the same as the number of elements of the corresponding
%   subarray.  If WS is a matrix, the first K entries in each column, where
%   K is the number of elements in the corresponding subarray, specifies
%   the weights for the elements in the corresponding subarray. If WS is a
%   cell array, each cell in the array is a column vector specifying the
%   weights for the elements in the corresponding subarray. 
% 
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   SteeringVector methods:
%
%   step     - Calculate the steering vector (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a steering vector object with same property values
%   isLocked - Locked status (logical)
%
%   SteeringVector properties:
%
%   SensorArray            - Sensor array 
%   PropagationSpeed       - Signal propagation speed 
%   IncludeElementResponse - Include element response in steering vector
%   NumPhaseShifterBits    - Number of bits in phase shifters
%   EnablePolarization     - Enable polarization simulation
%
%   % Example:
%   %   Calculate the steering vector for a 4-element uniform linear array 
%   %   at the direction of 30 degrees azimuth and 20 degrees elevation. 
%   %   Assume the array's operating frequency is 300 MHz. Compare the beam
%   %   pattern before and after the steering.
%
%   array = phased.ULA(4);
%   steervector = phased.SteeringVector('SensorArray',array);
%   sv = steervector(3e8,[30; 20])
%   c = steervector.PropagationSpeed;
%   subplot(211)
%   pattern(array,3e8,-180:180,0,'PropagationSpeed',c,...
%       'CoordinateSystem','rectangular'); 
%   title('Before steering');
%   subplot(212)
%   pattern(array,3e8,-180:180,0,'PropagationSpeed',c,'Weights',sv,...
%       'CoordinateSystem','rectangular'); 
%   title('After steering');
%
%   See also phased, phased.ElementDelay, phased.ArrayResponse,
%   phased.ArrayGain.

%   Copyright 2018-2019 The MathWorks, Inc.

%   Reference
%   [1] Van Trees, Optimum Array Processing, Wiley, 2002


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen    

    properties (Nontunable)
        %OperatingFrequency Operating frequency (Hz)
        %   Specify the operating frequency (in Hz) as a positive scalar.
        %   The default value is 3e8.
        OperatingFrequency = 3e8;    
    end
        
    methods

        function obj = SimulinkSteeringVector(varargin)
            %SteeringVector   Construct the SteeringVector class.
            
            obj@phased.internal.AbstractSteeringVector(varargin{:});

        end
        
    end
    
    methods
        function set.OperatingFrequency(obj,value)
            validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'OperatingFrequency');
            obj.OperatingFrequency = value;
        end    
    end
    
    methods (Access = protected)
        
        function setupImpl(obj,ang,stang) %#ok<INUSD>
            
            setupImpl@phased.internal.AbstractSteeringVector(obj);

            sz_freq = size(obj.OperatingFrequency);
            sz_angle = size (ang);
            obj.pNumInputAngles = sz_angle(2);
            obj.pIsFreqScalar = (sz_freq(2)==1);
        end
        
        function validateInputsImpl(obj,varargin)
            validateAngleInputs(obj,varargin{:});
        end
        
        function loadObjectImpl(obj,s,wasLocked) %#ok<INUSD>
            s = loadSubObjects(obj,s);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end
        
        function sv = stepImpl(obj,varargin)
            freq = obj.OperatingFrequency;
            sv = computeSteeringVector(obj,freq,varargin{:});
        end
        
        function num = getNumInputsImpl(obj) 
            num = getNumInputsImpl@phased.internal.AbstractSteeringVector(obj);
            if isa(obj.SensorArray,'phased.internal.AbstractSubarray') && ...
                    obj.IncludeElementResponse && ...
                    ~strncmp(obj.SensorArray.SubarraySteering,'None',1)
                num = num+2;
            
            else
                num = num+1;
            end
        end
        
    end

    methods (Access = protected)
        function varargout = getOutputSizeImpl(obj)
            szAng = propagatedInputSize(obj,1);
            numDOF = getDOF(obj.SensorArray);
            varargout{1} = [numDOF szAng(2)];
        end
        function varargout = isOutputFixedSizeImpl(obj)
            varargout{1} = propagatedInputFixedSize(obj,1);
        end
    end
    
    methods (Static,Hidden,Access=protected)
        function groups = getPropertyGroupsImpl
            groups = getPropertyGroupsImpl@phased.internal.AbstractSteeringVector;
            dEnablePolarization = ...
                matlab.system.display.internal.Property('EnablePolarization', ...
                'IsGraphical', false);
            props = {...
                'OperatingFrequency',...
                'IncludeElementResponse',...
                'NumPhaseShifterBits',...
                dEnablePolarization};
            groups(1).PropertyList = [groups(1).PropertyList props];
        end
        function header = getHeaderImpl
            header = matlab.system.display.Header(...
                'Title',getString(message('phased:library:block:SteeringVectorTitle')),...
                'Text',getString(message('phased:library:block:SteeringVectorDesc')));
        end
    end
    methods (Access = protected)
        function varargout = getInputNamesImpl(obj)  %#ok<MANU>
            varargout = {'Ang','Steer'};
        end
        
        function varargout = getOutputNamesImpl(obj) %#ok<MANU>
            varargout = {'SV'};
        end
        function str = getIconImpl(obj) %#ok<MANU>
            str = sprintf('Steering\nVector');
        end
        
    end    
end

