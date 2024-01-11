classdef (Sealed,Hidden,StrictDefaults) Directivity < phased.internal.AbstractDirectivity
%This class is for internal use only. It may be removed in the future.

%Directivity   Sensor directivity
%   H = phased.internal.Directivity creates a directivity System object, H.
%   This object calculates the directivity for the specified directions. By
%   default a 2-element uniform linear array (ULA) is used.
%
%   H = phased.internal.Directivity(Name,Value) creates a directivity
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   D = step(H,FREQ,ANG) returns the directivity (in dBi), D, of the sensor
%   for the directions specified in ANG (in degrees), at the given
%   operating frequency FREQ (in Hz). FREQ is a row vector of length L and
%   ANG can be either a row vector of length M or a 2xM matrix. If you set
%   the EnablePolarization property of H to false, D is an MxL matrix whose
%   elements contain the directivity of the sensor array at angles
%   specified in ANG and frequencies specified in FREQ. If you set the
%   EnablePolarization property of H to true, D is a struct with H and V
%   fields. Both H and V fields are MxL matrices whose elements contain the
%   directivity of the sensor array in horizontal and vertical
%   polarization, respectively, at angles specified in ANG and frequencies
%   specified in FREQ.
%
%   When ANG is a 2xM matrix, each column of the matrix specifies the
%   direction in the space in the [azimuth; elevation] form. The azimuth
%   angle should be between [-180 180] degrees and the elevation angle
%   should be between [-90 90] degrees. If ANG is a length M row vector,
%   each element specifies a direction's azimuth angle and the
%   corresponding elevation angle is assumed to be 0.
%
%   D = step(H,FREQ,ANG,W) uses W as the weights applied on the sensor
%   array when you set the WeightsInputPort property to true. W can be a
%   length-N column vector or an NxL matrix where N is the number of
%   subarrays if SensorArray contains subarrays, or the number of elements
%   otherwise. L is the number of frequencies specified in FREQ. If W is a
%   vector, the weights are applied at all frequencies. If W is a matrix,
%   each column of W represents the weights used at the corresponding
%   frequency specified in FREQ.
%
%   D = step(H,FREQ,ANG,STEER) uses STEER as the subarray steering angle
%   (in degrees). STEER can be a scalar or a length-2 column vector. If
%   STEER is a vector, it is in the form of [AzimuthAngle; ElevationAngle].
%   If STEER is a scalar, it represents the azimuth angle and the elevation
%   angle is assumed to be 0. This syntax is only applicable when you use
%   subarrays in the SensorArray property and set the SubarraySteering
%   property in the SensorArray to either 'Phase' or 'Time'.
%
%   D = step(H,FREQ,ANG,WS) ses WS as the weights applied to each element
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
%   You can combine optional input arguments when their enabling properties
%   are set. Optional inputs must be listed in the same order as the order
%   of the enabling properties. For example,
%
%   D = step(H,FREQ,ANG,W,STEER)
%
%   or
%
%   D = step(H,FREQ,ANG,W,WS)
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(H, x) and y = H(x) are
%   equivalent.
%
%   Directivity methods:
%
%   step     - Calculate the directivity of the sensor (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create a directivity object with same property values
%   isLocked - Locked status (logical)
%
%   Directivity properties:
%
%   Sensor             - Sensor 
%   PropagationSpeed   - Signal propagation speed 
%   WeightsInputPort   - Add input to specify the weights
%   EnablePolarization - Enable polarization simulation
%
%   % Examples:
%
%   % Example 1:
%   %   Calculate and plot the directivity for a 4-element ULA between -90 
%   %   and 90 degrees in azimuth.  
%   %   Assume the array's operating frequency is 300 MHz.
%
%   ha = phased.ULA(4);
%   hag = phased.internal.Directivity('Sensor',ha);
%   fc = 300e6; ang = -90:90;
%   g = step(hag,fc,ang);
%   plot(ang,g); title('Array Directivity'); 
%   xlabel('Angle (degrees)'); ylabel('Array Directivity (dB)');
%
%   % Example 2:
%   %   Calculate the directivity for the above array at the direction of
%   %   30 degrees azimuth and 20 degrees elevation. Assuming a Hamming 
%   %   taper is used.
%
%   ha = phased.ULA(4);
%   hag = phased.internal.Directivity('Sensor',ha,...
%               'WeightsInputPort',true);
%   fc = 300e6; ang = [30;20]; w = hamming(4);
%   g = step(hag,fc,ang,w)
%
%   See also phased, phased.ElementDelay, phased.SteeringVector,
%   phased.ArrayResponse, phased.ArrayGain.

%   Copyright 2013-2017 The MathWorks, Inc.

%   Reference
%   [1] Constantine A. Balanis, Antenna Theory, 3rd ed. Wiley-Interscience,
%       2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    methods

        function obj = Directivity(varargin)
            obj@phased.internal.AbstractDirectivity(varargin{:});
        end
    end
    
    methods (Access = protected)
        
        function g = stepImpl(obj,freq,ang, varargin)
            if ~obj.pUseArray && ~isElementFromAntenna(obj.Sensor) && ...
                    isDirectivityKnown(obj.Sensor)
                g = directivity(obj.Sensor,freq,ang);
            else
                g = stepImpl@phased.internal.AbstractDirectivity(obj,freq,ang, varargin{:});
            end
        end
        
        function resp = extractResponse(obj,respIn)
            if isPolarizationCapable(obj.Sensor)
                if obj.EnablePolarization
                    resp.H = abs(respIn.H).^2;
                    resp.V = abs(respIn.V).^2;
                else
                    if isstruct(respIn)
                        resp = hypot(respIn.H,respIn.V).^2;
                    else
                        resp = abs(respIn).^2;
                    end
                end
            else
                resp = abs(respIn).^2;
            end
        end
        
        function g = extractOutput(obj,resp,intresp)
            if obj.EnablePolarization
                g.H = phased.internal.normalizeIntegratedPower(resp.H,intresp.H,true);
                g.H = pow2db(g.H); % return dB

                g.V = phased.internal.normalizeIntegratedPower(resp.V,intresp.V,true);
                g.V = pow2db(g.V); % return dB
            else
                g = phased.internal.normalizeIntegratedPower(resp,intresp,true);
                g = pow2db(g); % return dB
            end
        end
        
    end
    
end

