classdef (Hidden) AbstractRadarTarget < phased.internal.AbstractTarget
%This class is for internal use only. It may be removed in the future.

%ABSTRACTRADARTARGET Define the ABSTRACTRADARTARGET class
% This is an abstract class in support of radar target functionality.

%   Copyright 2015-2017 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [2] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed., 
%       McGraw-Hill, 2001
%   [3] Harold Mott, Antenna for Radar and Communications: A polarimetric
%       Approach, Wiley-Interscience, 1992
%   [4] Eugene Knott, et al. Radar Cross Section, Scitech Publishing, 2004


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Logical, Nontunable)
        %EnablePolarization  Enable polarization
        %   Set this property to true to enable polarization. Set this
        %   property to false to ignore polarization. The default value of
        %   this property is false. This property applies when the sensor
        %   specified in the Sensor property is capable of simulating
        %   polarization.
        EnablePolarization = false
    end
    
    properties (Nontunable)
        %PropagationSpeed   Propagation speed (m/s)
        %   Specify the propagation speed (in m/s) of the signal as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('LightSpeed')
        %OperatingFrequency     Operating frequency (Hz)
        %   Specify the operating frequency (in Hz) of the radar target as
        %   a scalar. The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 3e8
    end
    
    properties(Access = protected)
        % Private flag for whether pRCS needs to be initialized
        pNeedRCSInit;
    end
    
    properties(Access = protected, Nontunable)        
        pLambda; % Private wavelength
    end
     
    methods
        function set.OperatingFrequency(obj,value)
            sigdatatypes.validateFrequency(value,'phased.RadarTarget',...
                'OperatingFrequency',{'scalar'});
            obj.OperatingFrequency = value;
        end
        function set.PropagationSpeed(obj,value)
            sigdatatypes.validateSpeed(value,...
                'phased.RadarTarget','PropagationSpeed',...
                {'scalar','positive'});
            obj.PropagationSpeed = value;
        end
    end

    methods (Access = protected)
        function obj = AbstractRadarTarget(varargin)
            obj = obj@phased.internal.AbstractTarget(varargin{:});
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractTarget(obj);
            if obj.pFluctuate
                obj.pNeedRCSInit = true;
            else
                obj.pNeedRCSInit = false;
            end
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractTarget(obj);
            if isLocked(obj)
                s.pNeedRCSInit = obj.pNeedRCSInit;
            end
            s.pLambda = obj.pLambda;
        end
        
        function s = loadSubObjects(obj,s,wasLocked)
            if isfield(s,'isLocked')  % for backwards compatibility
                if s.isLocked
                    obj.pFluctuate = s.pFluctuate;
                    if s.pFluctuate
                        s = rmfield(s,'pRandState');
                        obj.cNoiseSource = phased.internal.NoiseSource(...
                            'SeedSource','Property','Seed',s.pSeed);
                        obj.pFluctuateFunc = s.pFluctuateFunc;
                        s = rmfield(s,'pFluctuateFunc');
                    end
                    s = rmfield(s,'pSeed');
                    s = rmfield(s,'pFluctuate');
                end
                s = rmfield(s,'isLocked');
            else
                s = loadSubObjects@phased.internal.AbstractTarget(obj,s,wasLocked);
            end
        end
    end
end
