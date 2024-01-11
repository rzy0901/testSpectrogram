classdef (Sealed) BackscatterHumanTarget < phased.internal.AbstractSampleRateEngine & ...
        matlab.system.mixin.Propagates & matlab.system.mixin.CustomIcon
%This class is for internal use only. It may be removed in the future.

%BackscatterHumanTarget Backscatter radar target for human
%   H = phased.BackscatterHumanTarget creates a backscatter human target
%   System object, H, that computes the reflected signal from a human. A
%   backscatter target is designed to be used in a monostatic setting where
%   the incident and reflect angles of the signal at the target are the
%   same.
%
%   H = phased.BackscatterHumanTarget(Name,Value) creates a radar target
%   object, H, with the specified property Name set to the specified Value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax:
%
%   Y = step(H,X,ANG) returns the reflected signal Y due to the incident
%   signal X. 
%
%   The human model consists of 16 body parts: left and right feet, left
%   and right lower legs, left and right upper legs, left and right hip,
%   left and right lower arms, left and right upper arms, left and right
%   shoulders, neck, and head. Each body part is represented by a cylinder.
%   Therefore, X is a 16-column matrix whose columns are incident signals
%   to each body part.
% 
%   ANG is a 2x16 matrix representing the signal's incident direction to
%   each body part. Each column of ANG specifies the incident direction of
%   the corresponding signal in the form of an [AzimuthAngle;
%   ElevationAngle] pair (in degrees).
%
%   Y is a column vector containing the combined reflected signal from all
%   body parts. The number of rows in Y is the same as the number of rows
%   in X.
%
%   BackscatterHumanTarget methods:
%
%   step     - Reflect the incoming signal (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create backscatter human target object with same property 
%              values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of the radar target object
%   
%   BackscatterHumanTarget properties:
%
%   Height              - Pedestrian height 
%   PropagationSpeed    - Propagation speed
%   OperatingFrequency  - Operating frequency 
%
%   % Example
%   %   Compute the reflected signal from a human, assume the signal
%   %   incidents from far field so the directions are the same for all
%   %   body parts. Assume the operating frequency is 24 GHz.
%
%   humantarget = phased.BackscatterHumanTarget('OperatingFrequency',24e9);
%   x = ones(10,16);
%   ang = repmat([30;0],1,16);
%   y = humantarget(x,ang);
%
%   See also phased, phased.BackscatterRadarTarget.
    
%   Copyright 2018-2019 The MathWorks, Inc.

%   References
%   [1] Joaquim Fortuny-Guasch and Jean-Marc Chareau, Radar Cross Section
%   Measurements of Pedestrian Dummies and Humans in the 24/77 GHz
%   Frequency Bands, 2013
%
%   [2] Naoyuiki Yamada, Radar Cross Section for Pedestrian in 76 GHz Band,
%   R&D Review of Toyota CRDL, Vol 39, No. 4, 2004

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Height     Pedestrian height (m)
        %   Specify the height (in meter) of the pedestrian as a positive
        %   scalar. The default value is 1.65.
        Height = 1.65
        %PropagationSpeed Propagation speed (m/s)
        %   Specify the wave propagation speed (in m/s) in free space as a
        %   scalar. The default value of this property is the speed of
        %   light.
        PropagationSpeed = physconst('lightspeed')
        %OperatingFrequency     Operating frequency (Hz)
        %   Specify the operating frequency (in Hz) of the radar target as
        %   a scalar. The default value of this property is 3e8, i.e., 300
        %   MHz.
        OperatingFrequency = 300e6
    end
    
    properties (Access = private)
        pNumBodyParts
    end
    
    properties (Access = private, Nontunable)
        cBackscatterTarget
    end
    
    methods
        function obj = BackscatterHumanTarget(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj,x)
            setupImpl@phased.internal.AbstractSampleRateEngine(obj);
            [~,~,bodypartssz,bodypartsax] = ...
                phased.internal.getHumanBodyParts(obj.Height);
            [rcsmat,az,el] = coder.const(@feval,...
                'phased.internal.BackscatterHumanTarget.assembleBodyPartsRCS',bodypartssz,bodypartsax,...
                obj.PropagationSpeed,obj.OperatingFrequency);
            
            obj.pNumBodyParts = size(rcsmat,3);
            obj.cBackscatterTarget = phased.BackscatterRadarTarget(...
                'AzimuthAngles',az,'ElevationAngles',el,'RCSPattern',rcsmat,...
                'PropagationSpeed',obj.PropagationSpeed,'OperatingFrequency',...
                obj.OperatingFrequency);
            
            obj.pNumInputChannels = getNumChannels(obj,x);
            obj.pValidatedNumInputChannels = getNumChannels(obj,x);
            
            if isInputDataSizePropagated(obj)
                sz = propagatedInputSize(obj,1);
                obj.cBackscatterTarget.pOutputSizeBound = sz(1);
                obj.cBackscatterTarget.pNumInputColumns = obj.pValidatedNumInputChannels;
            end
        end
        
        function y = stepImpl(obj,x,ang)
            % validate angle
            cond = any(ang(1,:)<-180) || any(ang(1,:)>180) || ...
                    any(ang(2,:)<-90) || any(ang(2,:)>90);
            if cond
                coder.internal.errorIf(cond,'phased:phased:invalidAzElAngle','ANG',...
                    '-180','180','-90','90');
            end
            
            y = sum(step(obj.cBackscatterTarget,x,ang),2);
        end
        
        function validateInputsImpl(obj,x,ang)
            xsize = size(x);
            cond =  ~isa(x,'double');
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:invalidInputDataType','X','double');
            end
            cond =  ~ismatrix(x) || isempty(x);
            if cond
                coder.internal.errorIf(cond, ...
                     'MATLAB:system:inputMustBeMatrix','X');
            end
            ncol = 16;
            cond = xsize(2)~=ncol;
            if cond
                coder.internal.errorIf(cond,'phased:phased:invalidColumnNumbers','X',sprintf('%d',ncol));
            end
            
            validateNumChannels(obj,x);
            
            angsize = size(ang);
            cond = xsize(2) ~= angsize(2);
            if cond
                coder.internal.errorIf(cond,'phased:phased:collector:AngleDimensionMismatch');
            end
            sigdatatypes.validateAzElAngle(ang,'','Ang');
        end

        function flag = isInputSizeMutableImpl(obj,index) %#ok<INUSL>
            % Return false if input size cannot change
            % between calls to the System object
            if index == 1
                flag = true;
            else
                flag = false;
            end
        end

        function flag = isInputComplexityMutableImpl(obj,index) %#ok<INUSL>
            % Return false if input complexity cannot change
            % between calls to the System object
            if index == 1
                flag = true;
            else
                flag = false;
            end
        end

        function flag = isInputDataTypeMutableImpl(obj,index) %#ok<INUSD>
            % Return false if input data type cannot change
            % between calls to the System object
            flag = false;
        end

        function num = getNumInputsImpl(obj) %#ok<MANU>
            % Define total number of inputs for system with optional inputs
            num = 2;
        end

        function s = loadSubObjects(obj,s,wasLocked)
            if wasLocked
                obj.cBackscatterTarget = phased.BackscatterRadarTarget.loadobj(s.cBackscatterTarget);
                s = rmfield(s,'cBackscatterTarget');
            end
        end
    
        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s
            s = loadSubObjects(obj,s,wasLocked);
            fn = fieldnames(s);
            for m = 1:numel(fn)
                obj.(fn{m}) = s.(fn{m});
            end
        end

        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@phased.internal.AbstractSampleRateEngine(obj);
            
            if isLocked(obj)
                s.pNumBodyParts = obj.pNumBodyParts;
                s.cBackscatterTarget = saveobj(obj.cBackscatterTarget);
            end

        end

        function icon = getIconImpl(obj) %#ok<MANU>
            % Define icon for System block
            icon = "Backscatter Human"; % Example: text icon
            % icon = ["My","System"]; % Example: multi-line text icon
            % icon = matlab.system.display.Icon("myicon.jpg"); % Example: image file icon
        end

        function [name,name2] = getInputNamesImpl(obj) %#ok<MANU>
            % Return input port names for System block
            name = 'X';
            name2 = 'ang';
        end

        function name = getOutputNamesImpl(obj) %#ok<MANU>
            % Return output port names for System block
            name = 'Y';
        end

        function out = getOutputSizeImpl(obj)
            % Return size for each output port
            out = propagatedInputSize(obj,1);
            out(2) = 1;

            % Example: inherit size from first input port
            % out = propagatedInputSize(obj,1);
        end

        function out = getOutputDataTypeImpl(obj) %#ok<MANU>
            % Return data type for each output port
            out = "double";

            % Example: inherit data type from first input port
            % out = propagatedInputDataType(obj,1);
        end

        function out = isOutputComplexImpl(obj)
            % Return true for each output port with complex data
            out = propagatedInputComplexity(obj,1);

            % Example: inherit complexity from first input port
            % out = propagatedInputComplexity(obj,1);
        end

        function out = isOutputFixedSizeImpl(obj)
            % Return true for each output port with fixed size
            out = propagatedInputFixedSize(obj,1);

            % Example: inherit fixed-size status from first input port
            % out = propagatedInputFixedSize(obj,1);
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractSampleRateEngine(obj);
            reset(obj.cBackscatterTarget);
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractSampleRateEngine(obj);
            release(obj.cBackscatterTarget);
        end
    end
    
    methods (Static,Hidden)
        function [rcsmat,az,el] = assembleBodyPartsRCS(bodypartssz,bodypartsax,c,fc)
        % Use ellipsoid to model body segments and set up orientation of
        % these body segments.

            az = -180:180;
            el = -90:90;

            % foot
            foot_a = bodypartssz(1,1);
            foot_c = bodypartssz(2,1);
            footax = bodypartsax(:,:,1);
            footrcs = rotpat(phased.internal.ellipsoidrcs(foot_a,foot_a,foot_c,c,fc,az,el),az,el,footax);

            % lower leg
            lowerleg_a = bodypartssz(1,2);
            lowerleg_c = bodypartssz(2,2);
            lowerlegrcs = phased.internal.ellipsoidrcs(lowerleg_a,lowerleg_a,lowerleg_c,c,fc,az,el);

            % upper leg
            upperleg_a = bodypartssz(1,3);
            upperleg_c = bodypartssz(2,3);
            upperlegrcs = phased.internal.ellipsoidrcs(upperleg_a,upperleg_a,upperleg_c,c,fc,az,el);

            % left and right hip
            hip_a = bodypartssz(1,4);
            hip_c = bodypartssz(2,4);
            hipax = bodypartsax(:,:,4);
            hiprcs = rotpat(phased.internal.ellipsoidrcs(hip_a,hip_a,hip_c,c,fc,az,el),az,el,hipax);

            % lower arm
            lowerarm_a = bodypartssz(1,5);
            lowerarm_c = bodypartssz(2,5);
            lowerarmrcs = phased.internal.ellipsoidrcs(lowerarm_a,lowerarm_a,lowerarm_c,c,fc,az,el);

            % upper arm
            upperarm_a = bodypartssz(1,6);
            upperarm_c = bodypartssz(2,6);
            upperarmrcs = phased.internal.ellipsoidrcs(upperarm_a,upperarm_a,upperarm_c,c,fc,az,el);

            % left and right shoulder
            shoulder_a = bodypartssz(1,7);
            shoulder_c = bodypartssz(2,7);
            shoulderax = bodypartsax(:,:,7);
            shoulderrcs = rotpat(phased.internal.ellipsoidrcs(shoulder_a,shoulder_a,shoulder_c,c,fc,az,el),az,el,shoulderax);

            % head
            head_a = bodypartssz(1,8);
            head_c = bodypartssz(2,8);
            headrcs = phased.internal.ellipsoidrcs(head_a,head_a,head_c,c,fc,az,el);

            % torso
            torso_a = bodypartssz(1,9);
            torso_c = bodypartssz(2,9);
            torsorcs = phased.internal.ellipsoidrcs(torso_a,torso_a,torso_c,c,fc,az,el);

            % Build rcs matrix (performed this way for codegen; addresses g1932019) 
            rcsmat = zeros(181,361,16);
            rcsmat(:,:,1)  = footrcs;
            rcsmat(:,:,2)  = footrcs;
            rcsmat(:,:,3)  = lowerlegrcs;
            rcsmat(:,:,4)  = lowerlegrcs;
            rcsmat(:,:,5)  = upperlegrcs;
            rcsmat(:,:,6)  = upperlegrcs;
            rcsmat(:,:,7)  = hiprcs;
            rcsmat(:,:,8)  = hiprcs;
            rcsmat(:,:,9)  = lowerarmrcs;
            rcsmat(:,:,10) = lowerarmrcs;
            rcsmat(:,:,11) = upperarmrcs;
            rcsmat(:,:,12) = upperarmrcs;
            rcsmat(:,:,13) = shoulderrcs;
            rcsmat(:,:,14) = shoulderrcs;
            rcsmat(:,:,15) = headrcs;
            rcsmat(:,:,16) = torsorcs;

            % normalize, according to [1] and [2], the average RCS in azimuth is about
            % -11.5 dBsm at 24 GHz range and -8.1 dBsm at 77 GHz range.
            rcs_az_avg_db = interp1([24 77],[-11.5 -8.1],fc/1e9,'linear','extrap');
            rcsmat_sum = (sum(sqrt(rcsmat(91,:,:)),3)).^2;
            k = db2pow(rcs_az_avg_db)/mean(rcsmat_sum);
            rcsmat = k*rcsmat;
            
        end
    end
end

