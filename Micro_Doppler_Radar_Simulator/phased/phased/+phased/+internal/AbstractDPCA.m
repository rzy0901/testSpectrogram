classdef (Hidden) AbstractDPCA < phased.internal.AbstractSTAP
%This class is for internal use only. It may be removed in the future.

%AbstractDPCA   Define the AbstractDPCA class
%   This is an abstract class in support of DPCA functionality.

%   Copyright 2008-2017 The MathWorks, Inc.


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

    properties (Nontunable)
        %Direction    Receiving mainlobe direction (deg)
        %   Specify the receiving mainlobe direction of the receiving
        %   sensor array as a length 2 column vector. The direction is
        %   specified in the format of [AzimuthAngle; ElevationAngle] (in
        %   degrees). Azimuth angle should be between -180 and 180.
        %   Elevation angle should be between -90 and 90. This property
        %   applies when you set the DirectionSource property to
        %   'Property'. The default value of this property is [0; 0].
        Direction = [0; 0];
    end

    properties (Nontunable, Logical) 
        %PreDopplerOutput    Output pre-Doppler result
        %   Set this property to true to output the processing result
        %   before applying the Doppler filtering. Set this property to
        %   false to output the processing result after the Doppler
        %   filtering. The default value of this property is false.
        PreDopplerOutput = false;
    end

    properties (Access = protected, Nontunable)
        %NumDPCAPulses  Number of pulses used for DPCA based processing
        %   Specify the number of pulses used for DPCA based processing as
        %   an integer.  Because DPCA based processing results in the
        %   reduced aperture size, NumDPCAPulses is normally 2 or 3, and
        %   must be larger than 1.
        NumDPCAPulses = 2;
    end
    
    properties (Access = protected, Nontunable)
        pArrayNumElements
        pReducedDimension
        pNumPulses
        cSteeringVector;
    end
    
    methods(Access = protected, Abstract)
        w = calcPreDopplerWeights(obj,angle,tempdata,CUTIdx)
    end
    
    methods
        function set.Direction(obj,val)
            sigdatatypes.validateAzElAngle(val,'phased.DPCA',...
                'Direction',{'size',[2 1]});
            obj.Direction = val;
        end
        
    end
    
    methods (Access = protected)
        function obj = AbstractDPCA(varargin)
            obj@phased.internal.AbstractSTAP(varargin{:});
        end
    end
    
    methods (Access = protected)
        
        function privValidateSensorArray(obj,val) 
            privValidateSensorArray@phased.internal.AbstractSTAP(obj,val);
            validateattributes(val,{'phased.ULA','phased.HeterogeneousULA'},...
                {'scalar'},'','SensorArray');
        end     
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = false;
            if obj.PreDopplerOutput
                if strcmp(prop,'DopplerSource') || strcmp(prop,'Doppler')
                    flag = true;
                end
            else
                if (obj.DopplerSource(1) == 'I') && strcmp(prop,'Doppler')
                    flag = true;
                end
            end
            if (obj.DirectionSource(1) == 'I') && strcmp(prop,'Direction')
                flag = true;
            end
            if(obj.PRFSource(1) == 'I'&& strcmp(prop, 'PRF'))
                flag = true;
            end
        end
        
        function num = getNumInputsImpl(obj)
            if obj.PreDopplerOutput
                if (obj.DirectionSource(1) == 'P')
                    num = 2;
                else
                    num = 3;
                end
            else
                if (obj.DirectionSource(1) == 'P')
                    if (obj.DopplerSource(1) == 'P')
                        num = 2;
                    else
                        num = 3;
                    end
                else
                    if (obj.DopplerSource(1) == 'P')
                        num = 3;
                    else
                        num = 4;
                    end
                end
            end
            if (obj.PRFSource(1) == 'I')
                num = num+1 ; %5
            end
        end
                
        function validateInputsImpl(obj,x,cutidx,varargin)
            validateInputsImpl@phased.internal.AbstractSTAP(obj,x,cutidx);
            if (obj.PRFSource(1) == 'P')
                prf = [];
                if (obj.DirectionSource(1) == 'I')
                    adir = varargin{1};
                else
                    adir = [];
                end
                if ~obj.PreDopplerOutput
                    if (obj.DirectionSource(1) == 'P')
                        if (obj.DopplerSource(1) == 'I')
                            dop = varargin{1};
                        else
                            dop = [];
                        end
                    else
                        if (obj.DopplerSource(1) == 'I')
                            dop = varargin{2};
                        else
                            dop = [];
                        end
                    end
                else
                    dop = [];
                end
            else
                prf = varargin{1};
                if (obj.DirectionSource(1) == 'I')
                    adir = varargin{2};
                else
                    adir = [];
                end
                if ~obj.PreDopplerOutput
                    if (obj.DirectionSource(1) == 'P')
                        if (obj.DopplerSource(1) == 'I')
                            dop = varargin{2};
                        else
                            dop = [];
                        end
                    else
                        if (obj.DopplerSource(1) == 'I')
                            dop = varargin{3};
                        else
                            dop = [];
                        end
                    end
                else
                    dop = [];
                end
            end
            if ~isempty(prf)
                validateInputPRFSpec(obj,prf);
            end
            if ~isempty(adir)
                validateInputAngleSpec(obj,adir);
            end
            if ~isempty(dop)
                validateInputDopplerSpec(obj,dop);
            end
        end
        
        function setupImpl(obj,x,CUTIdx,varargin) %#ok<INUSD>
            setupImpl@phased.internal.AbstractSTAP(obj,x);
            obj.pArrayNumElements = getNumElements(obj.SensorArray);
            sz_x = size(x);
            obj.pNumPulses = sz_x(3);
            obj.pReducedDimension = reducedNumPulses(obj,obj.pNumPulses);            
            obj.cSteeringVector = phased.SteeringVector(...
                'SensorArray',obj.SensorArray,...
                'PropagationSpeed',obj.PropagationSpeed,...
                'NumPhaseShifterBits',obj.NumPhaseShifterBits);
        end
        
        function flag = isInputComplexityLockedImpl(obj,index)
            flag = isInputComplexityLockedImpl@phased.internal.AbstractSTAP(obj,index);
            if (obj.PRFSource(1) == 'P')
                if (obj.DirectionSource(1) == 'I') && (index == 3)
                    flag = true;
                end
                if ~obj.PreDopplerOutput
                    if (obj.DirectionSource(1) == 'P')
                        if (obj.DopplerSource(1) == 'I') && (index == 3)
                            flag = true;
                        end
                    else
                        if (obj.DopplerSource(1) == 'I') && (index == 4)
                            flag = true;
                        end
                    end
                end
            else
                if index == 3
                    flag = true;
                end
                if (obj.DirectionSource(1) == 'I') && (index == 4)
                    flag = true;
                end
                if ~obj.PreDopplerOutput
                    if (obj.DirectionSource(1) == 'P')
                        if (obj.DopplerSource(1) == 'I') && (index == 4)
                            flag = true;
                        end
                    else
                        if (obj.DopplerSource(1) == 'I') && (index == 5)
                            flag = true;
                        end
                    end
                end
            end
        end
        
        function releaseImpl(obj)
            release(obj.cSteeringVector);
        end
        
        function resetImpl(obj)
            reset(obj.cSteeringVector);
        end
                
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractSTAP(obj);
            s.NumDPCAPulses = obj.NumDPCAPulses;
            if isLocked(obj)
                s.cSteeringVector = saveobj(obj.cSteeringVector);
                s.pArrayNumElements = obj.pArrayNumElements;
                s.pReducedDimension = obj.pReducedDimension;
                s.pNumPulses = obj.pNumPulses;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractSTAP(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cSteeringVector = phased.SteeringVector.loadobj(s.cSteeringVector);
                    s = rmfield(s,'cSteeringVector');
                end
            end
        end
        
        function [yOut, wOut] = stepImpl(obj,x,CUTIdx,varargin)
            
            if (obj.PRFSource(1) == 'P')
                prf = obj.PRF;
                if obj.PreDopplerOutput
                    if (obj.DirectionSource(1) == 'P')
                        angle = obj.Direction;
                    else
                        angle = varargin{1};
                    end
                else
                    if (obj.DirectionSource(1) == 'P')
                        angle = obj.Direction;
                        if (obj.DopplerSource(1) == 'P')
                            dop = obj.Doppler;
                        else
                            dop = varargin{1};
                        end
                    else
                        angle = varargin{1};
                        if (obj.DopplerSource(1) == 'P')
                            dop = obj.Doppler;
                        else
                            dop = varargin{2};
                        end
                    end
                end
            else
                prf = varargin{1};
                if obj.PreDopplerOutput
                    if (obj.DirectionSource(1) == 'P')
                        angle = obj.Direction;
                    else
                        angle = varargin{2};
                    end
                else
                    if (obj.DirectionSource(1) == 'P')
                        angle = obj.Direction;
                        if (obj.DopplerSource(1) == 'P')
                            dop = obj.Doppler;
                        else
                            dop = varargin{2};
                        end
                    else
                        angle = varargin{2};
                        if (obj.DopplerSource(1) == 'P')
                            dop = obj.Doppler;
                        else
                            dop = varargin{3};
                        end
                    end
                end
            end
            
            prf = validateInputPRF(obj,prf);
            angle = validateInputAngle(obj,angle);
            CUTIdx = validateInputCUTIdx(obj,CUTIdx);
            
            
            [y, w] = preDopplerData(obj,x,CUTIdx,angle);
            
            if ~obj.PreDopplerOutput
                
                dop = validateInputDoppler(obj,dop,prf); 
                dpstv = dopsteeringvec(dop,size(y,2),prf);
                ytemp = dpstv'*y.';
                wOut = calcPostDopplerWeights(obj,w,dop,prf);
                yOut = ytemp(:);
            else
                wOut = w;
                yOut = y;
            end
            
        end

        function wdpca = formDPCAWeights(obj,angle)
        %formDPCAWeights Form DPCA weights
        %   W = formDPCAWeights(Hs,ANG) returns the DPCA weights
        %   of the DPCA processor Hs, whose look direction is specified by
        %   ANG. ANG is a 2x1 vector in the form of [AzimuthAngle; 
        %   ElevationAngle] (in degrees). 
        
        % [1] James Ward, "Space-Time Adaptive Processing for Airborne
        % Radar". MIT Lincoln Lab tech report 1015, 1994.
        % [2] J. R. Guerci, "Space-Time Adaptive Processing for Radar".
        % Artech House, 2003.
        
            % Angle validation is already done by calling function

            angleSteeringVec = step(obj.cSteeringVector,...
                obj.OperatingFrequency,angle);
            
            % preallocate
            N = obj.pArrayNumElements;
            Np = obj.NumDPCAPulses;
            wdpca = complex(zeros(Np*N,1));
            
            % generate DPCA weights
            NumZeroEntries = Np - 1;
            ReducedAperture = N - NumZeroEntries;
            for m = 1:Np
                temp = NumZeroEntries - (m-1);
                % DPCA weights are zeros for non-overlapping elements
                % across pulses and the corresponding steering vector
                % elements for overlapped elements across pulses. The
                % weights put on each pulses are based on binomial
                % coefficients.
                wdpca((m-1)*N+1:m*N) = ...
                    [zeros(temp,1);...
                    (-1)^(m-1)*nchoosek(Np-1,m-1)*...
                    angleSteeringVec(temp+1:temp+ReducedAperture);...
                    zeros(m-1,1)];
            end
        end
               
        function dataShft = prepareDPCAData(obj,x)
        %prepareDPCAData Prepare data cube for DPCA processing
        %   Y = prepareDPCAData(H,X) converts the data cube X into a new
        %   cube Y which is more convenient for the DPCA processor Hs to
        %   work with.  X is a 3-dimensional numeric array whose dimensions
        %   are (range, channels, pulses). Output Y is also a 3-dimensional
        %   numeric array whose dimensions are (NumDPCAPulses*channels,
        %   reducedNumPulses, range).
        
            ftdim = size(x,1);
            N = obj.pArrayNumElements;
            NDPCAPulses = obj.NumDPCAPulses;
            data = complex(zeros(ftdim,NDPCAPulses*N,...
                obj.pReducedDimension));
            for m = 1:NDPCAPulses
                data(:,(m-1)*N+1:m*N,:) = x(:,:,m:end-NDPCAPulses+m);
            end
            dataShft = shiftdim(data,1);
            % make the N pulse space time snapshot the first dimension. the
            % output of the ADPCA is multiple pulses, needs to pass through
            % the doppler filtering.
        end
        
        function rdim = reducedNumPulses(obj,NPulse)
        %reducedNumPulses Reduced pulse dimensionality after DPCA
        %   RDIM = reducedNumPulses(H,X) returns the reduced pulse
        %   dimensionality after the spatial DPCA specified in the DPCA
        %   processor Hs is performed on the data cube X. X must be a
        %   3-dimensional numeric array and its third dimension specifies
        %   the original number of pulses.

            rdim = NPulse-(obj.NumDPCAPulses-1);
        end

        function [rdata, w] = preDopplerData(obj,x,CUTIdx,angle)
        %preDopplerData Return pre-Doppler processed data 
        %   RDATA = preDopplerData(H,X,IDX,ANG) returns the reduced
        %   dimension data after performing spatial domain DPCA based
        %   processing defined in Hs. The reduced dimension data can then
        %   be further processed in the doppler domain.  X must be a
        %   3-dimensional data cube whose dimensions are [range, channels,
        %   pulses]. The processing focused on the cell specified by the
        %   index IDX. The look direction is specified in ANG. ANG
        %   is a 2x1 vector in the form of [AzimuthAngle; ElevationAngle]
        %   (in degrees).
        %
        %   [RDATA, W] = preDopplerData(...) also returns the pre-Doppler
        %   weights.
        
            % DPCA processing
            tempx = prepareDPCAData(obj,x);
                       
            % Calculate pre-doppler weights
            w = calcPreDopplerWeights(obj,angle,tempx,CUTIdx);
            
            rdim = obj.pReducedDimension;
            rdata = complex(zeros(size(x,1),rdim));   
            
            % Apply pre-doppler weights, the output is a matrix with
            % reduced dimension, range x reduced pulse dim.
            for m = 1:rdim
                rdata(:,m) = (w(:,m)'*squeeze(tempx(:,m,:))).'; % match dimension
            end
        end
        
        function w = calcPostDopplerWeights(obj,wdpca,dp,prf)
            % Pre-allocate
            spatialDim = obj.pArrayNumElements;
            w = complex(zeros(spatialDim*obj.pNumPulses,1));
            % Calculate Doppler steering vector
            reducedDim = obj.pReducedDimension;
            dpSteeringVec = dopsteeringvec(dp,reducedDim,prf);
            % Use Doppler steering vector and pre-Doppler weights to
            % calculate total weights
            for m = 1:reducedDim
                tempIdx = (m-1)*spatialDim+1:(m+obj.NumDPCAPulses-1)*spatialDim;
                w(tempIdx) = w(tempIdx) + wdpca(:,m)*dpSteeringVec(m);
            end
        end
    end
    
    methods (Static,Hidden,Access=protected)
      function groups = getPropertyGroupsImpl
        groups = getPropertyGroupsImpl@phased.internal.AbstractSTAP('ula');
        props = 'PreDopplerOutput';
        groups(1).PropertyList = [groups(1).PropertyList props];
      end
    end
    methods (Access = protected) %for Simulink    
        function varargout = getOutputSizeImpl(obj)
            if obj.PreDopplerOutput
                szX = propagatedInputSize(obj,1);
                varargout{1} = [szX(1) szX(3)-1];
                if obj.WeightsOutputPort
                    varargout{2} = [szX(2)*2 szX(3)-1];
                end            
            else
                [varargout{1:nargout}] = getOutputSizeImpl@phased.internal.AbstractSTAP(obj);
            end
        end
    end
end





