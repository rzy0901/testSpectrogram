classdef (Hidden) AbstractTimeDomainSMIBeamformer < phased.internal.AbstractTimeDomainBeamformer
%This class is for internal use only. It may be removed in the future.

%AbstractTimeDomainSMIBeamformer   Abstract class for time domain SMI beamformers

%   Copyright 2010-2014 The MathWorks, Inc.   


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen
    
    properties
        %DiagonalLoadingFactor  Diagonal loading factor
        %   Specify the diagonal loading factor as a positive scalar. The
        %   default value of this property is 0. Diagonal loading is a
        %   technique used to achieve robust beamforming performance,
        %   especially when the sample support is small. This property is
        %   tunable.
        DiagonalLoadingFactor = 0;
    end

    properties (Nontunable, PositiveInteger) 
        %FilterLength     FIR filter length
        %   Specify the length of FIR filter behind each sensor elements in
        %   the array as a positive integer. The default value of this
        %   property is 1.
        FilterLength = 1;
    end

    properties (Nontunable, Logical) 
        %TrainingInputPort  Enable training data input
        %   Set this property to true to add input to specify the
        %   additional training data. Set this property to false to not add
        %   input to specify the training data. The default value of this
        %   property is false. When this property is false, the input
        %   signal itself is used as the training data.
        TrainingInputPort = false;
    end
    
    properties (Nontunable, Logical) 
        %WeightsOutputPort    Enable weights output
        %   Set this property to true to output the weights used in the
        %   beamformer. Set this property to false to not output the
        %   weights. The default value of this property is false.
        WeightsOutputPort = false;
    end
    
    properties(Access = protected, Nontunable)
        cWeightsEstimator
        pSpaceTimeSnapshotDimension
    end
    
    properties(Access = protected)
        pDataBuffer
    end
    
    methods
        function set.DiagonalLoadingFactor(obj,val)
            validateattributes( val, { 'double','single' }, { 'scalar', 'nonnegative', 'real', 'finite', 'nonempty' }, '', 'DiagonalLoadingFactor');
            obj.DiagonalLoadingFactor = val;
        end
    end

    methods (Access = protected)
        function obj = AbstractTimeDomainSMIBeamformer(varargin)
            obj@phased.internal.AbstractTimeDomainBeamformer(varargin{:});

        end
    end
    
    methods (Access = protected)
        function validatePropertiesImpl(obj)
            validatePropertiesImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
        end
        
        function setupImpl(obj,x,varargin)

            setupImpl@phased.internal.AbstractTimeDomainBeamformer(obj,x);
            obj.pSpaceTimeSnapshotDimension = ...
                obj.pDOF*obj.FilterLength;
            obj.cWeightsEstimator = phased.internal.BlockSMIWeightsEstimator;
            processTunedPropertiesImpl(obj);
            if isreal(x)
                obj.pDataBuffer = zeros(obj.FilterLength-1,obj.pDOF,class(x));
            else
                obj.pDataBuffer = complex(zeros(obj.FilterLength-1,obj.pDOF,class(x)));
            end
        end
        
        function flag = isInputComplexityLockedImpl(obj,index) 
            flag = isInputComplexityLockedImpl@phased.internal.AbstractTimeDomainBeamformer(obj,index);
            if obj.TrainingInputPort
                if index == 2
                    flag = false;
                end
                if strncmpi(obj.DirectionSource,'i',1) && (index == 3)
                    flag = true;
                end
            else
                if strncmpi(obj.DirectionSource,'i',1) && (index == 2)
                    flag = true;
                end
            end
        end
        
        function processTunedPropertiesImpl(obj)
            obj.cWeightsEstimator.DiagonalLoadingFactor = ...
                obj.DiagonalLoadingFactor;
        end
        
        function releaseImpl(obj)
            releaseImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            release(obj.cWeightsEstimator);
        end
        
        function resetImpl(obj)
            resetImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            reset(obj.cWeightsEstimator);
            obj.pDataBuffer(:) = 0;
        end

        function num = getNumInputsImpl(obj)
            num = getNumInputsImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            if obj.TrainingInputPort
                num = num + 1;
            end
        end
        
        function num = getNumOutputsImpl(obj)
            num = getNumOutputsImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            if obj.WeightsOutputPort
                num = 2;
            end
        end
        
        function validateInputsImpl(obj,x,varargin)
            validateInputsImpl@phased.internal.AbstractTimeDomainBeamformer(obj,x);
                     
            if obj.TrainingInputPort
                xt = varargin{1};
                if strncmpi(obj.DirectionSource,'i',1)
                    ang = varargin{2};
                else
                    ang = [];
                end
            else
                xt = [];
                if strncmpi(obj.DirectionSource,'i',1)
                   ang = varargin{1};
                else
                   ang = [];
                end
            end
            
            if ~isempty(ang)
                validateInputAngle(obj,ang);
            end
            if ~isempty(xt)
                validateInputSignal(obj,xt,'XT');
                sz_xt = size(xt);
                validateTimeDomainTrainingSignalSize(obj,sz_xt(1),'XT');
            end
            
            sz_x = size(x);
            validateTimeDomainTrainingSignalSize(obj,sz_x(1),'X');
                    
        end
        
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            if isLocked(obj)
                s.cWeightsEstimator = saveobj(obj.cWeightsEstimator);
                s.pSpaceTimeSnapshotDimension = obj.pSpaceTimeSnapshotDimension;
                s.pDataBuffer = obj.pDataBuffer;
            end
        end
        
        function s = loadSubObjects(obj,s)
            s = loadSubObjects@phased.internal.AbstractTimeDomainBeamformer(obj,s);
            if isfield(s,'isLocked')
                if s.isLocked
                    obj.cWeightsEstimator = phased.internal.BlockSMIWeightsEstimator.loadobj(s.cWeightsEstimator);
                    s = rmfield(s,'cWeightsEstimator');
                end
                s = rmfield(s,'isLocked');
            end
        end
        
        function [y, w] = stepImpl(obj,x,varargin)

            classtouse=class(x);
            anginputflag = strncmpi(obj.DirectionSource,'i',1);
            
            if obj.TrainingInputPort
                xt = cast(varargin{1},classtouse);
                if anginputflag
                    ang = varargin{2};
                else
                    ang = obj.Direction;
                end
                x_presteered_temp = steer(obj,x,ang);
                xt_presteered_temp = steer(obj,xt,ang);
                if obj.FilterLength > 1
                    x_presteered = prepareSpaceTimeSnapshot(obj,x_presteered_temp,true);
                    xt_presteered = prepareSpaceTimeSnapshot(obj,xt_presteered_temp,false);
                else
                    x_presteered = x_presteered_temp;
                    xt_presteered = xt_presteered_temp;
                end
            else
                if anginputflag
                    ang = varargin{1};
                else
                    ang = obj.Direction;
                end
                x_presteered_temp = steer(obj,x,ang);
                if obj.FilterLength > 1         
                    x_presteered = prepareSpaceTimeSnapshot(obj,x_presteered_temp,true);
                else
                    x_presteered = x_presteered_temp;
                end
                xt_presteered = x_presteered;
            end
            
            w = step(obj.cWeightsEstimator,xt_presteered);
            y = x_presteered*conj(w);

        end
        
        function validateTimeDomainTrainingSignalSize(obj,NumTrainingSamples,str_signal)
            cond = NumTrainingSamples < obj.FilterLength;
            if cond
                coder.internal.errorIf(cond,'phased:TimeDomainBeamformer:NotEnoughSamples', str_signal, obj.FilterLength);
            end
        end
        
        function spaceTimeSnap = prepareSpaceTimeSnapshot(obj,x,streamflag)
            
            classtouse=class(x);
            % convert MxN signal to MxSpaceTimeSnapDimension
            % Each row, the organization can be considered as
            % [FilterTap1ForElements FilterTap2ForElements ...]
            % Each row corresponds to a time instance. For example, at
            % first row, signal only at first tap, nothing in other taps
            
            [M,N] = size(x);
            if isreal(x)
                spaceTimeSnap = zeros(M,obj.pSpaceTimeSnapshotDimension,classtouse);
            else
                spaceTimeSnap = complex(zeros(M,obj.pSpaceTimeSnapshotDimension,classtouse));
            end
            filtlen = obj.FilterLength;
            if streamflag
                x_augmented = [obj.pDataBuffer; x];
                for m = 1:filtlen
                    startIdx = filtlen-m+1;
                    spaceTimeSnap(:,(m-1)*N+1:m*N) = ...
                        x_augmented(startIdx:startIdx+M-1,:);   % equation (6.762) in Van Trees, 2002
                end
                obj.pDataBuffer = x_augmented(M+1:end,:);
            else
                for m = 1:filtlen
                    spaceTimeSnap(m:end,(m-1)*N+1:m*N) = ...
                        x(1:end-m+1,:);
                end
            end
        end
            

    end
    methods (Access = protected) %for Simulink
        function varargout = getInputNamesImpl(obj)
            %Insert XT
            if obj.TrainingInputPort
                [varargout{1:nargout-1}] = getInputNamesImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
                varargout = {varargout{1} 'XT' varargout{2:end}};
            else
                [varargout{1:nargout}] = getInputNamesImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            end
        end
        function varargout = getOutputSizeImpl(obj)
            [varargout{1:nargout}] = getOutputSizeImpl@phased.internal.AbstractTimeDomainBeamformer(obj);
            if obj.WeightsOutputPort
                szW = varargout{2};
                varargout{2} = [szW(1)*obj.FilterLength szW(2)];
            end
        end
        function aIdx = getAngInputIdx(obj)
            if obj.TrainingInputPort
                aIdx = 3;
            else
                aIdx = 2;
            end
        end
    end

end


