classdef (Hidden) AbstractSteppedFMWaveform < phased.internal.AbstractContinuousPhasePulseWaveform
%This class is for internal use only. It may be removed in the future.
       
%   Copyright 2017 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %FrequencyStep Frequency step (Hz)
    %   Specify the linear frequency step size (in Hz) as a positive
    %   scalar. The default value of this property is 2e4 (20 kHz).
    FrequencyStep = 2e4;
end

properties (Nontunable, PositiveInteger) 
    %NumSteps   Number of frequency steps
    %   Specify the number of frequency steps as a positive integer. The
    %   default value of this property is 5. When NumSteps is 1, the
    %   stepped FM waveform reduces to a rectangular waveform.
    NumSteps = 5;
end
    
methods
    function set.FrequencyStep(obj, value)
        validateattributes( value, { 'double' }, { 'scalar', 'positive', 'finite' }, '', 'FrequencyStep');
        obj.FrequencyStep = value;
    end
end

methods (Access = public)
    % Constructor
    function obj = AbstractSteppedFMWaveform(varargin)
        obj@phased.internal.AbstractContinuousPhasePulseWaveform(varargin{:});
    end
    
    function bw = bandwidth(obj,prfidx)
    %bandwidth Bandwidth of the waveform
    %   BW = bandwidth(H) returns the bandwidth (in Hz) of the pulses for
    %   the stepped FM pulse waveform H. If there are N frequency steps,
    %   the bandwidth is N times the frequency step. If there is no
    %   frequency stepping, the bandwidth is the reciprocal of the
    %   pulse width.
    %
    %   BW = bandwidth(H,PRFIDX) specifies the index of the PRF of choice,
    %   PRFIDX, as a positive scalar. The PRF is selected from a predefined
    %   list specified in the PRF property.
    %
    %   % Example:
    %   %   Determine the bandwidth of a stepped FM waveform H.
    %
    %   H = phased.SteppedFMWaveform;
    %   bw = bandwidth(H)
    %
    %   See also phased, phased.SteppedFMWaveform, dutycycle.
    
        phased.internal.narginchk(1,2,nargin);
        validateProperties(obj);
        
        if nargin < 2
            prfidx = 1;
        end
        
        sigdatatypes.validateIndex(prfidx,'bandwidth','PRFIDX',...
            {'scalar','<=',numel(obj.PRF)});
        
        % validate prfidx
        Nsteps = obj.NumSteps;
        if Nsteps < 2
            bw = 1/getPulseWidth(obj);
        else
            bw = obj.FrequencyStep*Nsteps;
        end

    end
        
end
    
methods (Access = protected)
    
    function s = getMatchingWaveform(obj,pidx)
        fs = obj.SampleRate;
        pw = getPulseWidth(obj);
        if obj.PRFSelectionInputPort
            [stepidx,prfidx] = ind2sub([obj.NumSteps,obj.pNumDistinctPRF],...
                pidx);
            nl = phased.internal.val2ind(...
                pw(prfidx),1/fs,0)-1; %non-zero length
            t = (0:nl-1)'/fs;
            s = exp(1i*2*pi*obj.FrequencyStep*t*(stepidx-1));
        else
            nl = phased.internal.val2ind(...
                pw(obj.getCircularIndex(pidx,numel(pw))),1/fs,0)-1; %non-zero length
            t = (0:nl(1)-1)'/fs;
            stepidx = obj.getCircularIndex(pidx,obj.NumSteps);
            s = exp(1i*2*pi*obj.FrequencyStep*t*(stepidx-1));
        end
    end
    
    function calcDistinctPulseParameters(obj)
        if obj.PRFSelectionInputPort
            % certain combination doesn't exist if sequencing through
            % property.
            obj.pNumDistinctPulses = obj.NumSteps*obj.pNumDistinctPRF;
            [~,prfgrid] = ndgrid(1:obj.NumSteps,1:obj.pNumDistinctPRF);
            obj.pDistinctPulseLength = obj.pPulseLength(prfgrid(:).');
            obj.pWaveformIndexForDistinctPulses = 1:obj.pNumDistinctPulses;
            pw = getPulseWidth(obj);
            nzl = phased.internal.val2ind(pw(prfgrid),1/obj.SampleRate,0)-1;
            obj.pNonZeroLength = nzl;
        else
            obj.pNumDistinctPulses = lcm(obj.pNumDistinctPRF,obj.NumSteps);
            obj.pDistinctPulseLength = repmat(obj.pPulseLength,1,...
                obj.pNumDistinctPulses/obj.pNumDistinctPRF);
            obj.pWaveformIndexForDistinctPulses = 1:obj.pNumDistinctPulses;
            pw = getPulseWidth(obj);
            nzl = phased.internal.val2ind(pw(obj.getCircularIndex(...
                    1:obj.pNumDistinctPulses,obj.pNumDistinctPRF)),...
                    1/obj.SampleRate,0)-1;
            obj.pNonZeroLength = nzl;
        end
    end
    
    function pidx = getDefaultMatchedFilterPulseIndex(obj) 
        pidx = 1:obj.NumSteps;
    end
    
    function resetImpl(obj)
        resetImpl@phased.internal.AbstractContinuousPhasePulseWaveform(obj);
        obj.pOutputStartPulseIndex = 1; % track step index
    end
    
    function sidx = getOutputSampleIndex(obj,prfidx)
        if obj.PRFSelectionInputPort
            nsamps = obj.NumSamples;
            sidx = zeros(1,nsamps);
            remainingSamples = nsamps;
            while remainingSamples > 0
                if obj.pOutputStartSampleIndex > 0
                    startidx = obj.pOutputStartSampleIndex;
                    tempidx = find(startidx >= obj.pDistinctNonZeroStart,1,'last');
                    if ~isempty(tempidx)
                        currentPulseIdx = tempidx(1);
                    else
                        currentPulseIdx = 1;
                    end
                    % currentStepIdx = obj.pWaveformIndexForDistinctPulses(currentPulseIdx);
                    if startidx ~= obj.pDistinctNonZeroStart(currentPulseIdx)
                        % just finish a pulse from the previous time
                        currentPulseStartIdx = obj.pDistinctNonZeroStart(currentPulseIdx);
                        currentPulseRemainingSamples = obj.pDistinctPulseLength(currentPulseIdx)-...
                            (startidx-currentPulseStartIdx);
                        currentPulseNumSamples = min(remainingSamples,currentPulseRemainingSamples);
                        currentPulseNumSamples = min(currentPulseNumSamples,nsamps); % codegen
                        sidx(nsamps-remainingSamples+1+(0:currentPulseNumSamples-1)) = ...
                            startidx+(0:currentPulseNumSamples-1);
                        remainingSamples = remainingSamples - currentPulseNumSamples;
                        if currentPulseRemainingSamples == currentPulseNumSamples
                            % all remaining samples for a pulse are out
                            % should point to the starting of the next
                            % pulse step. But since if it's the beginning of
                            % the pulse, it always recompute the starting
                            % index in the other branch, simply set to 0.
                            obj.pOutputStartSampleIndex = 0;
                            % use pOutputStartPulseIndex to track next step
                            obj.pOutputStartPulseIndex = obj.getCircularIndex(obj.pOutputStartPulseIndex+1,obj.NumSteps);
                        else % currentPulseRemainingSamples > currentPulseNumSamples
                            % in this case remainingSamples must be 0
                            obj.pOutputStartSampleIndex = sidx(end)+1;
                        end
                    end
                else
                    DistinctPulseIdx = sub2ind([obj.NumSteps obj.pNumDistinctPRF],...
                        obj.pOutputStartPulseIndex,prfidx);
                    startidx = obj.pDistinctNonZeroStart(DistinctPulseIdx);
                    currentPulseRemainingSamples = obj.pDistinctPulseLength(DistinctPulseIdx);
                    currentPulseNumSamples = min(remainingSamples,currentPulseRemainingSamples);
                    currentPulseNumSamples = min(currentPulseNumSamples,nsamps); % codegen
                    sidx(nsamps-remainingSamples+1+(0:currentPulseNumSamples-1)) = ...
                        startidx+(0:currentPulseNumSamples-1);
                    remainingSamples = remainingSamples - currentPulseNumSamples;
                    if currentPulseRemainingSamples == currentPulseNumSamples
                        obj.pOutputStartSampleIndex = 0;
                        obj.pOutputStartPulseIndex = obj.getCircularIndex(obj.pOutputStartPulseIndex+1,obj.NumSteps);
                    else % currentPulseRemainingSamples > currentPulseNumSamples
                        % in this case remainingSamples must be 0
                        obj.pOutputStartSampleIndex = sidx(end)+1;
                    end
                end
            end
        else
            sidx = obj.pOutputStartSampleIndex+obj.pOutputSampleInterval;
            sidx = obj.getCircularIndex(sidx,obj.pNumDistinctSamples);
            % reset start output sample index for next simulation
            obj.pOutputStartSampleIndex = sidx(end);
        end
            
    end
    
    function pidx = getOutputPulseIndex(obj,prfidx)
        if obj.PRFSelectionInputPort
            % use start pulse index to track the step index
            stepstartidx = obj.pOutputStartPulseIndex;
            stepidx = obj.getCircularIndex(stepstartidx+obj.pOutputPulseInterval,obj.NumSteps);
            DistinctPulseIdx = sub2ind([obj.NumSteps obj.pNumDistinctPRF],...
                stepidx,repmat(prfidx,1,numel(obj.pOutputPulseInterval)));
            obj.pOutputStartPulseIndex = stepidx(end); 
            pidx = DistinctPulseIdx(1:end-1);           
            
        else
            pidx = obj.pOutputStartPulseIndex+obj.pOutputPulseInterval;
            pidx = obj.getCircularIndex(pidx,obj.pNumDistinctPulses);
            % reset start pulse index for next simulation step
            obj.pOutputStartPulseIndex = pidx(end);
        end
    end
    
    function wname = getWaveformName(obj) %#ok<MANU>
        wname = 'Stepped FM pulse waveform';
    end
   
end

end

