classdef (Hidden) AbstractRectangularWaveform < phased.internal.AbstractContinuousPhasePulseWaveform
%This class is for internal use only. It may be removed in the future.

%   Copyright 2017 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

methods (Access = public)
    function obj = AbstractRectangularWaveform(varargin)
        obj@phased.internal.AbstractContinuousPhasePulseWaveform(varargin{:});
    end
    
    function bw = bandwidth(obj,prfidx)
    %bandwidth Bandwidth of the waveform
    %   BW = bandwidth(H) returns the bandwidth (in Hz) of the pulses for
    %   the rectangular pulse waveform H. The bandwidth is the reciprocal
    %   of the pulse width.
    %
    %   BW = bandwidth(H,PRFIDX) specifies the index of the PRF of choice,
    %   PRFIDX, as a positive scalar. The PRF is selected from a predefined
    %   list specified in the PRF property.
    %
    %   % Example:
    %   %   Determine the bandwidth of a rectangular pulse waveform H.
    %
    %   H = phased.RectangularWaveform;
    %   bw = bandwidth(H)
    %
    %   See also phased, phased.RectangularWaveform, dutycycle.
        phased.internal.narginchk(1,2,nargin);
        validateProperties(obj);
        
        if nargin < 2
            prfidx = 1;
        end
        
        sigdatatypes.validateIndex(prfidx,'bandwidth','PRFIDX',...
            {'scalar','<=',numel(obj.PRF)});
        
        % validate prfidx
        temp = 1./getPulseWidth(obj);
        bw = temp(prfidx);
    end
    
end

methods (Access = protected)    
    
    function s = getMatchingWaveform(obj,pidx)
        fs = obj.SampleRate;
        pw = getPulseWidth(obj);
        nl = phased.internal.val2ind(...
            pw(obj.getCircularIndex(pidx,numel(pw))),1/fs,0)-1; %non-zero length
        s = ones(nl,1);
    end
    
    function calcDistinctPulseParameters(obj)
        obj.pNumDistinctPulses = obj.pNumDistinctPRF;
        obj.pDistinctPulseLength = obj.pPulseLength;
        obj.pWaveformIndexForDistinctPulses = 1:obj.pNumDistinctPulses;
        pw = obj.pPulseWidth;
        nzl = phased.internal.val2ind(pw,1/obj.SampleRate,0)-1;
        obj.pNonZeroLength = nzl;
    end
    
    function pidx = getDefaultMatchedFilterPulseIndex(obj) %#ok<MANU>
        pidx = 1;
    end
    
    function sidx = getOutputSampleIndex(obj,prfidx)
        if obj.PRFSelectionInputPort
            nsamps = obj.NumSamples;
            sidx = zeros(1,nsamps);
            remainingSamples = nsamps;
            while remainingSamples > 0
                if obj.pOutputStartSampleIndex>0
                    startidx = obj.pOutputStartSampleIndex;
                    tempidx = find(startidx >= obj.pDistinctNonZeroStart,1,'last');
                    if ~isempty(tempidx)
                        currentPulseIdx = tempidx(1);
                    else
                        currentPulseIdx = 1;
                    end
                    if startidx ~= obj.pDistinctNonZeroStart(currentPulseIdx)
                        % just finish a pulse the previous time
                        currentPulseStartIdx = obj.pDistinctNonZeroStart(currentPulseIdx);
                        currentPulseRemainingSamples = ...
                            obj.pDistinctPulseLength(currentPulseIdx)-(startidx-currentPulseStartIdx);
                        currentPulseNumSamples = min(remainingSamples,...
                            currentPulseRemainingSamples);
                        currentPulseNumSamples = min(currentPulseNumSamples,nsamps); % codegen
                        sidx(nsamps-remainingSamples+1+(0:currentPulseNumSamples-1)) = ...
                            startidx+(0:currentPulseNumSamples-1);
                        remainingSamples = remainingSamples - currentPulseNumSamples;
                        if currentPulseRemainingSamples == currentPulseNumSamples 
                            % all remaining samples for a pulse are out
                            % should point to the starting of the current
                            % pulse. But since if it's the beginning of the
                            % pulse, it always recompute the starting index
                            % in the other branch, simply set to 0.
                            obj.pOutputStartSampleIndex = 0;
                        else % currentPulseRemainingSamples > currentPulseNumSamples
                            % in this case remainingSamples must be 0
                            obj.pOutputStartSampleIndex = sidx(end)+1;
                        end
                    end

                else
                    startidx = obj.pDistinctNonZeroStart(prfidx);
                    currentPulseRemainingSamples = obj.pDistinctPulseLength(prfidx);
                    currentPulseNumSamples = min(remainingSamples,...
                        currentPulseRemainingSamples);
                    currentPulseNumSamples = min(currentPulseNumSamples,nsamps); % codegen
                    sidx(nsamps-remainingSamples+1+(0:currentPulseNumSamples-1)) = ...
                        startidx+(0:currentPulseNumSamples-1);
                    remainingSamples = remainingSamples - currentPulseNumSamples;
                    if currentPulseRemainingSamples == currentPulseNumSamples 
                        % all remaining samples for a pulse are out
                        obj.pOutputStartSampleIndex = 0;
                    else % currentPulseRemainingSamples > currentPulseNumSamples
                        % in this case remainingSamples must be 0
                        obj.pOutputStartSampleIndex = sidx(end)+1;
                    end

                end
            end
        else  % no input port
            sidx = obj.pOutputStartSampleIndex+obj.pOutputSampleInterval;
            sidx = obj.getCircularIndex(sidx,obj.pNumDistinctSamples);
            % reset start output sample index for next simulation
            obj.pOutputStartSampleIndex = sidx(end);
        end
            
    end
    
    function pidx = getOutputPulseIndex(obj,prfidx)
        if obj.PRFSelectionInputPort
            pidx = prfidx*ones(1,obj.NumPulses);
        else
            pidx = obj.pOutputStartPulseIndex+obj.pOutputPulseInterval;
            pidx = obj.getCircularIndex(pidx,obj.pNumDistinctPulses);
            % reset start pulse index for next simulation step
            obj.pOutputStartPulseIndex = pidx(end);
        end
    end
    
    function wname = getWaveformName(obj) %#ok<MANU>
        wname = 'Rectangular pulse waveform';
    end
   
end

end 
