classdef (Hidden) AbstractLinearFMWaveform < phased.internal.AbstractContinuousPhasePulseWaveform
%This class is for internal use only. It may be removed in the future.
   
%   Copyright 2017 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005


%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %SweepBandwidth Sweep bandwidth (Hz)
    %   Specify the bandwidth of the linear FM sweeping (in Hz) as a
    %   positive scalar. The default value is 1e5 (100 kHz). 
    SweepBandwidth = 1e5;
    %SweepDirection Sweep direction
    %   Specify the direction of the linear FM sweep as one of 'Up' |
    %   'Down'. The default value is 'Up'.
    SweepDirection = 'Up';
    %SweepInterval  Sweep interval
    %   Specify the linear FM sweeping interval using one of 'Positive' |
    %   'Symmetric', where the default is 'Positive'. If SweepInterval is
    %   'Positive', the waveform sweeps in the interval between 0 and B
    %   where B is the sweeping bandwidth. If SweepInterval is 'Symmetric',
    %   the waveform sweeps in the interval between -B/2 and B/2.
    SweepInterval = 'Positive'
    %Envelope Envelope function
    %   Specify the envelope function as one of 'Rectangular' | 'Gaussian'.
    %   The default value is 'Rectangular'.
    Envelope = 'Rectangular';
end

methods
    function obj = set.SweepBandwidth(obj, value)
        validateattributes(value,{'double'},...
            {'scalar','positive','finite'},...
            '','SweepBandwidth');
        obj.SweepBandwidth = value;
    end
end

properties(Constant, Hidden)
    SweepDirectionSet = matlab.system.StringSet({'Up','Down'});
    SweepIntervalSet = matlab.system.StringSet({'Positive','Symmetric'});
    EnvelopeSet = matlab.system.StringSet({'Rectangular','Gaussian'});
end

methods 

    function obj = AbstractLinearFMWaveform(varargin)
        obj@phased.internal.AbstractContinuousPhasePulseWaveform(varargin{:});
    end   
    
    function h = getStretchProcessor(obj,refrng,rngspan,c,prfidx)
    %getStretchProcessor Create stretch processor for the waveform
    %   HS = getStretchProcessor(H) returns the stretch processor, HS, for
    %   the waveform. HS is a phased.StretchProcessor object. HS is set up
    %   so the reference range corresponds to 1/4 of the maximum
    %   unambiguous range of a pulse. The range span corresponds to 1/10 of
    %   the distance traveled by the wave within the pulse width. The
    %   propagation speed is the speed of light. When the PRF property has
    %   multiple entries, the first PRF is used.
    %
    %   HS = getStretchProcessor(H,REFRNG) specifies the reference range
    %   (in meters), REFRNG, as a positive scalar. 
    %
    %   HS = getStretchProcessor(H,REFRNG,RNGSPAN) specifies the range span
    %   (in meters), RNGSPAN, as a positive scalar. RNGSPAN is centered at
    %   REFRNG.
    %
    %   HS = getStretchProcessor(H,REFRNG,RNGSPAN,V) specifies the
    %   propagation speed (in m/s), V, as a positive scalar.
    %
    %   HS = getStretchProcessor(H,REFRNG,RNGSPAN,V,PRFIDX) specifies the
    %   index of the PRF of choice, PRFIDX, as a positive scalar. The PRF
    %   is selected from a predefined list specified in the PRF property.
    %
    %   % Example:
    %   %   Create a stretch processor for a linear FM waveform with a 
    %   %   reference range of 5000 meters and a range span of 500 meters.
    %
    %   hwav = phased.LinearFMWaveform;
    %   hs = getStretchProcessor(hwav,5000,500)
    %
    %   See also phased.LinearFMWaveform, phased.StretchProcessor.
    
        phased.internal.narginchk(1,5,nargin);
        validateProperties(obj);
        prf_vec = obj.PRF;
        if nargin < 5
            prfidx = 1;
        else
            sigdatatypes.validateIndex(prfidx,'','PRFIDX',...
                {'scalar','<=',numel(prf_vec)});
        end
        prf = prf_vec(prfidx);
        pw_vec = getPulseWidth(obj);
        pw = pw_vec(prfidx);
        if nargin < 4
            c = physconst('lightspeed');
        end
        if nargin < 3
            rngspan = c*pw/20;
        end
        if nargin < 2
            rmax = c/(2*prf);
            refrng = rmax/4;
        end
        sigdatatypes.validateSpeed(c,'','V',{'scalar','positive'});
        sigdatatypes.validateDistance(refrng,'','REFRNG',...
            {'scalar','positive'});
        sigdatatypes.validateDistance(rngspan,'','RNGSPAN',...
            {'scalar','positive'});
        phased.StretchProcessor.validateStretchRangeSpan(...
            c,prf,pw,refrng,rngspan);
        
        slope = obj.SweepBandwidth/pw;
        if (obj.SweepDirection(1) == 'D') %Down
            sSlope = -slope;
        else
            sSlope = slope;
        end
        h = phased.StretchProcessor(...
            'SampleRate',obj.SampleRate,...
            'PRF',prf,...
            'PulseWidth',pw,...
            'SweepSlope',sSlope,...
            'SweepInterval',obj.SweepInterval,...
            'PropagationSpeed',c,...
            'ReferenceRange',refrng,...
            'RangeSpan',rngspan);
    end
    
    function bw = bandwidth(obj,prfidx) 
    %bandwidth Bandwidth of the waveform
    %   BW = bandwidth(H) returns the bandwidth (in Hz) of the pulses for
    %   the linear FM pulse waveform H. The bandwidth is the same as the
    %   sweep bandwidth.
    %
    %   BW = bandwidth(H,PRFIDX) specifies the index of the PRF of choice,
    %   PRFIDX, as a positive scalar. The PRF is selected from a predefined
    %   list specified in the PRF property.
    %
    %   % Example:
    %   %   Determine the bandwidth of a linear FM pulse waveform H.
    %
    %   H = phased.LinearFMWaveform;
    %   bw = bandwidth(H)
    %
    %   See also phased, phased.LinearFMWaveform, dutycycle.

        phased.internal.narginchk(1,2,nargin);
        validateProperties(obj);
        
        if nargin < 2
            prfidx = 1;
        end
        
        sigdatatypes.validateIndex(prfidx,'bandwidth','PRFIDX',...
            {'scalar','<=',numel(obj.PRF)});
        
        % validate prfidx
        bw = obj.SweepBandwidth;
    end

end
    
methods (Access = protected)
    
    function wav = getMatchingWaveform(obj,pidx)
        tau_in = getPulseWidth(obj);
        tau = tau_in(obj.getCircularIndex(pidx,numel(tau_in)));
        fs = obj.SampleRate;
        phaseslope = obj.SweepBandwidth/tau;
        nl = phased.internal.val2ind(tau,1/fs,0)-1; %non-zero length
        t = (0:nl-1)'/fs;
        if obj.SweepDirection(1) == 'U' %Up
            if (obj.SweepInterval(1) == 'P') %Positive
                s = exp(1i*pi*phaseslope*t.^2);
            else
                s = exp(1i*pi*phaseslope*t.*(t-tau));
            end
        else
            if (obj.SweepInterval(1) == 'P')%Positive
                s =  exp(1i*pi*phaseslope*t.*(2*tau-t));
            else
                s = exp(1i*pi*phaseslope*t.*(tau-t));
            end
        end
        
        if (obj.Envelope(1) == 'G') %Gaussian
            wav =  exp(-t.^2/tau^2).*s;
        else
            wav = s;   %needed by codegen to constant fold
        end
    end
    
    function calcDistinctPulseParameters(obj)
        obj.pNumDistinctPulses = obj.pNumDistinctPRF;
        obj.pDistinctPulseLength = obj.pPulseLength;
        obj.pWaveformIndexForDistinctPulses = 1:obj.pNumDistinctPulses;
        pw = getPulseWidth(obj);
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
                if obj.pOutputStartSampleIndex > 0
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
        wname = 'Linear FM pulse waveform';
    end
   
end

end
