classdef (Hidden) AbstractPhaseCodedWaveform < phased.internal.AbstractPulseWaveform
%This class is for internal use only. It may be removed in the future.

%   Copyright 2017 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

properties (Nontunable)
    %Code   Phase code 
    %   Specify type of code used in phase modulation as one of 'Frank' |
    %   'P1' | 'P2' | 'Px' | 'Zadoff-Chu' | 'P3' | 'P4' | 'Barker', where
    %   the default is 'Frank'.
    Code = 'Frank';
    %ChipWidth     Chip width (s)
    %   Specify the duration of each chip (in seconds) in a phase-coded
    %   waveform as a positive scalar. The default value of this property
    %   is 1e-5.
    ChipWidth = 1e-5;
end

properties (Dependent, Hidden)
    %Type   Phase code 
    %   Specify type of code used in phase modulation as one of 'Frank' |
    %   'P1' | 'P2' | 'Px' | 'Zadoff-Chu' | 'P3' | 'P4' | 'Barker', where
    %   the default is 'Frank'.
    Type = 'Frank';
end

properties (PositiveInteger,Nontunable)
    %NumChips   Number of chips
    %   Specify the number of chips in a phase-coded waveform as a positive
    %   integer. The default value of this property is 4.
    %
    %   When you set the Code property to 'Frank', 'P1', 'P2', or 'Px', the
    %   value of this property must be a perfect square. In addition, when
    %   you set the Code property to 'P2', the value of this property must
    %   also be even.
    %
    %   When you set the Code property to 'Barker', the value of this
    %   property must be one of the following: 2, 3, 4, 5, 7, 11, and 13.
    NumChips = 4;
    %SequenceIndex  Zadoff-Chu sequence index
    %   Specify the sequence index used in Zadoff-Chu code as a positive
    %   integer. This property is applicable when you set the Code property
    %   to 'Zadoff-Chu'. The value of this property must be relatively
    %   prime to the value of the NumChips property. The default value of
    %   this property is 1.
    SequenceIndex = 1;
end

properties (Access = private, PositiveInteger, Nontunable)
    % Oversample ratio between SampleRate and ChipWidth, i.e., how many
    % samples are used for each chip.
    pOversampleRatio;
end

properties(Constant, Hidden)
    CodeSet = matlab.system.StringSet({'Barker','Frank','P1','P2',...
        'P3','P4','Px','Zadoff-Chu'});
end

methods
    function set.ChipWidth(obj,value)
        sigdatatypes.validateDuration(value,'','ChipWidth',{'scalar'});
        obj.ChipWidth = value;
    end
    
    function set.Type(obj,value)
        coder.internal.warning('phased:Waveform:ReplaceTypeInPhaseCoded',...
            'Type','phased.PhasedCodedWaveform','Code','Type');
        obj.Code = value;
    end
    
    function value = get.Type(obj)
        coder.internal.warning('phased:Waveform:ReplaceTypeInPhaseCoded',...
            'Type','phased.PhasedCodedWaveform','Code','Type');
        value = obj.Code;
    end
end

methods

    function obj = AbstractPhaseCodedWaveform(varargin)
        obj@phased.internal.AbstractPulseWaveform(varargin{:});
    end

    function bw = bandwidth(obj,prfidx)
    %bandwidth Bandwidth of the waveform
    %   BW = bandwidth(H) returns the bandwidth (in Hz) of the pulses for
    %   the phase-coded pulse waveform H. The bandwidth is the reciprocal
    %   of the chip width.
    %
    %   BW = bandwidth(H,PRFIDX) specifies the index of the PRF of choice,
    %   PRFIDX, as a positive scalar. The PRF is selected from a predefined
    %   list specified in the PRF property.
    %
    %   % Example:
    %   %   Determine the bandwidth of a Frank code waveform H.
    %
    %   H = phased.PhaseCodedWaveform;
    %   bw = bandwidth(H)
    %
    %   See also phased, phased.PhaseCodedWaveform, dutycycle.

        phased.internal.narginchk(1,2,nargin);
        validateProperties(obj);
        
        if nargin < 2
            prfidx = 1;
        end
               
        sigdatatypes.validateIndex(prfidx,'bandwidth','PRFIDX',...
            {'scalar','<=',numel(obj.PRF)});
        
        % validate prfidx
        bw = 1/obj.ChipWidth;
    end
end

methods (Access = protected)
    function validatePropertiesImpl(obj)
        validatePropertiesImpl@phased.internal.AbstractPulseWaveform(obj);
        
        temp = obj.SampleRate*obj.ChipWidth;
        cond =  abs(temp-round(temp)) > eps(temp);
        if cond
            coder.internal.errorIf(cond, ...
                                   'phased:Waveform:NeedProductInteger', 'SampleRate', 'ChipWidth'); 
        end
        cond =   any(obj.ChipWidth.*obj.NumChips.*obj.PRF > 1) ;
        if cond
            coder.internal.errorIf(cond, ...
                                   'phased:Waveform:NotLessThanOrEqualTo', 'ChipWidth', sprintf( '%5.2e', 1/(max(obj.PRF)*obj.NumChips)));
        end
        cond = (strcmp(obj.Code,'Frank') || strcmp(obj.Code,'P1') || ...
                strcmp(obj.Code,'P2') || strcmp(obj.Code,'Px')) && ...
                rem(sqrt(obj.NumChips),1);
        if cond
            coder.internal.errorIf( cond, ...
                 'phased:Waveform:NeedSquare','NumChips');  
        end
        cond =  strcmp(obj.Code,'P2') && rem(obj.NumChips,2) ;
        if cond
            coder.internal.errorIf(cond, ...
                 'phased:Waveform:NeedEven','NumChips');  
        end
        cond =  strcmp(obj.Code,'Zadoff-Chu') && gcd(obj.NumChips,obj.SequenceIndex) ~= 1 ;
        if cond
            coder.internal.errorIf(cond, ...
                 'phased:Waveform:NeedPrime','NumChips','SequenceIndex'); 
        end
        cond =  strcmp(obj.Code,'Barker') && all(obj.NumChips ~= [2 3 4 5 7 11 13]) ;
        if cond
            coder.internal.errorIf(cond, ...
                 'phased:Waveform:InvalidBarkerLength','NumChips');
        end
    end
    
    function setupImpl(obj)
        setupImpl@phased.internal.AbstractPulseWaveform(obj);
        obj.pOversampleRatio = round(obj.SampleRate*obj.ChipWidth);
    end
        
    function flag = isInactivePropertyImpl(obj, prop)
        flag = isInactivePropertyImpl@phased.internal.AbstractPulseWaveform(obj, prop);
        if strcmp(prop,'SequenceIndex') && ...
                ~strcmp(obj.Code,'Zadoff-Chu')
            flag = true;
        end
    end
    
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@phased.internal.AbstractPulseWaveform(obj);
        if isLocked(obj)
            s.pOversampleRatio = obj.pOversampleRatio;
        end
    end
    
    function loadObjectImpl(obj,s,~)
        s = loadSubObjects(obj,s);
        if isfield(s,'Type')
            obj.Code = s.Type;
            s = rmfield(s,'Type');
        end
        fn = fieldnames(s);
        for m = 1:numel(fn)
            obj.(fn{m}) = s.(fn{m});
        end
    end
    
    function value = getPulseWidth(obj)
        value = repmat(obj.ChipWidth*obj.NumChips,1,obj.pNumDistinctPRF);
    end

end

methods (Access = protected)
    
    function s = getMatchingWaveform(obj,pidx) %#ok<INUSD>
        switch obj.Code
            case 'Frank'
                s = phased.internal.frankcode(obj.NumChips);
            case 'P1'
                s = phased.internal.p1code(obj.NumChips);
            case 'P2'
                s = phased.internal.p2code(obj.NumChips);
            case 'Px'
                s = phased.internal.pxcode(obj.NumChips);
            case 'Zadoff-Chu'
                s = phased.internal.zadoffchucode(obj.NumChips,obj.SequenceIndex);
            case 'P3'
                s = phased.internal.p3code(obj.NumChips);
            case 'P4'
                s = phased.internal.p4code(obj.NumChips);
            case 'Barker'
                s = phased.internal.barkercode(obj.NumChips);
        end
        osratio = round(obj.SampleRate*obj.ChipWidth);
        s = repmat(s,1,osratio);
        s = s.';
        s = s(:);
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
        wname = 'Phase-coded pulse waveform';
    end
   
end

end

% [EOF]
