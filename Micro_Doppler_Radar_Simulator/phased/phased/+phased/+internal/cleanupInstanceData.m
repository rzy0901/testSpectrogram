function outData = cleanupInstanceData(inData)
%CLEANUPINSTANCEDATA Cleanup obsolete properties from the instance data
%   This function is a transformation function for phased specified in the
%   block forwarding table of phasedlib

%   Copyright 2017 The MathWorks, Inc.

outData.NewBlockPath = ''; % No change in the library block path
outData.NewInstanceData = [];

if ~isempty(inData) && isstruct(inData) && isfield(inData,'InstanceData') && ...
        isfield(inData.InstanceData,'Name')
    instanceData = inData.InstanceData;
    names = {instanceData.Name};
    for pndx = 1:length(names)
        if any(strcmp(names{pndx},{'Sensor','SensorArray'}))
            sensorExpr = instanceData(pndx).Value;
            sensorExpr = updateCustomAntennaRadiationPattern(sensorExpr);
            instanceData(pndx).Value = sensorExpr;
        end
    end
    outData.NewInstanceData = instanceData;
else
    outData = inData;
end

end

function custStr = updateCustomAntennaRadiationPattern(custStr)
RadiationPatternIdx = strfind(custStr,'RadiationPattern');
if ~isempty(RadiationPatternIdx)
    CommaIdx = strfind(custStr,',');
    LeftBracketIdx = strfind(custStr,'(');
    RightBracketIdx = strfind(custStr,')');
    RadiationPatternValueStartIdx = CommaIdx(find(CommaIdx>RadiationPatternIdx,1))+1;
    IdxPool = sort([CommaIdx,LeftBracketIdx,RightBracketIdx]);
    bcount = 0;
    for m = 1:numel(IdxPool)
        % There are only two possibilities. Either the value
        % ends with a comma, or the value with an right bracket
        % which matches the bracket of the entire constructor.
        % So we'll search through the string starting from the
        % starting point and determine when to stop. If we
        % encounter a left bracket, bcount adds 1, if we
        % encounter a right bracket, bcount subtracts 1. For
        % the first case, we'll get a comma when bcount is 0.
        % For the second case, we'll get a bcount of -1.
        current_idx = IdxPool(m);
        if current_idx > RadiationPatternValueStartIdx
            if custStr(current_idx) == ','
                if bcount == 0
                    RadiationPatternValueEndIdx = current_idx-1;
                    break;
                end
            elseif custStr(current_idx) == '('
                bcount = bcount+1;
            else % custStr(current_idx) == ')'
                bcount = bcount-1;
                if bcount == -1
                    RadiationPatternValueEndIdx = current_idx-1;
                    break;
                end
            end
        end
    end
    RadiationPatternValueStr = custStr(RadiationPatternValueStartIdx:RadiationPatternValueEndIdx);
    custStr = [custStr(1:RadiationPatternValueEndIdx) ',''PhasePattern'',zeros(size(' RadiationPatternValueStr '))' custStr(RadiationPatternValueEndIdx+1:end)];
    custStr = replace(custStr,"RadiationPattern","MagnitudePattern");
end
end