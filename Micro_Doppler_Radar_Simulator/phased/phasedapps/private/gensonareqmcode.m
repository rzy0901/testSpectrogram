function gensonareqmcode(sparms)
% This function is for internal use only. It may be removed in the future.

% GENSONAREQMCODE Generate MATLAB code for Sonar Equation Calculator

% Copyright 2017 The MathWorks, Inc.

% Get date/time and version numbers
date_time_str = datestr(now);
ml = ver('matlab');
pat = ver('phased');

modeconfig = {'Active', ...
    'Passive'};

% mcodebuffer to holds the text for matlab code
mcode = sigcodegen.mcodebuffer;
mcode.addcr('% MATLAB Code');
mcode.addcr(['% Generated by ' ml.Name ' ' ml.Version ' and ' pat.Name ' ' pat.Version]);
mcode.addcr(['% Generated on ' date_time_str]);
mcode.addcr();
mcode.addcr(['% Mode Of Operation: ' modeconfig{sparms.mode_type}]);
mcode.addcr();
mcode.addcr(['% All quantities are in standard units' sprintf('\n')]);
mcode.addcr(['NoiseLevel = ' num2str(sparms.NoiseLevel) ';' sprintf('\t\t\t\t\t\t\t\t\t') '% Noise Level (dB//1uPa)']);
mcode.addcr(['DirectivityIndex = ' num2str(sparms.directivity_index) ';' sprintf('\t\t\t\t\t\t\t\t') '% Directivity Index (dB)']);

cond = (sparms.mode_type == 1);
if (cond)  % If Active write Target Strength
    mcode.addcr(['TargetStrength = ' num2str(sparms.target_strength) ';' sprintf('\t\t\t\t\t\t\t\t') '% Target Strength (dB)']);
end

if (sparms.calc_type == 1) % If Range write Freq and depth
    mcode.addcr(['Frequency = ' num2str(sparms.rngfreq) ';' sprintf('\t\t\t\t\t\t\t\t\t') '% Frequency (Hz)']);
    mcode.addcr(['Depth = ' num2str(sparms.rngdepth) ';' sprintf('\t\t\t\t\t\t\t\t\t\t') '% Depth(m)']);
end

if ~((sparms.calc_type == 1)||(sparms.calc_type ==2 ))
    if sparms.tldetection == 0
        mcode.addcr(['TransmissionLoss = ' num2str(sparms.TL) ';' sprintf('\t\t\t\t\t\t\t\t') '% Transmission Loss (dB)']);
    elseif sparms.tldetection == 1
        mcode.addcr(['TransmissionLoss = range2tl(' num2str(sparms.range) ',' ...
            num2str(sparms.freq) ',' num2str(sparms.depth)  ...
            ');' sprintf('\t\t') '% Transmission Loss(dB)']);
    end
end

if ~(sparms.calc_type == 3)
    mcode.addcr(['SourceLevel = ' num2str(sparms.SL) ';' ...
        sprintf('\t\t\t\t\t\t\t\t\t') '% Source Level(dB//1uPa)']);
end

if ~(sparms.calc_type == 4)
    if sparms.detection == 0
        mcode.addcr(['SNR = ' num2str(sparms.SNR) ';' sprintf('\t\t\t\t\t\t\t\t\t\t\t') '% SNR (dB)']);
    elseif sparms.detection == 1
        mcode.addcr(['SNR = shnidman(' num2str(sparms.Pd) ',' ...
            num2str(sparms.Pfa) ',' num2str(sparms.numpulse) ',' ...
            num2str(sparms.swerlingcase) ');' sprintf('\t\t\t\t\t') '% SNR (dB)']);
    end
end
mcode.addcr();

% Code for actual calculation
switch sparms.mode_type
    case 1
        calc_str = {'TargetRange = tl2range(TransmissionLoss,Frequency,Depth);', ...
            'TransmissionLoss = sonareqtl(SourceLevel, SNR, NoiseLevel,DirectivityIndex,TargetStrength );', ...
            'SourceLevel = sonareqsl(SNR, NoiseLevel, DirectivityIndex,TransmissionLoss, TargetStrength);', ...
            'SNR = sonareqsnr(SourceLevel, NoiseLevel, DirectivityIndex,TransmissionLoss, TargetStrength);'};
    case 2
        calc_str = {'TargetRange = tl2range(TransmissionLoss,Frequency,Depth);', ...
            'TransmissionLoss = sonareqtl(SourceLevel, SNR, NoiseLevel,DirectivityIndex);', ...
            'SourceLevel = sonareqsl(SNR, NoiseLevel, DirectivityIndex,TransmissionLoss);', ...
            'SNR = sonareqsnr(SourceLevel, NoiseLevel, DirectivityIndex,TransmissionLoss);'};
end

if sparms.calc_type == 1
    mcode.addcr(calc_str{[2,1]});
else
    mcode.addcr(calc_str{sparms.calc_type});
end

% Launch the editor with the MATLAB code
editorDoc = matlab.desktop.editor.newDocument;

editorDoc.Text = mcode.string;
editorDoc.smartIndentContents();
% Scroll document to line 1
editorDoc.goToLine(1);
end

% [EOF]