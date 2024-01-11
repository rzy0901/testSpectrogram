function gensonareqreport(sparms)
% This function is for internal use only. It may be removed in the future.

% genradareqreport Print a text report of parameters of Sonar Equation App

% Copyright 2017 The MathWorks, Inc.

% Get date/time and version numbers
date_time_str = datestr(now);
ml = ver('matlab');
pat = ver('phased');

% String dictionaries
modeconfig = {'Active', ...
    'Passive'};

calcstr = {'Range', ...
    'Transmission Loss',...
    'Source Level', ...
    'SNR'};

% Store the report content as StringWriter
rstring = StringWriter;
rstring.addcr('% Sonar Equation Calculator Report');
rstring.addcr('%');
rstring.addcr(['% Generated by ' ml.Name ' ' ml.Version ' and ' pat.Name ' ' pat.Version]);
rstring.addcr(['% Generated on ' date_time_str]);
rstring.addcr('%');
rstring.addcr(['% Calculation: ' calcstr{sparms.calc_type}]);
rstring.addcr(['% Sonar Mode: ' modeconfig{sparms.mode_type}]);
rstring.addcr('%');
rstring.addcr(['% Noise Level (dB//1uPa)................... ' num2str(sparms.NoiseLevel)]);
rstring.addcr(['% Directivity Index (dB)................... ' num2str(sparms.directivity_index)]);
cond = (sparms.mode_type == 1);
if (cond)
    rstring.addcr(['% Target Strength (dB)..................... ' num2str(sparms.target_strength)]);
end

%Range Calculation
if(sparms.calc_type == 1)
    rstring.addcr(['% Frequency (Hz)........................... ' num2str(sparms.rngfreq)]);
    rstring.addcr(['% Depth (m)................................ ' num2str(sparms.rngdepth)]);
end
% Write Transmission Loss if not to be calculated
if~((sparms.calc_type == 1)||(sparms.calc_type == 2))
    rstring.addcr(['% Transmission Loss (dB)................... ' num2str(sparms.TL)]);
end

% Write source level if not to be calculated
if ~(sparms.calc_type == 3)
    rstring.addcr(['% Source Level (dB//1uPa).................. ' num2str(sparms.SL)]);
end

% Write snr if not to be calculated
if ~(sparms.calc_type == 4)
    rstring.addcr(['% SNR (dB)................................. ' num2str(sparms.SNR)]);
end

result_str = {'%  Range (m)............................... ', ...
    '%  Transmission Loss (dB).................. ', ...
    '% Source Level (dB//1uPa).................. ', ...
    '% SNR (dB)................................. '};

% Add the radar equation string
rstring.addcr();
rstring.addcr([result_str{sparms.calc_type}  num2str(sparms.result,4)]);
% Open untitled file in editor and show the report text there
editorDoc = matlab.desktop.editor.newDocument;
editorDoc.Text = rstring.string;
% Scroll document to line 1
editorDoc.goToLine(1);
end % Gen report func
