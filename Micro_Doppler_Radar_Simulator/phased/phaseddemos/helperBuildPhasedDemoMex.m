function helperBuildPhasedDemoMex(inputFileName)
% This function helperBuildPhasedDemoMex is only in support of
% RadarStreamExample. It may be removed in a future release.

tempDir = tempname;
if ~exist(tempDir,'dir')
    mkdir(tempDir);
end

currentDir = cd(tempDir);
codegen(inputFileName)
dotidx = strfind(inputFileName,'.');
if isempty(dotidx)
    dotidx = length(inputFileName)+1;
end
mexFileName = sprintf('%s_mex.*',inputFileName(1:dotidx-1));
copyfile(mexFileName,currentDir);
cd(currentDir);

