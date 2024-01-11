function idx = decomposeCube(numax1, numax2, numax3)
%This function is for internal use only. It may be removed in the future.

% idx = DECOMPOSECUBE(numax1, numax2, numax3) Takes the number of elements
% in the three axes of a cube (for example rows, columns and pages) and
% returns a cell array of indices idx that divide the cube into nearly
% equal chunks along the third axis -numax3. The maximum chunk size (in
% terms of elements) are controlled by a constant below.


%Constants
CHUNKMAXSIZE = 10e6;

sliceSize = numax1 * numax2;


slicesPerChunk = max( floor(CHUNKMAXSIZE /sliceSize), 1);

minChunks = ceil(numax3/slicesPerChunk);

chunkSizes = zeros(1,minChunks);
chunkSizes(:) = floor(numax3/minChunks);
chunkSizes(1:rem(numax3,minChunks)) = ceil(numax3/minChunks);

%Build the start and end indexes using prefix sums
endIdx = cumsum(chunkSizes);
startIdx = endIdx + 1 - chunkSizes;

%Build the cell array of indices using the colon function
idx = arrayfun(@colon, startIdx, endIdx, 'UniformOutput', false);