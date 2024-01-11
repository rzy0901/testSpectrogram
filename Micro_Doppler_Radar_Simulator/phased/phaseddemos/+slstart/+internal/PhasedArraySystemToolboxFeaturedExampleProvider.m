classdef PhasedArraySystemToolboxFeaturedExampleProvider < slstart.internal.FeaturedExampleProvider
    % Stub implementation of slstart.internal.FeaturedExampleProvider for
    % a product with no feature examples.

    % Copyright 2015 The MathWorks, Inc.

    properties (GetAccess = public, SetAccess = private)
        % The customer visible product name this example ships with
        Product = 'Phased Array System Toolbox';

        % The short name for this product as used by the Help Browser.
        ProductShortName = 'phased';

        % Names of featured examples in this product.
        FeaturedExamples = {'slexMonostaticRadarExampleExample',...
            'slexBistaticExampleExample',...
            'slexFMCWExampleExample',...
            'slexMicrophoneBeamformerExampleExample',...
            'slexBeamformerExampleExample',...
            'slexBeamscanMVDRDOAExampleExample',...
            'slexSTAPExampleExample'};
    end
end
