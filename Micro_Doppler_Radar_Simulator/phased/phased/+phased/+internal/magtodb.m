function out = magtodb(in)
%This function is used temporarily until constant timeouts are fixed

%   Copyright 2011-2012 The MathWorks, Inc.

%#codegen

coder.extrinsic('mag2db');
out = coder.internal.const(mag2db(in));
