function K = cel2kel(C)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2015 The MathWorks, Inc.

%cel2kel    Convert Celsius to Kelvin 

%#codegen
%#ok<*EMCA>

K = C+273.15;