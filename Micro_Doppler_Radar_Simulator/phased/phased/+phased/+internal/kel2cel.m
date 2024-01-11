function C = kel2cel(K)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2015 The MathWorks, Inc.

%kel2cel    Convert Kelvin to Celsius

%#codegen
%#ok<*EMCA>

C = K-273.15;