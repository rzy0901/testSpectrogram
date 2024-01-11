function [A,B,C,D] = getKinematicModel(modelidx,n,t)
%This function is for internal use only. It may be removed in the future.

%getKinematicModel  State space model for constant velocity and
%acceleration kinematic models
%   [A,B,C,D] = phased.internal.getKinematicModel(IDX,N,T) returns the
%   state space model matrices, A, B, C, and D, for either constant
%   velocity or constant acceleration kinematic models. If IDX is 1, the
%   model is for constant velocity. If IDX is 2, the model is for constant
%   acceleration. N is an integer between 1 and 3 indicating the dimensions
%   for kinematic models and T is the discrete sample time of the model.
%
%   For a one dimensional constant velocity kinematic model, A = [1 T;0 1];
%   B = [T^2/2;T]; C = [1 0]; D = 0.
%
%   For a one dimensional constant acceleration kinematic model, 
%   A = [1 T T^2/2;0 1 T;0 0 1]; B = [T^2/2;T;1]; C = [1 0 0]; D = 0.

%   Copyright 2015 The MathWorks, Inc.

%   Reference
%   [1] Blackman, Multiple-Target Tracking with Radar Applications, Artech
%   House, 1986
%   [2] Bar-Shalom, et al. Estimation With Applications to Tracking and
%   Navigation, John Wiley & Sons, 2001

%#codegen

if modelidx == 1  % ConstantVelocity, Discrete white noise acceleration (DWNA) [2] 6.3.2
    A1d = [1 t;0 1];
    B1d = [t^2/2;t];
    C1d = [1 0];
else % ConstantAcceleration, Discrete Wiener process acceleration (DWPA), [2], 6.3.3
    A1d = [1 t t^2/2;0 1 t;0 0 1];
    B1d = [t^2/2;t;1];
    C1d = [1 0 0];
end
switch n
    case 1
        A = A1d;
        B = B1d;
        C = C1d;
    case 2
        A = blkdiag(A1d,A1d);
        B = blkdiag(B1d,B1d);
        C = blkdiag(C1d,C1d);
    otherwise % case 3
        A = blkdiag(A1d,A1d,A1d);
        B = blkdiag(B1d,B1d,B1d);
        C = blkdiag(C1d,C1d,C1d);
end
D = 0;

