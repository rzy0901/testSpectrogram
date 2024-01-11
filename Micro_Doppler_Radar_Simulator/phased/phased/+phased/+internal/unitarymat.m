function  Q = unitarymat(M)
%This function is for internal use only. It may be removed in the future.

%unitarymat  Unitary transformation matrix Q
%   Q = unitarymat(Hdoa,M) returns a unitary transformation matrix Q
%   that transforms a signal covariance matrix Sx with dimension M x M to a
%   Real Symmetric Matrix Sq.

%   Copyright 2009 The MathWorks, Inc.

%#codegen
    
    half = floor(M/2);
    I = eye(half);
    J = fliplr(I);

    if rem(M,2),
        % M odd
        O = zeros(half,1);
        Qmid = [O.' sqrt(2) O.'];
    else
        % M even
        O = []; Qmid = [];
    end
    Q = 1/sqrt(2)*[I O 1i*I;Qmid;J O -1i*J]; % Eq (7.58) and (7.59) in [1]

end
