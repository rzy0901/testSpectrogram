function w = lcmvweights(xArg,C,G,dl)
%This function is for internal use only. It may be removed in the future.

%LCMVWEIGHTS  LCMV weights 
%   W = phased.internal.lcmvweights(X,C,G,DL) calculates the linear
%   constraint minimum variance (LCMV) weights using sample matrix
%   inversion (SMI) method. X is an MxN data matrix whose rows are
%   snapshots in time. W is a length N column vector representing the LCMV
%   weights. C is an NxL matrix whose columns are constraints and L stands
%   for the number of constraints. Note that L cannot exceed N. G is a
%   length L column vector whose elements are the desired responses
%   corresponding to the constraints specified in C. DL is the diagonal
%   loading factor for the data covariance matrix. DL must be a
%   non-negative scalar.
%
%   The calculation in this function avoids forming the covariance matrix
%   of the data to reduce the requirement on the dynamic range. The
%   calculation is performed in the data domain.
%
%   Example:
%   % Calculate MVDR weights for a ULA toward 45 degrees azimuth.
%
%   % signal simulation
%   t = (0:1000)';
%   x = sin(2*pi*0.01*t);
%   c = 3e8; Fc = 3e8;
%   incidentAngle = [45; 0];
%   ha = phased.ULA('NumElements',5);
%   x = collectPlaneWave(ha,x,incidentAngle,Fc,c);
%   noise = 0.1*(randn(size(x)) + 1j*randn(size(x)));
%   rx = x+noise;
%
%   % weights calculation
%   hstv = phased.SteeringVector('SensorArray',ha,'PropagationSpeed',c);
%   vd = step(hstv,Fc,incidentAngle);
%   rd = 1;
%   w = phased.internal.lcmvweights(x,vd,rd,0.01);
%   pattern(ha,Fc,-180:180,0,'Weights',w,...
%       'PropagationSpeed',c,'Type','powerdb');

%   Copyright 2010-2011 The MathWorks, Inc.

%   Reference
%   [1] Guerci, Space-Time Adaptive Processing for Radar, 2003, pp173
%   [2] Van Trees, Optimum Array Processing, 2002

% This is an internal method. For performance purpose, we omit the
% validation here.

% calculate weights for adaptive methods, MVDR and LCMV [2]. Limit 
% operation on data matrix to improve performance [1].

%#codegen

[m,n] = size(xArg);
xArg = xArg/sqrt(m);

% add diagonal loading+
delta = sqrt(dl);
if delta ~= 0
    x = [xArg;delta*eye(n)];
else
    x = xArg;
end

% the covariance matrix is defined as E{x.'*conj(x)}
if size(C,2) == 1
    % MVDR
    temp = qrlinsolve(x.',C);
    w = G*temp/(C'*temp);
else
    % LCMV
    if m >= n
        [temp,F] = qrlinsolve(x.',C);         
        w = temp*qrlinsolve(F',G);
    else
        % when matrix is fat, F is no longer square and we cannot play the
        % trick of thin matrix. Therefore, we have to form R2 and use LU.
        temp = qrlinsolve(x.',C);
        R2 = C'*temp;   % R2 = C'*R^(-1)*C
        [L2, U2] = lu(R2);
        temp2 = U2\(L2\G);
        w = temp*temp2;
    end
end

end

function [Xout,Fout] = qrlinsolve(A,B)
%qrlinsolve Solve linear equations with (A*A')X=B using QR decomposition

%   Reference
%   [1] Guerci, Space-Time Adaptive Processing for Radar, 2003, pp173

[~,AR,PIVOTPERM] = qr(A',0);         % AR is upper triangular, see [1]
% Forward substitution after backward substitution 
F = AR'\B(PIVOTPERM,:);
X = AR\F;
%X(PIVOTPERM,:) = X;
[~,PIVOTPERMSRT] = sort(PIVOTPERM);
Xout = X(PIVOTPERMSRT,:);
if nargout > 1
    Fout = F(PIVOTPERMSRT,:);
else
    Fout = F;
end

end


% [EOF]
