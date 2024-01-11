function [pp,dpp] = motionhermite(t,x,dx)
%This function is for internal use only. It may be removed in the future.

%MOTIONHERMITE Compute Hermite polynomial for motion
%   [PP,DPP] = motionhermite(T,X,DX) computes the Hermite polynomial PP and
%   DPP based on the position, X, and velocity, DX, obtained at time
%   instances specified in T. PP is the polynomial for X and DPP is the
%   polynomial for DX. To obtain the result in other time instances, Tq,
%   use
%
%   Xq = ppval(PP,Tq) and DXq = ppval(DPP,Tq)
%
%   The polynomials ensures that at known time instances T, Xq and DXq
%   match X and DX, respectively. Note X and DX are one dimensional. To
%   have it work in 3D, define a polynomial for each coordinate.
%
%   % Example:
%   %   Interpolate path based on given coordinates.
%   
%   t = [0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4 2*pi];
%   y = cos(t);
%   yd = -sin(t);
%   [pp,dpp] = phased.internal.motionhermite(t,y,yd);
%   tq = 0:pi/10:2*pi;
%   yq = ppval(pp,tq);
%   ydq = ppval(dpp,tq);
%   xlabel('Time (s)');
%   plot(t,y,'b-',tq,yq,'b--',t,yd,'r-',tq,ydq,'r--');
%   legend('X','Xq','DX','DXq')
%
%   See also phased, phased.Platform.

%   Copyright 2018 The MathWorks, Inc.

%#ok<*EMCA>
%#codegen

narginchk(2,3);

if nargin < 3  % pos
    % obtain piecewise cubic Hermite
    dx = zeros(size(x)); % real
    dx(1) = (x(2)-x(1))./(t(2)-t(1));
    dx(end) = (x(end)-x(end-1))./(t(end)-t(end-1));
    dx(2:end-1) = (x(3:end)-x(1:end-2))./(t(3:end)-t(1:end-2));
    pp = pwch(t(:).',x(:).',dx(:).');
    
    % do cubic Hermite for velocity, using piecewise slope
    dpp = pchip(t(:).',dx(:).');
else % pos & vel
    % obtain piecewise cubic Hermite
    pp = pwch(t(:).',x(:).',dx(:).');
    
    % do cubic Hermite for velocity
    dpp = pchip(t(:).',dx(:).');
end

% % extract the coefficients
% [breaks,coefs,npieces,order,dim] = unmkpp(pp);
% 
% % take the derivative of each polynomial
% dpp = mkpp(breaks,repmat(order-1:-1:1,dim*npieces,1).*coefs(:,1:order-1),dim);
