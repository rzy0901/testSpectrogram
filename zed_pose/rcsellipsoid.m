function rcs = rcsellipsoid(a,b,c,phi_i,theta_i,phi_s,theta_s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates the bistatic (or monostatic) RCS of an ellipsoid
% with semi-axis lengths of a, b, and c.
% phi_i, theta_i: the aspect and azimuth angles of incident waves.
% phi_s, theta_s: the aspect and azimuth angles of scattered waves.
% Seeing details for ellipsoid rcs settings, please refer to 
% https://rzy0901.github.io/post/rcs/#bistatic-rcs-estimation-of-an-ellipsoid5-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~=5 && nargin ~= 7
    nargin
    error("Please use correct input parameters.")
end
if nargin == 5 % monostatic rcs
    rcs = (pi*a^2*b^2*c^2)/(a^2*(sin(theta_i))^2*(cos(phi_i))^2+b^2*(sin(theta_i))^2*(sin(phi_i))^2+c^2*(cos(theta_i))^2)^2;
end
if nargin == 7 % bistatic rcs
    rcs = (4*pi*a^2*b^2*c^2)*((1+cos(theta_i)*cos(theta_s))*cos(phi_s-phi_i)+sin(theta_i)*sin(theta_s))^2/ ...
        (a^2*(sin(theta_i)*cos(phi_i)+sin(theta_s)*cos(phi_s))^2+ ...
        b^2*(sin(theta_i)*sin(phi_i)+sin(theta_s)*sin(phi_s))^2+ ...
        c^2*(cos(theta_i)+cos(theta_s))^2 ...
        )^2;
end