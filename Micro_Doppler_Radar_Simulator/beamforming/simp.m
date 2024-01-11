% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     simp.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       September 2015 created
%
%  *************************************************************************************
%    Description:
% 
%    Performs integration using Simpson method
%
%    out = pattern(el, az, W)
%
%    Inputs:
%
%       1. a, b - interval for integration
%       2. el   - elevation angle, [rad]
%       3. W    - weight coefficients
%       4. e    - precision
%
%    Outputs:
%
%       1. out  - array of output amplitudes weighted by antenna gain coefficients
%
%  *************************************************************************************/
function y = simp(a, b, el, W, e, paa)
% Number of points between [a, b]
N = 100;
Int2 = -3; % integral value for N/2 points
Int1 = 0; % integral value for N points
% Check precision
while (abs(Int1-Int2)>e)
    % Update values
    Int2 = Int1;

    % Increase precision
    N = 2*N;
    % Recalculate step
    h = (b-a)/(N-1);
    % Initial value for integral in Simpson method
    Int1 = pattern(el, a, W, paa) * sin(el) + pattern(el, b, W, paa) * sin(el);
    for i=1:N-2
        if ( mod(i,2)~=0 )
            Int1 = Int1 + 4 * ( pattern(el, a + h*i, W, paa) * sin(el) );
        else Int1 = Int1 + 2 * ( pattern(el, a + h*i, W, paa) * sin(el) );
        end
    end
    Int1 = Int1*h/3;
end
y = Int1;
end