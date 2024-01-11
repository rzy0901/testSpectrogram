% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_atten_coef.m
%    Authors:       A. Lomayev, R. Maslennikov
%    Version:       1.0
%    History:       May 2010 created
%
%  *************************************************************************************
%    Description:
% 
%    function generates attenuation coefficients for cluster due to
%    human blockage
%
%    [atten_coef] = cr_atten_coef(ap_sp)
%
%    Inputs:
%
%       1. ap_sp - parameter selects subscenario: 0 - STA-STA, 1 - STA-AP
%
%    Outputs:
%
%       1. atten_coef - attenuation coefficients array in linear scale
% 
%  *************************************************************************************/
function [atten_coef] = cr_atten_coef(ap_sp)

switch(ap_sp)
    case 0, % STA - STA
        
        % 1st order ceiling cluster
        atten_coef_db(1) = 0;
        % 1st order wall clusters
        atten_coef_db(2) = atten_fun(0);
        atten_coef_db(3) = atten_fun(0);
        atten_coef_db(4) = atten_fun(0);
        atten_coef_db(5) = atten_fun(0);
        % 2nd order wall-ceiling (ceiling-wall) clusters
        atten_coef_db(6) = atten_fun(1);
        atten_coef_db(7) = atten_fun(1);
        atten_coef_db(8) = atten_fun(1);
        atten_coef_db(9) = atten_fun(1);
        % 2nd order wall clusters
        atten_coef_db(10) = atten_fun(0);
        atten_coef_db(11) = atten_fun(0);
        atten_coef_db(12) = atten_fun(0);
        atten_coef_db(13) = atten_fun(0);
        atten_coef_db(14) = atten_fun(0);
        atten_coef_db(15) = atten_fun(0);
        atten_coef_db(16) = atten_fun(0);
        atten_coef_db(17) = atten_fun(0);
        
    case 1, % STA-AP
        
        % 1st order clusters
        atten_coef_db(1) = atten_fun(0);
        atten_coef_db(2) = atten_fun(0);
        atten_coef_db(3) = atten_fun(0);
        atten_coef_db(4) = atten_fun(0);
        % 2nd order clusters
        atten_coef_db(5) = atten_fun(0);
        atten_coef_db(6) = atten_fun(0);
        atten_coef_db(7) = atten_fun(0);
        atten_coef_db(8) = atten_fun(0);
        atten_coef_db(9) = atten_fun(0);
        atten_coef_db(10) = atten_fun(0);
        atten_coef_db(11) = atten_fun(0);
        atten_coef_db(12) = atten_fun(0);        
end

atten_coef = 10.^(atten_coef_db./20);        
        

function [A] = atten_fun(type)

switch(type)
    case 0, % 1st/2nd order reflections from walls
        
        A = 1;
        while(A>0) % truncation level is equal to 0 dB
            
            % GMM (Gaussian mixture model) in log-scale
            p = rand(1,1);
            
            if (p>0.83)
                mv = -47.2; % mean value
                std = 10;   % standard deviation
                A = std.*randn(1,1) + mv;
            else
                mv = -18.2; % mean value
                std = 8.3;  % standard deviation
                A = std.*randn(1,1) + mv;
            end
        end
                
    case 1, % 2nd order reflections from wall and then ceiling (or from ceiling and then wall)
        
        A = 1;
        while(A>0) % truncation level is equal to 0 dB
            
            % Gaussian distribution in log-scale
            mv = -18.4; % mean value
            std = 8.8;  % standard deviation            
            A = std.*randn(1,1) + mv;
        end
end
