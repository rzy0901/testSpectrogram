% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_ref_loss_conf5.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       November 2015 created
%
%  *************************************************************************************
%    Description:
% 
%    function returns clusters reflection loss according to normal distribution based
%    on the experimental data
%
%    [ref_loss11, ref_loss12] = cr_ref_loss_conf5()
%
%    Outputs:
%
%       1. ref_loss11 - reflection loss coefficients array for pol(1) -> pol(1) polarization
%       2. ref_loss12 - reflection loss coefficients array for pol(1) -> pol(2) polarization
%
%    Inputs:
%
%       1. ap_sp    - parameter selects subscenario: 0 - STA-STA, 1 - STA-AP
%       2. pol      - antennas polarization type for PAA on Tx side
%
%    Row dimension in ref_loss array for ap_sp = 0 :
%
%    1                       - 1st order ceiling cluster coefficient
%    2,3,4,5                 - 1st order walls clusters  coefficients
%
%    6,7,8,9                 - 2nd order wall-ceiling (ceiling-wall) clusters  coefficients
%    10,11,12,13,14,15,16,17 - 2nd order walls clusters coefficients
%
%    Row dimension in ref_loss array for ap_sp = 1 :
%
%    1,2,3,4                 - 1st order walls clusters coefficients
%
%    5,6,7,8,9,10,11,12      - 2nd order walls clusters coefficients
%
%  *************************************************************************************/
function [ ref_loss11, ref_loss12 ] = cr_ref_loss_conf5( ap_sp, pol )

% polarization vectors for TX and RX antennas
tx_pol_1 = polarization(pol(1));
                
rx_pol_1 = polarization(pol(1));
rx_pol_2 = polarization(double(~pol(1)));

% choose subscenario
switch (ap_sp)
   case 0, % subscenario STA-STA : 17 clusters = 5 (1st order) + 12 (2nd order)
                
   % 1st order ceiling cluster
   % reflection matrix
        R(1,1) = ref_coef(0,1);
        R(2,2) = ref_coef(0,0);
        R(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
        R(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
           
        ref_loss11(1,:) = rx_pol_1'*R*tx_pol_1;
        ref_loss12(1,:) = rx_pol_2'*R*tx_pol_1;
                
        % 1st order walls clusters
        for i=2:5
            % reflection matrix
            R(1,1) = ref_coef(0,0);
            R(2,2) = ref_coef(0,1);
            R(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            R(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            
            ref_loss11(i,:) = rx_pol_1'*R*tx_pol_1;
            ref_loss12(i,:) = rx_pol_2'*R*tx_pol_1;
        end
                
        % 2nd order wall-ceiling (ceiling-wall) clusters
        for i=6:9
            % reflection matrix
            R(1,1) = ref_coef(1,0);
            R(2,2) = ref_coef(1,0);
            R(1,2) = 0.2.*rand(1,1) - 0.1;
            R(2,1) = 0.2.*rand(1,1) - 0.1;
            
            ref_loss11(i,:) = rx_pol_1'*R*tx_pol_1;
            ref_loss12(i,:) = rx_pol_2'*R*tx_pol_1;
        end
                
        % 2nd order walls clusters
        for i=10:17
            % reflection matrix
            R(1,1) = ref_coef(1,1);
            R(2,2) = ref_coef(1,0.87);
            R(1,2) = 0.2.*rand(1,1) - 0.1;
            R(2,1) = 0.2.*rand(1,1) - 0.1;
            
            ref_loss11(i,:) = rx_pol_1'*R*tx_pol_1;
            ref_loss12(i,:) = rx_pol_2'*R*tx_pol_1;
        end
    case 1, % subscenario STA-AP : 12 clusters = 4 (1st order) + 8 (2nd order)
        
        % 1st order ceiling clusters
        for i=1:4
            % reflection matrix
            R(1,1) = ref_coef(0,0);
            R(2,2) = ref_coef(0,1);
            R(1,2) = 0.4.*rand(1,1) - 0.2;
            R(2,1) = 0.4.*rand(1,1) - 0.2;
            
            ref_loss11(i,:) = rx_pol_1'*R*tx_pol_1;
            ref_loss12(i,:) = rx_pol_2'*R*tx_pol_1;
        end
                
        % 2nd order walls clusters
        for i=5:12
            % reflection matrix
            R(1,1) = ref_coef(1,1);
            R(2,2) = ref_coef(1,0.73);
            R(1,2) = 0.3.*rand(1,1) - 0.15;
            R(2,1) = 0.3.*rand(1,1) - 0.15;
           
            ref_loss11(i,:) = rx_pol_1'*R*tx_pol_1;
            ref_loss12(i,:) = rx_pol_2'*R*tx_pol_1;
        end        
            otherwise,
                error('Prohibited value of "ap_sp" parameter');
end

% generate reflection coefficient in dB
function y = ref_coef_db(type,size)

switch (type)
    case 0, % 1st order
        mean_value = -10; % dB
        sigma   =   4; % dB
    case 1,% 2nd order
        mean_value = -16; % dB
        sigma   =   5; % dB
end

y = mean_value + randn(size,1)*sigma;

idx = find(y > -2);

while ~isempty(idx)
    
    y(idx) = mean_value + randn(length(idx),1)*sigma;
    idx = find(y > -2);
end


% generate reflection coefficient
function y = ref_coef(type,p)

switch (type)
    case 0, % 1st order reflections
        mean_value = -10; % db
        sigma = 4; % db
    case 1, % 2nd order reflections
        mean_value = -16; % db
        sigma = 5; % db        
end

y_db = mean_value + randn(1,1).*sigma;

idx = find(y_db > -2);

while ~isempty(idx)
    
    y_db(idx) = mean_value + randn(1,1).*sigma;
    idx = find(y_db > -2);
end

y = 10.^(y_db./20);

prob = rand(1,1);

index = find(prob>p);

y(index) = -y(index);
