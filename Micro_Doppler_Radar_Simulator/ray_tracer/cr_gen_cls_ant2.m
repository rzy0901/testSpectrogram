% /*************************************************************************************
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_gen_cls_ant2.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       April 2016 created
%  *************************************************************************************
%    Description:
% 
%    generates NLOS clusters for SU-MIMO configurations with double PAA on both sides
%    (conf.#3, 4)
% 
%  *************************************************************************************/
function [cls_11, cls_12, cls_21, cls_22] = cr_gen_cls_ant2(scenario, dist_Tx, dist_Rx, phi_Tx, phi_Rx)
    
global panels_list;

% Load scene for conference room
[panels_list] = cr_build_scene();

% Tx & Rx positions:
switch ( scenario )
    case 0, % STA-STA
        panels_list(5) = []; % remove reflections from table
        tx_z = 100;
        tx_x = 1 + dist_Tx + rand(1,1).*(100-2*dist_Tx);
        tx_y = 150 + dist_Tx + rand(1,1).*(250-2*dist_Tx);
        TX1 = [tx_x;tx_y;tx_z];
        
        TX2 = TX1;
        TX2(1) = TX2(1) + dist_Tx*cos(phi_Tx*pi/180);
        TX2(2) = TX2(2) + dist_Tx*sin(phi_Tx*pi/180);
        
        rx_z = 100;
        rx_x = 1 + dist_Rx + rand(1,1).*(100 - 2*dist_Rx);
        rx_y = 150 + dist_Rx + rand(1,1).*(250 - 2*dist_Rx);        
        RX1 = [rx_x;rx_y;rx_z];
        
        RX2 = RX1;
        RX2(1) = RX2(1) + dist_Rx*cos(phi_Rx*pi/180);
        RX2(2) = RX2(2) + dist_Rx*sin(phi_Rx*pi/180);

    case 1, % STA-AP
        panels_list(1) = [];
        panels_list(4) = [];
        tx_z = 290;
        tx_x = 150;
        tx_y = 50;
        TX1 = [tx_x;tx_y;tx_z];
        
        TX2 = TX1;
        TX2(1) = TX2(1) + dist_Tx;
        
        rx_z = 100;
        rx_x = 100 + dist_Rx + rand(1,1).*(100 - 2*dist_Rx);
        rx_y = 100 + dist_Rx + rand(1,1).*(250 - 2*dist_Rx);
        RX1 = [rx_x;rx_y;rx_z];
        
        RX2 = RX1;
        RX2(1) = RX2(1) + dist_Rx*cos(phi_Rx*pi/180);
        RX2(2) = RX2(2) + dist_Rx*sin(phi_Rx*pi/180);
end

toa_los11 = norm(RX1 - TX1)/ 3e+1;
toa_los12 = norm(RX2 - TX1)/ 3e+1;
toa_los21 = norm(RX1 - TX2)/ 3e+1;
toa_los22 = norm(RX2 - TX2)/ 3e+1;
    
toa_los = min(min(min(toa_los11, toa_los12), toa_los21), toa_los22);

% Run ray tracer between 1st PAA (TX) and 1st PAA (RX)
[shift_tx11, shift_rx11, str_id11, tx_az11, tx_el11, rx_az11, rx_el11, dist11, str_id12, tx_az12, tx_el12, rx_az12, rx_el12, dist12] = ray_tracer_ant2(TX1, RX1, RX2);
cls_11 = cr_reorder_clusters(scenario, str_id11, tx_az11, tx_el11, rx_az11, rx_el11, dist11, toa_los);
cls_12 = cr_reorder_clusters(scenario, str_id12, tx_az12, tx_el12, rx_az12, rx_el12, dist12, toa_los);

[shift_tx22, shift_rx22, str_id22, tx_az22, tx_el22, rx_az22, rx_el22, dist22, str_id21, tx_az21, tx_el21, rx_az21, rx_el21, dist21] = ray_tracer_ant2(TX2, RX2, RX1);
cls_22 = cr_reorder_clusters(scenario, str_id22, tx_az22, tx_el22, rx_az22, rx_el22, dist22, toa_los);
cls_21 = cr_reorder_clusters(scenario, str_id21, tx_az21, tx_el21, rx_az21, rx_el21, dist21, toa_los);

size = length(cls_11.tx_az) + 1;
cls_11.tx_az(size) = shift_tx11;
cls_11.rx_az(size) = shift_rx11;

cls_22.tx_az(size) = shift_tx22;
cls_22.rx_az(size) = shift_rx22;

end

