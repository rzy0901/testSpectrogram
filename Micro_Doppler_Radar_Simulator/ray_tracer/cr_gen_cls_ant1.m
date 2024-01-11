% /*************************************************************************************
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_gen_cls_ant1.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       April 2016 created
%  *************************************************************************************
%    Description:
% 
%    generates NLOS clusters for SU-MIMO configurations with single PAA on both sides
%    (conf.#1, 2, 5)
% 
%  *************************************************************************************/
function [cls_11] = cr_gen_cls_ant1(scenario,Radar_x,Radar_y,Radar_z)

% Room dimensions:
    
global panels_list;

% Load scene with conference room
[panels_list] = cr_build_scene();

% Tx & Rx positions:
switch ( scenario )
    case 0, % STA-STA
        panels_list(5) = []; % remove reflections from table
        tx_z = Radar_z*100;
        tx_x = Radar_x*100 + rand(1,1).*0;
        tx_y = Radar_y*100 + rand(1,1).*0;
        
        TX1 = [tx_x;tx_y;tx_z];
        
        rx_z = Radar_z*100;
        rx_x = Radar_x*100 + rand(1,1).*10;
        rx_y = Radar_y*100 + rand(1,1).*10;
        
        RX1 = [rx_x;rx_y;rx_z];
    case 1, % STA-AP
        panels_list(1) = []; % remove reflections from ceiling
        panels_list(4) = []; % remove reflections from table
        
        % AP position is fixed
        tx_z = 290;
        tx_x = 150;
        tx_y = 50;
        
        TX1 = [tx_x;tx_y;tx_z];
        
        rx_z = 100;
        rx_x = 100 + rand(1,1).*100;
        rx_y = 100 + rand(1,1).*250;
        
        RX1 = [rx_x;rx_y;rx_z];
end

toa_los = norm(RX1 - TX1)/ 3e+1;

% Run ray tracer between TX and RX
[str_id, tx_az, tx_el, rx_az, rx_el, dist] = ray_tracer_ant1(TX1, RX1);

cls_11 = cr_reorder_clusters(scenario, str_id, tx_az, tx_el, rx_az, rx_el, dist, toa_los);


end
