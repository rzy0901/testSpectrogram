% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     max_power_ray_conf3.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       December 2015 created
%
%  *************************************************************************************
%    Description:
%
%    Executes max power ray alg for configuration #3
%
%    [imp_res] = max_power_ray_conf3(cfg, ch)
%
%    Inputs:
%
%       1. cfg      - part of configuration structure defining beamforming related parameters
%       2. ch       - channel structure
%
%    Outputs:
%
%       1. imp_res - channel impulse response structure
%
%  *************************************************************************************/
function [imp_res] = max_power_ray_conf3(cfg, ch)

% conversion of angles to new TX/RX coordinates
% (x-axis along LOS, z-axis normal to horizontal plane, y-axis is added for right tern)
% azimuth angle thi: [0:360], elevation angle theta: [0:180]
shift_tx11 = ch.tx_az_11(end);
shift_rx11 = ch.rx_az_11(end);
shift_tx22 = ch.tx_az_22(end);
shift_rx22 = ch.rx_az_22(end);

ch.tx_az_11(end) = [];
ch.rx_az_11(end) = [];
ch.tx_az_22(end) = [];
ch.rx_az_22(end) = [];

ch.tx_az_11 = mod(360 + ch.tx_az_11, 360);
ch.rx_az_11 = mod(360 + ch.rx_az_11, 360);
ch.tx_el_11 = 90 - ch.tx_el_11;
ch.rx_el_11 = 90 - ch.rx_el_11;

ch.tx_az_12 = mod(360 + ch.tx_az_12, 360);
ch.rx_az_12 = mod(360 + ch.rx_az_12, 360);
ch.tx_el_12 = 90 - ch.tx_el_12;
ch.rx_el_12 = 90 - ch.rx_el_12;

ch.tx_az_21 = mod(360 + ch.tx_az_21, 360);
ch.rx_az_21 = mod(360 + ch.rx_az_21, 360);
ch.tx_el_21 = 90 - ch.tx_el_21;
ch.rx_el_21 = 90 - ch.rx_el_21;

ch.tx_az_22 = mod(360 + ch.tx_az_22, 360);
ch.rx_az_22 = mod(360 + ch.rx_az_22, 360);
ch.tx_el_22 = 90 - ch.tx_el_22;
ch.rx_el_22 = 90 - ch.rx_el_22;

[ch] = filter_channel(cfg.ant_type, ch, cfg.paa.phi_tx, cfg.paa.phi_rx, shift_tx11, shift_rx11, shift_tx22, shift_rx22);

% find the ray with maximum power for V -> V mode
[power_max_11, idx_max_11] = max(abs(ch.am_11).^2);
tx_az_max1 = ch.tx_az_11(idx_max_11);
tx_el_max1 = ch.tx_el_11(idx_max_11);
rx_az_max1 = ch.rx_az_11(idx_max_11);
rx_el_max1 = ch.rx_el_11(idx_max_11);

% find the ray with maximum power for H -> H mode
[power_max_22, idx_max_22] = max(abs(ch.am_22).^2);
tx_az_max2 = ch.tx_az_22(idx_max_22);
tx_el_max2 = ch.tx_el_22(idx_max_22);
rx_az_max2 = ch.rx_az_22(idx_max_22);
rx_el_max2 = ch.rx_el_22(idx_max_22);

% Processing of V -> V channel
% TX antenna gain
am_tx_11 = ant_gain(ch.am_11, ch.tx_az_11, ch.tx_el_11, tx_az_max1, tx_el_max1, cfg.paa);

% RX antenna gain
am_rx_11 = ant_gain(am_tx_11, ch.rx_az_11, ch.rx_el_11, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h11 = am_rx_11;

% Processing of H -> H channel
% TX antenna gain
am_tx_22 = ant_gain(ch.am_22, ch.tx_az_22, ch.tx_el_22, tx_az_max2, tx_el_max2, cfg.paa);

% RX antenna gain
am_rx_22 = ant_gain(am_tx_22, ch.rx_az_22, ch.rx_el_22, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h22 = am_rx_22;

% Processing of V -> H channel
% TX antenna gain
am_tx_12 = ant_gain(ch.am_12, ch.tx_az_12, ch.tx_el_12, tx_az_max1, tx_el_max1, cfg.paa);

% RX antenna gain
am_rx_12 = ant_gain(am_tx_12, ch.rx_az_12, ch.rx_el_12, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h12 = am_rx_12;

% Processing of H -> V channel
% TX antenna gain
am_tx_21 = ant_gain(ch.am_21, ch.tx_az_21, ch.tx_el_21, tx_az_max2, tx_el_max2, cfg.paa);

% RX antenna gain
am_rx_21 = ant_gain(am_tx_21, ch.rx_az_21, ch.rx_el_21, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h21 = am_rx_21;

% sort impulse response values
[toa_11, ind_sort] = sort(ch.toa_11);
imp_res.h11 = imp_res_h11(ind_sort);
imp_res.toa_11 = toa_11;

[toa_12, ind_sort] = sort(ch.toa_12);
imp_res.h12 = imp_res_h12(ind_sort);
imp_res.toa_12 = toa_12;

[toa_21, ind_sort] = sort(ch.toa_21);
imp_res.h21 = imp_res_h21(ind_sort);
imp_res.toa_21 = toa_21;

[toa_22, ind_sort] = sort(ch.toa_22);
imp_res.h22 = imp_res_h22(ind_sort);
imp_res.toa_22 = toa_22;
