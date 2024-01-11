% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     max_power_ray_conf2.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       December 2015 created
%
%  *************************************************************************************
%    Description:
%
%    Executes max power ray alg for configuration #2
%
%    [imp_res] = max_power_ray_conf2(cfg, ch)
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
function [imp_res] = max_power_ray_conf2(cfg, ch)

% conversion of angles to new TX/RX coordinates
% (x-axis along LOS, z-axis normal to horizontal plane, y-axis is added for right tern)
% azimuth angle thi: [0:360], elevation angle theta: [0:180]
ch.tx_az = mod(360 + ch.tx_az,360);
ch.rx_az = mod(360 + ch.rx_az,360);
ch.tx_el = 90 - ch.tx_el;
ch.rx_el = 90 - ch.rx_el;

[ch] = filter_channel(cfg.ant_type, ch, cfg.paa.phi_tx, cfg.paa.phi_rx);

% find the ray with maximum power for V -> V mode
[~, idx_max_h11] = max(abs(ch.am_h11).^2);
tx_az_max1 = ch.tx_az(idx_max_h11);
tx_el_max1 = ch.tx_el(idx_max_h11);
rx_az_max1 = ch.rx_az(idx_max_h11);
rx_el_max1 = ch.rx_el(idx_max_h11);

% find the ray with maximum power for H -> H mode
[~, idx_max_h22] = max(abs(ch.am_h22).^2);
tx_az_max2 = ch.tx_az(idx_max_h22);
tx_el_max2 = ch.tx_el(idx_max_h22);
rx_az_max2 = ch.rx_az(idx_max_h22);
rx_el_max2 = ch.rx_el(idx_max_h22);

% Processing of V -> V channel
% TX antenna gain
am_tx_h11 = ant_gain(ch.am_h11, ch.tx_az, ch.tx_el, tx_az_max1, tx_el_max1, cfg.paa);

% RX antenna gain
am_rx_h11 = ant_gain(am_tx_h11, ch.rx_az, ch.rx_el, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h11 = am_rx_h11;

% Processing of H -> H channel
% TX antenna gain
am_tx_h22 = ant_gain(ch.am_h22, ch.tx_az, ch.tx_el, tx_az_max2, tx_el_max2, cfg.paa);

% RX antenna gain
am_rx_h22 = ant_gain(am_tx_h22, ch.rx_az, ch.rx_el, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h22 = am_rx_h22;

% Processing of V -> H channel
% TX antenna gain
am_tx_h12 = ant_gain(ch.am_h12, ch.tx_az, ch.tx_el, tx_az_max1, tx_el_max1, cfg.paa);

% RX antenna gain
am_rx_h12 = ant_gain(am_tx_h12, ch.rx_az, ch.rx_el, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h12 = am_rx_h12;

% Processing of H -> V channel
% TX antenna gain
am_tx_h21 = ant_gain(ch.am_h21, ch.tx_az, ch.tx_el, tx_az_max2, tx_el_max2, cfg.paa);

% RX antenna gain
am_rx_h21 = ant_gain(am_tx_h21, ch.rx_az, ch.rx_el, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h21 = am_rx_h21;

% sort impulse response values
[toa, ind_sort] = sort(ch.toa);
imp_res.toa_11 = toa;
imp_res.h11 = imp_res_h11(ind_sort);
imp_res.h12 = imp_res_h12(ind_sort);
imp_res.h21 = imp_res_h21(ind_sort);
imp_res.h22 = imp_res_h22(ind_sort);
