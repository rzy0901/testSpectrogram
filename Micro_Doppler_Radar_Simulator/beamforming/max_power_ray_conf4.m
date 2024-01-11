% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     max_power_ray_conf4.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       December 2015 created
%
%  *************************************************************************************
%    Description:
%
%    Executes max power ray alg for configuration #4
%
%    [imp_res] = max_power_ray_conf4(cfg, ch)
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
function [imp_res] = max_power_ray_conf4(cfg, ch)

shift_tx11 = ch.tx_az_11(end);
shift_rx11 = ch.rx_az_11(end);
shift_tx22 = ch.tx_az_22(end);
shift_rx22 = ch.rx_az_22(end);

ch.tx_az_11(end) = [];
ch.rx_az_11(end) = [];
ch.tx_az_22(end) = [];
ch.rx_az_22(end) = [];

% conversion of angles to new TX/RX coordinates
% (x-axis along LOS, z-axis normal to horizontal plane, y-axis is added for right tern)
% azimuth angle thi: [0:360], elevation angle theta: [0:180]
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

%% Find max power rays in sub-channels
[power_max_11vv, idx_max_11vv] = max(abs(ch.am_11vv).^2);
% Tx1
tx_az_max1 = ch.tx_az_11(idx_max_11vv);
tx_el_max1 = ch.tx_el_11(idx_max_11vv);
% Rx1
rx_az_max1 = ch.rx_az_11(idx_max_11vv);
rx_el_max1 = ch.rx_el_11(idx_max_11vv);

[power_max_11hh, idx_max_11hh] = max(abs(ch.am_11hh).^2);
% Tx1
tx_az_max2 = ch.tx_az_11(idx_max_11hh);
tx_el_max2 = ch.tx_el_11(idx_max_11hh);
% Rx1
rx_az_max2 = ch.rx_az_11(idx_max_11hh);
rx_el_max2 = ch.rx_el_11(idx_max_11hh);

[power_max_22vv, idx_max_22vv] = max(abs(ch.am_22vv).^2);
% Tx2
tx_az_max3 = ch.tx_az_22(idx_max_22vv);
tx_el_max3 = ch.tx_el_22(idx_max_22vv);
% Rx2
rx_az_max3 = ch.rx_az_22(idx_max_22vv);
rx_el_max3 = ch.rx_el_22(idx_max_22vv);

[power_max_22hh, idx_max_22hh] = max(abs(ch.am_22hh).^2);
% Tx2
tx_az_max4 = ch.tx_az_22(idx_max_22hh);
tx_el_max4 = ch.tx_el_22(idx_max_22hh);
% Rx2
rx_az_max4 = ch.rx_az_22(idx_max_22hh);
rx_el_max4 = ch.rx_el_22(idx_max_22hh);

%% From Tx1 -> Rx1
% V -> V
am_tx_11vv = ant_gain(ch.am_11vv, ch.tx_az_11, ch.tx_el_11, tx_az_max1, tx_el_max1, cfg.paa);
am_rx_11vv = ant_gain(am_tx_11vv, ch.rx_az_11, ch.rx_el_11, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h11 = am_rx_11vv;

% V -> H
am_tx_11vh = ant_gain(ch.am_11vh, ch.tx_az_11, ch.tx_el_11, tx_az_max1, tx_el_max1, cfg.paa);
am_rx_11vh = ant_gain(am_tx_11vh, ch.rx_az_11, ch.rx_el_11, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h12 = am_rx_11vh;

% H -> V
am_tx_11hv = ant_gain(ch.am_11hv, ch.tx_az_11, ch.tx_el_11, tx_az_max2, tx_el_max2, cfg.paa);
am_rx_11hv = ant_gain(am_tx_11hv, ch.rx_az_11, ch.rx_el_11, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h21 = am_rx_11hv;

% H -> H
am_tx_11hh = ant_gain(ch.am_11hh, ch.tx_az_11, ch.tx_el_11, tx_az_max2, tx_el_max2, cfg.paa);
am_rx_11hh = ant_gain(am_tx_11hh, ch.rx_az_11, ch.rx_el_11, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h22 = am_rx_11hh;

%% From Tx2 -> Rx2
% V -> V
am_tx_22vv = ant_gain(ch.am_22vv, ch.tx_az_22, ch.tx_el_22, tx_az_max3, tx_el_max3, cfg.paa);
am_rx_22vv = ant_gain(am_tx_22vv, ch.rx_az_22, ch.rx_el_22, rx_az_max3, rx_el_max3, cfg.paa);
imp_res_h33 = am_rx_22vv;

% V -> H
am_tx_22vh = ant_gain(ch.am_22vh, ch.tx_az_22, ch.tx_el_22, tx_az_max3, tx_el_max3, cfg.paa);
am_rx_22vh = ant_gain(am_tx_22vh, ch.rx_az_22, ch.rx_el_22, rx_az_max4, rx_el_max4, cfg.paa);
imp_res_h34 = am_rx_22vh;

% H -> V
am_tx_22hv = ant_gain(ch.am_22hv, ch.tx_az_22, ch.tx_el_22, tx_az_max4, tx_el_max4, cfg.paa);
am_rx_22hv = ant_gain(am_tx_22hv, ch.rx_az_22, ch.rx_el_22, rx_az_max3, rx_el_max3, cfg.paa);
imp_res_h43 = am_rx_22hv;

% H -> H
am_tx_22hh = ant_gain(ch.am_22hh, ch.tx_az_22, ch.tx_el_22, tx_az_max4, tx_el_max4, cfg.paa);
am_rx_22hh = ant_gain(am_tx_22hh, ch.rx_az_22, ch.rx_el_22, rx_az_max4, rx_el_max4, cfg.paa);
imp_res_h44 = am_rx_22hh;

%% From Tx1 -> Rx2
% V -> V
am_tx_12vv = ant_gain(ch.am_12vv, ch.tx_az_12, ch.tx_el_12, tx_az_max1, tx_el_max1, cfg.paa);
am_rx_12vv = ant_gain(am_tx_12vv, ch.rx_az_12, ch.rx_el_12, rx_az_max3, rx_el_max3, cfg.paa);
imp_res_h13 = am_rx_12vv;

% V -> H
am_tx_12vh = ant_gain(ch.am_12vh, ch.tx_az_12, ch.tx_el_12, tx_az_max1, tx_el_max1, cfg.paa);
am_rx_12vh = ant_gain(am_tx_12vh, ch.rx_az_12, ch.rx_el_12, rx_az_max4, rx_el_max4, cfg.paa);
imp_res_h14 = am_rx_12vh;

% H -> V
am_tx_12hv = ant_gain(ch.am_12hv, ch.tx_az_12, ch.tx_el_12, tx_az_max2, tx_el_max2, cfg.paa);
am_rx_12hv = ant_gain(am_tx_12hv, ch.rx_az_12, ch.rx_el_12, rx_az_max3, rx_el_max3, cfg.paa);
imp_res_h23 = am_rx_12hv;

% H -> H
am_tx_12hh = ant_gain(ch.am_12hh, ch.tx_az_12, ch.tx_el_12, tx_az_max2, tx_el_max2, cfg.paa);
am_rx_12hh = ant_gain(am_tx_12hh, ch.rx_az_12, ch.rx_el_12, rx_az_max4, rx_el_max4, cfg.paa);
imp_res_h24 = am_rx_12hh;

%% From Tx2 -> Rx1
% V -> V
am_tx_21vv = ant_gain(ch.am_21vv, ch.tx_az_21, ch.tx_el_21, tx_az_max3, tx_el_max3, cfg.paa);
am_rx_21vv = ant_gain(am_tx_21vv, ch.rx_az_21, ch.rx_el_21, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h31 = am_rx_21vv;

% V -> H
am_tx_21vh = ant_gain(ch.am_21vh, ch.tx_az_21, ch.tx_el_21, tx_az_max3, tx_el_max3, cfg.paa);
am_rx_21vh = ant_gain(am_tx_21vh, ch.rx_az_21, ch.rx_el_21, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h32 = am_rx_21vh;

% H -> V
am_tx_21hv = ant_gain(ch.am_21hv, ch.tx_az_21, ch.tx_el_21, tx_az_max4, tx_el_max4, cfg.paa);
am_rx_21hv = ant_gain(am_tx_21hv, ch.rx_az_21, ch.rx_el_21, rx_az_max1, rx_el_max1, cfg.paa);
imp_res_h41 = am_rx_21hv;

% H -> H
am_tx_21hh = ant_gain(ch.am_21hh, ch.tx_az_21, ch.tx_el_21, tx_az_max4, tx_el_max4, cfg.paa);
am_rx_21hh = ant_gain(am_tx_21hh, ch.rx_az_21, ch.rx_el_21, rx_az_max2, rx_el_max2, cfg.paa);
imp_res_h42 = am_rx_21hh;

%% sort impulse response values
[toa_11, ind_sort] = sort(ch.toa_11);
imp_res.h11 = imp_res_h11(ind_sort);
imp_res.h12 = imp_res_h12(ind_sort);
imp_res.h21 = imp_res_h21(ind_sort);
imp_res.h22 = imp_res_h22(ind_sort);
imp_res.toa_11 = toa_11;

[toa_12, ind_sort] = sort(ch.toa_12);
imp_res.h13 = imp_res_h13(ind_sort);
imp_res.h14 = imp_res_h14(ind_sort);
imp_res.h23 = imp_res_h23(ind_sort);
imp_res.h24 = imp_res_h24(ind_sort);
imp_res.toa_12 = toa_12;

[toa_21, ind_sort] = sort(ch.toa_21);
imp_res.h31 = imp_res_h31(ind_sort);
imp_res.h32 = imp_res_h32(ind_sort);
imp_res.h41 = imp_res_h41(ind_sort);
imp_res.h42 = imp_res_h42(ind_sort);
imp_res.toa_21 = toa_21;

[toa_22, ind_sort] = sort(ch.toa_22);
imp_res.h33 = imp_res_h33(ind_sort);
imp_res.h34 = imp_res_h34(ind_sort);
imp_res.h43 = imp_res_h43(ind_sort);
imp_res.h44 = imp_res_h44(ind_sort);
imp_res.toa_22 = toa_22;

end
