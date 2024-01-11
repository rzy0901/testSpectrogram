% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_ch_model.m
%    Authors:       A. Lomayev, R. Maslennikov, Y.Gagiev
%    Version:       1.0
%    History:       May 2010 created
%                   April 2016 updated
%
%  *************************************************************************************
%    Description:
%
%    function returns channel impulse response for Conference Room (CR) environment
%
%    [imp_res] = cr_ch_model()
%
%    Inputs:
%
%       no inputs, parameters are set in cr_ch_cfg.m configuration file
%
%    Outputs:
%
%       1. imp_res - channel impulse response
%
%    Update: changed function interfaces
%  *************************************************************************************/
function [imp_res] = cr_ch_model(Radar_x,Radar_y,Radar_z)

% load configuration structure <- cr_ch_cfg.m
cfg = cr_ch_cfg;

% generate space-time channel impulse response realization
[ch] = gen_cr_ch(cfg,Radar_x,Radar_y,Radar_z);

% apply beamforming algorithm
[imp_res] = beamforming(cfg.bf, ch);

% continuous time to descrete time conversion
[imp_res] = digitize(cfg.bf.ant_type, cfg.sample_rate, imp_res);

% normalization according to Pnorm parameter
[imp_res] = normalize(cfg.bf.ant_type, cfg.Pnorm, imp_res);
