% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     beamforming.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       December 2015 created
%
%  *************************************************************************************
%    Description:
%
%    function calls max power ray alg for specified ant_type
%
%    [ imp_res ] = max_power_ray(cfg, ch)
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
function [ imp_res ] = max_power_ray(cfg, ch)

switch(cfg.ant_type)
    case 1, % phase antenna array with single polarization
       imp_res = max_power_ray_conf1(cfg, ch);
    case 2, % phased antenna array with dual polarization
       imp_res = max_power_ray_conf2(cfg, ch);
    case 3, % double phased antenna array
       imp_res = max_power_ray_conf3(cfg, ch);
    case 4, % single phased antenna array on Tx, Rx; Rx receives w/ dual polarization
       imp_res = max_power_ray_conf4(cfg, ch);
    case 5,
       imp_res = max_power_ray_conf5(cfg, ch);
end
