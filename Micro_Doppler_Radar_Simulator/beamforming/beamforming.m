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
%    function calls beamforming Max Power Ray algorithm
%
%    [imp_res,toa] = beamforming(cfg,ch)
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
function [imp_res] = beamforming(cfg, ch)

% Calling max power ray algorithm
% Others beamforming algorithms can be supported
imp_res = max_power_ray(cfg, ch);

end

