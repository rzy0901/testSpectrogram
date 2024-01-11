% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     digitize.m
%    Authors:       A. Lomayev, R. Maslennikov, Y. Gagiev
%    Version:       1.0
%    History:       May 2010 created
%                   April 2016 updated
%
%  *************************************************************************************
%    Description:
% 
%    function returns channel response due to target sample rate
%
%    [imp_res] = digitize(ant_type, sample_rate, imp_res)
%
%    Inputs:
%
%    1. ant_type    - antenna type that defines internal structure of imp_res
%    2. sample_rate - sample rate in [GHz]
%    3. imp_res     - structure with impulse responses after beamforming
% 
%    Outputs:
%
%    1. imp_res - structure with digitized impulse responses
%
%    Update: Conversion to discrete time is made with aligning on
%    channel tap with min time
%  *************************************************************************************/
function [samp_imp_res] = digitize(ant_type, sample_rate, imp_res)

% Continuous time to descrete time conversion
switch(ant_type)        
    case 1,
        min_time = min(imp_res.toa_11);
        samp_imp_res.h11 = ct2dt(imp_res.h11, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h12 = ct2dt(imp_res.h12, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h21 = ct2dt(imp_res.h21, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h22 = ct2dt(imp_res.h22, imp_res.toa_11, sample_rate, min_time);
        
    case 2,
        min_time = min(imp_res.toa_11);
        samp_imp_res.h11 = ct2dt(imp_res.h11, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h12 = ct2dt(imp_res.h12, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h21 = ct2dt(imp_res.h21, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h22 = ct2dt(imp_res.h22, imp_res.toa_11, sample_rate, min_time);
        
    case 3,
        min_toa_11 = min(imp_res.toa_11);
        min_toa_12 = min(imp_res.toa_12);
        min_toa_21 = min(imp_res.toa_21);
        min_toa_22 = min(imp_res.toa_22);
        
        min_time = min(min(min(min_toa_11, min_toa_12), min_toa_21), min_toa_22);
        samp_imp_res.h11 = ct2dt(imp_res.h11, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h12 = ct2dt(imp_res.h12, imp_res.toa_12, sample_rate, min_time);
        samp_imp_res.h21 = ct2dt(imp_res.h21, imp_res.toa_21, sample_rate, min_time);
        samp_imp_res.h22 = ct2dt(imp_res.h22, imp_res.toa_22, sample_rate, min_time);
        
    case 4,
        min_toa_11 = min(imp_res.toa_11);
        min_toa_12 = min(imp_res.toa_12);
        min_toa_21 = min(imp_res.toa_21);
        min_toa_22 = min(imp_res.toa_22);
        
        min_time = min(min(min(min_toa_11, min_toa_12), min_toa_21), min_toa_22);
        samp_imp_res.h11 = ct2dt(imp_res.h11, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h12 = ct2dt(imp_res.h12, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h21 = ct2dt(imp_res.h21, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h22 = ct2dt(imp_res.h22, imp_res.toa_11, sample_rate, min_time);

        samp_imp_res.h13 = ct2dt(imp_res.h13, imp_res.toa_12, sample_rate, min_time);
        samp_imp_res.h14 = ct2dt(imp_res.h14, imp_res.toa_12, sample_rate, min_time);
        samp_imp_res.h23 = ct2dt(imp_res.h23, imp_res.toa_12, sample_rate, min_time);
        samp_imp_res.h24 = ct2dt(imp_res.h24, imp_res.toa_12, sample_rate, min_time);

        samp_imp_res.h31 = ct2dt(imp_res.h31, imp_res.toa_21, sample_rate, min_time);
        samp_imp_res.h32 = ct2dt(imp_res.h32, imp_res.toa_21, sample_rate, min_time);
        samp_imp_res.h41 = ct2dt(imp_res.h41, imp_res.toa_21, sample_rate, min_time);
        samp_imp_res.h42 = ct2dt(imp_res.h42, imp_res.toa_21, sample_rate, min_time);
        
        samp_imp_res.h33 = ct2dt(imp_res.h33, imp_res.toa_22, sample_rate, min_time);
        samp_imp_res.h34 = ct2dt(imp_res.h34, imp_res.toa_22, sample_rate, min_time);
        samp_imp_res.h43 = ct2dt(imp_res.h43, imp_res.toa_22, sample_rate, min_time);
        samp_imp_res.h44 = ct2dt(imp_res.h44, imp_res.toa_22, sample_rate, min_time);
        
    case 5,
        min_time = min(imp_res.toa_11);
        samp_imp_res.h11 = ct2dt(imp_res.h11, imp_res.toa_11, sample_rate, min_time);
        samp_imp_res.h12 = ct2dt(imp_res.h12, imp_res.toa_11, sample_rate, min_time);
        
    otherwise,
       error('Prohibited value of "ant_type" parameter');
end


function [h_dt] = ct2dt(h_ct, t, sample_rate, min_time)

t_s = 1./sample_rate;
t_0 = t - min_time;
N = round(max(t_0)./t_s) + 1;
h_dt = zeros(1,N);

for ray_ix = 1:length(t_0)
    time_bin = round(t_0(ray_ix)./t_s) + 1;
    h_dt(time_bin) = h_dt(time_bin) + h_ct(ray_ix);
end
