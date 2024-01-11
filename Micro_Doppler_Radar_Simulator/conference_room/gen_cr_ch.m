% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     gen_cr_ch.m
%    Authors:       Y.Gagiev, A. Lomayev, R. Maslennikov
%    Version:       1.0
%    History:       May 2010 created
%                   April 2016 updated
%
%  *************************************************************************************
%    Description:
%
%    Generates channel regarding type of SU-MIMO configuration
%
%    [ch] = gen_cr_ch(cfg,ps,pol)
%
%    Inputs:
%
%      1.  cfg.ap_sp       - parameter selects subscenario: 0 - STA-STA, 1 - STA-AP
%      1.  cfg.D           - distance in [meters] between TX and RX
%      2.  cfg.Plos        - LOS parameter: 0 - LOS between TX and RX is blocked, 1 - non-blocked
%      3.  cfg.Psta_1st_c  - probability of 1st order reflections from ceiling in STA-STA subscenario
%      4.  cfg.Psta_1st_w  - probability of 1st order reflections from walls in STA-STA subscenario
%      5.  cfg.Psta_2nd_wc - probability of 2nd order wall-ceiling (ceiling-wall) reflections in STA-STA subscenario
%      6.  cfg.Psta_2nd_w  - probability of 2nd order reflections from walls in STA-STA subscenario
%      7.  cfg.Pap_1st     - probability of 1st order reflections from walls in STA-AP subscenario
%      8.  cfg.Pap_2nd     - probability of 2nd order reflections from walls in STA-AP subscenario
%      9.  ps              - polarization support parameter: 0 - TX/RX polarization vectors are not applied, 1 - polarization is applied
%      10. pol             - antennas polarization 1x2 vector: pol(1) - polarization type for TX antenna, pol(2) - polarization type for RX antenna
%
%    Outputs:
%
%       1. ch - structure with the channel
%
%    Update: Added switch to select required function for channel
%    generation
%  *************************************************************************************/
function [ch] = gen_cr_ch( cfg , Radar_x , Radar_y , Radar_z )

switch (cfg.bf.ant_type)
    case 1
        ch = gen_cr_ch_conf1(cfg.cr, cfg.bf.pol, cfg.bf.paa.lyam,Radar_x,Radar_y,Radar_z);
    case 2,
        ch = gen_cr_ch_conf2(cfg.cr, cfg.bf.paa.lyam);
    case 3,
        ch = gen_cr_ch_conf3(cfg.cr, cfg.bf.tx_pol, cfg.bf.rx_pol, cfg.bf.dist_tx, cfg.bf.dist_rx, cfg.bf.paa.phi_tx, cfg.bf.paa.phi_rx, cfg.bf.paa.lyam);
    case 4,
        ch = gen_cr_ch_conf4(cfg.cr, cfg.bf.dist_tx, cfg.bf.dist_rx, cfg.bf.paa.phi_tx, cfg.bf.paa.phi_rx, cfg.bf.paa.lyam);
    case 5,
        ch = gen_cr_ch_conf5(cfg.cr, cfg.bf.pol, cfg.bf.paa.lyam);
    otherwise,
        error('Prohibited value of "cfg.bf.ant_type" parameter');

end