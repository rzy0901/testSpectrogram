% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_ch_cfg.m
%    Authors:       A. Lomayev, R. Maslennikov, Y.Gagiev
%    Version:       1.0
%    History:       May 2010 created
%                   April 2016 updated
%
%  *************************************************************************************
%    Description:
%
%    configuration function returns structure contained main parameters for
%    channel function cr_ch_model.m generating channel impulse response for
%    Conference Room (CR) environment
%
%    [cfg] = cr_ch_cfg()
%
%    Inputs: no inputs
%
%    Outputs: configuration structure
%
%    cfg.field    - common parameters
%    cfg.cr.field - space temporal clusters distribution related parameters
%    cfg.bf.field - beamforming parameters
%    cfg.bf.paa.field - parameters of phased antenna array
%
%    Update: added new fields to support new SU-MIMO configurations
%  *************************************************************************************/
function [cfg] = cr_ch_cfg()

% COMMON PARAMETERS

cfg.Pnorm       = 0;        % normalization parameter: 0 - w/o normalization, 1 - apply normalization for output channel impulse response
cfg.sample_rate = 0.2;     % sample rate in [GHz] applied for continuous time to discrete time channel impulse response conversion

% CFG.CR SUBSTRUCTURE   
cfg.cr.ap_sp = 0;           % parameter selects subscenario: 0 - STA-STA subscenario, 1 - STA-AP subscenario

cfg.cr.D    = 0;            % distance in [meters] between TX and RX, note that when the AP is placed near the ceiling distance D between TX and RX is set in horizontal plane
cfg.cr.Plos = 0;            % LOS (Line-of-Sight) parameter, permitted values: 0 - corresponds to NLOS scenario, 1 - corresponds to LOS scenario

% probabilities for STA-STA subscenario
cfg.cr.Psta_1st_c  = 1;     % probability that the cluster is present (i.e. not blocked) for the 1st order reflections from ceiling
cfg.cr.Psta_1st_w  = 0.76;  % probability that the cluster is present (i.e. not blocked) for the 1st order reflections from walls
cfg.cr.Psta_2nd_wc = 0.963; % probability that the cluster is present (i.e. not blocked) for the 2nd order wall-ceiling (ceiling-wall) reflections
cfg.cr.Psta_2nd_w  = 0.825; % probability that the cluster is present (i.e. not blocked) for the 2nd order reflections from walls

% probabilities for STA-AP subscenario
cfg.cr.Pap_1st = 0.874;     % probability that the cluster is present (i.e. not blocked) for the 1st order reflections from walls
cfg.cr.Pap_2nd = 0.93;      % probability that the cluster is present (i.e. not blocked) for the 2nd order reflections from walls

% CFG.BF SUBSTRUCTURE (BEAMFORMING & ANTENNAS RELATED PARAMETERS)

% antenna type parameter selects antenna type in beamforming search procedure
cfg.bf.ant_type = 1;     % Antenna type: 1 - conf#1, 2 - conf#2, 3 - conf#3, 4 - conf#4, 5 - conf#5

cfg.bf.ps     = 1;          % polarization support parameter: 0 - TX/RX polarization vectors are not applied, 1 - polarization is applied
cfg.bf.pol(1) = 1;          % antenna polarization type on TX side: 0 - linear in theta direction, 1 - linear in thi direction, 2 - LHCP, 3 - RHCP
cfg.bf.pol(2) = 0;          % antenna polarization type on RX side: 0 - linear in theta direction, 1 - linear in thi direction, 2 - LHCP, 3 - RHCP

% configuration #3
cfg.bf.tx_pol(1) = 1;       % antenna polarization type for PAA#1 on TX side: 0 - linear in theta direction, 1 - linear in thi direction
cfg.bf.tx_pol(2) = 0;       % antenna polarization type for PAA#2 on TX side: 0 - linear in theta direction, 1 - linear in thi direction
cfg.bf.dist_tx = 0.3;       % distance between PAAs on TX side, [meters]

cfg.bf.rx_pol(1) = 1;       % antenna polarization type for PAA#1 on RX side: 0 - linear in theta direction, 1 - linear in thi direction
cfg.bf.rx_pol(2) = 0;       % antenna polarization type for PAA#2 on RX side: 0 - linear in theta direction, 1 - linear in thi direction
cfg.bf.dist_rx = 0.3;       % distance between PAAs on X side, [meters]

% Number of elements through X-axis
cfg.bf.paa.Nx = 1;
% Number of elements through Y-axis
cfg.bf.paa.Ny = 1;
% Wavelength
% cfg.bf.paa.lyam = 15e-3;
cfg.bf.paa.lyam = 0.0857;%3.5GHz
% Distance between elements on X-axis
cfg.bf.paa.dx = cfg.bf.paa.lyam / 2;
% Distance between elements on Y-axis
cfg.bf.paa.dy = cfg.bf.paa.lyam / 2;
% Angle of Tx rotation in azimuth plane along LOS direction
cfg.bf.paa.phi_tx = rand() * 360 * 1;
% Angle of Rx rotation in azimuth plane along LOS direction
cfg.bf.paa.phi_rx = rand() * 360 * 1;