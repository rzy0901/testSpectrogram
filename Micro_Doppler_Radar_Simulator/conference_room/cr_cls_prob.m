% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_cls_prob.m
%    Authors:       A. Lomayev, R. Maslennikov
%    Version:       1.0
%    History:       May 2010 created
%
%  *************************************************************************************
%    Description:
% 
%    function generates clusters probabilities
%
%    [P] = cr_cls_prob(cfg)
%
%    Inputs:
%
%       1. cfg.ap_sp       - parameter selects subscenario: 0 - STA-STA, 1 - STA-AP
%       2. cfg.Psta_1st_c  - probability of 1st order reflections from ceiling in STA-STA subscenario
%       3. cfg.Psta_1st_w  - probability of 1st order reflections from walls in STA-STA subscenario
%       4. cfg.Psta_2nd_wc - probability of 2nd order wall-ceiling (ceiling-wall) reflections in STA-STA subscenario
%       5. cfg.Psta_2nd_w  - probability of 2nd order reflections from walls in STA-STA subscenario
%       6. cfg.Pap_1st     - probability of 1st order reflections from walls in STA-AP subscenario
%       7. cfg.Pap_2nd     - probability of 2nd order reflections from walls in STA-AP subscenario
%
%    Outputs:
%
%       1. P - array of clusters probabilities
% 
%  *************************************************************************************/
function [P] = cr_cls_prob(cfg)

switch(cfg.ap_sp)
    case 0, % STA - STA
        % 1st order ceiling cluster
        P(1) = 1;
        % 1st order wall clusters
        prob = 4.*rand(1,1);
        num = ceil(prob);
        P(2:5) = 1;
        P(1+num) = cfg.Psta_1st_w;
        % 2nd order wall-ceiling (ceiling-wall) clusters
        prob = 4.*rand(1,1);
        num = ceil(prob);
        P(6:9) = 1;
        P(5+num) = cfg.Psta_2nd_wc;
        % 2nd order wall clusters
        P(10:17) = cfg.Psta_2nd_w;
        
    case 1, % STA - AP
        % 1st order clusters
        prob = 4.*rand(1,1);
        num = ceil(prob);
        P(1:4) = 1;
        P(num) = cfg.Pap_1st;
        % 2nd order clusters
        P(5:12) = cfg.Pap_2nd;        
end