% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  60 GHz Conference Room Channel Model
%    File Name:     gen_cr_ch_mimo.m
%    Authors:       Y. Gagiev
%    Version:      
%    History:       April 2016 created
%
%  *************************************************************************************
%    Description:
%
%    Generates channel for conf#3
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
%      9.  tx_pol          - polarization type on TX
%      10. rx_pol          - polarization type on RX
%      11. dist_tx         - distance between PAAs on Tx side
%      12. dist_rx         - distance between PAAs on Rx side
%      13. lambda          - wavelength
%
%    Outputs:
%
%      1.  ch.am_11        - amplitudes array for link from 1st PAA (TX) to 1st PAA (RX)
%      2.  ch.am_12        - amplitudes array for link from 1st PAA (TX) to 2nd PAA (RX)
%      3.  ch.am_21        - amplitudes array for link from 2nd PAA (TX) to 1st PAA (RX)
%      4.  ch.am_22        - amplitudes array for link from 2nd PAA (TX) to 2nd PAA (RX)
%
%      5.  ch.toa_11       - times of arrival array for link from 1st PAA (TX) to 1st PAA (RX)
%      6.  ch.toa_12       - times of arrival array for link from 1st PAA (TX) to 2nd PAA (RX)
%      7.  ch.toa_21       - times of arrival array for link from 2nd PAA (TX) to 1st PAA (RX)
%      8.  ch.toa_22       - times of arrival array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      
%      9.  ch.tx_az_11     - TX azimuths array for link from 1st PAA (TX) to 1st PAA (RX)
%      10. ch.tx_el_11     - TX elevations array for link from 1st PAA (TX) to 1st PAA (RX)
%      11. ch.rx_az_11     - RX azimuths array for link from 1st PAA (TX) to 1st PAA (RX)
%      12. ch.rx_el_11     - RX elevations array for link from 1st PAA (TX) to 1st PAA (RX)
%      
%      13. ch.tx_az_12     - TX azimuths array for link from 1st PAA (TX) to 2nd PAA (RX)
%      14. ch.tx_el_12     - TX elevations array for link from 1st PAA (TX) to 2nd PAA (RX)
%      15. ch.rx_az_12     - RX azimuths array for link from 1st PAA (TX) to 2nd PAA (RX)
%      16. ch.rx_el_12     - RX elevations array for link from 1st PAA (TX) to 2nd PAA (RX)
%      
%      17. ch.tx_az_21     - TX azimuths array for link from 2nd PAA (TX) to 1st PAA (RX)
%      18. ch.tx_el_21     - TX elevations array for link from 2nd PAA (TX) to 1st PAA (RX)
%      19. ch.rx_az_21     - RX azimuths array for link from 2nd PAA (TX) to 1st PAA (RX)
%      20. ch.rx_el_21     - RX elevations array for link from 2nd PAA (TX) to 1st PAA (RX)
%      
%      21. ch.tx_az_22     - TX azimuths array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      22. ch.tx_el_22     - TX elevations array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      23. ch.rx_az_22     - RX azimuths array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      24. ch.rx_el_22     - RX elevations array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      
%  *************************************************************************************/
function [ch] = gen_cr_ch_conf3(cfg, tx_pol, rx_pol, dist_tx, dist_rx, phi_tx, phi_rx, lambda)

% NLOS attenuation coefficients (due to reflection)
[Gr_11, Gr_12, Gr_21, Gr_22] = cr_ref_loss_conf3(cfg.ap_sp, tx_pol, rx_pol);

% NLOS clusters probabilities
P = cr_cls_prob(cfg);

% NLOS attenuation coefficients (due to human blockage)
Gb = cr_atten_coef(cfg.ap_sp);

% choose subscenario
switch (cfg.ap_sp)
    case 0, % STA - STA
        
        % distance between TX and RX
        D = cfg.D;
        
        % generate NLOS clusters
        [cls_11, cls_12, cls_21, cls_22] = cr_gen_cls_ant2(0, dist_tx * 100, dist_rx * 100, phi_tx, phi_rx);
            
        % clusters distances
        dist_11 = (cls_11.toa.*1e-9).*(3e8) + D;
        dist_12 = (cls_12.toa.*1e-9).*(3e8) + sqrt(D^2 + dist_tx * dist_rx);
        dist_21 = (cls_21.toa.*1e-9).*(3e8) + sqrt(D^2 + dist_tx * dist_rx);
        dist_22 = (cls_22.toa.*1e-9).*(3e8) + D;
        
        % calculate attenuation constants for clusters
        % attenuation due to propagation
        Gp_11 = sqrt((lambda.^2)./((4.*pi.*(dist_11)).^2));
        Gp_12 = sqrt((lambda.^2)./((4.*pi.*(dist_12)).^2));
        Gp_21 = sqrt((lambda.^2)./((4.*pi.*(dist_21)).^2));
        Gp_22 = sqrt((lambda.^2)./((4.*pi.*(dist_22)).^2));
        
        % attenuation constants
        G_11 = Gp_11.*Gr_11;
        G_12 = Gp_12.*Gr_12;
        G_21 = Gp_21.*Gr_21;
        G_22 = Gp_22.*Gr_22;
        
        % generate intra-clusters structure
        toa_11 = [];
        am_11 = [];
        tx_az_11 = [];
        tx_el_11 = [];
        rx_az_11 = [];
        rx_el_11 = [];
        
        toa_12 = [];
        am_12 = [];
        tx_az_12 = [];
        tx_el_12 = [];
        rx_az_12 = [];
        rx_el_12 = [];

        toa_21 = [];
        am_21 = [];
        tx_az_21 = [];
        tx_el_21 = [];
        rx_az_21 = [];
        rx_el_21 = [];

        toa_22 = [];
        am_22 = [];
        tx_az_22 = [];
        tx_el_22 = [];
        rx_az_22 = [];
        rx_el_22 = [];

        % LOS cluster
        if (cfg.Plos)
            toa_11 = 0;
            toa_12 = 0;
            toa_21 = 0;
            toa_22 = 0;
            am_11 = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            am_12 = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_21 = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_22 = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            
            % Apply polarization
            tx_pol_1 = polarization(tx_pol(1));
            tx_pol_2 = polarization(tx_pol(2));
                
            rx_pol_1 = polarization(rx_pol(1));
            rx_pol_2 = polarization(rx_pol(2));
            
            H = eye(2);
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            pol_coef_11 = rx_pol_1'*H*tx_pol_1;

            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            pol_coef_12 = rx_pol_1'*H*tx_pol_2;
            
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            pol_coef_21 = rx_pol_2'*H*tx_pol_1;

            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            pol_coef_22 = rx_pol_2'*H*tx_pol_2;
               
            am_11 = am_11.*pol_coef_11;
            am_12 = am_12.*pol_coef_12;
            am_21 = am_21.*pol_coef_21;
            am_22 = am_22.*pol_coef_22;
            
            tx_az_11 = 0;
            tx_el_11 = 0;
            rx_az_11 = 0;
            rx_el_11 = 0;
            
            tx_az_12 = 0;
            tx_el_12 = 0;
            rx_az_12 = 0;
            rx_el_12 = 0;

            tx_az_21 = 0;
            tx_el_21 = 0;
            rx_az_21 = 0;
            rx_el_21 = 0;

            tx_az_22 = 0;
            tx_el_22 = 0;
            rx_az_22 = 0;
            rx_el_22 = 0;

        end
        
        % NLOS clusters
        while (isempty(toa_11) | (toa_11 == 0))
            for i=1:17
                
                incls_11 = cr_gen_intra_cls(cls_11.toa(i));
                incls_12 = cr_gen_intra_cls(cls_12.toa(i));
                incls_21 = cr_gen_intra_cls(cls_21.toa(i));
                incls_22 = cr_gen_intra_cls(cls_22.toa(i));
                
                toa_11 = [toa_11; cls_11.toa(i) + incls_11.toa];
                toa_12 = [toa_12; cls_12.toa(i) + incls_12.toa];
                toa_21 = [toa_21; cls_21.toa(i) + incls_21.toa];
                toa_22 = [toa_22; cls_22.toa(i) + incls_22.toa];
                
                if rand(1,1) <= P(i)                    
                    am_11 = [am_11; incls_11.am.*G_11(i)];
                else
                    am_11 = [am_11; incls_11.am.*G_11(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_12 = [am_12; incls_12.am.*G_12(i)];
                else
                    am_12 = [am_12; incls_12.am.*G_12(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_21 = [am_21; incls_21.am.*G_21(i)];
                else
                    am_21 = [am_21; incls_21.am.*G_21(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_22 = [am_22; incls_22.am.*G_22(i)];
                else
                    am_22 = [am_22; incls_22.am.*G_22(i).*Gb(i)];
                end
                
                tx_az_11 = [tx_az_11; cls_11.tx_az(i) + incls_11.tx_az];
                tx_el_11 = [tx_el_11; cls_11.tx_el(i) + incls_11.tx_el];
                rx_az_11 = [rx_az_11; cls_11.rx_az(i) + incls_11.rx_az];
                rx_el_11 = [rx_el_11; cls_11.rx_el(i) + incls_11.rx_el];
                
                tx_az_12 = [tx_az_12; cls_12.tx_az(i) + incls_12.tx_az];
                tx_el_12 = [tx_el_12; cls_12.tx_el(i) + incls_12.tx_el];
                rx_az_12 = [rx_az_12; cls_12.rx_az(i) + incls_12.rx_az];
                rx_el_12 = [rx_el_12; cls_12.rx_el(i) + incls_12.rx_el];
                
                tx_az_21 = [tx_az_21; cls_21.tx_az(i) + incls_21.tx_az];
                tx_el_21 = [tx_el_21; cls_21.tx_el(i) + incls_21.tx_el];
                rx_az_21 = [rx_az_21; cls_21.rx_az(i) + incls_21.rx_az];
                rx_el_21 = [rx_el_21; cls_21.rx_el(i) + incls_21.rx_el];
                
                tx_az_22 = [tx_az_22; cls_22.tx_az(i) + incls_22.tx_az];
                tx_el_22 = [tx_el_22; cls_22.tx_el(i) + incls_22.tx_el];
                rx_az_22 = [rx_az_22; cls_22.rx_az(i) + incls_22.rx_az];
                rx_el_22 = [rx_el_22; cls_22.rx_el(i) + incls_22.rx_el];
            end
                tx_az_11 = [tx_az_11; cls_11.tx_az(18)];
                rx_az_11 = [rx_az_11; cls_11.rx_az(18)];
                
                tx_az_22 = [tx_az_22; cls_22.tx_az(18)];
                rx_az_22 = [rx_az_22; cls_22.rx_az(18)];
        end
        
    case 1, % STA - AP
        
        % distance between TX and RX
        D = sqrt( (cfg.D).^2 + (1.9).^2 );
        
        % generate NLOS clusters
        [cls_11, cls_12, cls_21, cls_22] = cr_gen_cls_ant2(1, dist_tx * 100, dist_rx * 100, phi_tx, phi_rx);
        
        % clusters distances
        dist_11 = (cls_11.toa.*1e-9).*(3e8) + D;
        dist_12 = (cls_12.toa.*1e-9).*(3e8) + D;
        dist_21 = (cls_21.toa.*1e-9).*(3e8) + D;
        dist_22 = (cls_22.toa.*1e-9).*(3e8) + D;
        
        % calculate attenuation constants for clusters
        % attenuation due to propagation
        Gp_11 = sqrt((lambda.^2)./((4.*pi.*(dist_11)).^2));
        Gp_12 = sqrt((lambda.^2)./((4.*pi.*(dist_12)).^2));
        Gp_21 = sqrt((lambda.^2)./((4.*pi.*(dist_21)).^2));
        Gp_22 = sqrt((lambda.^2)./((4.*pi.*(dist_22)).^2));
        
        % attenuation constants
        G_11 = Gp_11.*Gr_11;
        G_12 = Gp_12.*Gr_12;
        G_21 = Gp_21.*Gr_21;
        G_22 = Gp_22.*Gr_22;
        
        % generate intra-clusters structure
        toa_11 = [];
        toa_12 = [];
        toa_21 = [];
        toa_22 = [];
        am_11 = [];
        am_12 = [];
        am_21 = [];
        am_22 = [];
        
        tx_az_11 = [];
        tx_el_11 = [];
        rx_az_11 = [];
        rx_el_11 = [];
        
        tx_az_12 = [];
        tx_el_12 = [];
        rx_az_12 = [];
        rx_el_12 = [];
        
        tx_az_21 = [];
        tx_el_21 = [];
        rx_az_21 = [];
        rx_el_21 = [];
        
        tx_az_22 = [];
        tx_el_22 = [];
        rx_az_22 = [];
        rx_el_22 = [];
        
        % LOS cluster
        if (cfg.Plos)
            toa_11 = 0;
            toa_12 = 0;
            toa_21 = 0;
            toa_22 = 0;
            
            am_11 = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            am_12 = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_21 = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_22 = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            
            if (ps)
                tx_pol_1 = polarization(tx_pol(1));
                tx_pol_2 = polarization(tx_pol(2));
                
                rx_pol_1 = polarization(rx_pol(1));
                rx_pol_2 = polarization(tx_pol(2));
                
                H = eye(2);
                H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
                H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
                pol_coef_11 = rx_pol_1'*H*tx_pol_1;

                H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
                H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
                pol_coef_12 = rx_pol_1'*H*tx_pol_2;

                H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
                H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
                pol_coef_21 = rx_pol_2'*H*tx_pol_1;

                H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
                H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
                pol_coef_22 = rx_pol_2'*H*tx_pol_2;
                
                am_11 = am_11.*pol_coef_11;
                am_12 = am_12.*pol_coef_12;
                am_21 = am_21.*pol_coef_21;
                am_22 = am_22.*pol_coef_22;
            end
            
            tx_az_11 = 0;
            tx_el_11 = -(asin(1.9./D).*180./pi);
            rx_az_11 = 0;
            rx_el_11 = -tx_el_11;
            
            tx_az_12 = 0;
            tx_el_12 = -(asin(1.9./sqrt(D^2 + dist_tx * dist_rx)).*180./pi);
            rx_az_12 = 0;
            rx_el_12 = -tx_el_12;

            tx_az_21 = 0;
            tx_el_21 = -(asin(1.9./sqrt(D^2 + dist_tx * dist_rx)).*180./pi);
            rx_az_21 = 0;
            rx_el_21 = -tx_el_21;
            
            tx_az_22 = 0;
            tx_el_22 = -(asin(1.9./D).*180./pi);
            rx_az_22 = 0;
            rx_el_22 = -tx_el_22;
            
        end
        
        % NLOS clusters
        while (isempty(toa_11) | (toa_11 == 0))
            for i=1:12
                
                incls_11 = cr_gen_intra_cls(cls_11.toa(i));
                incls_12 = cr_gen_intra_cls(cls_12.toa(i));
                incls_21 = cr_gen_intra_cls(cls_21.toa(i));
                incls_22 = cr_gen_intra_cls(cls_22.toa(i));
                
                toa_11 = [toa_11; cls_11.toa(i) + incls_11.toa];
                toa_12 = [toa_12; cls_12.toa(i) + incls_12.toa];
                toa_21 = [toa_21; cls_21.toa(i) + incls_21.toa];
                toa_22 = [toa_22; cls_22.toa(i) + incls_22.toa];
                
                if rand(1,1) <= P(i)                    
                    am_11 = [am_11; incls_11.am.*G_11(i)];
                else
                    am_11 = [am_11; incls_11.am.*G_11(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_12 = [am_12; incls_12.am.*G_12(i)];
                else
                    am_12 = [am_12; incls_12.am.*G_12(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_21 = [am_21; incls_21.am.*G_21(i)];
                else
                    am_21 = [am_21; incls_21.am.*G_21(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_22 = [am_22; incls_22.am.*G_22(i)];
                else
                    am_22 = [am_22; incls_22.am.*G_22(i).*Gb(i)];
                end
                
                tx_az_11 = [tx_az_11; cls_11.tx_az(i) + incls_11.tx_az];
                tx_el_11 = [tx_el_11; cls_11.tx_el(i) + incls_11.tx_el];
                rx_az_11 = [rx_az_11; cls_11.rx_az(i) + incls_11.rx_az];
                rx_el_11 = [rx_el_11; cls_11.rx_el(i) + incls_11.rx_el];

                tx_az_12 = [tx_az_12; cls_12.tx_az(i) + incls_12.tx_az];
                tx_el_12 = [tx_el_12; cls_12.tx_el(i) + incls_12.tx_el];
                rx_az_12 = [rx_az_12; cls_12.rx_az(i) + incls_12.rx_az];
                rx_el_12 = [rx_el_12; cls_12.rx_el(i) + incls_12.rx_el];

                tx_az_21 = [tx_az_21; cls_21.tx_az(i) + incls_21.tx_az];
                tx_el_21 = [tx_el_21; cls_21.tx_el(i) + incls_21.tx_el];
                rx_az_21 = [rx_az_21; cls_21.rx_az(i) + incls_21.rx_az];
                rx_el_21 = [rx_el_21; cls_21.rx_el(i) + incls_21.rx_el];

                tx_az_22 = [tx_az_22; cls_22.tx_az(i) + incls_22.tx_az];
                tx_el_22 = [tx_el_22; cls_22.tx_el(i) + incls_22.tx_el];
                rx_az_22 = [rx_az_22; cls_22.rx_az(i) + incls_22.rx_az];
                rx_el_22 = [rx_el_22; cls_22.rx_el(i) + incls_22.rx_el];
            end
                tx_az_11 = [tx_az_11; cls_11.tx_az(13)];
                rx_az_11 = [rx_az_11; cls_11.rx_az(13)];
                
                tx_az_22 = [tx_az_22; cls_22.tx_az(13)];
                rx_az_22 = [rx_az_22; cls_22.rx_az(13)];
        end
    otherwise,
        error('Prohibited value of "cfg.cr.ap_sp" parameter');
end        

% check azimuth overflow

%11
ind_tx_az = find(tx_az_11 > 180);
tx_az_11(ind_tx_az) = tx_az_11(ind_tx_az) - 360;
ind_tx_az = find(tx_az_11 < -180);
tx_az_11(ind_tx_az) = tx_az_11(ind_tx_az) + 360;

ind_rx_az = find(rx_az_11 > 180);
rx_az_11(ind_rx_az) = rx_az_11(ind_rx_az) - 360;
ind_rx_az = find(rx_az_11 < -180);
rx_az_11(ind_rx_az) = rx_az_11(ind_rx_az) + 360;

% 12
ind_tx_az = find(tx_az_12 > 180);
tx_az_12(ind_tx_az) = tx_az_12(ind_tx_az) - 360;
ind_tx_az = find(tx_az_12 < -180);
tx_az_12(ind_tx_az) = tx_az_12(ind_tx_az) + 360;

ind_rx_az = find(rx_az_12 > 180);
rx_az_12(ind_rx_az) = rx_az_12(ind_rx_az) - 360;
ind_rx_az = find(rx_az_12 < -180);
rx_az_12(ind_rx_az) = rx_az_12(ind_rx_az) + 360;

% 21
ind_tx_az = find(tx_az_21 > 180);
tx_az_21(ind_tx_az) = tx_az_21(ind_tx_az) - 360;
ind_tx_az = find(tx_az_21 < -180);
tx_az_21(ind_tx_az) = tx_az_21(ind_tx_az) + 360;

ind_rx_az = find(rx_az_21 > 180);
rx_az_21(ind_rx_az) = rx_az_21(ind_rx_az) - 360;
ind_rx_az = find(rx_az_21 < -180);
rx_az_21(ind_rx_az) = rx_az_21(ind_rx_az) + 360;

% 22
ind_tx_az = find(tx_az_22 > 180);
tx_az_22(ind_tx_az) = tx_az_22(ind_tx_az) - 360;
ind_tx_az = find(tx_az_22 < -180);
tx_az_22(ind_tx_az) = tx_az_22(ind_tx_az) + 360;

ind_rx_az = find(rx_az_22 > 180);
rx_az_22(ind_rx_az) = rx_az_22(ind_rx_az) - 360;
ind_rx_az = find(rx_az_22 < -180);
rx_az_22(ind_rx_az) = rx_az_22(ind_rx_az) + 360;

% check elevation overflow

% 11
ind_tx_el = find(tx_el_11 > 90);
tx_el_11(ind_tx_el) = 180 - tx_el_11(ind_tx_el);
tx_az_11(ind_tx_el) = tx_az_11(ind_tx_el) + (-180).*(tx_az_11(ind_tx_el)>0) + 180.*(tx_az_11(ind_tx_el)<=0);

ind_tx_el = find(tx_el_11 < -90);
tx_el_11(ind_tx_el) = -(180 + tx_el_11(ind_tx_el));
tx_az_11(ind_tx_el) = tx_az_11(ind_tx_el) + (-180).*(tx_az_11(ind_tx_el)>0) + 180.*(tx_az_11(ind_tx_el)<=0);

ind_rx_el = find(rx_el_11 > 90);
rx_el_11(ind_rx_el) = 180 - rx_el_11(ind_rx_el);
rx_az_11(ind_rx_el) = rx_az_11(ind_rx_el) + (-180).*(rx_az_11(ind_rx_el)>0) + 180.*(rx_az_11(ind_rx_el)<=0);

ind_rx_el = find(rx_el_11 < -90);
rx_el_11(ind_rx_el) = -(180 + rx_el_11(ind_rx_el));
rx_az_11(ind_rx_el) = rx_az_11(ind_rx_el) + (-180).*(rx_az_11(ind_rx_el)>0) + 180.*(rx_az_11(ind_rx_el)<=0);

% 12
ind_tx_el = find(tx_el_12 > 90);
tx_el_12(ind_tx_el) = 180 - tx_el_12(ind_tx_el);
tx_az_12(ind_tx_el) = tx_az_12(ind_tx_el) + (-180).*(tx_az_12(ind_tx_el)>0) + 180.*(tx_az_12(ind_tx_el)<=0);

ind_tx_el = find(tx_el_12 < -90);
tx_el_12(ind_tx_el) = -(180 + tx_el_12(ind_tx_el));
tx_az_12(ind_tx_el) = tx_az_12(ind_tx_el) + (-180).*(tx_az_12(ind_tx_el)>0) + 180.*(tx_az_12(ind_tx_el)<=0);

ind_rx_el = find(rx_el_12 > 90);
rx_el_12(ind_rx_el) = 180 - rx_el_12(ind_rx_el);
rx_az_12(ind_rx_el) = rx_az_12(ind_rx_el) + (-180).*(rx_az_12(ind_rx_el)>0) + 180.*(rx_az_12(ind_rx_el)<=0);

ind_rx_el = find(rx_el_12 < -90);
rx_el_12(ind_rx_el) = -(180 + rx_el_12(ind_rx_el));
rx_az_12(ind_rx_el) = rx_az_12(ind_rx_el) + (-180).*(rx_az_12(ind_rx_el)>0) + 180.*(rx_az_12(ind_rx_el)<=0);

% 21
ind_tx_el = find(tx_el_21 > 90);
tx_el_21(ind_tx_el) = 180 - tx_el_21(ind_tx_el);
tx_az_21(ind_tx_el) = tx_az_21(ind_tx_el) + (-180).*(tx_az_21(ind_tx_el)>0) + 180.*(tx_az_21(ind_tx_el)<=0);

ind_tx_el = find(tx_el_21 < -90);
tx_el_21(ind_tx_el) = -(180 + tx_el_21(ind_tx_el));
tx_az_21(ind_tx_el) = tx_az_21(ind_tx_el) + (-180).*(tx_az_21(ind_tx_el)>0) + 180.*(tx_az_21(ind_tx_el)<=0);

ind_rx_el = find(rx_el_21 > 90);
rx_el_21(ind_rx_el) = 180 - rx_el_21(ind_rx_el);
rx_az_21(ind_rx_el) = rx_az_21(ind_rx_el) + (-180).*(rx_az_21(ind_rx_el)>0) + 180.*(rx_az_21(ind_rx_el)<=0);

ind_rx_el = find(rx_el_21 < -90);
rx_el_21(ind_rx_el) = -(180 + rx_el_21(ind_rx_el));
rx_az_21(ind_rx_el) = rx_az_21(ind_rx_el) + (-180).*(rx_az_21(ind_rx_el)>0) + 180.*(rx_az_21(ind_rx_el)<=0);

% 22
ind_tx_el = find(tx_el_22 > 90);
tx_el_22(ind_tx_el) = 180 - tx_el_22(ind_tx_el);
tx_az_22(ind_tx_el) = tx_az_22(ind_tx_el) + (-180).*(tx_az_22(ind_tx_el)>0) + 180.*(tx_az_22(ind_tx_el)<=0);

ind_tx_el = find(tx_el_22 < -90);
tx_el_22(ind_tx_el) = -(180 + tx_el_22(ind_tx_el));
tx_az_22(ind_tx_el) = tx_az_22(ind_tx_el) + (-180).*(tx_az_22(ind_tx_el)>0) + 180.*(tx_az_22(ind_tx_el)<=0);

ind_rx_el = find(rx_el_22 > 90);
rx_el_22(ind_rx_el) = 180 - rx_el_22(ind_rx_el);
rx_az_22(ind_rx_el) = rx_az_22(ind_rx_el) + (-180).*(rx_az_22(ind_rx_el)>0) + 180.*(rx_az_22(ind_rx_el)<=0);

ind_rx_el = find(rx_el_22 < -90);
rx_el_22(ind_rx_el) = -(180 + rx_el_22(ind_rx_el));
rx_az_22(ind_rx_el) = rx_az_22(ind_rx_el) + (-180).*(rx_az_22(ind_rx_el)>0) + 180.*(rx_az_22(ind_rx_el)<=0);

% output channel structure
ch.am_11 = am_11;
ch.am_12 = am_12;
ch.am_21 = am_21;
ch.am_22 = am_22;

ch.toa_11 = toa_11;
ch.toa_12 = toa_12;
ch.toa_21 = toa_21;
ch.toa_22 = toa_22;

ch.tx_az_11 = tx_az_11;
ch.tx_el_11 = tx_el_11;
ch.rx_az_11 = rx_az_11;
ch.rx_el_11 = rx_el_11;

ch.tx_az_12 = tx_az_12;
ch.tx_el_12 = tx_el_12;
ch.rx_az_12 = rx_az_12;
ch.rx_el_12 = rx_el_12;

ch.tx_az_21 = tx_az_21;
ch.tx_el_21 = tx_el_21;
ch.rx_az_21 = rx_az_21;
ch.rx_el_21 = rx_el_21;

ch.tx_az_22 = tx_az_22;
ch.tx_el_22 = tx_el_22;
ch.rx_az_22 = rx_az_22;
ch.rx_el_22 = rx_el_22;
