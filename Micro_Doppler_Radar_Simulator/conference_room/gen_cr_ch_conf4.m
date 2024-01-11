% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     gen_cr_ch_mimo.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       November 2015 created
%
%  *************************************************************************************
%    Description:
%
%    Generates channel for conf#4
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
%      9.  dist_tx         - distance between PAAs on Tx side
%      10. dist_rx         - distance between PAAs on Rx side
%      11. lambda          - wavelength
%
%    Outputs:
%
%      1.  ch.am_11vv      - amplitudes array for link from 1st PAA (TX) to 1st PAA (RX) for direct polarization
%      2.  ch.am_11vh      - amplitudes array for link from 1st PAA (TX) to 1st PAA (RX) for cross polarization
%      3.  ch.am_11hv      - amplitudes array for link from 1st PAA (TX) to 1st PAA (RX) for cross polarization
%      4.  ch.am_11hh      - amplitudes array for link from 1st PAA (TX) to 1st PAA (RX) for direct polarization
%
%      5.  ch.am_12vv      - amplitudes array for link from 1st PAA (TX) to 2nd PAA (RX) for direct polarization
%      6.  ch.am_12vh      - amplitudes array for link from 1st PAA (TX) to 2nd PAA (RX) for cross polarization
%      7.  ch.am_12hv      - amplitudes array for link from 1st PAA (TX) to 2nd PAA (RX) for cross polarization
%      8.  ch.am_12hh      - amplitudes array for link from 1st PAA (TX) to 2nd PAA (RX) for direct polarization
%
%      9.  ch.am_21vv      - amplitudes array for link from 2nd PAA (TX) to 1st PAA (RX) for direct polarization
%      10. ch.am_21vh      - amplitudes array for link from 2nd PAA (TX) to 1st PAA (RX) for cross polarization
%      11. ch.am_21hv      - amplitudes array for link from 2nd PAA (TX) to 1st PAA (RX) for cross polarization
%      12. ch.am_21hh      - amplitudes array for link from 2nd PAA (TX) to 1st PAA (RX) for direct polarization
%
%      13. ch.am_22vv      - amplitudes array for link from 2nd PAA (TX) to 2nd PAA (RX) for direct polarization
%      14. ch.am_22vh      - amplitudes array for link from 2nd PAA (TX) to 2nd PAA (RX) for cross polarization
%      15. ch.am_22hv      - amplitudes array for link from 2nd PAA (TX) to 2nd PAA (RX) for cross polarization
%      16. ch.am_22hh      - amplitudes array for link from 2nd PAA (TX) to 2nd PAA (RX) for direct polarization
%
%      17. ch.toa_11       - times of arrival array for link from 1st PAA (TX) to 1st PAA (RX)
%      18. ch.toa_12       - times of arrival array for link from 1st PAA (TX) to 2nd PAA (RX)
%      19. ch.toa_21       - times of arrival array for link from 2nd PAA (TX) to 1st PAA (RX)
%      20. ch.toa_22       - times of arrival array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      
%      21. ch.tx_az_11     - TX azimuths array for link from 1st PAA (TX) to 1st PAA (RX)
%      22. ch.tx_el_11     - TX elevations array for link from 1st PAA (TX) to 1st PAA (RX)
%      23. ch.rx_az_11     - RX azimuths array for link from 1st PAA (TX) to 1st PAA (RX)
%      24. ch.rx_el_11     - RX elevations array for link from 1st PAA (TX) to 1st PAA (RX)
%      
%      25. ch.tx_az_12     - TX azimuths array for link from 1st PAA (TX) to 2nd PAA (RX)
%      26. ch.tx_el_12     - TX elevations array for link from 1st PAA (TX) to 2nd PAA (RX)
%      27. ch.rx_az_12     - RX azimuths array for link from 1st PAA (TX) to 2nd PAA (RX)
%      28. ch.rx_el_12     - RX elevations array for link from 1st PAA (TX) to 2nd PAA (RX)
%      
%      29. ch.tx_az_21     - TX azimuths array for link from 2nd PAA (TX) to 1st PAA (RX)
%      30. ch.tx_el_21     - TX elevations array for link from 2nd PAA (TX) to 1st PAA (RX)
%      31. ch.rx_az_21     - RX azimuths array for link from 2nd PAA (TX) to 1st PAA (RX)
%      32. ch.rx_el_21     - RX elevations array for link from 2nd PAA (TX) to 1st PAA (RX)
%      
%      33. ch.tx_az_22     - TX azimuths array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      34. ch.tx_el_22     - TX elevations array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      35. ch.rx_az_22     - RX azimuths array for link from 2nd PAA (TX) to 2nd PAA (RX)
%      36. ch.rx_el_22     - RX elevations array for link from 2nd PAA (TX) to 2nd PAA (RX)
%
%  *************************************************************************************/
function [ch] = gen_cr_ch_conf4(cfg, dist_tx, dist_rx, phi_tx, phi_rx, lambda)

% NLOS attenuation coefficients (due to reflection)
[Gr_vv, Gr_vh, Gr_hv, Gr_hh] = cr_ref_loss_conf2(cfg.ap_sp);

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
        G_11vv = Gp_11.*Gr_vv;
        G_11vh = Gp_11.*Gr_vh;
        G_11hv = Gp_11.*Gr_hv;
        G_11hh = Gp_11.*Gr_hh;
        
        G_12vv = Gp_12.*Gr_vv;
        G_12vh = Gp_12.*Gr_vh;
        G_12hv = Gp_12.*Gr_hv;
        G_12hh = Gp_12.*Gr_hh;
        
        G_21vv = Gp_21.*Gr_vv;
        G_21vh = Gp_21.*Gr_vh;
        G_21hv = Gp_21.*Gr_hv;
        G_21hh = Gp_21.*Gr_hh;
        
        G_22vv = Gp_22.*Gr_vv;
        G_22vh = Gp_22.*Gr_vh;
        G_22hv = Gp_22.*Gr_hv;
        G_22hh = Gp_22.*Gr_hh;
        
        % generate intra-clusters structure
        % From Tx1 to Rx1
        toa_11 = [];
        am_11vv = [];
        am_11vh = [];
        am_11hv = [];
        am_11hh = [];
        tx_az_11 = [];
        tx_el_11 = [];
        rx_az_11 = [];
        rx_el_11 = [];
        
        % From Tx1 to Rx2
        toa_12 = [];
        am_12vv = [];
        am_12vh = [];
        am_12hv = [];
        am_12hh = [];
        tx_az_12 = [];
        tx_el_12 = [];
        rx_az_12 = [];
        rx_el_12 = [];

        % From Tx2 to Rx1
        toa_21 = [];
        am_21vv = [];
        am_21vh = [];
        am_21hv = [];
        am_21hh = [];
        tx_az_21 = [];
        tx_el_21 = [];
        rx_az_21 = [];
        rx_el_21 = [];
        
        % From Tx2 to Rx2
        toa_22 = [];
        am_22vv = [];
        am_22vh = [];
        am_22hv = [];
        am_22hh = [];
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
            am_11vv = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            am_11vh = am_11vv;
            am_11hv = am_11vv;
            am_11hh = am_11vv;
            
            am_12vv = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_12vh = am_12vv;
            am_12hv = am_12vv;
            am_12hh = am_12vv;

            am_21vv = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_21vh = am_21vv;
            am_21hv = am_21vv;
            am_21hh = am_21vv;

            am_22vv = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            am_22vh = am_22vv;
            am_22hv = am_22vv;
            am_22hh = am_22vv;
            
            % Apply polarization
            tx_pol_v = polarization(0);
            tx_pol_h = polarization(1);
                
            rx_pol_v = polarization(0);
            rx_pol_h = polarization(1);
            
            H = eye(2);
            % For Tx1 to Rx1
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
               
            pol_coef_vv = rx_pol_v'*H*tx_pol_v;
            pol_coef_vh = rx_pol_v'*H*tx_pol_h;
            pol_coef_hv = rx_pol_h'*H*tx_pol_v;
            pol_coef_hh = rx_pol_h'*H*tx_pol_h;
               
            am_11vv = am_11vv.*pol_coef_vv;
            am_11vh = am_11vh.*pol_coef_vh;
            am_11hv = am_11hv.*pol_coef_hv;
            am_11hh = am_11hh.*pol_coef_hh;
            
            % For Tx1 to Rx2
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            am_12vv = am_12vv.*pol_coef_vv;
            am_12vh = am_12vh.*pol_coef_vh;
            am_12hv = am_12hv.*pol_coef_hv;
            am_12hh = am_12hh.*pol_coef_hh;

            % For Tx2 to Rx1
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            am_21vv = am_21vv.*pol_coef_vv;
            am_21vh = am_21vh.*pol_coef_vh;
            am_21hv = am_21hv.*pol_coef_hv;
            am_21hh = am_21hh.*pol_coef_hh;

            % For Tx2 to Rx2
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            am_22vv = am_22vv.*pol_coef_vv;
            am_22vh = am_22vh.*pol_coef_vh;
            am_22hv = am_22hv.*pol_coef_hv;
            am_22hh = am_22hh.*pol_coef_hh;
            
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
                    am_11vv = [am_11vv; incls_11.am.*G_11vv(i)];
                    am_11vh = [am_11vh; incls_11.am.*G_11vh(i)];
                    am_11hv = [am_11hv; incls_11.am.*G_11hv(i)];
                    am_11hh = [am_11hh; incls_11.am.*G_11hh(i)];
                else
                    am_11vv = [am_11vv; incls_11.am.*G_11vv(i).*Gb(i)];
                    am_11vh = [am_11vh; incls_11.am.*G_11vh(i).*Gb(i)];
                    am_11hv = [am_11hv; incls_11.am.*G_11hv(i).*Gb(i)];
                    am_11hh = [am_11hh; incls_11.am.*G_11hh(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_12vv = [am_12vv; incls_12.am.*G_12vv(i)];
                    am_12vh = [am_12vh; incls_12.am.*G_12vh(i)];
                    am_12hv = [am_12hv; incls_12.am.*G_12hv(i)];
                    am_12hh = [am_12hh; incls_12.am.*G_12hh(i)];
                else
                    am_12vv = [am_12vv; incls_12.am.*G_12vv(i).*Gb(i)];
                    am_12vh = [am_12vh; incls_12.am.*G_12vh(i).*Gb(i)];
                    am_12hv = [am_12hv; incls_12.am.*G_12hv(i).*Gb(i)];
                    am_12hh = [am_12hh; incls_12.am.*G_12hh(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_21vv = [am_21vv; incls_21.am.*G_21vv(i)];
                    am_21vh = [am_21vh; incls_21.am.*G_21vh(i)];
                    am_21hv = [am_21hv; incls_21.am.*G_21hv(i)];
                    am_21hh = [am_21hh; incls_21.am.*G_21hh(i)];
                else
                    am_21vv = [am_21vv; incls_21.am.*G_21vv(i).*Gb(i)];
                    am_21vh = [am_21vh; incls_21.am.*G_21vh(i).*Gb(i)];
                    am_21hv = [am_21hv; incls_21.am.*G_21hv(i).*Gb(i)];
                    am_21hh = [am_21hh; incls_21.am.*G_21hh(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_22vv = [am_22vv; incls_22.am.*G_22vv(i)];
                    am_22vh = [am_22vh; incls_22.am.*G_22vh(i)];
                    am_22hv = [am_22hv; incls_22.am.*G_22hv(i)];
                    am_22hh = [am_22hh; incls_22.am.*G_22hh(i)];
                else
                    am_22vv = [am_22vv; incls_22.am.*G_22vv(i).*Gb(i)];
                    am_22vh = [am_22vh; incls_22.am.*G_22vh(i).*Gb(i)];
                    am_22hv = [am_22hv; incls_22.am.*G_22hv(i).*Gb(i)];
                    am_22hh = [am_22hh; incls_22.am.*G_22hh(i).*Gb(i)];
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
        G_11vv = Gp_11.*Gr_vv;
        G_11vh = Gp_11.*Gr_vh;
        G_11hv = Gp_11.*Gr_hv;
        G_11hh = Gp_11.*Gr_hh;
        
        G_12vv = Gp_12.*Gr_vv;
        G_12vh = Gp_12.*Gr_vh;
        G_12hv = Gp_12.*Gr_hv;
        G_12hh = Gp_12.*Gr_hh;
        
        G_21vv = Gp_21.*Gr_vv;
        G_21vh = Gp_21.*Gr_vh;
        G_21hv = Gp_21.*Gr_hv;
        G_21hh = Gp_21.*Gr_hh;
        
        G_22vv = Gp_22.*Gr_vv;
        G_22vh = Gp_22.*Gr_vh;
        G_22hv = Gp_22.*Gr_hv;
        G_22hh = Gp_22.*Gr_hh;
        
        % generate intra-clusters structure
        toa_11 = [];
        toa_12 = [];
        toa_21 = [];
        toa_22 = [];
        am_11vv = [];
        am_11vh = [];
        am_11hv = [];
        am_11hh = [];
        
        am_12vv = [];
        am_12vh = [];
        am_12hv = [];
        am_12hh = [];
        
        am_21vv = [];
        am_21vh = [];
        am_21hv = [];
        am_21hh = [];
        
        am_22vv = [];
        am_22vh = [];
        am_22hv = [];
        am_22hh = [];
        
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
            
            am_11vv = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            am_11vh = am_11vv;
            am_11hv = am_11vv;
            am_11hh = am_11vv;
            
            am_12vv = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_12vh = am_12vv;
            am_12hv = am_12vv;
            am_12hh = am_12vv;
            
            am_21vv = sqrt((lambda.^2)./((4.*pi.*sqrt(D^2 + dist_tx * dist_rx)).^2)).*exp(2j.*pi.*(sqrt(D^2 + dist_tx * dist_rx)./lambda));
            am_21vh = am_21vv;
            am_21hv = am_21vv;
            am_21hh = am_21vv;
            
            am_22vv = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            am_22vh = am_22vv;
            am_22hv = am_22vv;
            am_22hh = am_22vv;
            
            tx_pol_v = polarization(0);
            tx_pol_h = polarization(1);
                
            rx_pol_v = polarization(0);
            rx_pol_h = polarization(1);
                
            H = eye(2);
            % From Tx1 -> Rx1
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
                
            pol_coef_vv = rx_pol_v'*H*tx_pol_v;
            pol_coef_vh = rx_pol_v'*H*tx_pol_h;
            pol_coef_hv = rx_pol_h'*H*tx_pol_v;
            pol_coef_hh = rx_pol_h'*H*tx_pol_h;
                
            am_11vv = am_11vv.*pol_coef_vv;
            am_11vh = am_11vh.*pol_coef_vh;
            am_11hv = am_11hv.*pol_coef_hv;
            am_11hh = am_11hh.*pol_coef_hh;

            % From Tx1 -> Rx2
            am_12vv = am_12vv.*pol_coef_vv;
            am_12vh = am_12vh.*pol_coef_vh;
            am_12hv = am_12hv.*pol_coef_hv;
            am_12hh = am_12hh.*pol_coef_hh;

            % From Tx2 -> Rx1
            am_21vv = am_21vv.*pol_coef_vv;
            am_21vh = am_21vh.*pol_coef_vh;
            am_21hv = am_21hv.*pol_coef_hv;
            am_21hh = am_21hh.*pol_coef_hh;
            
            % From Tx2 -> Rx2
            am_22vv = am_22vv.*pol_coef_vv;
            am_22vh = am_22vh.*pol_coef_vh;
            am_22hv = am_22hv.*pol_coef_hv;
            am_22hh = am_22hh.*pol_coef_hh;
            
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
                    am_11vv = [am_11vv; incls_11.am.*G_11vv(i)];
                    am_11vh = [am_11vh; incls_11.am.*G_11vh(i)];
                    am_11hv = [am_11hv; incls_11.am.*G_11hv(i)];
                    am_11hh = [am_11hh; incls_11.am.*G_11hh(i)];
                else
                    am_11vv = [am_11vv; incls_11.am.*G_11vv(i).*Gb(i)];
                    am_11vh = [am_11vh; incls_11.am.*G_11vh(i).*Gb(i)];
                    am_11hv = [am_11hv; incls_11.am.*G_11hv(i).*Gb(i)];
                    am_11hh = [am_11hh; incls_11.am.*G_11hh(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_12vv = [am_12vv; incls_12.am.*G_12vv(i)];
                    am_12vh = [am_12vh; incls_12.am.*G_12vh(i)];
                    am_12hv = [am_12hv; incls_12.am.*G_12hv(i)];
                    am_12hh = [am_12hh; incls_12.am.*G_12hh(i)];
                else
                    am_12vv = [am_12vv; incls_12.am.*G_12vv(i).*Gb(i)];
                    am_12vh = [am_12vh; incls_12.am.*G_12vh(i).*Gb(i)];
                    am_12hv = [am_12hv; incls_12.am.*G_12hv(i).*Gb(i)];
                    am_12hh = [am_12hh; incls_12.am.*G_12hh(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_21vv = [am_21vv; incls_21.am.*G_21vv(i)];
                    am_21vh = [am_21vh; incls_21.am.*G_21vh(i)];
                    am_21hv = [am_21hv; incls_21.am.*G_21hv(i)];
                    am_21hh = [am_21hh; incls_21.am.*G_21hh(i)];
                else
                    am_21vv = [am_21vv; incls_21.am.*G_21vv(i).*Gb(i)];
                    am_21vh = [am_21vh; incls_21.am.*G_21vh(i).*Gb(i)];
                    am_21hv = [am_21hv; incls_21.am.*G_21hv(i).*Gb(i)];
                    am_21hh = [am_21hh; incls_21.am.*G_21hh(i).*Gb(i)];
                end
                
                if rand(1,1) <= P(i)                    
                    am_22vv = [am_22vv; incls_22.am.*G_22vv(i)];
                    am_22vh = [am_22vh; incls_22.am.*G_22vh(i)];
                    am_22hv = [am_22hv; incls_22.am.*G_22hv(i)];
                    am_22hh = [am_22hh; incls_22.am.*G_22hh(i)];
                else
                    am_22vv = [am_22vv; incls_22.am.*G_22vv(i).*Gb(i)];
                    am_22vh = [am_22vh; incls_22.am.*G_22vh(i).*Gb(i)];
                    am_22hv = [am_22hv; incls_22.am.*G_22hv(i).*Gb(i)];
                    am_22hh = [am_22hh; incls_22.am.*G_22hh(i).*Gb(i)];
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
ch.am_11vv = am_11vv;
ch.am_11vh = am_11vh;
ch.am_11hv = am_11hv;
ch.am_11hh = am_11hh;

ch.am_12vv = am_12vv;
ch.am_12vh = am_12vh;
ch.am_12hv = am_12hv;
ch.am_12hh = am_12hh;

ch.am_21vv = am_21vv;
ch.am_21vh = am_21vh;
ch.am_21hv = am_21hv;
ch.am_21hh = am_21hh;

ch.am_22vv = am_22vv;
ch.am_22vh = am_22vh;
ch.am_22hv = am_22hv;
ch.am_22hh = am_22hh;

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
