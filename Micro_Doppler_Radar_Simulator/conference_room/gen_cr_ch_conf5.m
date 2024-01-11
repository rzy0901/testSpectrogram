% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     gen_cr_ch_mimo.m
%    Authors:       Y. Gagiev
%    Version:       1.0
%    History:       April 2016 created
%
%  *************************************************************************************
%    Description:
%
%    Generates channel for conf#5
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
%      9.  pol             - polarization type
%     10.  lambda          - wavelength
%
%    Outputs:
%
%       1. ch.am_h11       - amplitudes array for direct link (with the same transmit and receive polarizations)
%       2. ch.am_h11       - amplitudes array for cross link (with the different transmit and receive polarizations)
%       2. ch.toa          - times of arrival array
%       3. ch.tx_az        - TX azimuths array
%       4. ch.tx_el        - TX elevations array
%       5. ch.rx_az        - RX azimuths array
%       6. ch.rx_el        - RX elevations array
%
%  *************************************************************************************/
function [ch] = gen_cr_ch_conf5(cfg, pol, lambda)

% vv(11)  vh(12)
% hv(21)  hh(22)

% NLOS attenuation coefficients (due to reflection)
[Gr_h11, Gr_h12] = cr_ref_loss_conf5(cfg.ap_sp, pol);

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
        cls = cr_gen_cls_ant1(0);
            
        % clusters distances
        dist = (cls.toa.*1e-9).*(3e8) + D;
        
        % calculate attenuation constants for clusters
        % attenuation due to propagation
        Gp = sqrt((lambda.^2)./((4.*pi.*(dist)).^2));        
        
        % attenuation constants
        G_h11 = Gp.*Gr_h11;
        G_h12 = Gp.*Gr_h12;
        
        % generate intra-clusters structure
        toa = [];
        am_h11 = [];
        am_h12 = [];
        tx_az = [];
        tx_el = [];
        rx_az = [];
        rx_el = [];
        
        % LOS cluster
        if (cfg.Plos)
            toa = 0;
            am = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            
            tx_pol_1 = polarization(pol(1));
                
            rx_pol_1 = polarization(pol(1));
            rx_pol_2 = polarization(double(~pol(1)));
                
            H = eye(2);
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
            
            pol_coef_h11 = rx_pol_1'*H*tx_pol_1;
            pol_coef_h12 = rx_pol_2'*H*tx_pol_1;
                
            am_h11 = am.*pol_coef_h11;
            am_h12 = am.*pol_coef_h12;
            
            tx_az = 0;
            tx_el = 0;
            rx_az = 0;
            rx_el = 0;
        end
        
        % NLOS clusters
        while (isempty(toa) | (toa == 0))
            for i=1:17
                
                incls = cr_gen_intra_cls(cls.toa(i));
                toa = [toa; cls.toa(i) + incls.toa];
                
                if rand(1,1) <= P(i)                    
                    am_h11 = [am_h11; incls.am.*G_h11(i)];
                    am_h12 = [am_h12; incls.am.*G_h12(i)];
                else
                    am_h11 = [am_h11; incls.am.*G_h11(i).*Gb(i)];
                    am_h12 = [am_h12; incls.am.*G_h12(i).*Gb(i)];
                end
                
                tx_az = [tx_az; cls.tx_az(i) + incls.tx_az];
                tx_el = [tx_el; cls.tx_el(i) + incls.tx_el];
                rx_az = [rx_az; cls.rx_az(i) + incls.rx_az];
                rx_el = [rx_el; cls.rx_el(i) + incls.rx_el];                               
            end
        end
        
    case 1, % STA - AP
        
        % distance between TX and RX
        D = sqrt( (cfg.D).^2 + (1.9).^2 );
        
        % generate NLOS clusters
        cls = cr_gen_cls_ant1(1);
        
        % clusters distances
        dist = (cls.toa.*1e-9).*(3e8) + D;
        
        % calculate attenuation constants for clusters
        % attenuation due to propagation
        Gp = sqrt((lambda.^2)./((4.*pi.*(dist)).^2));
        
        % attenuation constants
        G_h11 = Gp .* Gr_h11;
        G_h12 = Gp .* Gr_h12;
        
        % generate intra-clusters structure
        toa = [];
        am_h11 = [];
        am_h12 = [];
        tx_az = [];
        tx_el = [];
        rx_az = [];
        rx_el = [];
        
        % LOS cluster
        if (cfg.Plos)
            toa = 0;
            am = sqrt((lambda.^2)./((4.*pi.*D).^2)).*exp(2j.*pi.*(D./lambda));
            
            tx_pol_1 = polarization(pol(1));
            
            rx_pol_1 = polarization(pol(1));
            rx_pol_2 = polarization(double(~pol(1)));
                               
            H = eye(2);
            H(1,2) = 0.1.*(2.*(randn(1,1)>0)-1);
            H(2,1) = 0.1.*(2.*(randn(1,1)>0)-1);
                
            pol_coef_h11 = rx_pol_1'*H*tx_pol_1;
            pol_coef_h12 = rx_pol_2'*H*tx_pol_1;
        
            am_h11 = am.*pol_coef_h11;
            am_h12 = am.*pol_coef_h12;
            
            tx_az = 0;
            tx_el = -(asin(1.9./D).*180./pi);
            rx_az = 0;
            rx_el = -tx_el;
        end
        
        % NLOS clusters
        while (isempty(toa) | (toa == 0))
            for i=1:12
                
                incls = cr_gen_intra_cls(cls.toa(i));
                toa = [toa; cls.toa(i) + incls.toa];
                
                if rand(1,1) <= P(i)                    
                    am_h11 = [am_h11; incls.am.*G_h11(i)];
                    am_h12 = [am_h12; incls.am.*G_h12(i)];
                else
                    am_h11 = [am_h11; incls.am.*G_h11(i).*Gb(i)];
                    am_h12 = [am_h12; incls.am.*G_h12(i).*Gb(i)];
                end  
                
                tx_az = [tx_az; cls.tx_az(i) + incls.tx_az];
                tx_el = [tx_el; cls.tx_el(i) + incls.tx_el];
                rx_az = [rx_az; cls.rx_az(i) + incls.rx_az];
                rx_el = [rx_el; cls.rx_el(i) + incls.rx_el];                              
            end
        end
    otherwise,
        error('Prohibited value of "cfg.cr.ap_sp" parameter');
end        

% check azimuth overflow
ind_tx_az = find(tx_az > 180);
tx_az(ind_tx_az) = tx_az(ind_tx_az) - 360;
ind_tx_az = find(tx_az < -180);
tx_az(ind_tx_az) = tx_az(ind_tx_az) + 360;

ind_rx_az = find(rx_az > 180);
rx_az(ind_rx_az) = rx_az(ind_rx_az) - 360;
ind_rx_az = find(rx_az < -180);
rx_az(ind_rx_az) = rx_az(ind_rx_az) + 360;

% check elevation overflow
ind_tx_el = find(tx_el > 90);
tx_el(ind_tx_el) = 180 - tx_el(ind_tx_el);
tx_az(ind_tx_el) = tx_az(ind_tx_el) + (-180).*(tx_az(ind_tx_el)>0) + 180.*(tx_az(ind_tx_el)<=0);

ind_tx_el = find(tx_el < -90);
tx_el(ind_tx_el) = -(180 + tx_el(ind_tx_el));
tx_az(ind_tx_el) = tx_az(ind_tx_el) + (-180).*(tx_az(ind_tx_el)>0) + 180.*(tx_az(ind_tx_el)<=0);

ind_rx_el = find(rx_el > 90);
rx_el(ind_rx_el) = 180 - rx_el(ind_rx_el);
rx_az(ind_rx_el) = rx_az(ind_rx_el) + (-180).*(rx_az(ind_rx_el)>0) + 180.*(rx_az(ind_rx_el)<=0);

ind_rx_el = find(rx_el < -90);
rx_el(ind_rx_el) = -(180 + rx_el(ind_rx_el));
rx_az(ind_rx_el) = rx_az(ind_rx_el) + (-180).*(rx_az(ind_rx_el)>0) + 180.*(rx_az(ind_rx_el)<=0);

% output channel structure
ch.am_h11 = am_h11;
ch.am_h12 = am_h12;
ch.toa = toa;
ch.tx_az = tx_az;
ch.tx_el = tx_el;
ch.rx_az = rx_az;
ch.rx_el = rx_el;
