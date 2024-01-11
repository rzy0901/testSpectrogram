function Simulator(v,ped_height,rho,Pedestrian_heading,Pedestrian_pos,Radar_x,Radar_y,Radar_z,Radar_pos)
        % add path
        rootdir=pwd;
        if strfind(rootdir,'\')
            subdir=strcat(rootdir,'\beamforming'); addpath(subdir);
            subdir=strcat(rootdir,'\common'); addpath(subdir);
            subdir=strcat(rootdir,'\conference_room'); addpath(subdir);
            subdir=strcat(rootdir,'\ray_tracer'); addpath(subdir);
            subdir=strcat(rootdir,'\work'); addpath(subdir);
            subdir=strcat(rootdir,'\phased\phased'); addpath(subdir);
        else %unix
            subdir=strcat(rootdir,'/beamforming'); addpath(subdir);
            subdir=strcat(rootdir,'/common'); addpath(subdir);
            subdir=strcat(rootdir,'/conference_room'); addpath(subdir);
            subdir=strcat(rootdir,'/ray_tracer'); addpath(subdir);
            subdir=strcat(rootdir,'/work'); addpath(subdir);
        end  

        % set parameters
        bw = 10e6; % bw for FMCW, 10 MHZ
        fs = 1*bw; % sampling rate 10 MHZ
        fc = 3.5e9; % carrier frequency
        tm = 10e-6; % 10 us
        c = 3e8; % light speed

        wav = phased.FMCWWaveform('SampleRate',fs,'SweepTime',tm,...
            'SweepBandwidth',bw); % transmit signal
        wave = step(wav);

        tx_pos = Radar_pos; % on the desk
        tx_vel = [0;0;0]; % fixed
        radar_tx = phased.Platform('InitialPosition',tx_pos,'Velocity',tx_vel,...
            'OrientationAxesOutputPort',true); % generate tx class
        
        ped = backscatterPedestrian('InitialPosition',Pedestrian_pos,'InitialHeading',Pedestrian_heading,...
            'PropagationSpeed',c,'OperatingFrequency',fc,'Height',ped_height,'WalkingSpeed',v); % generate the person

        chan_ped = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,...
            'TwoWayPropagation',true,'SampleRate',fs); % generate free space channel

        tx = phased.Transmitter('PeakPower',1,'Gain',25); % generate tx
        rx = phased.ReceiverPreamp('Gain',25,'NoiseFigure',10); % generate rx

        Tsamp = 0.001; % duration for each pulse
        npulse = 3000; % 5number of pluses, total time 3000 * 0.001 = 3.0 seconds
        xr_ped = complex(zeros(round(fs*tm),npulse));
        xr_ped_perfect = complex(zeros(round(fs*tm),npulse));

%         % generate DAHC model
%         num = 1;
%         freq = 1; % 1ms, the changing speed of static channel
%         for m = 1:npulse+1
%             if mod(m, freq) == 0
%                 [imp_res] = cr_ch_model(Radar_x,Radar_y,Radar_z); % see cr_ch_cfg.m for detailed configuration
%                 channels(:, m) = imp_res.h11(1:10);
%                 num = num + 1;
%             end
%         end
% 
%         corr_channels = zeros(10, npulse+1);
%         corr_channels(:,1) = channels(:,1);
%         for m = 1:npulse
%             if mod(m, freq) == 0
%                 corr_channels(:, m+1) = rho*corr_channels(:, m) + (1-rho)*channels(:, m+1);
%             else
%                 corr_channels(:, m+1) = corr_channels(:, m);
%             end
%         end
%         save(['Channel\','Un_Ch_','_',num2str(rho),'.mat'],'corr_channels');

        % load generated DAHC model
        load(['./Channel/','Un_Ch_',num2str(rho),'.mat']);
        
        temp = reshape(corr_channels, [10*(npulse+1), 1]);
        coin = randi([2, 10]);
        shift_corr_channels = reshape(temp(1:10*npulse), [10,npulse]);
        for m = 1:npulse
                
            [pos_tx, vel_tx, ~] = radar_tx(Tsamp); % tx position, velocity, and angle
            [pos_ped, vel_ped, ax_ped] = move(ped,Tsamp,Pedestrian_heading); % the person moves

            plot(ped);
            xlabel('X(m)','Fontsize',14);
            ylabel('Y(m)','Fontsize',14);
            zlabel('Z(m)','Fontsize',14);
            title('Pedestrian Trajectory','Fontsize',14);
            
            angrt_ped = zeros(2,size(pos_ped,2));
            for k = 1:size(pos_ped,2)
                [range,angle] = rangeangle(pos_tx,pos_ped(:,k),ax_ped(:,:,k)); % compute the reflection angle
                angrt_ped(:,k) = angle;
            end
            x = tx(wav()); % generate the tx signal

            % person
            xt_ped = chan_ped(repmat(x,1,size(pos_ped,2)),pos_tx,pos_ped,vel_tx,vel_ped); % channel to the person
            xt_ped = reflect(ped,xt_ped,angrt_ped); % reflection from the person

            xt_wall = conv(x, shift_corr_channels(:,m)); % the received signal due to the conference room  
            xr_ped(:,m) = rx(xt_ped + xt_wall(1:round(fs*tm))); % superposition signal
            xr_ped_perfect(:,m) = rx(xt_ped); % perfect case without the environment
        end
        
        xd_ped = conj(dechirp(xr_ped,x)); % compute the correlation between the received and transmit signals
        xd_ped_perfect = conj(dechirp(xr_ped_perfect,x)); % perfect case

        % figure 1: perfect case
        clf;
        spectrogram(sum(xd_ped_perfect),kaiser(128,10),120,256,1/Tsamp,'centered','yaxis');
        colormap('jet');
        clim = get(gca,'CLim');
        set(gca,'CLim',clim(2)+[-50 0]);
        set(gca,'YLim',[-500, 500]);
       
        % figure 2: imperfect case
        clf;
        spectrogram(sum(xd_ped),kaiser(128,10),120,256,1/Tsamp,'centered','yaxis');
        colormap('jet');
        clim = get(gca,'CLim');
        set(gca,'CLim',clim(2)+[-50 0])

        % figure 3: principal component analysis
        [uxd,sxd,vxd] = svd(xr_ped);
        clf
        plot(10*log10(diag(sxd)));
        xlabel('Rank','Fontsize',14);
        ylabel('Singular Values','Fontsize',14);
        set(gca,'XLim',[0 50]);
        hold on;
        plot([2 2],[-40 10],'r--');
        plot([18 18],[-40 10],'r--');
        text(0.5,-10,'A');
        text(10,-10,'B');
        text(30,-10,'C');

        % figure 4: time-frequency after removing the impact from conference room
        rk = 2:18;
        xdr = uxd(:,rk)*sxd(rk,:)*vxd';
        clf
        [S,F,T,P] = spectrogram(sum(xdr),kaiser(128,10),120,256,1/Tsamp,'centered','yaxis');
        xlabel('Time (secs)','Fontsize',14);
        ylabel('Frequency (Hz)','Fontsize',14);
        colormap('jet');
        clim = get(gca,'CLim');
        set(gca,'CLim',clim(2)+[-50 0])
end