clear; clc; close all;
load('./data_18/lqr/squat/squat_1.mat')
% 舍弃第一帧和最后一帧.
keypoints([1,end-1],:,:)=[]; % keypoints in XYZ
timestampList([1,end-1])=[];
timestampList = timestampList - timestampList(1);
Njoints = size(keypoints,2);
joints_18 = ["NOSE","NECK","RIGHT_SHOULDER","RIGHT_ELBOW","RIGHT_WRIST","LEFT_SHOULDER","LEFT_ELBOW","LEFT_WRIST","RIGHT_HIP","RIGHT_KNEE","RIGHT_ANKLE","LEFT_HIP","LEFT_KNEE","LEFT_ANKLE","RIGHT_EYE","LEFT_EYE","RIGHT_EAR","LEFT_EAR"];
joints_34 = ["PELVIS","NAVAL_SPINE","CHEST_SPINE","NECK","LEFT_CLAVICLE","LEFT_SHOULDER","LEFT_ELBOW","LEFT_WRIST","LEFT_HAND","LEFT_HANDTIP","LEFT_THUMB","RIGHT_CLAVICLE","RIGHT_SHOULDER","RIGHT_ELBOW","RIGHT_WRIST","RIGHT_HAND","RIGHT_HANDTIP","RIGHT_THUMB","LEFT_HIP","LEFT_KNEE","LEFT_ANKLE","LEFT_FOOT","RIGHT_HIP","RIGHT_KNEE","RIGHT_ANKLE","RIGHT_FOOT","HEAD","NOSE","LEFT_EYE","LEFT_EAR","RIGHT_EYE","RIGHT_EAR","LEFT_HEEL","RIGHT_HEEL"];
joints_70 = ["PELVIS","SPINE_1","SPINE_2","SPINE_3","NECK","NOSE","LEFT_EYE","RIGHT_EYE","LEFT_EAR","RIGHT_EAR","LEFT_CLAVICLE","RIGHT_CLAVICLE","LEFT_SHOULDER","RIGHT_SHOULDER","LEFT_ELBOW","RIGHT_ELBOW","LEFT_WRIST","RIGHT_WRIST","LEFT_HIP","RIGHT_HIP","LEFT_KNEE","RIGHT_KNEE","LEFT_ANKLE","RIGHT_ANKLE","LEFT_BIG_TOE","RIGHT_BIG_TOE","LEFT_SMALL_TOE","RIGHT_SMALL_TOE","LEFT_HEEL","RIGHT_HEEL","LEFT_HAND_THUMB_1","LEFT_HAND_THUMB_2","LEFT_HAND_THUMB_3","LEFT_HAND_THUMB_4","LEFT_HAND_INDEX_1","LEFT_HAND_INDEX_2","LEFT_HAND_INDEX_3","LEFT_HAND_INDEX_4","LEFT_HAND_MIDDLE_1","LEFT_HAND_MIDDLE_2","LEFT_HAND_MIDDLE_3","LEFT_HAND_MIDDLE_4","LEFT_HAND_RING_1","LEFT_HAND_RING_2","LEFT_HAND_RING_3","LEFT_HAND_RING_4","LEFT_HAND_PINKY_1","LEFT_HAND_PINKY_2","LEFT_HAND_PINKY_3","LEFT_HAND_PINKY_4","RIGHT_HAND_THUMB_1","RIGHT_HAND_THUMB_2","RIGHT_HAND_THUMB_3","RIGHT_HAND_THUMB_4","RIGHT_HAND_INDEX_1","RIGHT_HAND_INDEX_2","RIGHT_HAND_INDEX_3","RIGHT_HAND_INDEX_4","RIGHT_HAND_MIDDLE_1","RIGHT_HAND_MIDDLE_2","RIGHT_HAND_MIDDLE_3","RIGHT_HAND_MIDDLE_4","RIGHT_HAND_RING_1","RIGHT_HAND_RING_2","RIGHT_HAND_RING_3","RIGHT_HAND_RING_4","RIGHT_HAND_PINKY_1","RIGHT_HAND_PINKY_2","RIGHT_HAND_PINKY_3","RIGHT_HAND_PINKY_4"];
if Njoints ==34
    connections = [1 2; 2 3; 3 5; 5 6; 6 7; 7 8; 8 9;9 10;8 11;3 12;12 13;13 14;14 15;15 16;16 17;15 18;1 19;19 20;20 21;21 22;1 23;23 24;24 25;25 26;3 4;4 27;27 28;28 29;29 30;28 31;31 32;21 33;25 34;33 22;34 26];
elseif Njoints == 18
    connections = [1 2;2 3;3 4;4 5;2 6;6 7;7 8;3 9;9 10;10 11;6 12;12 13;13 14;3 6;9 12;1 15;15 17;1 16;16 18];
elseif Njoints == 70
    connections = [1 2;2 3;3 4;4 5;5 6;6 7;7 9;6 8;8 10;4 11;11 13;13 15;15 17;17 31;31 32;32 33;33 34;17 35;35 36;36 37;37 38;17 39;39 40;40 41;41 42;17 43;43 44;44 45;45 46;17 47;47 48;48 49;49 50;4 12;12 14;14 16;16 18;18 51;51 52;52 53;53 54;18 55;55 56;56 57;57 58;18 59;59 60;60 61;61 62;18 63;63 64;64 65;65 66;18 67;67 68;68 69;69 70;1 19;19 21;21 23;23 29;23 25;23 27;1 20;20 22;22 24;24 30;24 26;24 28];
end
Nframes = length(timestampList);
frameLength = 1/30; % fps =  30;
T = frameLength*Nframes;
Tx_pos = [0 -0.1 -1.5]; % XYZ
Rx_pos = [0 -0.1 0]; % XYZ
drawScenario = true;
%% plot
if drawScenario == true
    hf = figure;
    hf.Color = 'white';
    for ii = 1:1:length(timestampList)
        cla
        x = squeeze(keypoints(ii,:,1));
        y = squeeze(keypoints(ii,:,2));
        z = squeeze(keypoints(ii,:,3));
        % plot
        human = plot3(z,x,y,'.','markersize', 13,'Color',"blue");
        hold on;
        axis equal;
        xmin = min([0 Tx_pos(1) Rx_pos(1) min(keypoints(:,:,1),[],'all')]);
        xmax = max([0 Tx_pos(1) Rx_pos(1) max(keypoints(:,:,1),[],'all')]);
        ymin = min([0 Tx_pos(2) Rx_pos(2) min(keypoints(:,:,2),[],'all')]);
        ymax = max([0 Tx_pos(2) Rx_pos(2) max(keypoints(:,:,2),[],'all')]);
        zmin = min([0 Tx_pos(3) Rx_pos(3) min(keypoints(:,:,3),[],'all')]);
        zmax = max([0 Tx_pos(3) Rx_pos(3) max(keypoints(:,:,3),[],'all')]);
        xlim([zmin zmax]); % Z
        ylim([xmin xmax]); % X
        zlim([ymin ymax]); % Y
%         view(30,30)
%         view(2)
%         view(30,5)
        camera = scatter3(0,0,0,200,"red",'o');
        tx = scatter3(Tx_pos(3),Tx_pos(1),Tx_pos(2),200,"magenta",'*');
        rx = scatter3(Rx_pos(3),Rx_pos(1),Rx_pos(2),200,"black",'*');
        for jj = 1:1:size(connections,1)
            plot3(z(connections(jj,:)),x(connections(jj,:)),y(connections(jj,:)),'Color','b','LineWidth',0.05);
        end
        for nj=1:Njoints        
            line([Tx_pos(3) z(nj)],[Tx_pos(1) x(nj)],...
                           [Tx_pos(2) y(nj)],'LineStyle',':',...
                           'color',[0.5 0.5 0.5],'LineWidth',0.05)
            line([Rx_pos(3) z(nj)],[Rx_pos(1) x(nj)],...
                           [Rx_pos(2) y(nj)],'LineStyle',':',...
                           'color',[0.5 0.5 0.5],'LineWidth',0.05)
        end    
        xlabel('Z(m)'); ylabel('X(m)'); zlabel('Y(m)'); title(sprintf('Timestamp: %f (ms)',timestampList(ii)));
        grid on;
        legend([camera tx rx human] ,'camera','transmitter','receiver','human');
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%         set(gca,'ZDir','reverse'); %Y
        drawnow;
        Frame=getframe(gcf);
        Image=frame2im(Frame);
        [Image,map]=rgb2ind(Image,256);
        if ii == 1
            imwrite(Image,map,'test.gif','gif', 'Loopcount',inf,'DelayTime',0.033);
        else
            imwrite(Image,map,'test.gif','gif','WriteMode','append','DelayTime',0.033);
        end
    end
end
%% Spectrogram
% Interpolation of the data:
fs = 1000; % new frame rate
TimeSamples = linspace(0,T,Nframes);
NframesNew = round(T*fs); % Number of frame after interpolation
TimeSamplesNew = linspace(0,T,NframesNew);
keypointsNew = zeros(length(TimeSamplesNew),Njoints,3);
for j=1:Njoints
    for k=1:3
        keypointsNew(:,j,k) = interp1(TimeSamples, keypoints(:,j,k),...
            TimeSamplesNew,'spline','extrap');
    end
end
% Calculate Radar returns from target
% Radar parameters
c = 3e8; % m/s
fc = 7.5e9; % carrier frequency
lambda = c/fc; %(m) wavelength
for nf = 1:NframesNew
    rcs = 0;
    for jj = 1:1:size(connections,1)
        joint1(1:3) = keypointsNew(nf,connections(jj,1),1:3);
        joint2(1:3) = keypointsNew(nf,connections(jj,2),1:3);
        % origin of constructed ellipsoid
        mid = 0.5*(joint1+joint2);
        R_Tx = norm(mid-Tx_pos);
        R_Rx = norm(mid-Rx_pos);
        % aspect vector
        aspect = joint1 - joint2;
        % semi-axis length
        a = norm(aspect)/4;
        b = norm(aspect)/4;
        c = norm(aspect)/2;
        % Calculate theta
        Cos_Theta_i = dot(Tx_pos-mid,aspect)/norm(mid-Tx_pos)/norm(aspect);
        Theta_i = acos(Cos_Theta_i);
        Cos_Theta_s = dot(Rx_pos-mid,aspect)/norm(mid-Rx_pos)/norm(aspect);
        Theta_s = acos(Cos_Theta_s);
        % Calculate phi
        Sin_Phi_i = (Tx_pos(2) - mid(2))/sqrt((Tx_pos(1)-mid(1))^2+(Tx_pos(2)-mid(2))^2);
        Phi_i = asin(Sin_Phi_i);
        Sin_Phi_s = (Rx_pos(2) - mid(2))/sqrt((Rx_pos(1)-mid(1))^2+(Rx_pos(2)-mid(2))^2);
        Phi_s = asin(Sin_Phi_s);
        % rcsellipsoid/R^2 is based on monostatic radar range equation
        Amp = rcsellipsoid(a,b,c,Phi_i,Theta_i,Phi_s,Theta_s)/R_Rx/R_Tx;
        Phase = exp(-1i*4*pi*(R_Tx+R_Rx)/lambda); 
        rcs_joint = Amp*Phase;
        rcs = rcs + rcs_joint;
    end
    RCS(nf) = rcs;
end

% micro-Doppler signature
F = fs;
figure;% figure('Position',[500 200 900 600])
colormap(jet)
spectrogram(RCS,kaiser(256,15),250,512,F,'centered','yaxis');
clim = get(gca,'CLim');
set(gca,'CLim',clim(2) + [-60 0]);
title('Micro-Doppler Signature', 'Fontsize',12,'color','k')
drawnow