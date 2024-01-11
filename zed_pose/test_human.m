clear; clc; close all;
person = 'lqr';% rzy
Radar_pos = [-3 0 -4]; % XYZ
drawScenario = true;
input_paths = sprintf('./data_34/%s/',person);
output_path = sprintf('./monostatic_radar_%dm_%dm_%dm/%s/',Radar_pos(1),Radar_pos(2),Radar_pos(3),person);
mkdir(output_path);
files=dir(fullfile(input_paths, '/**/', '*.mat'));
for ii = 1:length(files)
    input_path = fullfile(files(ii).folder,files(ii).name);
    calSpectrogram(input_path,Radar_pos,drawScenario,output_path);
    close all;
end
%%
function [s,f,t] = calSpectrogram(input_path,Radar_pos,drawScenario,output_path)
load(input_path)
[~,name,~] = fileparts(input_path);
name = [output_path name];
% 舍弃第一帧和最后一帧.
keypoints([1,end-1],:,:)=[]; % keypoints in XYZ
timestampList([1,end-1])=[];
timestampList = timestampList - timestampList(1);
Njoints = size(keypoints,2);
joints_18 = ["NOSE","NECK","RIGHT_SHOULDER","RIGHT_ELBOW","RIGHT_WRIST","LEFT_SHOULDER","LEFT_ELBOW","LEFT_WRIST","RIGHT_HIP","RIGHT_KNEE","RIGHT_ANKLE","LEFT_HIP","LEFT_KNEE","LEFT_ANKLE","RIGHT_EYE","LEFT_EYE","RIGHT_EAR","LEFT_EAR"];
joints_34 = ["PELVIS","NAVAL_SPINE","CHEST_SPINE","NECK","LEFT_CLAVICLE","LEFT_SHOULDER","LEFT_ELBOW","LEFT_WRIST","LEFT_HAND","LEFT_HANDTIP","LEFT_THUMB","RIGHT_CLAVICLE","RIGHT_SHOULDER","RIGHT_ELBOW","RIGHT_WRIST","RIGHT_HAND","RIGHT_HANDTIP","RIGHT_THUMB","LEFT_HIP","LEFT_KNEE","LEFT_ANKLE","LEFT_FOOT","RIGHT_HIP","RIGHT_KNEE","RIGHT_ANKLE","RIGHT_FOOT","HEAD","NOSE","LEFT_EYE","LEFT_EAR","RIGHT_EYE","RIGHT_EAR","LEFT_HEEL","RIGHT_HEEL"];
if Njoints ==34
    connections = [1 2; 2 3; 3 5; 5 6; 6 7; 7 8; 8 9;9 10;8 11;3 12;12 13;13 14;14 15;15 16;16 17;15 18;1 19;19 20;20 21;21 22;1 23;23 24;24 25;25 26;3 4;4 27;27 28;28 29;29 30;28 31;31 32;21 33;25 34;33 22;34 26];
else % Njoints = 18
    connections = [1 2;2 3;3 4;4 5;2 6;6 7;7 8;3 9;9 10;10 11;6 12;12 13;13 14;3 6;9 12;1 15;15 17;1 16;16 18];
end
Nframes = length(timestampList);
frameLength = 1/30; % fps =  30;
T = frameLength*Nframes;
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
        human = plot3(z,x,y,'.','markersize', 13);
        hold on;
        axis equal;
        xlim([-8 1]); % Z
        ylim([-3 4]); % X
        zlim([-1.5 1.5]); % Y
        view(30,30)
        %         view(2)
        camera = scatter3(0,0,0,[],"red",'o','DisplayName','Camera');
        radar = scatter3(Radar_pos(3),Radar_pos(1),Radar_pos(2),[],"black",'*','DisplayName','Radar');
        for jj = 1:1:size(connections,1)
            plot3(z(connections(jj,:)),x(connections(jj,:)),y(connections(jj,:)),'Color','b','LineWidth',0.05);
        end
        for nj=1:Njoints
            line([Radar_pos(3) z(nj)],[Radar_pos(1) x(nj)],...
                [Radar_pos(2) y(nj)],'LineStyle',':',...
                'color',[0.5 0.5 0.5],'LineWidth',0.05)
        end
        xlabel('Z(m)'); ylabel('X(m)'); zlabel('Y(m)'); title(sprintf('Timestamp: %d (ms)',timestampList(ii)));
        grid on;
        legend([camera radar human] ,'camera','radar','human');
        drawnow;
        Frame=getframe(gcf);
        Image=frame2im(Frame);
        [Image,map]=rgb2ind(Image,256);
        if ii == 1
            imwrite(Image,map,fullfile(pwd,sprintf('%s.gif',name)),'gif', 'Loopcount',inf,'DelayTime',0.03);
        else
            imwrite(Image,map,fullfile(pwd,sprintf('%s.gif',name)),'gif','WriteMode','append','DelayTime',0.03);
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
        r1(1:3) = keypointsNew(nf,connections(jj,1),1:3);
        r2(1:3) = keypointsNew(nf,connections(jj,2),1:3);
        r1r2_mid=0.5*(r1+r2);
        R = norm(r1r2_mid-Radar_pos);
        R1 = Radar_pos - r1r2_mid;
        R2 = r2-r1;
        Cos_Theta  = ((R1(1)*R2(1))+(R1(2)*R2(2))+(R1(3)*R2(3)))/norm(R1)/norm(R2);
        Theta = acos(Cos_Theta);
        Phase = exp(-1i*4*pi*R/lambda)/R^2;
        height = norm(r2-r1);
        radius = height/4;
        Amp = sqrt(1/4*pi*radius^4*height^2/(radius^2*(sin(Theta))^2+1/4*height^2*(cos(Theta))^2));
        rcs_joint = Amp*Phase;
        rcs = rcs + rcs_joint;
    end
    RCS(nf) = rcs;
end

% micro-Doppler signature
F = fs;
figure;% figure('Position',[500 200 900 600])
colormap(jet)
[s,f,t] = spectrogram(RCS,kaiser(256,15),250,512,F,'centered','yaxis');
spectrogram(RCS,kaiser(256,15),250,512,F,'centered','yaxis');
clim = get(gca,'CLim');
set(gca,'CLim',clim(2) + [-60 0]);
title('Micro-Doppler Signature', 'Fontsize',12,'color','k')
saveas(gcf,sprintf('%s.jpg',name));
end