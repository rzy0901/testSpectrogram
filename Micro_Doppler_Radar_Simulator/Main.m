close all;
clc;
clear;

rho = 0.99; % time correlation of the static channel, 1 -> unchanged; 0 -> uncorrelated
Radar_x = 0.1;
Radar_y = 4.4;
Radar_z = 1;
Radar_pos = [Radar_x;Radar_y;Radar_z]; % the position of radar

ped_height = 1.75; % the height of pedestrian
v = 1; % moving speed
Pedestrian_x = 1.5;
Pedestrian_y = 3.5;
Pedestrian_z = 0;
Pedestrian_pos = [Pedestrian_x;Pedestrian_y;Pedestrian_z]; % the position of pedestrian
Pedestrian_heading = 340; % the moving direction of pedestrian

Simulator(v,ped_height,rho,Pedestrian_heading,Pedestrian_pos,Radar_x,Radar_y,Radar_z,Radar_pos);