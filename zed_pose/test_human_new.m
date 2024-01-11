clc; clear; close all;
addpath("../mediapipe_spectrogram/utils/");
if exist('./data0103_18/', 'dir') == 7
    rmdir('./data0103_18/', 's');
end
mkdir('./data0103_18/');
Tx_pos = [0 -0.1 0]; % XYZ
Rx_pos = [0 -0.1 1.5]; % XYZ
fc = 60.48e9;
fs = 2000;
AWGN_mean = 0;
AWGN_var = 0;
thres_A_TRD = -30;
drawScenario = false;
rcsRendering = false;
using_camera_coordinate = true;
connections18 = [1 2;2 3;3 4;4 5;2 6;6 7;7 8;3 9;9 10;10 11;6 12;12 13;13 14;3 6;9 12;1 15;15 17;1 16;16 18];
connections34 = [1 2; 2 3; 3 5; 5 6; 6 7; 7 8; 8 9;9 10;8 11;3 12;12 13;13 14;14 15;15 16;16 17;15 18;1 19;19 20;20 21;21 22;1 23;23 24;24 25;25 26;3 4;4 27;27 28;28 29;29 30;28 31;31 32;21 33;25 34;33 22;34 26];
connections = connections18;
folder = '../testZED/data_new_18';
subfolders = ["back" "boxing" "clapping" "hand_cycling" "jump1" "jump2" "leg" "run" "sorry" "swim" "walk"];
subfolders = ["clapping"];
for ii = 1:length(subfolders)
    subfolder = subfolders(ii)
    index = 2;
    filename_mat = sprintf("%s_%d.mat",subfolder,index)
    filename_jpg = sprintf("%s_%d.jpg",subfolder,index)
    filename_gif = sprintf("%s_%d.gif",subfolder,index)
    input_mat_path = absPath(fullfile(folder,subfolder,filename_mat))
    output_jpg_path = absPath(sprintf("./data0103_18/%s",filename_jpg))
    output_gif_path = absPath(sprintf("./data0103_18/%s",filename_gif))
    pic_save = true;
    tic
    simuSpectrogram(Tx_pos,Rx_pos,fc,fs,AWGN_mean,AWGN_var,thres_A_TRD, ...
        drawScenario,rcsRendering,input_mat_path,using_camera_coordinate, ...
        connections,output_jpg_path,output_gif_path,pic_save);
    toc
end
function absolutePath = absPath(relativePath)
currentPath = pwd;
absolutePath = fullfile(currentPath, relativePath);
end
