% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     cr_test.m
%    Authors:       A. Lomayev, R. Maslennikov, Y. Gagiev
%    Version:       1.0
%    History:       May 2010 created
%
%  *************************************************************************************
%    Description:
%
%    test
%
%  *************************************************************************************/
clear all
clc

seed = 1;
randn('state',seed);
rand('state',seed);
for iter = 1 : 10
[imp_res] = cr_ch_model;

figure;
stem(abs(imp_res.h11),'b')
grid on
title('example of channel impulse response h11 realization','FontSize',8,'FontWeight','bold')
xlabel('samples','FontSize',8,'FontWeight','bold');
set(gca,'FontSize',8,'FontWeight','bold');
end

