clc;
clear;
close all;
model =CreateModel();
tic;
plotmap(model);
global_chromosome =Muti_Uav_Ga(model);
toc;


