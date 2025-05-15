clear all;
close all;
step=3000; % number of timesteps

color4=[159 159 159]/255;
color8=[176,  92,  176]/255;
blue3=[86,160,148]/255;
color5=[141, 159, 191]/255;

nematic_order=importdata("nematic_order.txt");
disp('global:');
mean(mean(nematic_order,"omitnan"))

local=importdata("local_order.txt");
disp('local:');
mean(mean(local,"omitnan"))