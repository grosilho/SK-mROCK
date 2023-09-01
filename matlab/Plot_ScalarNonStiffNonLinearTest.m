% Plot solutions 
clear;
clc;
% close all;

folder = '../../Results/Tests/ScalarNonStiffNonLinearTest/';
sol_name = 'RKC1';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

figure;
plot(t,y);