% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../dist/Release/GNU-Linux/Tests/PopulationDynamics/';
sol_name = 'sol_rkc';
file_name = [folder sol_name '.m'];

run(file_name);

figure;
plot(t,y(:,1))
hold on;
plot(t,y(:,2));