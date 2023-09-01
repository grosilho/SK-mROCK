% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../../Results/Tests/Brusselator/';
sol_name = 'sol';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

y1 = y(:,1);
y2 = y(:,2);

figure;
subplot(2,1,1);
plot(t,y1);
subplot(2,1,2);
plot(t,y2);

