% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../dist/Release/GNU-Linux/Tests/VanDerPol/';
sol_name = 'sol_rkc_PPI_2';
file_name = [folder sol_name '.m'];

run(file_name);

figure;
plot(t,y(:,1))
hold on;
plot(t,y(:,2));

figure;
plot(y(:,1),y(:,2));