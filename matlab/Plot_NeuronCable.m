% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../../Results/Tests/NeuronCable/';
sol_name = 'RKU2';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

neqn = numel(y(1,:));

x = linspace(0,1,neqn);
[X,T]=meshgrid(x,t);

figure;
surf(X,T,y);
shading interp;

