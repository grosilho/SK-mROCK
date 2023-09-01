% Plot solutions Integro differential equaiton
clear;
clc;
% close all;

folder = '../../Results/Tests/IntegroDifferentialEquation/';
sol_name = 'RKC1';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

neqn = numel(y(1,:));

x = linspace(0,1,neqn+1);
x = x(2:end);
[X,T]=meshgrid(x,t);

figure;
surf(X,T,y);
shading interp;