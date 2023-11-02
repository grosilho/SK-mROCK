% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../results/StochasticReactionDiffusion/';
sol_name = 'sol';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

neqn = numel(y(1,:));

x = linspace(0,1,neqn);
[X,T]=meshgrid(x,t);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure;
surf(X,T,y);
shading interp;
xlabel('x')
ylabel('t')
set(gca,'FontSize',15);
title('Reaction-Diffusion Equation')

