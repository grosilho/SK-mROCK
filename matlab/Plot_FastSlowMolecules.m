% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../results/FastSlowMolecules/';
solname = 'sol';

n_states = 2;
fileID = fopen([folder solname '_evolution.bin']);
y = fread(fileID,'double');
y = reshape(y,[1+n_states,numel(y)/(1+n_states)]);
t = y(1,:);
y = y(2:end,:);

t_end = 100.;
t_end = min(t_end,t(end));
tol = 1e-6;
last = find(t>=t_end-tol,1);

t = t(1:last);
y = y(:,1:last);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure;
plot(t,y(1,:),'LineWidth',1.5);
hold on;
plot(t,y(2,:),'LineWidth',1.5);
legend('Fast','Slow','fontsize',14);
set(gca,'FontSize',15);
title('Population Dynamics Model')


