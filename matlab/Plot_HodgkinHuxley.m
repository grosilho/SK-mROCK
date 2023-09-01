% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../../Results/Tests/HodgkinHuxley/';
sol_name = 'sol';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

npoints = numel(y(1,:))/4;

V = y(:,1:npoints);
n = y(:,(npoints+1):3:end);
m = y(:,(npoints+2):3:end);
h = y(:,(npoints+3):3:end);

x = linspace(0,1,npoints);
[X,T]=meshgrid(x,t);

figure;
subplot(2,2,1);
surf(X,T,V);
shading interp;
subplot(2,2,2);
surf(X,T,n);
shading interp;
subplot(2,2,3);
surf(X,T,m);
shading interp;
subplot(2,2,4);
surf(X,T,h);
shading interp;



