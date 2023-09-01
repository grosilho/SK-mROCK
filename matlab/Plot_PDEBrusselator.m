
clear;
clc;
close all;

folder = '../dist/Release/GNU_Version_10-MacOSX/Tests/PDEBrusselator/';
folder = '../dist/Release_Linux/GNU-11-Linux/Tests/PDEBrusselator/';
sol_name = 'IE_AD_DF';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

neqn = numel(y(1,:));
Nu = (2*neqn+1)/3;
Nv = neqn-Nu;

U = y(:,1:Nu);
V = y(:,(Nu+1):end);


xu = linspace(0,1,Nu+2);
xu = xu(2:(end-1));
xv = linspace(0,1,Nv+2);
xv = xv(2:(end-1));
[Xu,Tu]=meshgrid(xu,t);
[Xv,Tv]=meshgrid(xv,t);

figure;
subplot(1,2,1);
surf(Xu,Tu,U);
shading interp;
subplot(1,2,2);
surf(Xv,Tv,V);
shading interp;