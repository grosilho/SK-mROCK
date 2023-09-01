% Plot solutions uBidomain
clear;
clc;
close all;

folder = '../dist/Release/GNU-11-MacOSX/Tests/uBidomain/';
sol_name = 'RKC1_BCF';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

neqn = numel(y(1,:));
Bnd_pts = neqn/4;

V = y(:,1:Bnd_pts);
v = y(:,(Bnd_pts+1):(2*Bnd_pts));
w = y(:,(2*Bnd_pts+1):(3*Bnd_pts));
s = y(:,(3*Bnd_pts+1):(4*Bnd_pts));

x = linspace(0,1,Bnd_pts);
[X,T]=meshgrid(x,t);

figure;
subplot(2,2,1);
surf(X,T,V);
shading interp;

subplot(2,2,2);
surf(X,T,v);
shading interp;

subplot(2,2,3);
surf(X,T,w);
shading interp;

subplot(2,2,4);
surf(X,T,s);
shading interp;
