clear;
clc;
close all;

folder = '../dist/Release/GNU_Version_10-MacOSX/Tests/CUSP/';
sol_name = 'spIE';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

npoints = numel(y(1,:))/3;

Y = y(:,1:npoints);
A = y(:,(npoints+1):(2*npoints));
B = y(:,(2*npoints+1):end);

x = linspace(0,1,npoints+1);
x = x(1:npoints);
[X,T]=meshgrid(x,t);

figure;
subplot(1,3,1);
surf(X,T,Y);
shading interp;
subplot(1,3,2);
surf(X,T,A);
shading interp;
subplot(1,3,3);
surf(X,T,B);
shading interp;