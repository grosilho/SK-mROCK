% Read multirate parameters

clear;
clc;
% close all;

folder = '../../Results/Tests/MonoDomain/';
sol_name = 'mRKC_dt_0p5';
file_name = [folder sol_name '_params.csv'];

T = readtable(file_name);
T.t = cumsum(T.dt);

subplot(1,2,1);
plot(T.t,T.rhoF,'r');
hold on;
plot(T.t,T.rhoS,'b');

subplot(1,2,2);
plot(T.t, T.s, 'b');
hold on;
plot(T.t,T.m,'r');
