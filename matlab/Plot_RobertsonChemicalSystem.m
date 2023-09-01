% Plot solutions Van Der Pol
clear;
clc;
close all;

folder = '../dist/Release/GNU-Linux/Tests/RobertsonChemicalSystem/';
sol_name = 'rkc';
file_name = [folder sol_name '.m'];

run(file_name);


y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);

figure;
subplot(3,1,1);
loglog(t,y1);
subplot(3,1,2);
loglog(t,y2);
subplot(3,1,3);
loglog(t,y3);


figure;
subplot(3,1,1);
plot(t,y1);
subplot(3,1,2);
plot(t,y2);
subplot(3,1,3);
plot(t,y3);

