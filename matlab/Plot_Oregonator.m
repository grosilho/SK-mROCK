% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../dist/Release/GNU-Linux/Tests/Oregonator/';
sol_name = 'sol_rkc1';
file_name = [folder sol_name '.m'];

run(file_name);


y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);
y4 = y(:,4);
y5 = y(:,5);

figure;
subplot(5,1,1);
semilogy(t,y1);
subplot(5,1,2);
semilogy(t,y2);
subplot(5,1,3);
semilogy(t,y3);
subplot(5,1,4);
semilogy(t,y4);
subplot(5,1,5);
semilogy(t,y5);

figure;
loglog(y2,y3);
figure;
loglog(y2,y5);

