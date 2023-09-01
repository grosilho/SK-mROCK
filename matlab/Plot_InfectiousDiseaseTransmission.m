% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../dist/Release/GNU-Linux/Tests/InfectiousDiseaseTransmission/';
sol_name = 'rkc';
file_name = [folder sol_name '.m'];

run(file_name);


S = y(:,1);
I = y(:,2);
V = y(:,3);

% figure;
% subplot(3,1,1);
% plot(t,S);
% subplot(3,1,2);
% plot(t,I);
% subplot(3,1,3);
% plot(t,V);

ind = find(t>10,1);

figure;
subplot(3,1,1);
plot(t(ind:end),S(ind:end));
subplot(3,1,2);
plot(t(ind:end),I(ind:end));
subplot(3,1,3);
plot(t(ind:end),V(ind:end));
% 
% figure;
% subplot(3,1,1);
% loglog(t(ind:end),S(ind:end));
% subplot(3,1,2);
% loglog(t(ind:end),I(ind:end));
% subplot(3,1,3);
% loglog(t(ind:end),V(ind:end));

