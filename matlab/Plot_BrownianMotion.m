% Plot brownian motions
clear;
clc;
% close all;

folder = '../dist/Release/GNU_Version_10-MacOSX/';
sol_name1 = 'brow1';
sol_name2 = 'brow2';
sol_name3 = 'brow3';
file_name1 = [folder sol_name1 '.csv'];
T1=csvread(file_name1);
file_name2 = [folder sol_name2 '.csv'];
T2=csvread(file_name2);
file_name3 = [folder sol_name3 '.csv'];
T3=csvread(file_name3);

m1 = size(T1,2)-1;
time1 = [0; T1(:,end)];
W1 = [zeros(1,m1);T1(:,1:(end-1))];
m2 = size(T2,2)-1;
time2 = [0; T2(:,end)];
W2 = [zeros(1,m2);T2(:,1:(end-1))];
m3 = size(T3,2)-1;
time3 = [0; T3(:,end)];
W3 = [zeros(1,m3);T3(:,1:(end-1))];

figure;
hold on;
for i=1:m1
   plot(time1,W1(:,i));
   hold on;
end
for i=1:m2
   plot(time2,W2(:,i));
   hold on;
end
for i=1:m3
   plot(time3,W3(:,i));
   hold on;
end
