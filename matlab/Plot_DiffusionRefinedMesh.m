% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../dist/Release/GNU_Version_10-MacOSX/Tests/DiffusionRefinedMesh/';
sol_name = 'sol_IE';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

N1 = 100;
N2 = 200;

x1 = linspace(0,1,N1+1);
x2 = linspace(1,2,N2+1);
x = [x1 x2(2:end)];

nt = numel(y(:,1));
I = round(linspace(0,1,10)*nt);
I(1)=1;

figure;
hold on;
for i=I
    plot(x,y(i,:));
end

