
clear;
clc;
close all;

% folder = '../dist/Release/GNU_Version_10-MacOSX/Tests/RadiationDiffusion/';
folder = '../../Results/Tests/RadiationDiffusion/';
sol_name = 'RKC1';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

neqn = numel(y(1,:));
Nsq = neqn/2;
N = sqrt(Nsq);
nt = numel(t);

ntout = 5;
tout_ind = round(linspace(1,nt,ntout));

% tout = [1,2,3];
% tout_ind = find(t>tout,1);
% ntout = numel(tout);

U = zeros(N,N,ntout,2);
for k=1:2
    for n=1:ntout
        U(:,:,n,k) = reshape(y(tout_ind(n),((k-1)*Nsq+1):(k*Nsq))',N,N);
    end
end

x = linspace(0,1,N);
[Y,X]=meshgrid(x,x);

for n=1:ntout
    figure;
    for k=1:2
        subplot(2,1,k);
        surf(X,Y,U(:,:,n,k));
        xlabel('x');
        ylabel('y');
        shading interp;        
    end
    sgtitle(['t = ' string t(tout_ind(n))]);
end





