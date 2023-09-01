clear;
clc;
close all;

folder = '../dist/Release/GNU_Version_10-MacOSX/Tests/mReactionDiffusion2DEquations/';
sol_name = 'sol_rkc';
% sol_name = 'sol_rkcsq';
file_name = [folder sol_name '_evolution.m'];

run(file_name);

m = 4;
neqn = numel(y(1,:));
Nsq = neqn/m;
N = sqrt(Nsq);
nt = numel(t);

ntout = 5;
tout_ind = round(linspace(1,nt,ntout));

% tout = 0.5;
% tout_ind = find(t>tout,1);
% ntout = numel(tout);

U = zeros(N,N,ntout,m);
for k=1:m
    for n=1:ntout
        U(:,:,n,k) = reshape(y(tout_ind(n),((k-1)*Nsq+1):(k*Nsq))',N,N);
    end
end



x = linspace(0,1,N+1);
x = x(1:(end-1));
[Y,X]=meshgrid(x,x);

for n=1:ntout
    figure;
    for k=1:m
        subplot(m,1,k);
        surf(X,Y,U(:,:,n,k));
        xlabel('x');
        ylabel('y');
        shading interp;        
    end
    sgtitle(['t = ' string t(tout_ind(n))]);
end


