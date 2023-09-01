% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../../Results/Tests/MonoDomain_in_1D/';
sol_name = 'RKC1';
file_name = [folder sol_name '_eigs.m'];
run(file_name);

[X,T]=meshgrid(1:size(eigvals,2),t);

fsa = 20;
figure;
subplot(2,1,1);
surf(T,X,real(eigvals));
shading interp;
ax = gca;
set(gca,'fontsize',fsa);
set(gca,'TickLabelInterpreter','latex')
xl=xlabel('$t$','fontsize',fsa,'interpreter','LaTeX');
yl=ylabel('$i$','fontsize',fsa,'interpreter','LaTeX');
zl=zlabel('$\mathcal{R}e(\lambda_i$)','fontsize',fsa,'interpreter','LaTeX');
subplot(2,1,2);
surf(T,X,imag(eigvals));
shading interp;
ax = gca;
set(gca,'fontsize',fsa);
set(gca,'TickLabelInterpreter','latex')
xl=xlabel('$t$','fontsize',fsa,'interpreter','LaTeX');
yl=ylabel('$i$','fontsize',fsa,'interpreter','LaTeX');
zl=zlabel('$\mathcal{I}m(\lambda_i$)','fontsize',fsa,'interpreter','LaTeX');

