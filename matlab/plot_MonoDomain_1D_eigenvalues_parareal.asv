% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../../Results/Tests/MonoDomain_in_1D/';
sol_name = 'RKC1_RKC1_damp_1_no_omp';

iters = 1:4;
threads = 1:8;

fsa = 20;

dt_c=100;
dt_f=0.1;

for k=iters
    figure(k+1-iters(1));
    hold on;
    for n=threads
        file_name = [folder sol_name '_in_sol_iter_' num2str(k) '_thread_' num2str(n-1) '_eigs.m'];
        run(file_name);
        [X,T]=meshgrid(1:size(eigvals,2),t);

        subplot(2,1,1);
        surf(T,X,real(eigvals));
        hold on;
        shading interp;
        ax = gca;
        set(gca,'fontsize',fsa);
        set(gca,'TickLabelInterpreter','latex')
        xl=xlabel('$t$','fontsize',fsa,'interpreter','LaTeX');
        yl=ylabel('$i$','fontsize',fsa,'interpreter','LaTeX');
        zl=zlabel('$\mathcal{R}e(\lambda_i$)','fontsize',fsa,'interpreter','LaTeX');
        subplot(2,1,2);
        surf(T,X,imag(eigvals));
        hold on;
        shading interp;
        ax = gca;
        set(gca,'fontsize',fsa);
        set(gca,'TickLabelInterpreter','latex')
        xl=xlabel('$t$','fontsize',fsa,'interpreter','LaTeX');
        yl=ylabel('$i$','fontsize',fsa,'interpreter','LaTeX');
        zl=zlabel('$\mathcal{I}m(\lambda_i$)','fontsize',fsa,'interpreter','LaTeX');
    end
end