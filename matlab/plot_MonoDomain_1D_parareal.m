% Plot solutions MonoDomain in 1D
clear;
clc;
%close all;

folder = '../../Results/Tests/MonoDomain_in_1D/';
refsol_name = 'RKU1';
sol_name = 'RKU1_RKU1';
ref_file_name = [folder refsol_name '_evolution.bin'];
N_gating_vars = 1;

plot_sol = 1;
save_png=0;
color_ref = [0.15 0.45 0.09];
sol_color = [0.15, 0.38, 0.61]; %RKC
%sol_color = [0.87, 0.19, 0.39]; %RKL
%sol_color = [1.00, 0.40, 0.0]; %RKU

iters = 1;
threads = 1:18;


fileID = fopen(ref_file_name);
A = fread(fileID,'double');
n_el_A = numel(A);
neqn = 501;
n_y_var = (1+N_gating_vars)*neqn;
n_time_steps = round(n_el_A/(n_y_var+1));
A = reshape(A,[n_y_var+1,n_time_steps]);
t_ref = A(1,:);
V_ref = A(2:(neqn+1),:);
g_ref = A((neqn+2):end,:);
clear A;

x_ref = linspace(0,5,neqn);
[T_ref,X_ref]=meshgrid(t_ref,x_ref);

%file_name = [folder sol_name '_out_sol_iter_' num2str(0) '_evolution.bin'];
%run(file_name);
%g_outer = zeros(N_gating_vars,numel(t),neqn);
%V_outer = y(:,space_indeces);
%for k=1:N_gating_vars
%    g_outer(k,:,:) = y(:,(neqn+k):N_gating_vars:end);
%end
%t_outer = t;
%clear y t;
%[X_outer,T_outer]=meshgrid(x_ref,t_outer);

fsa = 30;
fs = [800 450];
scrsz = get(0,'ScreenSize');
az = -170;
el = 15;
%figure;
%surf(X_ref,T_ref,V_ref,'FaceColor','r', 'FaceAlpha',0.8, 'EdgeColor','r');
%hold on;
%surf(X_outer,T_outer,V_outer,'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','b','LineStyle','-');

err = [];
for k=iters
    %file_name = [folder sol_name '_out_sol_iter_' num2str(k) '_evolution.m'];
    %run(file_name);
    %V_outer = y(:,space_indeces);
    %g_outer = zeros(N_gating_vars,numel(t),neqn);
    %for j=1:N_gating_vars
    %    g_outer(j,:,:) = y(:,(neqn+j):N_gating_vars:end);
    %end
    %t_outer = t;
    %clear y t;
    %[X_outer,T_outer]=meshgrid(x_ref,t_outer);

    if(plot_sol)
        figure('Position',[scrsz(3)/2 scrsz(4)/2 fs(1) fs(2)]);
        subplot(2,1,1);
        %surf(X_ref,T_ref,V_ref,'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','r');
        %hold on;
        %surf(X_outer,T_outer,V_outer,'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','b','LineStyle','-');
    
        surf(X_ref,T_ref,V_ref,'SpecularExponent',1,'SpecularStrength',1,...
            'DiffuseStrength',1,'AmbientStrength',0.4,'FaceAlpha',0.5,...
            'FaceColor',color_ref,'AlignVertexCenters','on',...
            'LineWidth',0.2,'EdgeAlpha',0);
        hold on;
        view(az,el);
        ax = gca;
        %ax.YAxisLocation = 'right';
        set(gca,'fontsize',fsa);
        set(gca,'TickLabelInterpreter','latex')
        xl=xlabel('$x$','fontsize',fsa,'interpreter','LaTeX');
        yl=ylabel('$t$','fontsize',fsa,'interpreter','LaTeX');
        zl=zlabel('$y$','fontsize',fsa,'interpreter','LaTeX');
        axis([min(x_ref) max(x_ref) min(t_ref) max(t_ref)]);

        subplot(2,1,2);
        surf(X_ref,T_ref,g_ref,'SpecularExponent',1,'SpecularStrength',1,...
            'DiffuseStrength',1,'AmbientStrength',0.4,'FaceAlpha',0.5,...
            'FaceColor',color_ref,'AlignVertexCenters','on',...
            'LineWidth',0.2,'EdgeAlpha',0);
        hold on;
        view(az,el);
        ax = gca;
        %ax.YAxisLocation = 'right';
        set(gca,'fontsize',fsa);
        set(gca,'TickLabelInterpreter','latex')
        xl=xlabel('$x$','fontsize',fsa,'interpreter','LaTeX');
        yl=ylabel('$t$','fontsize',fsa,'interpreter','LaTeX');
        zl=zlabel('$g$','fontsize',fsa,'interpreter','LaTeX');
        axis([min(x_ref) max(x_ref) min(t_ref) max(t_ref)]);
    end
    for n=threads
        file_name = [folder sol_name '_in_sol_iter_' num2str(k) '_thread_' num2str(n-1) '_evolution.bin'];
        fileID = fopen(file_name);
        A = fread(fileID,'double');
        n_el_A = numel(A);
        n_time_steps = round(n_el_A/(n_y_var+1));
        A = reshape(A,[n_y_var+1,n_time_steps]);
        t_inner = A(1,:);
        V_inner = A(2:(neqn+1),:);
        g_inner = A((neqn+2):end,:);
        clear A;

        [T_inner,X_inner]=meshgrid(t_inner,x_ref);

        %surf(X_inner,T_inner,V_inner,'FaceColor','b', 'FaceAlpha',0., 'EdgeColor','b','LineStyle','-');
        
        if(plot_sol)
            subplot(2,1,1);
            surf(X_inner,T_inner,V_inner,'FaceLighting','gouraud','MeshStyle','column',...
                'SpecularColorReflectance',0,'SpecularExponent',5,'SpecularStrength',0.2,...
                'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on',...
                'LineWidth',0.2,'FaceAlpha',0.5,'FaceColor',sol_color,...
                'EdgeAlpha',0.2,'EdgeColor','none');

            subplot(2,1,2);
            surf(X_inner,T_inner,g_inner,'FaceLighting','gouraud','MeshStyle','column',...
                'SpecularColorReflectance',0,'SpecularExponent',5,'SpecularStrength',0.2,...
                'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on',...
                'LineWidth',0.2,'FaceAlpha',0.5,'FaceColor',sol_color,...
                'EdgeAlpha',0.2,'EdgeColor','none');
        end
    end
    if(plot_sol)
        subplot(2,1,1);
        legend('$y$',['$y^{' num2str(k) '}$'],'fontsize',fsa,'interpreter','LaTeX');
        subplot(2,1,2);
        legend('$g$',['$g^{' num2str(k) '}$'],'fontsize',fsa,'interpreter','LaTeX');
        if(save_png==1)
            screen2eps(['~/' sol_name '_k_' num2str(k) '.png'],'png'); 
        end
    end

    err = [err, max(abs(V_inner(:,end)-V_ref(:,end))./(1+abs(V_ref(:,end))))];
    %fprintf('Error iter %i = %d\n',k,err);
end

figure;
fsa=18;
method = refsol_name(1:4);
if(method=='RKC1')
    method = [method ' $\varepsilon = ' refsol_name(end) '$'];
end
semilogy(iters,err,'DisplayName',method,'LineWidth',2);
hold on;
set(gca,'fontsize',fsa);
set(gca,'TickLabelInterpreter','latex')
xl=xlabel('iteration Parareal','fontsize',fsa,'interpreter','LaTeX');
yl=ylabel('$\ell_\infty$ error','fontsize',fsa,'interpreter','LaTeX');
%legend(pl,method,'fontsize',fsa,'interpreter','LaTeX');
legend('fontsize',fsa,'interpreter','LaTeX');
legend show;
