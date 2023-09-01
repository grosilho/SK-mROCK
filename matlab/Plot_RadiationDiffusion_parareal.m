
clear;
%clc;
close all;

n_threads = 12;
iter = 3;

% folder = '../dist/Release/GNU_Version_10-MacOSX/Tests/RadiationDiffusion/';
folder = '../../Results/Tests/RadiationDiffusion/';
refsol_name = 'mRKC';
reffile_name = [folder refsol_name '_evolution.m'];
sol_name = 'mRKC_mRKC';
%sol_name = 'RKC1';
if(n_threads>0)
    file_name = [folder sol_name '_in_sol_iter_' num2str(iter) ...
             '_thread_' num2str(n_threads-1) '_evolution.m'];
else
    file_name = [folder sol_name '_evolution.m'];
end

run(reffile_name);
neqn = numel(y(1,:));
Nsq = neqn/2;
N = sqrt(Nsq);
tend = t(end);

U_ref = zeros(N,N,2);
U_ref(:,:,1) = reshape(y(end,1:Nsq)',N,N);
U_ref(:,:,2) = reshape(y(end,(Nsq+1):(2*Nsq))',N,N);
clear y t;

x = linspace(0,1,N);
[Y,X]=meshgrid(x,x);


fsa = 30;
fs = [800 900];
scrsz = get(0,'ScreenSize');
az = -100;
el = 15;
color_ref = [0.15 0.45 0.09];
color_sol = [0.15, 0.38, 0.61]; %RKC
%sol_color = [0.87, 0.19, 0.39]; %RKL
%sol_color = [1.00, 0.40, 0.0]; %RKU


run(file_name);
U = zeros(N,N,2);
U(:,:,1) = reshape(y(end,1:Nsq)',N,N);
U(:,:,2) = reshape(y(end,(Nsq+1):(2*Nsq))',N,N);
clear y t;

figure('Position',[scrsz(3)/2 scrsz(4)/2 fs(1) fs(2)]);
for k=1:2
    subplot(2,1,k);
    hold on;
    surf(X,Y,U_ref(:,:,k),'SpecularExponent',1,'SpecularStrength',1,...
        'DiffuseStrength',1,'AmbientStrength',0.4,'FaceAlpha',0.5,...
        'FaceColor',color_ref,'AlignVertexCenters','on',...
        'LineWidth',0.2,'EdgeAlpha',0);
    surf(X,Y,U(:,:,k),'FaceLighting','gouraud','MeshStyle','column',...
            'SpecularColorReflectance',0,'SpecularExponent',5,'SpecularStrength',0.2,...
            'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on',...
            'LineWidth',0.2,'FaceAlpha',0.5,'FaceColor',color_sol,...
            'EdgeAlpha',0.2,'EdgeColor','none');
    view(az,el);
    ax = gca;
    %ax.YAxisLocation = 'right';
    set(gca,'fontsize',fsa);
    set(gca,'TickLabelInterpreter','latex')
    xl=xlabel('$x$','fontsize',fsa,'interpreter','LaTeX');
    yl=ylabel('$y$','fontsize',fsa,'interpreter','LaTeX');
    if(k==1)
        lab_var = 'E';
    else
        lab_var = 'T';
    end
    zl=zlabel(['$' lab_var '$'],'fontsize',fsa,'interpreter','LaTeX');
    err = norm(U_ref(:,:,k)-U(:,:,k),"fro")/N;
    fprintf(['Error ' lab_var ' = %f\n'], err );
    axis([min(x) max(x) min(x) max(x)]);      
end






