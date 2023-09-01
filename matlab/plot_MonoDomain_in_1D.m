% Plot solutions Van Der Pol
clear;
clc;
% close all;

folder = '../../Results/Tests/MonoDomain_in_1D/';
sol_name = 'RKC1_damp_1_dt_100';
file_name = [folder sol_name '_evolution'];
N_gating_vars = 1;

tic;
% Reading a .m file (slow)
%file_name = [file_name '.m'];
%run(file_name);
%neqn = numel(y(1,:))/(1+N_gating_vars);
%V = y(:,1:neqn)';
%n = y(:,(neqn+1):end)';

% Reading a .bin file (faster)
file_name = [file_name '.bin'];
fileID = fopen(file_name);
A = fread(fileID,'double');
n_el_A = numel(A);
neqn = 501;
n_y_var = (1+N_gating_vars)*neqn;
n_time_steps = round(n_el_A/(n_y_var+1));
A = reshape(A,[n_y_var+1,n_time_steps]);
t = A(1,:);
V = A(2:(neqn+1),:);
n = A((neqn+2):end,:);
% End of reading file
toc;

x = linspace(0,5,neqn);
[T,X]=meshgrid(t,x);

fsa = 30;
fs = [800 450];
scrsz = get(0,'ScreenSize');
az = -128;
el = 40;
figure('Position',[scrsz(3)/2 scrsz(4)/2 fs(1) fs(2)]);
hold on;
subplot(1+N_gating_vars,1,1);
surf(X,T,V,'SpecularExponent',1,...
        'SpecularStrength',1,...
         'DiffuseStrength',1,...
    'AmbientStrength',0.4,...
    'FaceAlpha',0.7,...
    'FaceColor',[0.15 0.45 0.09],...
    'AlignVertexCenters','on',...
    'LineWidth',0.2,...
    'EdgeAlpha',0);
view(az,el);
shading interp;
ax = gca;
%ax.YAxisLocation = 'right';
set(gca,'fontsize',fsa);
set(gca,'TickLabelInterpreter','latex')
xl=xlabel('$x$','fontsize',fsa,'interpreter','LaTeX');
yl=ylabel('$t$','fontsize',fsa,'interpreter','LaTeX');
zl=zlabel('$y$','fontsize',fsa,'interpreter','LaTeX');
axis([min(x) max(x) min(t) max(t)]);

for k=1:N_gating_vars
    subplot(1+N_gating_vars,1,1+k);
    surf(X,T,n);
    shading interp;
end

