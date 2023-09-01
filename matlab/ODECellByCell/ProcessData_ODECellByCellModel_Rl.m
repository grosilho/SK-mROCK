clear;
%clc;
%close all;

folder = '../../results/ODECellByCellModel/';
folder = '../../results_cluster/ODECellByCellModel/nx=30/';
subfolder = 'CVl_C_Rl_cl_100/';

Rl_list = readtable([folder subfolder 'Rl_list.csv']);

first_Rl = 1;
last_Rl= Rl_list.n(end);
Rl_list = Rl_list(first_Rl:last_Rl,:);

n_states = 21; % 8 for LuoRudy, 21 for Courtemanche, 41 for Ohara Rudy
cl = 100e-4;
cw = 10e-4;
nc = 30;
ps = zeros(2,2);
% longitudinal 
ps(1,:) = [9.5*cl,0];
ps(2,:) = [19.5*cl,0];
V_th = -20;

CV = zeros(size(Rl_list,1),1);
for Rli=1:size(Rl_list,1)
    solname = ['Rl_' num2str(Rl_list.n(Rli))];
    run([folder subfolder solname '_geo.m']);
    transmembrane_G = [];
    dom=mdom{1};
    for i=2:size(dom,1) 
        transmembrane_G = [transmembrane_G;dom{i}.G];
    end
    n_pts = size(transmembrane_G,1);
    %fileID = fopen([folder subfolder solname '_evolution.bin']);
    %precision = [num2str(1+n_pts) '*double'];
    %skip = n_pts*(n_states-1)*8;
    %y = fread(fileID,precision,skip);
    %y = reshape(y,[1+n_pts,numel(y)/(1+n_pts)]);
    %t = y(1,:);
    %y = y(2:end,:);
    fileID = fopen([folder subfolder solname '_V_evolution.bin']);
    y = fread(fileID,'double');
    y = reshape(y,[1+n_pts,numel(y)/(1+n_pts)]);
    t = y(1,:);
    y = y(2:end,:);
    hist_j=zeros(2,1);
    p = zeros(2,2);
    ps0=ps;
    
    t2_j=double.empty(1,0);
    while (isempty(t2_j) && ps0(2,1)>0)
        d = transmembrane_G-ps0(2,:);
        d = sqrt(d(:,1).^2+d(:,2).^2);
        [~,j] = min(d);
        hist_j(2) = j;
        p(2,:) = transmembrane_G(j,:);
        t2_j = find(y(hist_j(2),:)>=V_th,1);
        ps0(2,:)=ps0(2,:)-[cl,0];
    end

    t1_j=double.empty(1,0);
    ps0(1,1) = min(ps0(1,1),ps0(2,1)-cl);
    while (isempty(t1_j) && ps0(1,1)>0)
        d = transmembrane_G-ps0(1,:);
        d = sqrt(d(:,1).^2+d(:,2).^2);
        [~,j] = min(d);
        hist_j(1) = j;
        p(1,:) = transmembrane_G(j,:);
        t1_j = find(y(hist_j(1),:)>=V_th,1);
        ps0(1,:)=ps0(1,:)-[cl,0];
    end
    if(ps0(1,1)<0|| ps0(2,1)<0)
        fprintf('Break at %i because no propagation anymore\n',Rli);
        break;
    end
    dt2 = t(t2_j)-t(t1_j); % in ms
    d = norm(p(1,:)-p(2,:)); % in cm
    d = 10*d; % in mm
    CV(Rli) = d/dt2;
    fprintf(solname);
    fprintf(' CV  = %f mm/ms (m/s)\n',CV(Rli));
end


figure;
subplot(2,1,1);
semilogx(Rl_list.Rl,CV);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
xlabel('$R_g$','interpreter','latex')
ylabel('$CV$','interpreter','latex')
subplot(2,1,2);
semilogx(1./Rl_list.Rl,CV);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
xlabel('$\kappa$','interpreter','latex')
ylabel('$CV$','interpreter','latex')

T=table(Rl_list.Rl,1./Rl_list.Rl,CV,'VariableNames',{'Rl','k','CV'});
writetable(T,[subfolder(1:end-1) '.csv']);


