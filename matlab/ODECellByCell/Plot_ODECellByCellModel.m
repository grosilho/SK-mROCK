clear;
%clc;
%close all;

folder = '../../results/ODECellByCellModel/';
% folder = '../../results_cluster/ODECellByCellModel/';
% folder = '../../results_cluster/ODECellByCellModel/CVltd_C_flat_randomly_deactivated_gj_nx_20_ny_20/p_100/';
% folder = '../../results/ODECellByCellModel/CVtld_C_nx_3_ny_3/';
% folder = '../../results_shared/ODECellByCellModel/NUOVI/CVl_C_amp_smooth/';   
% folder = '../../results_shared/ODECellByCellModel/CVl_C_amp/';   
% solname = 'freq_2_amp_10';
% solname = 'squaredwave';
%solname = 'CVt_C_squaredwave_op/center/both/length_4_amplitude_4_p_1';
solname = 'sol';

run([folder solname '_geo.m']);

read_gating_vars = false;
t_max = 4200; 
V_th = -20;

cl = 100e-4;
cw = 10e-4;
% nc = 30;
ncx = 20;
ncy = 20;
direction = 'l';
old_file = 0;

transmembrane_G = [];
dom=mdom{1};
for i=2:size(dom,1) 
    transmembrane_G = [transmembrane_G;dom{i}.G];
end
n_pts = size(transmembrane_G,1);

if old_file
    fileID = fopen([folder solname '_evolution.bin']);
    fprintf([folder solname '_evolution.bin\n'])
    n_states=21;
    precision = [num2str(1+n_pts) '*double'];
    skip = n_pts*(n_states-1)*8;
    y = fread(fileID,precision,skip);
    y = reshape(y,[1+n_pts,numel(y)/(1+n_pts)]);
    t = y(1,:);
    y = y(2:end,:);
else
    fileID = fopen([folder solname '_V_evolution.bin']);
    y = fread(fileID,'double');
    y = reshape(y,[1+n_pts,numel(y)/(1+n_pts)]);
    t = y(1,:);
    y = y(2:end,:);
end

tol = 1e-6;
t_max = min(t_max,t(end));
t_max_i = find(t>=t_max-tol,1);
t = t(1:t_max_i);   
y = y(:,1:t_max_i);

if 1
    t_plot = 2.3;    
    t_i = find(t>=t_plot-tol,1);
    t_plot = t(t_i);
    
    figure;
    hold on;
    domj = 1:size(mdom,1);
    a=1;
    for j=domj
        dom=mdom{j};
        for i=1:size(dom,1) 
            G = dom{i}.G;
            plot(G(:,1),G(:,2));
            if(j==1 && i>1)
                b = a + size(G,1)-1;
                V_plot = y(a:b,t_i);
                plot3(G(:,1),G(:,2),V_plot);
                a = b+1;
            end
        end
    end
end

if 0
    if strcmp(direction,'l')
        ps = zeros(ncx,2);
        for i=1:ncx        
              ps(i,:) = [(i-0.5)*cl,0];                      
        end 
    elseif strcmp(direction,'t')
        ps = zeros(ncy,2);
        for i=1:ncy        
              ps(i,:) = [0,(i-0.5)*cw];
        end
    elseif strcmp(direction,'d')
        ps = zeros(2,2);
        ps(1,:) = [0,0];
        ps(2,:) = [ncx*cl,ncy*cw];
    end   
    V = zeros(size(y,2),size(ps,1));
    figure;
    hold on;
    for i=1:size(ps,1)
        d = transmembrane_G-ps(i,:);
        d = sqrt(d(:,1).^2+d(:,2).^2);
        [~,j] = min(d);
        p = transmembrane_G(j,:);
        V(:,i) = y(j,:);         
        plot(t,V(:,i),'LineWidth',1.5,'DisplayName',num2str(i));
    end
    title('$V(t)$','Interpreter','latex','FontSize',16);
    %legend('FontSize',16,'Interpreter','latex');
    %legend show;   
    %reduce precision to 0.1ms
    % n=round(0.1/t(2));
    % T_t=t(1:n:end);
    % T_V = V(1:n:end,:);
    % plot(T_t,T_V(:,1),'o');
    % hold on;
    % plot(T_t,T_V(:,2),'o');
    % T = table(T_t',T_V(:,1),T_V(:,2),'variablenames',{'t','Vclose','Vfar'});
    % writetable(T,[solname '_V.csv']);
end

if 0
    if strcmp(direction,'l')
%         pa = [0.5*cl,0]; %4.5 e 14.5
%         pb = [9.5*cl,0];
        pa = [0,cw/2];
        pb = [ncx*cl,cw/2];
        pa = [0,0];
        pb = [ncx*cl,0];
    elseif strcmp(direction,'t')
        pa = [0,1.5*cw];
        pb = [0,14.5*cw];
        pa = [0,0];
        pb = [0,ncy*cw];  
    elseif strcmp(direction,'d')
        pa = [0,0];
        pb = [ncx*cl,ncy*cw];    
    end

    CV=get_CV(folder,solname,pa,pb,V_th,old_file);
    %fprintf([folder solname]);
    fprintf(' CV = %f mm/ms (m/s)\n',CV);
end

