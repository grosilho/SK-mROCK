clear;
%clc;
%close all;

kappa = {'1e_4','2e_4','5e_4','1e_3','1e_2','1e_1'};
%kappa=kappa([2,5,10,15])

folder = '../results/ODECellByCellModel/';
n_states = 8; % 8 for LuoRudy, 13 for Fox, 21 for Courtemanche, 41 for Ohara Rudy
cl = 1.25e-2;

figure;
hold on;
for k=1:numel(kappa)
    solname = ['mRKC_dt_2e_2_dG_1e_3_kappa_' kappa{k}];
    run([folder solname '_geo.m']);
    
    transmembrane_G = [];
    dom=mdom{1};
    for i=2:size(dom,1) 
        transmembrane_G = [transmembrane_G;dom{i}.G];
    end

    n_pts = size(transmembrane_G,1);
    fileID = fopen([folder solname '_evolution.bin']);
    A = fread(fileID,'double');
    A = reshape(A,[1+n_pts*n_states,numel(A)/(1+n_pts*n_states)]);
    t = A(1,:);
    y = A(2:end,:);

    %for i=[11]
    i=11;
        ps = [(i-0.5)*cl,0];
        d = transmembrane_G-ps;
        d = sqrt(d(:,1).^2+d(:,2).^2);
        [~,j] = min(d);
        % p = transmembrane_G(j,:);
        V = y(j,:);         
        plot(t,V,'LineWidth',1.5,'DisplayName',['k = ' kappa{k}]);
    %end

    T = table(t',V','variablenames',{'t','V'});
    writetable(T,['V_for_k_' kappa{k} '.csv']);
end
title('$V(t)$','Interpreter','latex','FontSize',16);
legend show;
ax=gca;
