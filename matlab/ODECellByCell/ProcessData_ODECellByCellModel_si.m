clear;
%clc;
%close all;

folder = '../../results_cluster/ODECellByCellModel/nx=30/';
subfolder = 'CVl_C_si_cl_50/';

si_list = readtable([folder subfolder 'si_list.csv']);
Imax_list = readtable([folder subfolder 'Imax_list.csv']);

first_si = 1;
last_si= si_list.n(end);
si_list = si_list(first_si:last_si,:);

cl = 50e-4;
cw = 10e-4;
nc = 30;
np = 5;
pa = zeros(np,2);
pb = zeros(np,2);
for i=1:np
    pa(i,:) = [(7.5+i)*cl,0];
    pb(i,:) = [(17.5+i)*cl,0];
end
V_th = -20;

CV = zeros(size(si_list,1),1);
for i=1:size(si_list,1)
    solname = ['si_Imax_' num2str(si_list.n(i))];
    CV(i)=get_CV([folder subfolder],solname,pa,pb,V_th);
    fprintf(solname);
    fprintf(' si = %f, Imax = %f, CV  = %f mm/ms (m/s)\n',si_list.si(i),Imax_list.Imax(i),CV(i));
end


figure;
semilogx(si_list.si,CV);
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;
xlabel('$\sigma_i$','interpreter','latex')
ylabel('$CV$','interpreter','latex')

T=table(si_list.si,Imax_list.Imax,CV,'VariableNames',{'si','Imax','CV'});
writetable(T,[subfolder(1:end-1) '.csv']);


